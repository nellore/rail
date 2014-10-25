#!/usr/bin/env python
"""
Rail-RNA-intron_config

Follows Rail-RNA-intron
Precedes Rail-RNA-intron_fasta

Reduce step in MapReduce pipelines that outputs all possible configurations of
nonoverlapping introns on the same strand that a readlet of length readlet_size
(the maximum read length in the data being analyzed) + args.fudge can span,
given a minimum exon size min_exon_size.

Input (read from stdin)
----------------------------
Tab-delimited tuple columns:

1. Reference name (RNAME in SAM format) +
    '+' or '-' indicating which strand is the sense strand
2. Sample index
3. Intron start position (inclusive)
4. Intron end position (exclusive)

Input is partitioned by strand/sample index (fields 1-2) and sorted by the
remaining fields. INPUT COORDINATES ARE ASSUMED TO BE 1-INDEXED.

Hadoop output (written to stdout)
----------------------------
Tab-delimited tuple columns:
1. Reference name (RNAME in SAM format) + 
    '+' or '-' indicating which strand is the sense strand
2. Comma-separated list of intron start positions in configuration
3. Comma-separated list of intron end positions in configuration
4. left_extend_size: by how many bases on the left side of an intron the
    reference should extend
5. right_extend_size: by how many bases on the right side of an intron the
    reference should extend
6. By how many bases on the left side of an intron the reference COULD extend,
    or NA if beginning of strand
7. By how many bases on the right side of an intron the reference COULD extend,
    or NA if end of strand
"""

import sys
import argparse
import time
from collections import defaultdict
from collections import deque
import os
import site

base_path = os.path.abspath(
                    os.path.dirname(os.path.dirname(os.path.dirname(
                        os.path.realpath(__file__)))
                    )
                )
utils_path = os.path.join(base_path, 'rna', 'utils')
site.addsitedir(utils_path)
site.addsitedir(base_path)

from dooplicity.tools import xstream
_input_line_count, _output_line_count = 0, 0

def edges_from_input_stream(input_stream, readlet_size=20,
    min_overlap_exon_size=1):
    """ Generates edges of directed acyclic graph (DAG) of introns.

        A DAG is constructed for each strand. Each node of the DAG represents
        a unique intron and is labeled by the tuple (intron_start, intron_end),
        where intron_start is the (1-based) coordinate of the first base of the
        intron, and intron_end is the coordinate of the first base after the
        intron. An edge occurs between two introns A and B iff they do not
        overlap, and no intron C occurs between A and B such that A, B, and C 
        do not overlap. The intron with larger coordinates is the child of
        the intron with smaller coordinates. Weight each edge by the number of
        exonic bases between the introns it connects.

        The DAG has sources and sinks. Pad the DAG with new sources and sinks
        as follows. Extend an edge to each original source
        (source_start, source_end) from a new source labeled by the tuple
        (None, max(1, source_start - readlet_size + 1)). Extend an edge
        from each original sink (sink_start, sink_end) to a new sink labeled by
        the tuple (sink_end + source_start + readlet_size - 1, None).
        So each original source is assigned exactly one new source, and each
        original sink is assigned exactly one new sink.

        The paths through this DAG span all possible combinations of
        nonoverlapping introns on the strand. Finding all subpaths (sequences
        of introns), each of whose weights is <= readlet_size, redundantly
        enumerates all possible combinations of introns a readlet
        can overlap. Unfortunately, obtaining such combinations in
        "hot spots," where there are many alternative splicings and short
        exons, can become computationally intractable for large readlet sizes.
        To (optionally) control these combinatorial blow-ups, impose an
        effective minimum exon size min_overlap_exon_size by redefining overlap
        between two introns: two introns overlap if there are fewer than
        min_overlap_exon_size exonic bases between them.

        The algorithm for generating the edges of the DAG operates on an input
        stream whose lines are composed of the following tab-separated fields:

        1. Reference name (RNAME in SAM format) +
            '+' or '-' indicating which strand is the sense strand
        2. Sample index
        3. Intron start position (inclusive)
        4. Intron end position (exclusive)

        The input is partitioned by strand/sample index (fields 1-2) and sorted
        by the remaining fields. INPUT COORDINATES ARE ASSUMED TO BE 1-INDEXED.

        Introns are sorted by start position. To begin, the first set of
        mutually overlapping introns are connected to their new sources. Two
        data structures encode the graph as it is constructed [1]: the set
        unlinked_nodes, containing introns that do not yet have any children,
        and linked_nodes, a dictionary that, where possible, maps each intron A
        to its corresponding successive nonoverlapping intron B with the
        smallest end position read so far. With each new intron N read, the
        nodes in unlinked_nodes and linked_nodes are checked for whether N is
        their child, and these edges are yielded. (Note that by construction,
        nodes are streamed in topological order.) Nodes from unlinked_nodes may
        thus be promoted to linked_nodes. Each value node in linked_nodes is
        also replaced with N if N has a smaller end position. Then every node A
        in linked_nodes is checked for a path A -> B -> N. If such a path
        exists, an edge can NEVER connect A with any successive introns, and A
        is removed from linked_nodes. The algorithm continues until the end of
        the strand is reached, when the edges connecting remaining
        unlinked_nodes and their new sinks are yielded.

        [1] Technically, there are three data structures. unlinked_nodes and
        linked_nodes contain only indices of introns, and "introns" is
        a dictionary that maps indices to tuples (intron_start, intron_end).

        input_stream: where to find sorted introns of the form specified above.
        fudge: by how much a readlet_size should be extended.
        min_exon_size: minimum number of exonic bases between two introns
            for them to be considered nonoverlapping.

        Yield value: An edge tuple (strand,
                                    sample_index,
                                    (intron A start, intron A end),
                                    (intron B start, intron B end)) or None
                     at the beginning of a new partition.
    """
    global _input_line_count
    for key, xpartition in xstream(input_stream, 2, skip_duplicates=True):
        unlinked_nodes = set()
        for q, value in enumerate(xpartition):
            assert len(value) == 2
            _input_line_count += 1
            if not q:
                # Denote start of new partition
                yield None
                # Handle first intron from partition
                intron_start, intron_end = int(value[0]), int(value[1])
                # Create fake source before first intron
                fake_source = (
                        None,
                        max(intron_start - (readlet_size - 1), 1)
                    )
                introns = {
                        0 : fake_source,
                        1 : (intron_start, intron_end)
                    }
                linked_nodes = { 0 : 1 }
                unlinked_nodes = set([1])
                index = 2
                # Yield first edge for strand (before first intron)
                yield key + (fake_source,
                                (intron_start, intron_end))
                first_intron = False
            else:
                # Handle next introns from partition
                intron_start, intron_end = int(value[0]), int(value[1])
                introns[index] = (intron_start, intron_end)
                nodes_to_trash = []
                for node in unlinked_nodes:
                    if intron_start >= introns[node][1] + \
                        min_overlap_exon_size:
                        nodes_to_trash.append(node)
                for node in nodes_to_trash:
                    linked_nodes[node] = index
                    unlinked_nodes.remove(node)
                unlinked_nodes.add(index)
                nodes_to_trash = []
                for node in linked_nodes:
                    intermediate_node = linked_nodes[node]
                    if intermediate_node in linked_nodes:
                        nodes_to_trash.append(node)
                    else:
                        yield key + (introns[node],
                                            (intron_start, intron_end))
                        if introns[intermediate_node][1] > intron_end:
                            linked_nodes[node] = index
                for node in nodes_to_trash:
                    del linked_nodes[node]
                    del introns[node]
                index += 1
        # Yield final edges for strand
        for node in unlinked_nodes:
            current_intron = introns[node]
            yield key + (current_intron, (current_intron[1]
                                            + readlet_size - 1, None))

def paths(graph, source, in_node, readlet_size, last_node, edge_span=2,
    min_edge_span_size=25, can_yield=False):
    """ Generates intron combos spanning readlet_size exonic bases.

        The algorithm is nonrecursive to ensure Python doesn't choke. Consider
        all paths through "graph" that originate at "source" and pass through
        "in_node." This generator yields all possible maximal subpaths
        (possibly repetitively) (source, in_node, ... , next_to_last_node,
        last_node), each of which represents a sequence of introns beginning at
        in_node and ending at next_to_last_node that a readlet of length
        readlet_size can overlap. But there is a caveat: if an intron sequence
        itself has a subpath composed of "edge_span" number of edges whose
        total weight is not at least min_edge_span_size, that path is
        suppressed. This controls blowups around many short exons/alternative
        splicings. source and last_node are included in yielded paths to
        determine by how many bases on either side of the intron sequence the
        reference should be extended.

        graph: a dictionary. Each key is an out node (intron_start, intron_end)
            of the graph, and its corresponding value is a list of in nodes
            to which it connects.
        source, in_node: an edge that terminates on in_node originates at
            source.
        readlet_size: maximum readlet size
        edge_span, min_edge_span_size: if an intron sequence from a
            maximal subpath itself has a subpath composed of edge_span number
            of edges whose total weight is not at least min_edge_span_size,
            that path is suppressed

        Yield value: a list of node tuples representing a path or None if
            a path is found to violate the condition involving edge_span and
            min_edge_span_size described above.
    """
    assert isinstance(edge_span, int) and edge_span >= 1, \
        'Edge span must be integer >= 1; was %d' % edge_span
    if not can_yield:
        '''Ensure in_node has all its children; that is, ensure in_node has at
        least one grandchild.'''
        all_kids = False
        for child in graph[in_node]:
            if child in graph:
                all_kids = True
                break
        if not all_kids: return
        '''Ensure stream is past point where a given child of in_node could
        connect to another node that gives rise to a maximal path.'''
        for child in graph[in_node]:
            if child[0] - in_node[1] + last_node[0] - child[1] \
                < readlet_size - 1:
                return
    path, base_sum  = [source], 0
    path_queue = deque([(in_node, base_sum, path)])
    while path_queue:
        in_node, base_sum, path = path_queue.popleft()
        path = path + [in_node]
        node_count = len(path)
        yielded = False
        if node_count >= 3:
            '''When the path spans at least three nodes, the intron sequence
            spans at least one node. If the path weight is >= readlet_size - 1,
            the terminal node cannot possibly be overlapped by a readlet
            overlapping the start node.'''
            base_sum += path[-1][0] - path[-2][1]
            if base_sum >= readlet_size - 1:
                yield path
                yielded = True
            elif node_count >= edge_span + 2 and \
                sum([path[-i][0] - path[-(i+1)][1]
                        for i in xrange(1, edge_span+1)]) \
                < min_edge_span_size:
                '''If the last edge_span edges of the path being constructed
                have a weight less than min_edge_span_size, a maximal subpath
                is forbidden. Yield None.'''
                yield None
                yielded = True
        if not yielded:
            if in_node in graph:
                for node in graph[in_node]:
                    path_queue.append((node, base_sum, path))
            else:
                '''in_node is not in the graph, and the end of a path has been
                reached. Any nodes added to the path will only give rise to an
                extension to the right of in_node by readlet_size - 1 bases.
                Tack on a fake final node that will give rise to this
                extension.'''
                try:
                    yield path + [(path[-1][1] + readlet_size - 1, None)]
                except TypeError:
                    # Path terminates on final fake intron; don't do anything
                    assert path[-1][1] is None

def consume_graph_and_print_combos(DAG, reverse_DAG, readlet_size, strand,
    last_node, output_stream, edge_span=2, min_edge_span_size=25,
    full_graph=False):
    """ Consumes graph, printing intron combos that can be overlapped by reads.

        See edges_from_input_stream()'s docstring for a detailed description of
        the directed acylic graph (DAG).

        To enumerate all possible combinations of introns that can be
        overlapped by a readlet, consider each node separately, and find the
        ways the intron represented by that node can be the first intron on a
        read; that is, walk every path starting from that node edge by edge
        until its weight exceeds readlet_size, the maximal readlet length.

        It is not necessary to keep the entire graph in memory, however. The 
        DAG has one or more sources. A source S is removed from the graph
        when all paths originating at every child node C_i of S have been
        reviewed to find intron sequences that can be overlapped by
        readlet_size exonic bases; then S is no longer needed to find by how
        many exonic bases to the left of C_i the reference should extend.
        More specifically, an edge from S to C_i is removed iff:
        1) All edges that start at C_i have been generated. By construction,
        this occurs if C_i has at least one grandchild.
        2) Every path of maximal length originating at C_i whose weight is
        < readlet_size - 1 can be constructed. If {G_j} are the children of C,
        this occurs if for every j, (the difference between the genomic start
        position of the last node (intron) streamed and the genomic end
        position of G_j) >= readlet_size - 1 - (weight of edge between
        C_i and G_j).

        Once sources are removed, the graph has new sources that can also
        be "consumed." So the DAG can alternately be generated, making it
        expand towards the right end of a strand, and consumed, making it
        retreat from the left end of the strand. The code contained here also
        prints intron combinations for maximal paths as described below:

        Tab-delimited tuple columns:
        1. Reference name (RNAME in SAM format) + 
            '+' or '-' indicating which strand is the sense strand
        2. Comma-separated list of intron start positions in configuration
        3. Comma-separated list of intron end positions in configuration
        4. left_extend_size: by how many bases on the left side of an intron
            the reference should extend
        5. right_extend_size: by how many bases on the right side of an intron
            the reference should extend

        DAG: a dictionary. Each key is a parent node (intron_start, intron_end)
            of the graph, and its corresponding value is the set of its child
            nodes.
        reverse_DAG: a dictionary. Each key is a child node
            (intron_start, intron_end), and its corresponding value is the set
            of its parent nodes.
        readlet_size: maximum readlet size
        strand: current strand
        last_node: child node from last edge added to graph. Used to
            determine of a source can be trashed.
        output_stream: where to write output
        edge_span, min_edge_span_size: parameters used by paths() function.
            See its docstring for more information.
        full_graph: True iff there are no more nodes to stream on the graph.
            Used to determine if the rest of the graph can be consumed.

        No return value.

        NOTE: THIS FUNCTION HAS SIDE EFFECTS. DAG and reverse_DAG are altered.
    """
    global _output_line_count
    source_queue = deque([node for node in DAG if node not in reverse_DAG])
    while source_queue:
        source = source_queue.popleft()
        nodes_to_remove = []
        for m, node in enumerate(DAG[source]):
            i = -1
            source_node_weight = min(node[0] - source[1], readlet_size - 1)
            parents_to_remove = []
            for parent in reverse_DAG[node]:
                if parent == source: continue
                if parent in source_queue \
                    and min(node[0] - parent[1], readlet_size - 1) \
                    <= source_node_weight:
                    '''Optimization: trash an edge between node and a different
                    parent early if it would give rise to the same paths with
                    shorter leftward extensions.'''
                    parents_to_remove.append(parent)
                    DAG[parent].remove(node)
            for parent in parents_to_remove:
                reverse_DAG[node].remove(parent)
            for i, path in enumerate(
                                paths(DAG, source, node, readlet_size,
                                        last_node,
                                        edge_span=edge_span,
                                        min_edge_span_size=min_edge_span_size,
                                        can_yield=full_graph)
                            ):
                try:
                    node_count = len(path)
                    left_size = path[1][0] - path[0][1]
                    right_size = path[-1][0] - path[-2][1]
                    print >>output_stream, '%s\t%s\t%s\t%d' \
                        '\t%d\t%s\t%s' % (
                        strand,
                        ','.join([str(path[k][0])
                                  for k in xrange(1, node_count - 1)]),
                        ','.join([str(path[k][1])
                                  for k in xrange(1, node_count - 1)]),
                        min(readlet_size - 1, left_size),
                        min(readlet_size - 1, right_size),
                        str(left_size) if path[0][0] is not None else 'NA',
                        str(right_size) if path[-1][1] is not None else 'NA'
                    )
                    _output_line_count += 1
                    sys.stdout.flush()
                except TypeError:
                    # Path is verboten
                    assert path is None
            if i != -1:
                '''paths() yielded something, which means all of node's
                children in the graph have been probed, no new children of node
                will be obtained from the stream, and no new maximal paths
                from node are available. Trash the edge from source to node.'''
                reverse_DAG[node].remove(source)
                if not len(reverse_DAG[node]):
                    # node is now a source
                    source_queue.append(node)
                    del reverse_DAG[node]
                nodes_to_remove.append(node)
        for node in nodes_to_remove:
            DAG[source].remove(node)
        if not len(DAG[source]):
            # Edges from source are no longer needed to construct extensions
            del DAG[source]

def go(input_stream=sys.stdin, output_stream=sys.stdout, readlet_size=20,
        min_overlap_exon_size=1, edge_span=2, min_edge_span_size=25, 
        verbose=False, fudge=0, flush_base_count=10000000):
    """ Runs Rail-RNA-intron_config.

        Reduce step in MapReduce pipelines that outputs all possible
        configurations of nonoverlapping introns on the same strand that a
        readlet of length readlet_size (the maximum readlet length in the data
        being analyzed) + args.fudge can span, given a minimum exon size
        min_exon_size.

        The code contained here switches between generating and consuming
        the DAG/printing intron combos. The switch occurs every
        flush_base_count bases along a given strand. See docstrings for
        functions as well as comments for more information about the algorithm.

        Input (read from stdin)
        ----------------------------
        Tab-delimited input tuple columns:

        1. Reference name (RNAME in SAM format) + 
            '+' or '-' indicating which strand is the sense strand
        2. Sample index
        3. Intron start position (inclusive)
        4. Intron end position (exclusive)

        Input is partitioned by strand+sample index (fields 1-2) and sorted by
        the remaining fields. INPUT COORDINATES ARE ASSUMED TO BE 1-INDEXED.

        Hadoop output (written to stdout)
        ----------------------------
        Tab-delimited tuple columns:
        1. Reference name (RNAME in SAM format) + 
            '+' or '-' indicating which strand is the sense strand
        2. Comma-separated list of intron start positions in configuration
        3. Comma-separated list of intron end positions in configuration
        4. left_extend_size: by how many bases on the left side of an intron
            the reference should extend
        5. right_extend_size: by how many bases on the right side of an intron
            the reference should extend
        6. By how many bases on the left side of an intron the reference COULD
            extend, or NA if beginning of strand
        7. By how many bases on the right side of an intron the reference COULD
            extend, or NA if end of strand

        input_stream: where to get input
        output_stream: where to write output
        min_overlap_exon_size: if two introns are separated by
            min_overlap_exon_size bases, they are regarded as overlapping
        edge_span, min_edge_span_size: see paths() and
            consume_graph_and_print_combos() for information about these
            parameters.
        verbose: True iff extra debugging messages should be written to stderr
        fudge: a splice junction may be detected at any position along the read
            besides directly before or after it; thus, the sequences recorded
            in the new index should include max_read_size - 1 exonic bases
            before and after an intron in reference space (on the forward
            strand). These extensions are extended further by the number of
            bases fudge to accommodate possible small insertions
        flush_base_count: algorithm switches between generating and consuming
            the graph every flush_base_count bases along the strand.

        No return value.
    """
    effective_readlet_size = readlet_size + fudge
    for edge in edges_from_input_stream(
                        input_stream, 
                        readlet_size=effective_readlet_size,
                        min_overlap_exon_size=min_overlap_exon_size
                    ):
        try:
            strand, sample_index, start_node, end_node = edge
        except TypeError:
            try:
                if verbose:
                    print >>sys.stderr, \
                        ('Consuming rest of graph on strand ' + strand +
                         ' for sample ' + sample_index)
                    print >>sys.stderr, 'Before consumption, DAG has %d ' \
                        'nodes, and reverse DAG has %d nodes.' \
                        % (len(DAG), len(reverse_DAG))
                    consume_start_time = time.time()
                consume_graph_and_print_combos(
                        DAG, 
                        reverse_DAG, 
                        effective_readlet_size,
                        strand,
                        end_node,
                        output_stream,
                        edge_span=edge_span,
                        min_edge_span_size=min_edge_span_size,
                        full_graph=True
                    )
                if verbose:
                    print >>sys.stderr, 'After consumption, DAG has %d ' \
                    'nodes, and reverse DAG has %d nodes.' \
                        % (len(DAG), len(reverse_DAG))
                    print >>sys.stderr, 'Time taken: %0.3f s' \
                        % (time.time() - consume_start_time)
            except NameError: pass
            DAG, reverse_DAG = defaultdict(set), defaultdict(set)
            flush_threshold = flush_base_count
            continue
        DAG[start_node].add(end_node)
        reverse_DAG[end_node].add(start_node)
        if end_node[0] >= flush_threshold:
            '''Ideally, one would use a memlimit, but can't take control of
            memory management in Python. In any event, the memory footprint
            of tends to be so small (just a few MB); all bets are off, however,
            if there's blow-up.'''
            if verbose:
                print >>sys.stderr, \
                    'Consuming graph up to end node', end_node, 'on ' \
                    'strand', strand, 'for sample', sample_index
                print >>sys.stderr, 'Before consumption, DAG has %d ' \
                    'nodes, and reverse DAG has %d nodes.' \
                        % (len(DAG), len(reverse_DAG))
                consume_start_time = time.time()
            consume_graph_and_print_combos(
                        DAG, 
                        reverse_DAG,
                        effective_readlet_size,
                        strand,
                        end_node,
                        output_stream,
                        edge_span=edge_span,
                        min_edge_span_size=min_edge_span_size,
                    )
            if verbose:
                print >>sys.stderr, 'After consumption, DAG has %d ' \
                    'nodes, and reverse DAG has %d nodes.' \
                        % (len(DAG), len(reverse_DAG))
                print >>sys.stderr, 'Time taken: %f s' \
                        % (time.time() - consume_start_time)
            flush_threshold += flush_base_count
    # End of stream reached; consume rest of graph from last strand
    try:
        if verbose:
            print >>sys.stderr, \
                ('Consuming rest of graph on strand ' + strand +
                 ' for sample ' + sample_index)
            print >>sys.stderr, 'Before consumption, DAG has %d ' \
                'nodes, and reverse DAG has %d nodes.' \
                % (len(DAG), len(reverse_DAG))
        consume_graph_and_print_combos(
                DAG, 
                reverse_DAG, 
                effective_readlet_size, 
                strand,
                end_node,
                output_stream,
                edge_span=edge_span,
                min_edge_span_size=min_edge_span_size,
                full_graph=True
            )
        if verbose:
            print >>sys.stderr, 'After consumption, DAG has %d ' \
            'nodes, and reverse DAG has %d nodes.' \
                % (len(DAG), len(reverse_DAG))
    except NameError: pass

if __name__ == '__main__':
    # Print file's docstring if -h is invoked
    parser = argparse.ArgumentParser(description=__doc__, 
                formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(\
        '--verbose', action='store_const', const=True, default=False,
        help='Print out extra debugging statements')
    parser.add_argument('--test', action='store_const', const=True,
        default=False, help='Run unit tests')
    parser.add_argument('--min-overlap-exon-size', type=int, required=False,
        default=1, help='Minimum size of exons; two introns are considered '
                        'to overlap if the number of bases between them on '
                        'the reference is less than this value')
    parser.add_argument('--edge-span', type=int, required=False,
        default=1, help='--edge-span edges of path from DAG representing an '
                        'intron combination must have weight >= '
                         '--min-edge-span-size to be output')
    parser.add_argument('--min-edge-span-size', type=int, required=False,
        default=1, help='--edge-span edges of path from DAG representing an '
                         'intron combination must have weight >= '
                         '--min-edge-span-size to be output')
    parser.add_argument('--readlet-size', type=int, required=False,
        default=20, help='Size of readlets to be aligned to reference '
                         'including only exonic bases around introns')
    parser.add_argument('--fudge', type=int, required=False, default=0, 
        help='A splice junction may be detected at any position along the '
             'read besides directly before or after it; thus, the sequences '
             'recorded in the new index should include max_read_size - 1 '
             'exonic bases before and after an intron in reference space (on '
             'the forward strand). These extensions are extended further by '
             'the number of bases --fudge to accommodate possible small '
             'insertions.')
    
    args = parser.parse_args(sys.argv[1:])

if __name__ == '__main__' and not args.test:
    start_time = time.time()
    go(input_stream=sys.stdin, output_stream=sys.stdout,
        min_overlap_exon_size=args.min_overlap_exon_size,
        edge_span=args.edge_span,
        min_edge_span_size=args.min_edge_span_size,
        readlet_size=args.readlet_size,
        verbose=args.verbose,
        fudge=args.fudge)
    print >>sys.stderr, 'DONE with intron_config.py; in/out=%d/%d; ' \
                        'time=%0.3f s' % (_input_line_count, 
                                            _output_line_count,
                                            time.time() - start_time)
elif __name__ == '__main__':
    # Test units
    del sys.argv[1:] # Don't choke on extra command-line parameters
    import unittest
    import tempfile
    import shutil

    class TestGo(unittest.TestCase):
        """ Tests go(). """
        def setUp(self):
            # Set up temporary directory
            self.temp_dir_path = tempfile.mkdtemp()
            self.input_file = os.path.join(self.temp_dir_path, 'introns.temp')
            self.output_file = os.path.join(self.temp_dir_path, 'configs.temp')

        def test_overlapping_intron_configuration_1(self):
            """ Fails if intron configurations are not enumerated properly. """
            with open(self.input_file, 'w') as input_stream:
                '''Recall that input is partitioned by first column and 
                sorted by the next three columns.'''
                input_stream.write(
                        'chr1\t1\ti\t10\t50\n'
                        'chr1\t1\ti\t30\t70\n'
                        'chr1\t1\ti\t75\t101\n'
                        'chr1\t1\ti\t90\t1300\n'
                        'chr1\t1\ti\t91\t101\n'
                    )
            with open(self.output_file, 'w') as output_stream:
                with open(self.input_file) as input_stream:
                    go(input_stream=input_stream, output_stream=output_stream,
                        fudge=1, edge_span=1, min_edge_span_size=1,
                        readlet_size=20)
            '''Read output; store configurations as frozen sets so there are
            no duplicate configurations.'''
            intron_configs = set()
            with open(self.output_file) as result_stream:
                for line in result_stream:
                    tokens = line.strip().split('\t')
                    intron_configs.add(frozenset(zip(
                            [int(pos) for pos in tokens[2].split(',')],
                            [int(end_pos) for end_pos in tokens[3].split(',')]
                        )))
            self.assertEqual(
                    set([
                        frozenset([(10, 50)]),
                        frozenset([(30, 70), (75, 101)]),
                        frozenset([(30, 70)]),
                        frozenset([(75, 101)]),
                        frozenset([(90, 1300)]),
                        frozenset([(91, 101)]),
                    ]), intron_configs
                )

        def test_overlapping_intron_configuration_2(self):
            """ Fails if intron configurations are not enumerated properly. """
            with open(self.input_file, 'w') as input_stream:
                '''Recall that input is partitioned by first column and 
                sorted by the next three columns.'''
                input_stream.write(
                        'chr1\t1\ti\t11\t200\n'
                        'chr1\t1\ti\t31\t56\n'
                        'chr1\t1\ti\t75\t201\n'
                        'chr1\t1\ti\t91\t101\n'
                        'chr1\t1\ti\t205\t225\n'
                        'chr2\t1\ti\t21\t76\n'
                    )
            with open(self.output_file, 'w') as output_stream:
                with open(self.input_file) as input_stream:
                    go(input_stream=input_stream, output_stream=output_stream,
                        fudge=1, readlet_size=20)
            '''Read output; store configurations as frozen sets so there are
            no duplicate configurations.'''
            intron_configs = defaultdict(set)
            with open(self.output_file) as result_stream:
                for line in result_stream:
                    tokens = line.strip().split('\t')
                    intron_configs[tokens[1]].add(frozenset(zip(
                            [int(pos) for pos in tokens[2].split(',')],
                            [int(end_pos) for end_pos in tokens[3].split(',')]
                        )))
            self.assertEqual(
                    set([
                        frozenset([(11, 200), (205, 225)]),
                        frozenset([(31, 56), (75, 201)]),
                        frozenset([(31, 56)]),
                        frozenset([(75, 201), (205, 225)]),
                        frozenset([(91, 101)]),
                        frozenset([(205, 225)])
                    ]), intron_configs['chr1']
                )
            self.assertEqual(
                    set([
                        frozenset([(21, 76)])
                    ]), intron_configs['chr2']
                )

        def tearDown(self):
            # Kill temporary directory
            shutil.rmtree(self.temp_dir_path)
    unittest.main()
