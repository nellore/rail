"""
cojunction_enum_delegate.py 

Output of Bowtie 2 from cojunction_enum.py is streamed to this script to obtain
final output. See cojunction_enum.py for output format information.
"""

import sys
import os
import site
from collections import defaultdict
import random

base_path = os.path.abspath(
                    os.path.dirname(os.path.dirname(os.path.dirname(
                        os.path.realpath(__file__)))
                    )
                )
utils_path = os.path.join(base_path, 'rna', 'utils')
site.addsitedir(utils_path)
site.addsitedir(base_path)

from dooplicity.tools import xstream, register_cleanup
from dooplicity.counters import Counter
from alignment_handlers import multiread_with_junctions, \
    indels_junctions_exons_mismatches

counter = Counter('cojunction_enum_delegate')
register_cleanup(counter.flush)

import string
_reversed_complement_translation_table = string.maketrans('ATCG', 'TAGC')

def cojunction_length(cojunction):
    """ Computes number of exonic bases spanned by cojunction

        cojunction: list of junctions [(intron_start, intron_end, ...), ...]

        Return value: number of exonic bases spanned by cojunction
    """
    return sum([cojunction[i][0] - cojunction[i-1][1]
                    for i in xrange(1, len(cojunction))])

def paths_from_cojunctions(cojunctions, span=50):
    """ Finds junction combinations that can be overlapped by span bases

        Consider a directed acyclic graph (DAG) where each node is a different
        intron (intron_start, intron_end). There is an edge extending from a
        node A to a node B only if:
        1) their corresponding introns have no overlapping bases
        2) A's start and end positions are less than B's start and end
            positions
        3) there is no intron C between A and B such that C overlaps neither
            A nor B
        Such a graph is constructed from input cojunctions, or junction
        combinations that may be overlapped by a read. A cojunction is a
        linear subgraph of the DAG: if it has three junctions, it looks like
        this---

           A ---> B ---> C

        ---for its constituent node A, B, and C. It can be connected to the 
        rest of the DAG only via its source and/or its sink. All possible
        edges among sources and sinks of cojunctions are included in the graph.

        Assign a weight to each edge equal to the number of exonic bases
        between the introns it connects. Now use DFS to enumerate paths
        through the graph that span at least "span" nodes.

        cojunctions: a list of lists of tuples (intron start, intron end,
            left displacement, right displacement); each list of tuples
            represents a different cojunction. EACH COJUNCTION MUST BE 
            COORDINATE-SORTED. SIDE EFFECT: cojunctions ends
            up resorted
        span: the number of exonic bases (total weight) that can be spanned by
            a path through nodes to return

        Return value: list of junction combinations that can be spanned by
            span bases; junction combos with "edge" junctions removed are
            also included
    """
    counter.add('paths_from_cojunctions')
    cojunctions_count = len(cojunctions)
    counter.add('cojunctions', cojunctions_count)
    if cojunctions_count <= 1: return cojunctions
    '''Make each node of DAG a cojunction; cojunctions are linear subgraphs
    anyway'''
    DAG = defaultdict(set)
    reverse_DAG = defaultdict(set)
    cojunctions.sort(key=lambda x: (x[0][0], x[-1][1]))
    for i in xrange(cojunctions_count):
        source_cojunction_length = cojunction_length(cojunctions[i])
        for j in xrange(i+1, cojunctions_count):
            negative_separation = (
                    min(cojunctions[j][-1][1], cojunctions[i][-1][1])
                    - max(cojunctions[j][0][0], cojunctions[i][0][0])
                )
            if negative_separation >= 0:
                # Overlap between cojunctions, so continue
                continue
            if -negative_separation > span:
                # Too far away
                break
            if (-negative_separation
                + cojunction_length(cojunctions[j])
                + source_cojunction_length + 2) <= span:
                # span can span both cojunctions
                DAG[i].add(j)
                reverse_DAG[j].add(i)
            '''Add edges while there's overlap between jth cojunction and
            a successive cojunction'''
            for k in xrange(j+1, cojunctions_count):
                overlap = (
                        min(cojunctions[k][-1][1], cojunctions[j][-1][1])
                        - max(cojunctions[k][0][0], cojunctions[j][0][0])
                    )
                if overlap >= 0:
                    negative_separation = (
                        min(cojunctions[k][-1][1], cojunctions[i][-1][1])
                        - max(cojunctions[k][0][0], cojunctions[i][0][0])
                    )
                    if (-negative_separation
                        + cojunction_length(cojunctions[k])
                        + source_cojunction_length + 2) <= span:
                        counter.add('dag_edge')
                        DAG[i].add(k)
                        reverse_DAG[k].add(i)
                else:
                    break
    sources = set(DAG.keys()) - set(reverse_DAG.keys())
    sinks = set(reverse_DAG.keys()) - set(DAG.keys())
    # DFS to enumerate paths; start with isolated vertices
    paths = [[i] for i in xrange(cojunctions_count) if (
                                i not in sources and i not in sinks
                            )]
    for source in sources:
        for sink in sinks:
            stack = [(source, [source])]
            while stack:
                (vertex, path) = stack.pop()
                for next_vertex in DAG[vertex] - set(path):
                    if next_vertex == sink:
                        paths.append(path + [next_vertex])
                    else:
                        stack.append((next_vertex, path + [next_vertex]))
    counter.add('paths', len(paths))
    prereturn = [[j for i in [cojunctions[k] for k in path] for j in i]
                 for path in paths]
    to_return = set()
    # Add edge combos
    for junction_combo in prereturn:
        to_return.add(tuple(junction_combo[:]))
        to_return.add(tuple(junction_combo[1:-1]))
        to_return.add(tuple(junction_combo[1:]))
        to_return.add(tuple(junction_combo[:-1]))
    return [junction_combo for junction_combo in to_return if junction_combo]

def selected_cojunctions(cojunctions, max_refs=300,
                            seq='ATC', rname='chr1', sense='+'):
    """ Selects the max_ref cojunctions with the fewest junctions.

        Ties are broken at random.

        cojunctions: list of junction combinations
        max_refs: maximum number of cojunctions to return
        seq, rname, and sense are used for random seed if necessary

        Return value: truncated cojunctions list
    """
    if len(cojunctions) <= max_refs:
        return cojunctions
    cojunctions.sort(key=len)
    if len(cojunctions[max_refs]) == len(cojunctions[max_refs - 1]):
        # Shuffle cojunctions, then resort
        random.seed(seq + rname + sense)
        random.shuffle(cojunctions)
        cojunctions.sort(key=len)
    return cojunctions[:max_refs]

def go(input_stream=sys.stdin, output_stream=sys.stdout, fudge=5,
        stranded=False, verbose=False, max_refs=300, report_multiplier=1.2):
    """ Emits junction combinations associated with reads.

        Soft-clipped Bowtie 2 alignments of read sequences to the transcript
        fragment index are used infer which cojunctions could possibly be
        overlapped by reads. Then maximal cliques of the graph described in
        the maximal_cliques() function are enumerated to obtain which
        junction combinations could possibly be overlapped by reads.

        input_stream: where to retrieve Bowtie 2 output
        output_stream: where to emit exon and junction tuples; typically, this
            is sys.stdout.
        fudge: by how many bases to extend left and right extend sizes
            to accommodate potential indels
        stranded: True iff input reads are strand-specific; this affects
            whether an output partition has a terminal '+' or '-' indicating
            the sense strand. Further, if stranded is True, an alignment is
            returned only if its strand agrees with the junction's strand.
        verbose: True if alignments should occasionally be written to stderr.
        max_refs: maximum number of reference sequences to enumerate per read;
            if more are present, prioritize those sequences that overlap
            the fewest junctions
        report_multiplier: if verbose is True, the line number of an
            alignment written to stderr increases exponentially with base
            report_multiplier.
    """
    output_line_count, next_report_line, i = 0, 0, 0
    for (qname,), xpartition in xstream(input_stream, 1):
        '''While labeled multireadlet, this list may end up simply a
        unireadlet.'''
        multiread = []
        for tokens in xpartition:
            flag = int(tokens[0])
            if verbose and next_report_line == i:
                print >>sys.stderr, \
                    'SAM output record %d: rdname="%s", flag=%d' % (i,
                                                                    qname,
                                                                    flag)
                next_report_line = int((next_report_line + 1)
                                        * report_multiplier + 1) - 1
            i += 1
            multiread.append((qname,) + tokens)
        if flag & 4: continue
        cojunctions, all_junctions = defaultdict(set), {}
        for alignment in multiread_with_junctions(multiread, stranded):
            cigar = alignment[5]
            md = [field for field in alignment
                    if field[:5] == 'MD:Z:'][0][5:]
            pos = int(alignment[3])
            seq = alignment[9]
            reversed_complement_seq = seq[::-1].translate(
                    _reversed_complement_translation_table
                )
            if seq < reversed_complement_seq:
                seq_to_print = seq
            else:
                seq_to_print = reversed_complement_seq
            seq_size = len(seq)
            rname = alignment[2]
            sense = [field for field in alignment
                        if field[:5] == 'XS:A:'][0][5:]
            if (rname, sense) not in all_junctions:
                all_junctions[(rname, sense)] = defaultdict(list)
            _, _, junctions, _, _ = indels_junctions_exons_mismatches(
                                                cigar, md, pos, seq,
                                                junctions_only=True
                                            )
            cojunctions[(rname, sense)].add(
                    tuple([(junction[0], junction[1])
                                for junction in junctions])
                )
            for junction in junctions:
                if (junction[0], junction[1]) \
                    not in all_junctions[(rname, sense)]:
                    all_junctions[(rname, sense)][(junction[0], junction[1])] \
                        = [junction[2], junction[3]]
                else:
                    all_junctions[(rname, sense)][
                            (junction[0], junction[1])
                        ][0] = max(all_junctions[(rname, sense)][
                                (junction[0], junction[1])
                            ][0], junction[2])
                    all_junctions[(rname, sense)][
                            (junction[0], junction[1])
                        ][1] = max(all_junctions[(rname, sense)][
                                (junction[0], junction[1])
                            ][1], junction[3])
        for rname, sense in all_junctions:
            to_write = set()
            for cojunction in selected_cojunctions(paths_from_cojunctions(
                    list(cojunctions[(rname, sense)]), span=(seq_size + fudge)
                ), max_refs=max_refs, seq=seq, rname=rname, sense=sense):
                left_extend_size = all_junctions[(rname, sense)][
                                        cojunction[0]
                                    ][0]
                right_extend_size = all_junctions[(rname, sense)][
                                        cojunction[-1]
                                    ][1]
                to_write.add(('{rname}{sense}\t{starts}'
                       '\t{ends}\t{left_size}'
                       '\t{right_size}\t{seq}').format(
                            rname=rname,
                            sense=sense,
                            starts=','.join(
                                    [str(junction[0])
                                        for junction in cojunction]
                                ),
                            ends=','.join(
                                    [str(junction[1])
                                        for junction in cojunction]
                                ),
                            left_size=(left_extend_size
                                        + fudge),
                            right_size=(right_extend_size
                                        + fudge),
                            seq=seq_to_print
                       ))
            counter.add('paths_out', len(to_write))
            for line_to_write in to_write:
                print line_to_write
                output_line_count += 1
    output_stream.flush()
    print >>sys.stderr, ('cojunction_enum_delegate.py reports %d output lines.'
                            % output_line_count)

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__, 
                formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--verbose', action='store_const', const=True,
        default=False,
        help='Print out extra debugging statements')
    parser.add_argument('--fudge', type=int, required=False,
        default=5,
        help='Permits a sum of exonic bases for a junction combo to be '
             'within the specified number of bases of a read sequence\'s '
             'size; this allows for indels with respect to the reference')
    parser.add_argument('--max-refs', type=int, required=False,
        default=300,
        help='Hard limit on the number of reference sequences to emit '
             'per read per strand. Prioritizes reference sequences that '
             'overlap the fewest junctions')
    parser.add_argument(
        '--stranded', action='store_const', const=True, default=False,
        help='Assume input reads come from the sense strand; then partitions '
             'in output have terminal + and - indicating sense strand')
    parser.add_argument('--test', action='store_const', const=True,
        default=False,
        help='Run unit tests; DOES NOT NEED INPUT FROM STDIN')
    parser.add_argument('--report-multiplier', type=float, required=False,
        default=1.2,
        help='When --verbose is also invoked, the only lines of lengthy '
             'intermediate output written to stderr have line number that '
             'increases exponentially with this base')
    args = parser.parse_args()

if __name__ == '__main__' and not args.test:
    go(stranded=args.stranded, fudge=args.fudge,
        verbose=args.verbose, max_refs=args.max_refs,
        report_multiplier=args.report_multiplier)
elif __name__ == '__main__':
    # Test units
    del sys.argv[1:] # Don't choke on extra command-line parameters
    import unittest
    import shutil
    import tempfile

    class TestCojunctionLength(unittest.TestCase):
        """ Tests cojunction_length(). """

        def test_lengths(self):
            """ Fails if output of cojunction_length() is wrong. """
            self.assertEquals(
                    cojunction_length(
                            ((2, 5), (10, 100), (110, 150))
                        ), 15
                )
            self.assertEquals(
                    cojunction_length(
                            ((2, 5),)
                        ), 0
                )

    class PathsFromCojunctions(unittest.TestCase):
        """ Tests paths_from_cojunctions(). """

        def test_lengths(self):
            """ Fails if output of paths_from_cojunctions() is wrong. """
            # Test that cojunctions are merged for long-enough span
            paths = paths_from_cojunctions(
                            [((2, 5), (10, 100), (110, 150)),
                             ((160, 170), (190, 200))], span=50
                        )
            self.assertTrue(
                    ((2, 5), (10, 100), (110, 150), (160, 170), (190, 200))
                    in paths
                )
            # Test that cojunctions remain distinct for short span
            paths = paths_from_cojunctions(
                            [((2, 5), (10, 100), (110, 150)),
                             ((160, 170), (190, 200))], span=20
                        )
            self.assertTrue(
                    ((2, 5), (10, 100), (110, 150))
                    in paths
                )
            self.assertTrue(
                    ((160, 170), (190, 200))
                    in paths
                )
            # Test edge cases: cojunctions remain distinct for span 46
            paths = paths_from_cojunctions(
                            [((2, 5), (10, 100), (110, 150)),
                             ((160, 170), (190, 200))], span=46
                        )
            self.assertTrue(
                    ((2, 5), (10, 100), (110, 150))
                    in paths
                )
            self.assertTrue(
                    ((160, 170), (190, 200))
                    in paths
                )
            # ...but not for span 47
            paths = paths_from_cojunctions(
                            [((2, 5), (10, 100), (110, 150)),
                             ((160, 170), (190, 200))], span=47
                        )
            self.assertTrue(
                    ((2, 5), (10, 100), (110, 150), (160, 170), (190, 200))
                    in paths
                )
            # Test complicated cases
            paths = paths_from_cojunctions(
                            [((2, 5), (10, 100), (110, 150)),
                             ((10, 110), (123, 221)),
                             ((110, 150), (180, 210)),
                             ((220, 240),)], span=75
                        )
            for cojunction in [
                    ((2, 5), (10, 100), (110, 150)),
                    ((10, 100), (110, 150)),
                    ((10, 100),),
                    ((2, 5), (10, 100)),
                    ((110, 150), (180, 210), (220, 240)),
                    ((110, 150), (180, 210)),
                    ((180, 210), (220, 240)),
                    ((180, 210),),
                    ((10, 110), (123, 221)),
                    ((10, 110),),
                    ((123, 221),)
                ]:
                self.assertTrue(
                        cojunction in paths
                    )
                self.assertEquals(
                        len(paths), 11
                    )
            paths = paths_from_cojunctions(
                            [((2, 5), (10, 100), (110, 150)),
                             ((10, 110), (123, 221)),
                             ((110, 150), (180, 210)),
                             ((220, 240),)], span=200
                        )
            for cojunction in [
                    ((2, 5), (10, 100), (110, 150), (220, 240)),
                    ((10, 100), (110, 150), (220, 240)),
                    ((2, 5), (10, 100), (110, 150)),
                    ((10, 100), (110, 150)),
                    ((110, 150), (180, 210), (220, 240)),
                    ((110, 150), (180, 210)),
                    ((180, 210), (220, 240)),
                    ((180, 210),),
                    ((10, 110), (123, 221)),
                    ((10, 110),),
                    ((123, 221),)
                ]:
                self.assertTrue(
                        cojunction in paths
                    )
                self.assertEquals(
                        len(paths), 11
                    )

    class TestSelectedCojunctions(unittest.TestCase):
        """ Tests selected_cojunctions(). """

        def test_length(self):
            """ Fails if selected_cojunctions() returns wrong jx count """
            self.assertEquals(
                    len(selected_cojunctions(
                            [((2, 5), (10, 100), (110, 150)),
                             ((2, 5), (10, 30)),
                             ((2, 5), (10, 50))], max_refs=2)
                        ), 2
                )
            # Tie must be broken
            self.assertEquals(
                    len(selected_cojunctions(
                            [((2, 5), (10, 100), (110, 150)),
                             ((2, 5), (10, 30), (35, 39)),
                             ((2, 5), (10, 50))], max_refs=2)
                        ), 2
                )

    unittest.main()
