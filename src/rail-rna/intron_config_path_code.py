#!/usr/bin/env python
"""
Rail-RNA-intron_config
Follows Rail-RNA-intron
Precedes Rail-RNA-intron_fasta

Reduce step in MapReduce pipelines that outputs all possible configurations of
nonoverlapping introns on the same strand that a read of length max_len
(the maximum read length in the data being analyzed) + args.fudge can span.

Input (read from stdin)
----------------------------
Two kinds of input lines, one from Rail-RNA-align (max_len) and the other
from Rail-RNA-intron (intron):

(max_len)
1. Reference name (RNAME in SAM format) +
    '+' or '-' indicating which strand is the sense strand
2. Sample label
3. The character 'a'.
4. A maximum read length output by a Rail-RNA-align mapper.
5. The character '-'.

(intron)
1. Reference name (RNAME in SAM format) +
    '+' or '-' indicating which strand is the sense strand
2. Sample label
3. The character 'i', which will place the row after all 'a's (maximum read
    length rows) in lexicographic sort order.
4. Intron start position (inclusive)
5. Intron end position (exclusive)

Input is partitioned by strand/sample label (fields 1-2) and sorted by the
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
"""

import sys
import argparse
import time

_input_line_count, _output_line_count = 0, 0

def running_sum(an_iterable):
    """ Generates the running sum of the elements in an iterable.

        Input:  [a, b, c, d, ...]
        Output: [a, a+b, a+b+c, a+b+c+d, ...]

        Yield value: next item in output above.
    """
    total = 0
    for item in an_iterable:
        total += item
        yield total

def edges_and_extend_sizes_from_input_stream(input_stream, fudge=0):
    global _input_line_count
    partition, last_line, max_len_mode = None, None, False
    # Store pairs of successive nonoverlapping introns
    unlinked_nodes = set()
    linked_nodes = {}
    introns = {}
    index = 0
    while True:
        line = input_stream.readline()
        if not line:
            assert not (len(introns) and partition is None)
            # Yield final edges for strand
            for node in unlinked_nodes:
                current_intron = introns[node]
                yield partition + (current_intron,
                                    (current_intron[1] + extend_size, None))
            break
        else:
            if line == last_line: continue
            _input_line_count += 1
            tokens = line.rstrip().split('\t')
            assert len(tokens) == 5
            if tokens[2] == 'a':
                if not max_len_mode:
                    max_len = None
                    max_len_mode = True
                max_len = max(max_len, int(tokens[3]))
            elif tokens[2] == 'i':
                if max_len_mode:
                    # Yield final edges for strand
                    for node in unlinked_nodes:
                        current_intron = introns[node]
                        yield partition + (current_intron,
                                            (current_intron[1] + extend_size,
                                                None))
                    extend_size = max_len - 1 + fudge
                    # Yield extend_size, which denotes start of new partition
                    yield extend_size
                    unlinked_nodes = set()
                    linked_nodes = {}
                    introns = {}
                    index = 0
                    last_partition = partition
                    partition = (tokens[0], tokens[1])
                    assert partition != last_partition
                    intron_start, intron_end = int(tokens[3]), int(tokens[4])
                    introns[index] = (intron_start, intron_end)
                    unlinked_nodes.add(index)
                    # Yield first edge for strand (before first intron)
                    yield partition + ((None,
                                        max(intron_start - extend_size, 1)),
                                            (intron_start, intron_end))
                    index += 1
                    max_len_mode = False
                else:
                    partition = (tokens[0], tokens[1])
                    intron_start, intron_end = int(tokens[3]), int(tokens[4])
                    introns[index] = (intron_start, intron_end)
                    nodes_to_remove = []
                    for node in unlinked_nodes:
                        if intron_start > introns[node][1]:
                            nodes_to_remove.append(node)
                    for node in nodes_to_remove:
                        linked_nodes[node] = index
                        unlinked_nodes.remove(node)
                    unlinked_nodes.add(index)
                    nodes_to_remove = []
                    for node in linked_nodes:
                        intermediate_node = linked_nodes[node]
                        if intermediate_node in linked_nodes:
                            nodes_to_remove.append(node)
                        else:
                            yield partition + (introns[node],
                                                (intron_start, intron_end))
                            if introns[intermediate_node][1] > intron_end:
                                linked_nodes[node] = index
                    for node in nodes_to_remove:
                        del linked_nodes[node]
                        del introns[node]
                    index += 1
            else:
                raise RuntimeError('Invalid line: %s. A line should have '
                                   'either an "a" (denoting a max_len) or '
                                   'an "i" (denoting an intron) in its '
                                   'third field.' % line)
        last_line = line

def go(input_stream=sys.stdin, output_stream=sys.stdout, min_exon_size=25,
        verbose=False, fudge=1):
    """ Runs Rail-RNA-intron_config.

        Outputs all possible configurations of nonoverlapping introns on the
        same strand that a read of length max_len (the maximum read length in
        the data being analyzed) + fudge can span.

        Input (read from stdin)
        ----------------------------
        Two kinds of input lines, one from Rail-RNA-align (max_len) and the
        other from Rail-RNA-intron (intron):

        (max_len)
        1. Reference name (RNAME in SAM format) + 
            '+' or '-' indicating which strand is the sense strand
        2. The character 'a'.
        3. A maximum read length output by a Rail-RNA-align mapper.
        4. The character '-'.

        (intron)
        1. Reference name (RNAME in SAM format) + 
            '+' or '-' indicating which strand is the sense strand
        2. The character 'i', which will place the row after all 'a's (maximum
            read length rows) in lexicographic sort order.
        3. Intron start position (inclusive)
        4. Intron end position (exclusive)

        Input is partitioned by strand (field 1) and sorted by the remaining
        fields. INPUT COORDINATES ARE ASSUMED TO BE 1-INDEXED.

        Hadoop output (written to stdout)
        ----------------------------
        Tab-delimited tuple columns:
        1. Reference name (RNAME in SAM format) + 
            '+' or '-' indicating which strand is the sense strand
        2. Comma-separated list of intron start positions in configuration
        3. Comma-separated list of intron end positions in configuration
        4. extend_size: by how many bases on either side of an intron the
            reference should extend

        input_stream: where to get input
        output_stream: where to write output
        min_exon_size: if two introns are separated by min_exon_size bases,
            they are regarded as overlapping
        verbose: True iff extra debugging messages should be written to stderr
        fudge: a splice junction may be detected at any position along the read
            besides directly before or after it; thus, the sequences recorded
            in the new index should include max_read_size - 1 exonic bases
            before and after an intron in reference space (on the forward
            strand). These extensions are extended further by the number of
            bases fudge to accommodate possible small insertions

        No return value.
    """
    DAG = {}
    for edge_or_extend_size \
        in edges_and_extend_sizes_from_input_stream(input_stream, fudge=fudge):
        try:
            strand, _, start_node, end_node = edge_or_extend_size
            DAG[end_node] = start_node
            if not len(paths):
                paths.add((start_node, end_node))
                continue
            new_paths = set()
            for path in paths:
                if path[-1] == start_node:
                    new_paths.add(path + (end_node,))
                else:
                    new_paths.add(path)
            paths = new_paths
            print paths
        except TypeError:
            '''TypeError occurs when generator yields a new extend_size. Output
            last paths and start new strand.'''
            extend_size = edge_or_extend_size
            paths = set()

if __name__ == '__main__':
    # Print file's docstring if -h is invoked
    parser = argparse.ArgumentParser(description=__doc__, 
                formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(\
        '--verbose', action='store_const', const=True, default=False,
        help='Print out extra debugging statements')
    parser.add_argument('--test', action='store_const', const=True,
        default=False, help='Run unit tests')
    parser.add_argument('--min-exon-size', type=int, required=False,
        default=25, help='Minimum size of exons; two introns are considered '
                         'to overlap if the number of bases between them on '
                         'the reference is less than this value')
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
        min_exon_size=args.min_exon_size, verbose=args.verbose,
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
    from collections import defaultdict
    import os

    class TestGo(unittest.TestCase):
        """ Tests go(); needs no fixture. """
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
                        'chr1\t1\ta\t20\t-\n'
                        'chr1\t1\ti\t10\t50\n'
                        'chr1\t1\ti\t30\t70\n'
                        'chr1\t1\ti\t75\t101\n'
                        'chr1\t1\ti\t90\t1300\n'
                        'chr1\t1\ti\t91\t101\n'
                    )
                # Extend size is 20 above with fudge
            with open(self.output_file, 'w') as output_stream:
                with open(self.input_file) as input_stream:
                    go(input_stream=input_stream, output_stream=output_stream,
                        fudge=1)

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
            print intron_configs
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
                        'chr1\t1\ta\t15\t-\n'
                        'chr1\t1\ta\t20\t-\n'
                        'chr1\t1\ti\t11\t200\n'
                        'chr1\t1\ti\t31\t56\n'
                        'chr1\t1\ti\t75\t201\n'
                        'chr1\t1\ti\t91\t101\n'
                        'chr1\t1\ti\t205\t225\n'
                        'chr2\t1\ti\t21\t76\n'
                    )
                # Extend size is 20 above with fudge
            with open(self.output_file, 'w') as output_stream:
                with open(self.input_file) as input_stream:
                    go(input_stream=input_stream, output_stream=output_stream,
                        fudge=1)
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