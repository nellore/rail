#!/usr/bin/env python
"""
Rail-RNA-intron_filter

Follows Rail-RNA-intron_search
Precedes Rail-RNA-intron_config OR Rail-RNA-intron_collect

Reduce step in MapReduce pipelines that filters out introns if they are not
found in a specified percentage of samples or are not covered by a given number
of reads.

In the special mode --collect-introns, this step just collects and outputs
introns across samples, and no step should follow it.

Input (read from stdin)
----------------------------
Tab-delimited columns:
1. Reference name (RNAME in SAM format) +
    '+' or '-' indicating which strand is the sense strand
2. Intron start position (inclusive)
3. Intron end position (exclusive)
4. '\x1f'-separated list of sample indexes in which intron was found
5. '\x1f'-separated list of numbers of reads in which intron was found
    in respective sample specified by field 4

Input is partitioned by fields 1-3.

Hadoop output (written to stdout)
----------------------------
If --collect-introns is False:
Tab-delimited tuple columns:
1. Reference name (RNAME in SAM format) +
    '+' or '-' indicating which strand is the sense strand
2. Sample index
3. Intron start position (inclusive)
4. Intron end position (exclusive)

If --collect-introns is True:
Tab-delimited tuple columns:
1. Reference name (RNAME in SAM format) +
    '+' or '-' indicating which strand is the sense strand
2. Intron start position (inclusive)
3. Intron end position (exclusive)
4. comma-separated list of sample indexes in which intron was found
5. comma-separated list of numbers of reads in which intron was found
    in respective sample specified by field 4

ALL OUTPUT COORDINATES ARE 1-INDEXED.
"""
import os
import site
import sys
from collections import defaultdict
import time

base_path = os.path.abspath(
                    os.path.dirname(os.path.dirname(os.path.dirname(
                        os.path.realpath(__file__)))
                    )
                )
utils_path = os.path.join(base_path, 'rna', 'utils')
site.addsitedir(utils_path)
site.addsitedir(base_path)

from dooplicity.tools import xstream
import manifest

def go(manifest_object, input_stream=sys.stdin, output_stream=sys.stdout,
        sample_fraction=0.05, coverage_threshold=5, collect_introns=False,
        verbose=False):
    """ Runs Rail-RNA-intron_filter.

        Filters out every intron from input_stream that is not either:
          (1) in round(sample_fraction * (total number of samples)) samples OR
          (2) found in at least coverage_threshold reads in at least one
            sample.

        Input (read from stdin)
        ----------------------------
        Tab-delimited columns:
        1. Reference name (RNAME in SAM format) +
            '+' or '-' indicating which strand is the sense strand
        2. Intron start position (inclusive)
        3. Intron end position (exclusive)
        4. '\x1f'-separated list of sample indexes in which intron was found
        5. '\x1f'-separated list of numbers of reads in which intron was found
            in respective sample specified by field 4

        Input is partitioned by fields 1-3.

        Hadoop output (written to stdout)
        ----------------------------
        If --collect-introns is False:
        Tab-delimited tuple columns:
        1. Reference name (RNAME in SAM format) +
            '+' or '-' indicating which strand is the sense strand
        2. Sample index
        3. Intron start position (inclusive)
        4. Intron end position (exclusive)

        If --collect-introns is True:
        Tab-delimited tuple columns:
        1. Reference name (RNAME in SAM format) +
            '+' or '-' indicating which strand is the sense strand
        2. Intron start position (inclusive)
        3. Intron end position (exclusive)
        4. '\x1f'-separated list of sample indexes in which intron was found
        5. '\x1f'-separated list of numbers of reads in which intron was found
            in respective sample specified by field 4

        ALL OUTPUT COORDINATES ARE 1-INDEXED.

        input_stream: where to find input introns
        manifest_object: object of class LabelsAndIndices that maps indices
            to labels and back; used to count number of samples.
        sample_fraction: fraction of samples in which an intron must appear
            to pass filter if coverage_threshold criterion is not satisfied
        coverage_threshold: number of reads that must overlap intron in at
            least one sample to pass filter of sample_fraction criterion is not
            satisfied
        collect_introns: collects and outputs introns across samples; ignores
            sample_fraction and coverage_threshold
    """
    input_line_count, output_line_count = 0, 0
    min_sample_count = int(round(
            len(manifest_object.label_to_index) * sample_fraction
        ))
    for (rname_and_strand, pos, end_pos), xpartition in xstream(
                                                                input_stream, 3
                                                            ):
        sample_indexes = defaultdict(int)
        print_intron = False
        for current_sample_indexes, current_sample_counts in xpartition:
            input_line_count += 1
            current_sample_counts = current_sample_counts.split('\x1f')
            for i, sample_index in enumerate(
                                        current_sample_indexes.split('\x1f')
                                    ):
                sample_indexes[sample_index] += int(current_sample_counts[i])
        if collect_introns:
            samples_to_dump = sorted(sample_indexes.items(),
                                        key=lambda sample: sample[0])
            print >>output_stream, '%s\t%012d\t%012d\t%s\t%s' % (
                    rname_and_strand, int(pos), int(end_pos),
                    ','.join([sample[0] for sample in samples_to_dump]),
                    ','.join([str(sample[1]) for sample in samples_to_dump])
                )
            output_line_count += 1
        else:
            sample_count = len(sample_indexes)
            max_coverage = max(sample_indexes.values())
            if (sample_count >= min_sample_count
                or max_coverage >= coverage_threshold):
                for sample_index in sample_indexes:
                    print >>output_stream, '%s\t%s\t%012d\t%012d' % (
                            rname_and_strand, sample_index,
                            int(pos), int(end_pos)
                        )
                    output_line_count += 1
            elif verbose:
                print >>sys.stderr, (
                        'Intron (%s, %s, %s) filtered out; it appeared in %d '
                        'sample(s), and its coverage in any one sample did '
                        'not exceed %d.'
                    ) % (rname_and_strand, pos, end_pos,
                            sample_count, max_coverage)
    return input_line_count, output_line_count

if __name__ == '__main__':
    import argparse
    # Print file's docstring if -h is invoked
    parser = argparse.ArgumentParser(description=__doc__, 
                formatter_class=argparse.RawDescriptionHelpFormatter)
    # Add command-line arguments
    parser.add_argument('--verbose', action='store_const', const=True,
        default=False,
        help='Print out extra debugging statements')
    parser.add_argument('--test', action='store_const', const=True,
        default=False,
        help='Run unit tests; DOES NOT NEED INPUT FROM STDIN, AND DOES NOT '
             'WRITE TO STDOUT')
    parser.add_argument('--collect-introns', action='store_const', const=True,
        default=False,
        help=('Just collects and outputs unfiltered introns; overrides '
              '--sample-fraction and --coverage-threshold'))
    parser.add_argument('--sample-fraction', type=float, required=False,
        default=0.05,
        help=('An intron passes the filter if it is present in this fraction '
              'of samples; if the intron meets neither this criterion nor '
              'the --coverage-threshold criterion, it is filtered out'))
    parser.add_argument('--coverage-threshold', type=int, required=False,
        default=5,
        help=('An intron passes the filter if it is covered by at least this '
              'many reads in a single sample; if the intron meets neither '
              'this criterion nor the --sample-fraction criterion, it is '
              'filtered out'))

    # Add command-line arguments for dependencies
    manifest.add_args(parser)

    '''Now collect arguments. While the variable args declared below is
    global, properties of args are also arguments of the go() function so
    different command-line arguments can be passed to it for unit tests.'''
    args = parser.parse_args(sys.argv[1:])

if __name__ == '__main__' and not args.test:
    start_time = time.time()
    manifest_object = manifest.LabelsAndIndices(args.manifest)
    input_line_count, output_line_count = go(
            manifest_object=manifest_object,
            input_stream=sys.stdin,
            output_stream=sys.stdout,
            sample_fraction=args.sample_fraction,
            coverage_threshold=args.coverage_threshold,
            collect_introns=args.collect_introns,
            verbose=args.verbose
        )
    print >>sys.stderr, 'DONE with intron_filter.py; in/out=%d/%d; ' \
        'time=%0.3f s' % (input_line_count, output_line_count,
                            time.time() - start_time)
elif __name__ == '__main__':
    # Test units
    del sys.argv[1:] # Don't choke on extra command-line parameters
    import unittest
    import shutil
    import tempfile

    class TestGo(unittest.TestCase):
        """ Tests print_readletized_output(). """
        def setUp(self):
            # For storing output of print_readletized_output()
            self.temp_dir_path = tempfile.mkdtemp()
            self.input_file = os.path.join(self.temp_dir_path,
                                                    'sample_input.tsv')
            self.output_file = os.path.join(self.temp_dir_path,
                                                    'sample_output.tsv')
            self.manifest_file = os.path.join(self.temp_dir_path,
                                                    'sample.manifest')

        def test_output_1(self):
            """ Fails if output of go() is wrong. """
            manifest_content = (
                    'file1.fastq\t0\tfile2.fastq\t0\t0\n'
                    'file3.fastq\t0\tfile4.fastq\t0\t1\n'
                    'file5.fastq\t0\tfile6.fastq\t0\t2\n'
                )
            with open(self.manifest_file, 'w') as manifest_stream:
                manifest_stream.write(manifest_content)
            introns = (
                    'chr1+\t100\t140\t0\t6\n'
                    'chr2-\t171\t185\t1\x1f2\t2\x1f3\n'
                    'chr3+\t23\t85\t1\t2\n'
                )
            '''Above, for a coverage_threshold 5 and a sample_fraction 0.5,
            the introns on chr1+ and chr2- should pass while the intron on
            chr3+ should not.'''
            with open(self.input_file, 'w') as output_stream:
                output_stream.write(introns)
            manifest_object = manifest.LabelsAndIndices(self.manifest_file)
            with open(self.input_file) as input_stream, \
                open(self.output_file, 'w') as output_stream:
                input_line_count, output_line_count = go(
                    manifest_object=manifest_object,
                    input_stream=input_stream,
                    output_stream=output_stream,
                    sample_fraction=0.5,
                    coverage_threshold=5,
                    verbose=False
                )
            self.assertEquals(
                    input_line_count, 3
                )
            self.assertEquals(
                    output_line_count, 3
                )
            output_lines = []
            with open(self.output_file) as output_stream:
                for line in output_stream:
                    output_lines.append(line.strip())
            self.assertTrue(
                    'chr1+\t0\t%012d\t%012d' % (100, 140) in output_lines
                )
            self.assertTrue(
                    'chr2-\t1\t%012d\t%012d' % (171, 185) in output_lines
                )
            self.assertTrue(
                    'chr2-\t2\t%012d\t%012d' % (171, 185)in output_lines
                )
            self.assertEquals(
                    len(output_lines), 3
                )

        def test_output_2(self):
            """ Fails if output of go() is wrong. """
            manifest_content = (
                    'file1.fastq\t0\tfile2.fastq\t0\t0\n'
                    'file3.fastq\t0\tfile4.fastq\t0\t1\n'
                    'file5.fastq\t0\tfile6.fastq\t0\t2\n'
                    'file7.fastq\t0\tfile8.fastq\t0\t3\n'
                    'file9.fastq\t0\tfile10.fastq\t0\t4\n'
                )
            with open(self.manifest_file, 'w') as manifest_stream:
                manifest_stream.write(manifest_content)
            introns = (
                    'chr1+\t100\t140\t0\t6\n'
                    'chr2-\t171\t185\t1\x1f2\t2\x1f3\n'
                    'chr3+\t23\t85\t1\t2\n'
                    'chr3+\t23\t85\t1\t5\n'
                )
            '''Above, for a coverage_threshold 5 and a sample_fraction 0.5,
            the intron on chr2- should NOT pass while the rest should.'''
            with open(self.input_file, 'w') as output_stream:
                output_stream.write(introns)
            manifest_object = manifest.LabelsAndIndices(self.manifest_file)
            with open(self.input_file) as input_stream, \
                open(self.output_file, 'w') as output_stream:
                input_line_count, output_line_count = go(
                    manifest_object=manifest_object,
                    input_stream=input_stream,
                    output_stream=output_stream,
                    sample_fraction=0.5,
                    coverage_threshold=5,
                    verbose=False
                )
            self.assertEquals(
                    input_line_count, 4
                )
            self.assertEquals(
                    output_line_count, 2
                )
            output_lines = []
            with open(self.output_file) as output_stream:
                for line in output_stream:
                    output_lines.append(line.strip())
            self.assertTrue(
                    'chr1+\t0\t%012d\t%012d' % (100, 140) in output_lines
                )
            self.assertTrue(
                    'chr3+\t1\t%012d\t%012d' % (23, 85)in output_lines
                )
            self.assertEquals(
                    len(output_lines), 2
                )

        def tearDown(self):
            # Kill temporary directory
            shutil.rmtree(self.temp_dir_path)
    
    unittest.main()