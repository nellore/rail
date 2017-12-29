#!/usr/bin/env python
"""
Rail-RNA-bed_pre
Follows Rail-RNA-realign_reads
Precedes Rail-RNA-bed / Rail-RNA-tsv

Reducer for MapReduce pipelines that sums coverage of indels/junctions and
computes other statistics by sample from Rail-RNA-realign.

Input (read from stdin)
----------------------------
Tab-delimited input tuple columns:
1. 'I', 'D', or 'N' for insertion, deletion, or junction line
2. Number string representing RNAME
3. Start position (Last base before insertion, first base of deletion,
                    or first base of intron)
4. End position (Last base before insertion, last base of deletion (exclusive),
                    or last base of intron (exclusive))
5. '+' or '-' indicating which strand is the sense strand for junctions,
   inserted sequence for insertions, or deleted sequence for deletions
6. Sample index
----Next fields are for junctions only; they are '\x1c' for indels----
7. Number of nucleotides between 5' end of intron and 5' end of read from which
    it was inferred, ASSUMING THE SENSE STRAND IS THE FORWARD STRAND. That is,
    if the sense strand is the reverse strand, this is the distance between the
    3' end of the read and the 3' end of the intron.
8. Number of nucleotides between 3' end of intron and 3' end of read from which
    it was inferred, ASSUMING THE SENSE STRAND IS THE FORWARD STRAND.
--------------------------------------------------------------------
9. Number of instances of junction, insertion, or deletion in sample; this is
    always +1 before bed_pre combiner/reducer

Input is partitioned by fields 1-5 and sorted by field 6.

Hadoop output (written to stdout)
----------------------------
Tab-delimited output tuple columns (bed):
1. 'I', 'D', or 'N' for insertion, deletion, or junction line
2. Sample index
3. Number string representing RNAME (+ '+ or -' if junction; same as field 6)
4. Start position (Last base before insertion, first base of deletion,
                    or first base of intron)
5. End position (Last base before insertion, last base of deletion (exclusive),
                    or last base of intron (exclusive))
6. '+' or '-' indicating which strand is the sense strand for junctions,
   inserted sequence for insertions, or deleted sequence for deletions
----Next fields are for junctions only; they are '\x1c' for indels----
7. MAX number of nucleotides between 5' end of intron and 5' end of read from
    which it was inferred, ASSUMING THE SENSE STRAND IS THE FORWARD STRAND.
    That is, if the sense strand is the reverse strand, this is the distance
    between the 3' end of the read and the 3' end of the intron.
8. MAX number of nucleotides between 3' end of intron and 3' end of read from
    which it was inferred, ASSUMING THE SENSE STRAND IS THE FORWARD STRAND.
9. MAXIMIN (number of nucleotides between 5' end of intron and 5' end of read,
            number of nucleotides between 3' end of intron and 3' end of read);
   min is between the args above; max is across reads.

Tab-delimited output tuple columns (collect)
1. '0' if insertion, '1' if deletion, or '2' if junction line
2. Number string representing RNAME (+ '+ or -' if junction; same as field 6)
3. Start position (Last base before insertion, first base of deletion,
                    or first base of intron)
4. End position (Last base before insertion, last base of deletion (exclusive),
                    or last base of intron (INCLUSIVE HERE))
5. '+' or '-' indicating which strand is the sense strand for junctions,
   inserted sequence for insertions, or deleted sequence for deletions
6. Coverage of feature for sample with index N
...
N + 6. Coverage of feature in sample with index N
--------------------------------------------------------------------
10. SUMMED number of instances of junction, insertion, or deletion in sample

OUTPUT COORDINATES ARE 1-INDEXED.
"""
import os
import sys
import argparse
import site
import time
import itertools

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
import manifest

counter = Counter('bed_pre')
register_cleanup(counter.flush)

def go(manifest_object, input_stream=sys.stdin, output_stream=sys.stdout,
        sample_fraction=0.05, coverage_threshold=5, verbose=False):
    """ Runs Rail-RNA-bed_pre

        Writes indels and junctions for outputting BEDs by sample and
        TSVs across samples.

        Input (read from stdin)
        ----------------------------
        Tab-delimited input tuple columns:
        1. 'I', 'D', or 'N' for insertion, deletion, or junction line
        2. Number string representing RNAME
        3. Start position (Last base before insertion, first base of deletion,
                            or first base of intron)
        4. End position (Last base before insertion, last base of deletion
                            (exclusive), or last base of intron (exclusive))
        5. '+' or '-' indicating which strand is the sense strand for
            junctions, inserted sequence for insertions, or deleted sequence
            for deletions
        6. Sample index
        ----Next fields are for junctions only; they are '\x1c' for indels----
        7. Number of nucleotides between 5' end of intron and 5' end of read
            from which it was inferred, ASSUMING THE SENSE STRAND IS THE
            FORWARD STRAND. That is, if the sense strand is the reverse strand,
            this is the distance between the 3' end of the read and the 3' end
            of the intron.
        8. Number of nucleotides between 3' end of intron and 3' end of read
            from which it was inferred, ASSUMING THE SENSE STRAND IS THE
            FORWARD STRAND.
        --------------------------------------------------------------------
        9. Number of instances of junction, insertion, or deletion in sample;
            this is always +1 before bed_pre combiner/reducer

        Input is partitioned by fields 1-5 and sorted by field 6.

        Hadoop output (written to stdout)
        ----------------------------
        Tab-delimited output tuple columns (bed):
        1. 'I', 'D', or 'N' for insertion, deletion, or junction line
        2. Sample index
        3. Number string representing RNAME (+ '+ or -' if junction; same as
            field 6)
        4. Start position (Last base before insertion, first base of deletion,
                            or first base of intron)
        5. End position (Last base before insertion, last base of deletion
                            (exclusive), or last base of intron (exclusive))
        6. '+' or '-' indicating which strand is the sense strand for
            junctions, inserted sequence for insertions, or deleted sequence
            for deletions
        ----Next fields are for junctions only; they are '\x1c' for indels----
        7. MAX number of nucleotides between 5' end of intron and 5' end of
            read from which it was inferred, ASSUMING THE SENSE STRAND IS THE
            FORWARD STRAND. That is, if the sense strand is the reverse strand,
            this is the distance between the 3' end of the read and the 3' end
            of the intron.
        8. MAX number of nucleotides between 3' end of intron and 3' end of
            read from which it was inferred, ASSUMING THE SENSE STRAND IS THE
            FORWARD STRAND.
        9. MAXIMIN (number of nucleotides between 5' end of intron and 5' end
                    of read, number of nucleotides between 3' end of intron and
                    3' end of read);
           min is between the args above; max is across reads.

        Tab-delimited output tuple columns (collect)
        1. '0' if insertion, '1' if deletion, or '2' if junction line
        2. Number string representing RNAME (+ '+ or -' if junction; same as
                                                field 6)
        3. Start position (Last base before insertion, first base of deletion,
                            or first base of intron)
        4. End position (Last base before insertion, last base of deletion
                            (exclusive), or last base of intron (exclusive))
        5. '+' or '-' indicating which strand is the sense strand for
            junctions, inserted sequence for insertions, or deleted sequence
            for deletions
        6. Coverage of feature for sample with index N
        ...
        N + 6. Coverage of feature in sample with index N
        --------------------------------------------------------------------
        10. SUMMED number of instances of junction, insertion, or deletion in
            sample

        OUTPUT COORDINATES ARE 1-INDEXED.

        input_stream: where to find input indels/junctions
        output_stream: where to write output
        manifest_object: object of class LabelsAndIndices that maps indices
            to labels and back; used to count number of samples.
        sample_fraction: fraction of samples in which an indel must appear
            to pass filter if coverage_threshold criterion is not satisfied
        coverage_threshold: number of reads that must overlap indel in at
            least one sample to pass filter of sample_fraction criterion is not
            satisfied
        verbose: output extra debugging statements

        Return value: tuple (input line count, output line count)
    """
    input_line_count, output_line_count = 0, 0

    '''Compute minimum number of samples in which indel should appear to be
    output if coverage threshold not met.'''
    total_sample_count = len(manifest_object.label_to_index)
    min_sample_count = int(round(
                total_sample_count * sample_fraction
            ))

    for (line_type, rname, pos, end_pos, strand_or_seq), xpartition in xstream(
                input_stream, 5
            ):
        collect_specs = [rname, pos, end_pos if line_type != 'N'
                                             else str(int(end_pos) - 1),
                                     strand_or_seq]
        sample_indexes, coverages = [], []
        if line_type == 'N':
            counter.add('junction_line')
            for sample_index, data in itertools.groupby(
                                                    xpartition, 
                                                    key=lambda val: val[0]
                                                ):
                coverage_sum = 0
                max_left_displacement, max_right_displacement = None, None
                maximin_displacement = None
                for _, left_displacement, right_displacement, coverage in data:
                    input_line_count += 1
                    left_displacement = int(left_displacement)
                    right_displacement = int(right_displacement)
                    max_left_displacement = max(left_displacement,
                                                max_left_displacement)
                    max_right_displacement = max(right_displacement,
                                                 max_right_displacement)
                    maximin_displacement = max(
                            min(left_displacement, right_displacement),
                            maximin_displacement
                        )
                    coverage_sum += int(coverage)
                assert max_left_displacement is not None
                assert max_right_displacement is not None
                assert maximin_displacement is not None
                counter.add('bed_line')
                print >>output_stream, \
                    'bed\tN\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d' % (
                        sample_index, rname, pos, end_pos, strand_or_seq,
                        max_left_displacement, max_right_displacement,
                        maximin_displacement, coverage_sum
                    )
                sample_indexes.append(sample_index)
                coverages.append(coverage_sum)
                output_line_count += 1
            counter.add('collect_line')
            output_stream.write('collect\t2\t')
            print >>output_stream, '\t'.join(
                        collect_specs
                        + [','.join(sample_indexes),
                           ','.join(map(str, coverages))]
                    )
            output_line_count += 1
        else:
            counter.add('insertion_line' if line_type == 'I'
                            else 'deletion_line')
            assert line_type in 'ID'
            sample_count = 0
            for sample_index, data in itertools.groupby(
                                                    xpartition, 
                                                    key=lambda val: val[0]
                                                ):
                coverage_sum = 0
                for _, _, _, coverage in data:
                    input_line_count += 1
                    coverage_sum += int(coverage)
                counter.add('bed_line')
                print >>output_stream, \
                    'bed\t%s\t%s\t%s\t%s\t%s\t%s\t\x1c\t\x1c\t\x1c\t%d' % (
                        line_type, sample_index, rname, pos, end_pos,
                        strand_or_seq, coverage_sum
                    )
                sample_indexes.append(sample_index)
                coverages.append(coverage_sum)
                sample_count += 1
                output_line_count += 1
            max_coverage = max(coverages)
            if (sample_count >= min_sample_count
                or (max_coverage >= coverage_threshold
                    and coverage_threshold != -1)):
                counter.add('collect_line')
                if line_type == 'I':
                    output_stream.write('collect\t0\t')
                else:
                    output_stream.write('collect\t1\t')
                print >>output_stream, \
                    '\t'.join(
                        collect_specs 
                        + [','.join(sample_indexes),
                           ','.join(map(str, coverages))]
                    )
                output_line_count += 1
            else:
                counter.add('indel_filtered_out')
                if verbose:
                    print >>sys.stderr, (
                        'Indel (%s, %s, %s, %s) filtered out; it appeared in '
                        '%d sample(s), and its coverage in any one sample did '
                        'not exceed %d.'
                    ) % (rname, strand_or_seq, pos, end_pos, sample_count,
                            max_coverage)
    counter.flush()
    return input_line_count, output_line_count

if __name__ == '__main__':
    # Print file's docstring if -h is invoked
    parser = argparse.ArgumentParser(description=__doc__, 
                formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--sample-fraction', type=float, required=False,
        default=0.05,
        help=('An indel passes the filter if it is present in this fraction '
              'of samples; if the indel meets neither this criterion nor '
              'the --coverage-threshold criterion, it is filtered out'))
    parser.add_argument('--coverage-threshold', type=int, required=False,
        default=5,
        help=('An indel passes the filter if it is covered by at least this '
              'many reads in a single sample; if the indel meets neither '
              'this criterion nor the --sample-fraction criterion, it is '
              'filtered out; use -1 to disable'))
    parser.add_argument('--keep-alive', action='store_const', const=True,
        default=False,
        help='Periodically print Hadoop status messages to stderr to keep ' \
             'job alive')
    parser.add_argument('--verbose', action='store_const', const=True,
        default=False,
        help='Print out extra debugging statements')
    parser.add_argument('--test', action='store_const', const=True,
        default=False,
        help='Run unit tests; DOES NOT NEED INPUT FROM STDIN, AND DOES NOT '
             'WRITE TO STDOUT')

    # Add command-line arguments for dependencies
    manifest.add_args(parser)

    '''Now collect arguments. While the variable args declared below is
    global, properties of args are also arguments of the go() function so
    different command-line arguments can be passed to it for unit tests.'''
    args = parser.parse_args(sys.argv[1:])

    # Start keep_alive thread immediately
    if args.keep_alive:
        from dooplicity.tools import KeepAlive
        keep_alive_thread = KeepAlive(sys.stderr)
        keep_alive_thread.start()

if __name__ == '__main__' and not args.test:
    start_time = time.time()
    manifest_object = manifest.LabelsAndIndices(
                                os.path.expandvars(args.manifest)
                            )
    input_line_count, output_line_count = go(
            manifest_object=manifest_object,
            input_stream=sys.stdin,
            output_stream=sys.stdout,
            sample_fraction=args.sample_fraction,
            coverage_threshold=args.coverage_threshold,
            verbose=args.verbose
        )
    print >>sys.stderr, 'DONE with bed_pre.py; in/out =%d/%d; time=%0.3f s' \
                         % (input_line_count, output_line_count,
                            time.time() - start_time)
elif __name__ == '__main__':
    # Test units
    del sys.argv[1:] # Don't choke on extra command-line parameters
    import unittest
    import shutil
    import tempfile

    class TestGo(unittest.TestCase):
        """ Tests go(). """
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
            junctions_and_indels = (
                    'I\t000000000000\t140\t140\tATAC\t0\t\x1c\t\x1c\t6\n'
                    'I\t000000000000\t140\t140\tATAC\t0\t\x1c\t\x1c\t6\n'
                    'I\t000000000002\t170\t170\tCT\t1\t\x1c\t\x1c\t4\n'
                    'D\t000000000001\t150\t156\tAACCTT\t0\t\x1c\t\x1c\t3\n'
                    'D\t000000000001\t150\t156\tAACCTT\t1\t\x1c\t\x1c\t3\n'
                    'N\t000000000003\t3567\t3890\t+\t2\t3\t3\t2'
                )
            '''Above, for a coverage_threshold 5 and a sample_fraction 0.5,
            the first insertion should be kept, the second insertion
            (third line) should be filtered out, the deletion should be kept,
            and the junction must be kept; junctions are not filtered.'''
            with open(self.input_file, 'w') as output_stream:
                output_stream.write(junctions_and_indels)
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
                    input_line_count, 6
                )
            self.assertEquals(
                    output_line_count, 8
                )
            output_lines = []
            with open(self.output_file) as output_stream:
                for line in output_stream:
                    output_lines.append(line.strip())
            self.assertTrue(
                    ('bed\tI\t0\t000000000000\t140\t140\tATAC'
                     '\t\x1c\t\x1c\t\x1c\t12')
                    in output_lines
                )
            self.assertTrue(
                    ('bed\tI\t1\t000000000002\t170\t170\tCT'
                     '\t\x1c\t\x1c\t\x1c\t4')
                    in output_lines
                )
            self.assertTrue(
                    ('bed\tD\t0\t000000000001\t150\t156\tAACCTT'
                     '\t\x1c\t\x1c\t\x1c\t3')
                    in output_lines
                )
            self.assertTrue(
                    ('bed\tD\t1\t000000000001\t150\t156\tAACCTT'
                     '\t\x1c\t\x1c\t\x1c\t3')
                    in output_lines
                )
            self.assertTrue(
                    'bed\tN\t2\t000000000003\t3567\t3890\t+\t3\t3\t3\t2'
                    in output_lines
                )
            self.assertTrue(
                    'collect\t0\t000000000000\t140\t140\tATAC\t12\t0\t0'
                    in output_lines
                )
            self.assertTrue(
                    'collect\t1\t000000000001\t150\t156\tAACCTT\t3\t3\t0'
                    in output_lines
                )
            self.assertTrue(
                    'collect\t2\t000000000003\t3567\t3889\t+\t0\t0\t2'
                    in output_lines
                )
            self.assertEquals(
                    len(output_lines), 8
                )
        def tearDown(self):
            # Kill temporary directory
            shutil.rmtree(self.temp_dir_path)
    
    unittest.main()
