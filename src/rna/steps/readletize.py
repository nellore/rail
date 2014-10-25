#!/usr/bin/env python
"""
Rail-RNA-readletize

Follows Rail-RNA-align_reads (after an optional intermediate Rail-RNA-sum step)
Precedes Rail-RNA-align_readlets (with an optional intermediate Rail-RNA-sum
    step)

Reduce step for MapReduce pipelines that divides each input sequence into
several overlapping segments, or readlets. These readlets will later be
aligned and used to infer the introns that original read sequences without
end-to-end alignments overlap.

Input (read from stdin)
----------------------------
Tab-delimited input tuple columns (readletize)
1. SEQ or its reversed complement, whichever is first in alphabetical order
    (does not have to be unique)
2. (('+' if primary alignment is substantially soft-clipped, else '-')
    + the sample label if field 1 is the read sequence); if the list is empty,
    '\x1c' is filler
3. (('+' if primary alignment is substantially soft-clipped, else '-') + the
    sample label if field 1 is the read sequence's reversed complement); if
    the list is empty, '\x1c' is filler

Input is partitioned by field 1, the read sequence.

Hadoop output (written to stdout)
----------------------------
(readletized) Tab-delimited output tuple columns:
1. Readlet sequence or its reversed complement, whichever is first in
    alphabetical order
2. Read sequence ID + ('-' if readlet sequence is reverse-complemented; else
    '+') + '\x1e' + displacement of readlet's 5' end from read's 5' end +
    '\x1e' + displacement of readlet's 3' end from read's 3' end (+, for
    EXACTLY one readlet of a read sequence, '\x1e' + read sequence + '\x1e'
    + number of instances of read sequence + '\x1e' + number of instances of
    read sequence's reversed complement + '\x1e' + (an '\x1f'-separated set of
    unique sample labels with read sequences that match the original read
    sequence) + '\x1e' + (an '\x1f'-separated set of unique sample labels with
    read sequences that match the reversed complement of the original read
    sequence). Here, a read sequence ID takes the form X:Y, where X is the
    "mapred_task_partition" environment variable -- a unique index for a task
    within a job -- and Y is the index of the read sequence relative to the
    beginning of the input stream.

(unique) Single column:
1. A unique read sequence
"""
import sys
import os
import string
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

# Initialize global variables for tracking number of input/output lines
_input_line_count, _output_line_count = 0, 0
_reversed_complement_translation_table = string.maketrans('ATCG', 'TAGC')

def go(output_stream=sys.stdout, input_stream=sys.stdin, min_readlet_size=8,
        max_readlet_size=25, readlet_interval=5, capping_multiplier=1.5,
        verbose=False, report_multiplier=1.2):
    """ Runs Rail-RNA-readletize

        Input (read from stdin)
        ----------------------------
        Tab-delimited input tuple columns:
        1. Read sequence (does not have to be unique)
        2. '\x1d'-separated list of samples, one for each sample read whose
            sequence matches the sequence in field 1. If list is empty,
            '\x1c' is filler.
        3. '\x1d'-separated list of samples, one for each sample read whose
            sequence matches the reversed complement of the sequence in field
            1. If list is empty, '\x1c' is filler.

        Input is partitioned by field 1, the read sequence.

        Hadoop output (written to stdout)
        ----------------------------
        Tab-delimited output tuple columns

        1. Readlet sequence or its reversed complement, whichever is first in
            alphabetical order
        2. Read sequence ID + ('-' if readlet sequence is reverse-complemented;
            else '+') + '\x1e' + displacement of readlet's 5' end from read's
            5' end + '\x1e' + displacement of readlet's 3' end from read's 3'
            end (+, for EXACTLY one readlet of a read sequence, '\x1e' + read
            sequence + '\x1e' + number of instances of read sequence + '\x1e'
            + number of instances of read sequence's reversed complement +
            '\x1e' + (an '\x1f'-separated set of unique sample labels with read
            sequences that match the original read sequence) + '\x1e' + (an
            '\x1f'-separated set of unique sample labels with read sequences
            that match the reversed complement of the original read sequence).
            Here, a read sequence ID takes the form X:Y, where X is the
            "mapred_task_partition" environment variable -- a unique index for
            a task within a job -- and Y is the index of the read sequence
            relative to the beginning of the input stream.

        output_stream: where to write output, typically a file stream.
        input_stream: where to retrieve sequences, typically sys.stdin on first
            pass of Bowtie and a file stream on second pass.
        min_readlet_size: "capping" readlets (that is, readlets that terminate
            at a given end of the read) are never smaller than this value.
            Ignored if readletize=False.
        max_readlet_size: size of every noncapping readlet. Ignored if 
            readletize=False.
        readlet_interval: number of bases separating successive readlets along
            the read. Ignored if readletize=False.
        capping_multiplier: successive capping readlets on a given end of a read
            are increased in size exponentially with base
            capping_multiplier.
        verbose: True if reads/readlets should occasionally be written 
            to stderr.
        report_multiplier: if verbose is True, the line number of a read or its 
            first readlet written to stderr increases exponentially with base
            report_multiplier.

        No return value.
    """
    global _input_line_count, _output_line_count
    next_report_line = 0
    try:
        task_partition = os.environ['mapred_task_partition']
    except KeyError:
        # Hadoop 2.x?
        try:
            task_partition = os.environ['mapreduce_task_partition']
        except KeyError:
            # A unit test is being run
            task_partition = '0'
    # Build list of capping readlet sizes
    cap_sizes = []
    cap_size = min_readlet_size
    while cap_size <= max_readlet_size:
        cap_sizes.append(cap_size)
        cap_size = int(cap_size*capping_multiplier)
        if cap_size == cap_sizes[-1]:
          cap_size += 1
    if cap_size != max_readlet_size:
        # Always have a start or end read of length max_readlet_size
        cap_sizes.append(max_readlet_size)
    for seq_count, ((seq,), xpartition) in enumerate(xstream(input_stream, 1)):
        print >>output_stream, 'unique\t%s' % seq
        seq_id = task_partition + ':' + str(seq_count)
        if len(seq) < min_readlet_size:
            continue
        samples = set()
        reversed_complement_samples = set()
        sample_count, reversed_complement_sample_count = 0, 0
        readletize = False
        for line_samples, line_reversed_complement_samples in xpartition:
            _input_line_count += 1
            if line_samples != '\x1c':
                samples_to_add = [sample[1:] for sample
                                    in line_samples.split('\x1d')]
                if '+' in line_samples:
                    '''Only readletize if sequence is sourced by at least one
                    soft clip.'''
                    readletize = True
                sample_count += len(samples_to_add)
                samples.update(samples_to_add)
            if line_reversed_complement_samples != '\x1c':
                reversed_complement_samples_to_add = [sample[1:] for sample
                    in line_reversed_complement_samples.split('\x1d')]
                if '+' in line_reversed_complement_samples:
                    '''Only readletize if sequence is sourced by at least one
                    soft clip.'''
                    readletize = True
                reversed_complement_sample_count += len(
                        reversed_complement_samples_to_add
                    )
                reversed_complement_samples.update(
                        reversed_complement_samples_to_add
                    )
        if not readletize: continue
        '''Construct a readlet identifier as follows: read sequence ID + 
        ('-' if readlet sequence is reverse-complemented; else '+') + '\x1e' +
        displacement of readlet's 5' end from read's 5' end + '\x1e' +
        displacement of readlet's 3' end from read's 3' end (+, for EXACTLY one
        readlet of a read sequence, '\x1e' + read sequence + '\x1e' +
        number of instances of read sequence + '\x1e' + number of instances of
        read sequence's reversed complement + '\x1e' + (an '\x1f'-separated set
        of unique sample labels with read sequences that match the original
        read sequence) + '\x1e' + (an '\x1f'-separated set of unique sample
        labels with read sequences that match the reversed complement of the
        original read sequence). Here, a read sequence ID takes the form X:Y,
        where X is the "mapred_task_partition" environment variable -- a unique
        index for a task within a job -- and Y is the index of the read
        sequence relative to the beginning of the input stream.'''
        to_write = []
        seq_size = len(seq)
        # Add capping readlets
        for cap_size in cap_sizes:
            readlet_seq = seq[:cap_size]
            reversed_complement_readlet_seq = readlet_seq[::-1].translate(
                            _reversed_complement_translation_table
                        )
            if readlet_seq < reversed_complement_readlet_seq:
                to_write.append('%s\t%s\x1e%d\x1e%d' % (readlet_seq, seq_id 
                                                        + '+', 0,
                                                        seq_size - cap_size))
            else:
                to_write.append('%s\t%s\x1e%d\x1e%d' \
                                    % (reversed_complement_readlet_seq,
                                        seq_id + '-', 0, seq_size - cap_size))
            readlet_seq = seq[-cap_size:]
            reversed_complement_readlet_seq = readlet_seq[::-1].translate(
                            _reversed_complement_translation_table
                        )
            if readlet_seq < reversed_complement_readlet_seq:
                to_write.append('%s\t%s\x1e%d\x1e%d' % (readlet_seq,
                                                    seq_id + '+',
                                                    seq_size - cap_size, 0))
            else:
                to_write.append('%s\t%s\x1e%d\x1e%d' \
                                    % (reversed_complement_readlet_seq,
                                        seq_id + '-', seq_size - cap_size, 0))
        # Add noncapping readlets
        for j in xrange(readlet_interval, seq_size - max_readlet_size,
                            readlet_interval):
            readlet_seq = seq[j:j+max_readlet_size]
            reversed_complement_readlet_seq = readlet_seq[::-1].translate(
                            _reversed_complement_translation_table
                        )
            if readlet_seq < reversed_complement_readlet_seq:
                to_write.append('%s\t%s\x1e%d\x1e%d' % (readlet_seq,
                                                    seq_id + '+', j, 
                                                    seq_size - j
                                                        - max_readlet_size))
            else:
                to_write.append('%s\t%s\x1e%d\x1e%d' \
                                    % (reversed_complement_readlet_seq,
                                        seq_id + '-', j, 
                                        seq_size - j - max_readlet_size))
        # Add additional info to first readlet in to_write
        to_write[0] = '\x1e'.join([to_write[0], seq,
                                    str(sample_count),
                                    str(reversed_complement_sample_count),
                                    '\x1f'.join(samples),
                                    '\x1f'.join(reversed_complement_samples)])
        if verbose and next_report_line == i:
            print >>sys.stderr, 'First readlet from read %d: %s' \
                % (i + 1, to_write[0])
            next_report_line = int((next_report_line + 1)
                * report_multiplier + 1) - 1
        for readlet in to_write:
            print >>output_stream, 'readletized\t%s' % readlet
        _output_line_count += len(to_write)
    output_stream.flush()

if __name__ == '__main__':
    import argparse
    # Print file's docstring if -h is invoked
    parser = argparse.ArgumentParser(description=__doc__, 
                formatter_class=argparse.RawDescriptionHelpFormatter)
    # Add command-line arguments
    parser.add_argument('--min-readlet-size', type=int, required=False,
        default=15, 
        help='Capping readlets (that is, readlets that terminate '
             'at a given end of the read) are never smaller than this value')
    parser.add_argument('--max-readlet-size', type=int, required=False,
        default=25, 
        help='Size of every noncapping readlet')
    parser.add_argument('--readlet-interval', type=int, required=False,
        default=12, 
        help='Number of bases separating successive noncapping readlets along '
             'the read')
    parser.add_argument('--capping-multiplier', type=float, required=False,
        default=1.5, 
        help='Successive capping readlets on a given end of a read are '
             'increased in size exponentially with this base')
    parser.add_argument('--report_multiplier', type=float, required=False,
        default=1.2,
        help='When --verbose is also invoked, the only lines of lengthy '
             'intermediate output written to stderr have line number that '
             'increases exponentially with this base')
    parser.add_argument('--verbose', action='store_const', const=True,
        default=False,
        help='Print out extra debugging statements')
    parser.add_argument(\
        '--stranded', action='store_const', const=True, default=False,
        help='Assume input reads come from the sense strand; then partitions '
             'in output have terminal + and - indicating sense strand')
    parser.add_argument('--test', action='store_const', const=True,
        default=False,
        help='Run unit tests; DOES NOT NEED INPUT FROM STDIN, AND DOES NOT '
             'WRITE EXONS AND INTRONS TO STDOUT')

    '''Now collect other arguments. While the variable args declared below is
    global, properties of args are also arguments of the go() function so
    different command-line arguments can be passed to it for unit tests.'''
    args = parser.parse_args(sys.argv[1:])

if __name__ == '__main__' and not args.test:
    import time
    start_time = time.time()
    go(verbose=args.verbose, 
        min_readlet_size=args.min_readlet_size, 
        max_readlet_size=args.max_readlet_size,
        readlet_interval=args.readlet_interval,
        capping_multiplier=args.capping_multiplier,
        report_multiplier=args.report_multiplier)
    print >>sys.stderr, 'DONE with readletize.py; in/out=%d/%d; ' \
        'time=%0.3f s' % (_input_line_count, _output_line_count,
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
            '''Every sequence has 100 bases.'''
            input_seqs = \
                    'TTACATACCATACAGTGCGCTAGCGGGTGACAGATATAATGCAGATCCAT' \
                    'ACAGACCAGATGGCAGACATGTGTTGCAGSCTGCAAGTGCAACGCGGTGA' \
                    '\t1\x1d1\t1\x1d1\x1d1\n' \
                    'GCAGAGTGCCGCAATGACGTGCGCCAAAGCGGTGACAGGGTGACAGTGAA' \
                    'CCAAGTGACAAGTGAACAGGTGCCAGAGTGACCGAGTGACCAGTGGACCA' \
                    '\t1\x1d2\t1\x1d2\x1d3\n' \
                    'CAGAGTGCCGCAATGACGTGCGCCAAAGCGGACAAAGCACCATGACAAGT' \
                    'ACACAGGTGACAGTGACAAGACAGAGGTGACACAGAGAAAGtGGGTGTGA' \
                    '\t1\x1d3\t1\x1d3\x1d2\x1d3\x1d1\n' \
                    'ATCGATTAAGCTATAACAGATAACATAGACATTGCGCCCATAATAGATAA' \
                    'CTGACACCTGACCAGTGCCAGATGACCAGTGCCAGATGGACGACAGTAGC' \
                    '\t1\t1\n'
            self.temp_dir_path = tempfile.mkdtemp()
            self.input_file = os.path.join(self.temp_dir_path,
                                'sample_input.tsv')
            # Store reads in temporary file
            with open(self.input_file, 'w') as reads_stream:
                reads_stream.write(input_seqs)
            # For storing output of go()
            self.output_file = os.path.join(self.temp_dir_path,
                                'sample_output.tsv')

        def test_output(self):
            """ Fails if output of go() is not in the right form. """
            with open(self.input_file) as input_stream:
                with open(self.output_file, 'w') as output_stream:
                    go(output_stream=output_stream, input_stream=input_stream,
                                    capping_multiplier=2, min_readlet_size=25,
                                    readlet_interval=5, max_readlet_size=50)

            collected_readlets = []
            with open(self.output_file) as processed_stream:
                for readlet in processed_stream:
                    collected_readlets.append(readlet.rstrip().split('\t'))
            '''Each read from input_reads spans 100 bases, and from the
            arguments passed to go() above, noncapping readlets should span 50.
            The capping fraction should arrange for one readlet spanning 25
            bases on either end of a given read. There should be 13 readlets in
            total. Spot-check some readlets after checking for read info.'''
            read_info = [info.split('\x1e') for _, info in collected_readlets]
            read_info = [info[3:4] + [set(info[-2].split('\x1f')), 
                            set(info[-1].split('\x1f'))]
                            for info in read_info if len(info) > 3]
            self.assertTrue([
                    'TTACATACCATACAGTGCGCTAGCGGGTGACAGATATAATGCAGATCCAT'
                    'ACAGACCAGATGGCAGACATGTGTTGCAGSCTGCAAGTGCAACGCGGTGA',
                    set(['1']),
                    set(['1'])
                ] in read_info)
            self.assertTrue([
                    'GCAGAGTGCCGCAATGACGTGCGCCAAAGCGGTGACAGGGTGACAGTGAA'
                    'CCAAGTGACAAGTGAACAGGTGCCAGAGTGACCGAGTGACCAGTGGACCA',
                    set(['1','2']),
                    set(['1','2','3'])
                ] in read_info)
            self.assertTrue([
                    'CAGAGTGCCGCAATGACGTGCGCCAAAGCGGACAAAGCACCATGACAAGT'
                    'ACACAGGTGACAGTGACAAGACAGAGGTGACACAGAGAAAGtGGGTGTGA',
                    set(['1','3']),
                    set(['1','2','3'])
                ] in read_info)
            self.assertTrue([
                    'ATCGATTAAGCTATAACAGATAACATAGACATTGCGCCCATAATAGATAA'
                    'CTGACACCTGACCAGTGCCAGATGACCAGTGCCAGATGGACGACAGTAGC',
                    set(['1']),
                    set(['1'])
                ] in read_info)
            # Capping readlets
            self.assertTrue([       'ATCGATTAAGCTATAACAGATAACA'
                                    'TAGACATTGCGCCCATAATAGATAA',
                                    '0:3+\x1e0\x1e50'
                                ] in collected_readlets or
                                [   'ATCGATTAAGCTATAACAGATAACA'
                                    'TAGACATTGCGCCCATAATAGATAA',
                                    '0:3+\x1e0\x1e50\x1e'
                                    'ATCGATTAAGCTATAACAGATAACA'
                                    'TAGACATTGCGCCCATAATAGATAA'
                                    'CTGACACCTGACCAGTGCCAGATGA'
                                    'CCAGTGCCAGATGGACGACAGTAGC'
                                    '\x1e1\x1e1\x1e1\x1e1'
                                ] in collected_readlets
                            )
            self.assertTrue([       'CCAGTGCCAGATGGACGACAGTAGC',
                                    '0:3+\x1e75\x1e0'
                                ] in collected_readlets or
                                [   'CCAGTGCCAGATGGACGACAGTAGC',
                                    '0:3+\x1e75\x1e0\x1e'
                                    'ATCGATTAAGCTATAACAGATAACA'
                                    'TAGACATTGCGCCCATAATAGATAA'
                                    'CTGACACCTGACCAGTGCCAGATGA'
                                    'CCAGTGCCAGATGGACGACAGTAGC'
                                    '\x1e1\x1e1\x1e1\x1e1'
                                ] in collected_readlets
                            )
            # Noncapping readlets
            self.assertTrue([       'GTCAG'
                                    'TTATCTATTATGGGCGCAATGTCTA'
                                    'TGTTATCTGTTATAGCTTAA',
                                    '0:3-\x1e5\x1e45'
                                ] in collected_readlets or
                                [   'GTCAG'
                                    'TTATCTATTATGGGCGCAATGTCTA'
                                    'TGTTATCTGTTATAGCTTAA',
                                    '0:3-\x1e5\x1e45\x1e'
                                    'ATCGATTAAGCTATAACAGATAACA'
                                    'TAGACATTGCGCCCATAATAGATAA'
                                    'CTGACACCTGACCAGTGCCAGATGA'
                                    'CCAGTGCCAGATGGACGACAGTAGC'
                                    '\x1e1\x1e1\x1e1\x1e1'
                                ] in collected_readlets
                            )
            self.assertTrue([       'TAATAGATAA'
                                    'CTGACACCTGACCAGTGCCAGATGA'
                                    'CCAGTGCCAGATGGA',
                                    '0:3+\x1e40\x1e10'
                                ] in collected_readlets or
                                [   'TAATAGATAA'
                                    'CTGACACCTGACCAGTGCCAGATGA'
                                    'CCAGTGCCAGATGGA',
                                    '0:3+\x1e40\x1e10\x1e'
                                    'ATCGATTAAGCTATAACAGATAACA'
                                    'TAGACATTGCGCCCATAATAGATAA'
                                    'CTGACACCTGACCAGTGCCAGATGA'
                                    'CCAGTGCCAGATGGACGACAGTAGC'
                                    '\x1e1\x1e1\x1e1\x1e1'
                                ] in collected_readlets
                            )

        def tearDown(self):
            # Kill temporary directory
            shutil.rmtree(self.temp_dir_path)
   
    unittest.main()
