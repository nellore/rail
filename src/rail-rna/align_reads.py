#!/usr/bin/env python
"""
Rail-RNA-align_reads

Follows Rail-RNA-preprocess
Precedes Rail-RNA-readletize, Rail-RNA-coverage_pre (both optionally with an
    intermediate Rail-RNA-sum combine/reduce step), Rail-RNA-realign,
    and Rail-RNA-bam

Alignment script for MapReduce pipelines that wraps Bowtie. Obtains exonic 
chunks from end-to-end alignments. Outputs sequences that do not align for
readletizing and, later, intron inference.

Input (read from stdin)
----------------------------
Tab-delimited input tuple columns in a mix of any of the following three
formats:
 Format 1 (single-end, 3-column):
  1. Name
  2. Nucleotide sequence
  3. Quality sequence
 Format 2 (paired-end, 5-column):
  1. Name
  2. Nucleotide sequence for mate 1
  3. Quality sequence for mate 1
  4. Nucleotide sequence for mate 2
  5. Quality sequence for mate 2
 Format 3 (paired, 6-column):
  1. Name for mate 1
  2. Nucleotide sequence for mate 1
  3. Quality sequence for mate 1
  4. Name for mate 2
  5. Nucleotide sequence for mate 2
  6. Quality sequence for mate 2

Hadoop output (written to stdout)
----------------------------
A given RNAME sequence is partitioned into intervals ("bins") of some 
user-specified length (see partition.py).

Exonic chunks (aka ECs; three formats, any or all of which may be emitted):

Format 1 (exon_ival); tab-delimited output tuple columns:
1. Reference name (RNAME in SAM format) + ';' + bin number
2. Sample label
3. EC start (inclusive) on forward strand
4. EC end (exclusive) on forward strand

Format 2 (exon_diff); tab-delimited output tuple columns:
1. Reference name (RNAME in SAM format) + ';' + bin number
2. Sample label
3. max(EC start, bin start) (inclusive) on forward strand IFF diff is
    positive and EC end (exclusive) on forward strand IFF diff is negative
4. +1 or -1.

Note that only unique alignments are currently output as ivals and/or diffs.

Format 3 (end_to_end_sam); tab-delimited output tuple columns:
Standard SAM output except fields are in different order, and the first field
corresponds to sample label. (Fields are reordered to facilitate partitioning
by sample name/RNAME and sorting by POS.) Each line corresponds to a
spliced alignment. The order of the fields is as follows.
1. Sample label
2. Number string representing RNAME; see BowtieIndexReference class in
    bowtie_index for conversion information
3. POS
4. QNAME
5. FLAG
6. MAPQ
7. CIGAR
8. RNEXT
9. PNEXT
10. TLEN
11. SEQ
12. QUAL
... + optional fields

Reads with no end-to-end alignments

Tab-delimited output tuple columns (unmapped):
1. QNAME
2. SEQ
3. QUAL

Tab-delimited output tuple columns (readletize)
1. SEQ or its reversed complement, whichever is first in alphabetical order
2. The sample label if field 1 is the read sequence; else '\x1c'
3. The sample label if field 1 is the read sequence's reversed complement;
    else '\x1c'

ALL OUTPUT COORDINATES ARE 1-INDEXED.
"""
import sys
import os
import site
import threading
import subprocess
import string
import shutil
import tempfile
import atexit

base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for directory_name in ['bowtie', 'sample', 'interval', 'manifest']:
    site.addsitedir(os.path.join(base_path, directory_name))

import bowtie
import bowtie_index
import sample
import partition
import manifest

# Initialize global variables for tracking number of input/output lines
_input_line_count = 0
_output_line_count = 0

_reversed_complement_translation_table = string.maketrans('ATCG', 'TAGC')

def handle_temporary_directory(temp_dir_path):
    """ Deletes temporary directory.

        temp_dir_paths: path to temporary directory to delete

        No return value.
    """
    shutil.rmtree(temp_dir_path)

def write_reads(output_stream, input_stream=sys.stdin, verbose=False,
    report_multiplier=1.2):
    """ Writes input reads/readlets in tab-separated format parsable by Bowtie.

        Unpaired reads are marked 0 after original read ID, while paired-end
        reads are marked 1 and 2.

        Input formats:
        3-token: name TAB seq TAB qual
        5-token: name TAB seq1 TAB qual1 TAB seq2 TAB qual2
        6-token: name1 TAB seq1 TAB qual1 TAB name2 TAB seq2 TAB qual2

        output_stream: where to write reads, typically a file stream.
        input_stream: where to retrieve reads, typically sys.stdin on first
            pass of Bowtie and a file stream on second pass.
        verbose: True if reads/readlets should occasionally be written 
            to stderr.
        report_multiplier: if verbose is True, the line number of a read or its 
            first readlet written to stderr increases exponentially with base
            report_multiplier.

        No return value.
    """
    next_report_line = 0
    global _input_line_count
    for i, line in enumerate(input_stream):
        tokens = line.rstrip().split('\t')
        if len(tokens) not in (3, 5, 6):
            raise RuntimeError('The following line has an invalid number of '
                'tab-separated tokens:\n%sA valid line has 3, 5, or 6 such '
                'tokens.' % line)
        # Check that a properly formed label is embedded in the read name
        sample.hasLab(tokens[0], mustHave=True)
        if len(tokens) == 3:
            to_write = '%s\t%s\t%s' \
                % (tokens[0], tokens[1], tokens[2])
        else:
            to_write = '%s\t%s\t%s\n%s\t%s\t%s' \
                % (tokens[0], tokens[1], tokens[2],
                    (tokens[0] if len(tokens) == 5 else tokens[3]),
                        tokens[-2], tokens[-1])
        print >>output_stream, to_write
        if verbose and next_report_line == i:
            print >>sys.stderr, 'Read(s) %d: %s' % (i, to_write)
            next_report_line = int((next_report_line + 1)
                * report_multiplier + 1) - 1
        _input_line_count += 1
    output_stream.flush()

class BowtieOutputThread(threading.Thread):
    """ Processes Bowtie alignments, emitting tuples for exons and introns. """
    
    def __init__(self, input_stream, reference_index, manifest_object, 
        output_stream=sys.stdout, exon_differentials=True,
        exon_intervals=False, end_to_end_sam=True, verbose=False,
        bin_size=10000, report_multiplier=1.2):
        """ Constructor for BowtieOutputThread.

            input_stream: where to retrieve Bowtie's SAM output, typically a
                Bowtie process's stdout.
            reference_index: object of class bowtie_index.BowtieIndexReference
                that permits access to reference; used for realignment of
                unmapped regions between exonic chunks.
            manifest_object: object of class LabelsAndIndices that maps indices
                to labels and back; used to shorten intermediate output.
            output_stream: where to emit exon and intron tuples; typically,
                this is sys.stdout.
            exon_differentials: True iff EC differentials are to be emitted.
            exon_intervals: True iff EC intervals are to be emitted.
            end_to_end_sam: True iff SAM with end_to_end alignments should
                be output. See docstring for more information.
            verbose: True if alignments should occasionally be written 
                to stderr.
            bin_size: genome is partitioned in units of bin_size for later load
                balancing.
            report_multiplier: if verbose is True, the line number of an
                alignment written to stderr increases exponentially with base
                report_multiplier.
        """
        super(BowtieOutputThread, self).__init__()
        self.daemon = True
        self.input_stream = input_stream
        self.output_stream = output_stream
        self.reference_index = reference_index
        self.verbose = verbose
        self.bin_size = bin_size
        self.report_multiplier = report_multiplier
        self.exon_differentials = exon_differentials
        self.exon_intervals = exon_intervals
        self.end_to_end_sam = end_to_end_sam
        self.manifest_object = manifest_object

    def run(self):
        """ Prints exons for end-to-end alignments.

            Overrides default method containing thread activity.

            No return value.
        """
        global _output_line_count
        next_report_line = 0
        sample_labels = set()
        '''Next read must be known to tell if a read mapped to multiple
        locations, so always work with previous read.'''
        while True:
            line = self.input_stream.readline()
            if not line: return # Bowtie output nothing
            # Skip header line
            if line[0] == '@': continue
            last_tokens = line.rstrip().split('\t')
            (last_qname, last_flag, last_rname, last_pos, last_mapq,
                last_cigar, last_rnext, last_pnext,
                last_tlen, last_seq, last_qual) = last_tokens[:11]
            last_flag = int(last_flag)
            last_pos = int(last_pos)
            last_seq_size = len(last_seq)
            last_end_pos = last_pos + last_seq_size
            '''Find XM:i field, which is > 0 if read had several valid
            alignments, but all were suppressed because bowtie -m X was
            invoked for X an integer >= 1. See Bowtie documentation.'''
            last_multimapped = False
            for field in last_tokens[::-1]:
                if field[:5] == 'XM:i:':
                    if int(field[5:]) > 0 and (last_flag & 4):
                        '''If read is multimapped and all alignments were
                        suppressed.'''
                        last_multimapped = True
                    break
            break
        # Initialize counter
        i = 0
        # While labeled multiread, this list may end up simply a uniread
        multiread = []
        while True:
            line = self.input_stream.readline()
            if line:
                tokens = line.rstrip().split('\t')
                (qname, flag, rname, pos, mapq, cigar, rnext,
                    pnext, tlen, seq, qual) = tokens[:11]
                flag = int(flag)
                pos = int(pos)
                seq_size = len(seq)
                end_pos = pos + seq_size
                '''Find XM:i field, which is > 0 if read had several valid
                alignments, but all were suppressed because bowtie -m X was
                invoked for X an integer >= 1. See Bowtie documentation.'''
                multimapped = False
                for field in tokens[::-1]:
                    if field[:5] == 'XM:i:':
                        if int(field[5:]) > 0 and (flag & 4):
                            '''If read is multimapped and all alignments were
                            suppressed.'''
                            multimapped = True
                        break
            if self.verbose and next_report_line == i:
                print >>sys.stderr, \
                    'SAM output record %d: rdname="%s", flag=%d' % (i,
                                                                    last_qname,
                                                                    last_flag)
                next_report_line = int((next_report_line + 1)
                    * self.report_multiplier + 1) - 1
            # Use indices as sample labels rather than sample labels themselves
            last_sample_label = self.manifest_object.label_to_index[
                                        sample.parseLab(last_qname)
                                    ]
            sample_labels.add(last_sample_label)
            multiread.append(last_tokens)
            if not line or qname != last_qname:
                '''If the next qname doesn't match the last qname or there are
                no more lines, all of a multiread's alignments have been
                collected.'''
                if (last_flag & 4) \
                    and not (len(multiread) > 1 or last_multimapped):
                    '''Write unmapped reads for realignment in another map
                    step.'''
                    print >>self.output_stream, '%s\t%s\t%s\t%s' % ('unmapped',
                        last_qname, last_seq, last_qual)
                    _output_line_count += 1
                    '''Write unmapped sequences for readletizing and
                    inferring introns. To reduce alignment burden, find
                    reversed complement and only write sequence that comes
                    first in lexicographic order.'''
                    last_seq = last_seq.upper()
                    reversed_complement_last_seq = last_seq[::-1].translate(
                            _reversed_complement_translation_table
                        )
                    if last_seq < reversed_complement_last_seq:
                        self.output_stream.write('readletize\t%s\t%s\t\x1c\n' \
                            % (last_seq, last_sample_label))
                    else:
                        print >>self.output_stream, \
                            'readletize\t%s\t\x1c\t%s' \
                            % (reversed_complement_last_seq, last_sample_label)
                elif self.end_to_end_sam:
                    '''End-to-end SAM is output for every line with at least
                    one possible alignment.'''
                    for alignment_tokens in multiread:
                        print >>self.output_stream, (
                            ('%s\t%s\t%s\t%012d\t%s\t%s\t'
                                % ('end_to_end_sam', last_sample_label, 
                                self.reference_index.rname_to_string[
                                        alignment_tokens[2]
                                    ],
                                int(alignment_tokens[3]),
                                alignment_tokens[0][:-2],
                                alignment_tokens[1])) 
                            + '\t'.join(alignment_tokens[4:])
                            + ('\tNH:i:%d' % len(multiread)))
                        _output_line_count += 1
                if not (last_flag & 4) and len(multiread) == 1 \
                    and not last_multimapped:
                    '''Read maps uniquely; the full alignment is to be called
                    as an exonic chunk (EC).'''
                    partitions = partition.partition(last_rname, last_pos, 
                        last_end_pos, self.bin_size)
                    if self.exon_differentials:
                        for (partition_id, 
                            partition_start, partition_end) in partitions:
                            # Print increment at interval start
                            assert last_pos < partition_end
                            print >>self.output_stream, \
                                'exon_diff\t%s\t%s\t%012d\t1' \
                                % (partition_id,
                                    last_sample_label,
                                    max(partition_start, last_pos))
                            _output_line_count += 1
                            assert last_end_pos > partition_start
                            if last_end_pos < partition_end:
                                '''Print decrement at interval end iff exon
                                ends before partition ends.'''
                                print >>self.output_stream, \
                                    'exon_diff\t%s\t%s\t%012d\t-1' \
                                    % (partition_id,
                                        last_sample_label,
                                        last_end_pos)
                                _output_line_count += 1
                    if self.exon_intervals:
                        for partition_id, _, _ in partitions:
                            print >>self.output_stream, \
                                'exon_ival\t%s\t%012d\t%012d\t%s' \
                                % (partition_id, last_pos, 
                                    last_end_pos, last_sample_label)
                            _output_line_count += 1
                multiread = []
            if not line:
                break
            last_tokens = tokens
            (last_qname, last_flag, last_rname, last_pos, last_mapq,
                last_cigar, last_rnext, last_pnext, last_tlen, last_seq,
                last_qual) = (qname, flag, rname, pos, mapq, cigar,
                rnext, pnext, tlen, seq, qual)
            (last_seq_size, last_end_pos, last_multimapped) = (seq_size,
                end_pos, multimapped)
            i += 1

def go(input_stream=sys.stdin, output_stream=sys.stdout, bowtie_exe='bowtie',
    bowtie_index_base='genome', manifest_file='manifest', bowtie_args=None, 
    bin_size=10000, verbose=False, exon_differentials=True,
    exon_intervals=False, end_to_end_sam=True, report_multiplier=1.2):
    """ Runs Rail-RNA-align_reads.

        A single pass of Bowtie is run to find end-to-end alignments. Unmapped
        reads are saved for readletizing to determine introns in sucessive
        reduce steps as well as for realignment in a later map step.

        Input (read from stdin)
        ----------------------------
        Tab-delimited input tuple columns in a mix of any of the following
        three formats:
         Format 1 (single-end, 3-column):
          1. Name
          2. Nucleotide sequence
          3. Quality sequence
         Format 2 (paired-end, 5-column):
          1. Name
          2. Nucleotide sequence for mate 1
          3. Quality sequence for mate 1
          4. Nucleotide sequence for mate 2
          5. Quality sequence for mate 2
         Format 3 (paired, 6-column):
          1. Name for mate 1
          2. Nucleotide sequence for mate 1
          3. Quality sequence for mate 1
          4. Name for mate 2
          5. Nucleotide sequence for mate 2
          6. Quality sequence for mate 2

        Hadoop output (written to stdout)
        ----------------------------
        A given RNAME sequence is partitioned into intervals ("bins") of some 
        user-specified length (see partition.py).

        Exonic chunks (aka ECs; three formats, any or all of which may be
        emitted):

        Format 1 (exon_ival); tab-delimited output tuple columns:
        1. Reference name (RNAME in SAM format) + ';' + bin number
        2. Sample label
        3. EC start (inclusive) on forward strand
        4. EC end (exclusive) on forward strand

        Format 2 (exon_diff); tab-delimited output tuple columns:
        1. Reference name (RNAME in SAM format) + ';' + 
            max(EC start, bin start) (inclusive) on forward strand IFF diff is
            positive and EC end (exclusive) on forward strand IFF diff is
            negative
        2. Bin number
        3. Sample label
        4. +1 or -1.

        Note that only unique alignments are currently output as ivals and/or
        diffs.

        Format 3 (end_to_end_sam); tab-delimited output tuple columns:
        Standard SAM output except fields are in different order, and the first
        field corresponds to sample label. (Fields are reordered to facilitate
        partitioning by sample name/RNAME and sorting by POS.) Each line
        corresponds to a spliced alignment. The order of the fields is as
        follows.
        1. Sample label
        2. Number string representing RNAME; see BowtieIndexReference class in
            bowtie_index for conversion information
        3. POS
        4. QNAME
        5. FLAG
        6. MAPQ
        7. CIGAR
        8. RNEXT
        9. PNEXT
        10. TLEN
        11. SEQ
        12. QUAL
        ... + optional fields

        Reads with no end-to-end alignments

        Tab-delimited output tuple columns (unmapped):
        1. QNAME
        2. SEQ
        3. QUAL

        Tab-delimited output tuple columns (readletize)
        1. SEQ or its reversed complement, whichever is first in alphabetical
            order
        2. The sample label if field 1 is the read sequence; else '\x1c'
        3. The sample label if field 1 is the read sequence's reversed
            complement; else '\x1c'

        Maximum read lengths found

        Tab-delimited output tuple columns (max_len):
        1. RNAME + ('-' or '+'; all combinations are included)
        2. Sample label
        3. The character 'a', which places it before 'i' in 
            lexicograhic sort order for reading in Rail-RNA-intron_config
        4. Maximum read length found
        5. The character '-'.

        ALL OUTPUT COORDINATES ARE 1-INDEXED.

        input_stream: where to find input reads.
        output_stream: where to emit exonic chunks and introns.
        bowtie_exe: filename of Bowtie executable; include path if not in
            $PATH.
        bowtie_index_base: the basename of the Bowtie index files associated
            with the reference.
        manifest_file: filename of manifest
        bowtie_args: string containing precisely extra command-line arguments
            to pass to first-pass Bowtie, e.g., "--tryhard --best"; or None.
        bin_size: genome is partitioned in units of bin_size for later load
            balancing.
        verbose: True iff more informative messages should be written to
            stderr.
        exon_differentials: True iff EC differentials are to be emitted.
        exon_intervals: True iff EC intervals are to be emitted.
        end_to_end_sam: True iff SAM with end_to_end alignments should be
            output. See docstring for more information.
        report_multiplier: if verbose is True, the line number of an alignment
            or read written to stderr increases exponentially with base
            report_multiplier.

        No return value.
    """
    temp_dir = tempfile.mkdtemp()
    atexit.register(handle_temporary_directory, temp_dir)
    reads_file = os.path.join(temp_dir, 'reads.temp')
    reference_index = bowtie_index.BowtieIndexReference(bowtie_index_base)
    manifest_object = manifest.LabelsAndIndices(manifest_file)
    with open(reads_file, 'w') as read_stream:
        write_reads(
            read_stream, input_stream=input_stream, verbose=verbose,
            report_multiplier=report_multiplier
        )
    bowtie_process, bowtie_command, threads = bowtie.proc(
            bowtieExe=bowtie_exe, bowtieIdx=bowtie_index_base,
            readFn=reads_file, bowtieArgs=bowtie_args, sam=True,
            stdoutPipe=True, stdinPipe=False
        )
    output_thread = BowtieOutputThread(
            bowtie_process.stdout, reference_index, manifest_object,
            exon_differentials=exon_differentials, 
            exon_intervals=exon_intervals, 
            bin_size=bin_size,
            verbose=verbose, 
            output_stream=output_stream,
            end_to_end_sam=end_to_end_sam,
            report_multiplier=report_multiplier
        )
    threads.append(output_thread)
    output_thread.start()
    # Join threads to pause execution in main thread
    for thread in threads:
        if verbose: print >>sys.stderr, 'Joining thread...'
        thread.join()
    output_stream.flush()

if __name__ == '__main__':
    import argparse
    # Print file's docstring if -h is invoked
    parser = argparse.ArgumentParser(description=__doc__, 
                formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--report_multiplier', type=float, required=False,
        default=1.2,
        help='When --verbose is also invoked, the only lines of lengthy '
             'intermediate output written to stderr have line number that '
             'increases exponentially with this base')
    parser.add_argument('--verbose', action='store_const', const=True,
        default=False,
        help='Print out extra debugging statements')
    parser.add_argument('--exon-differentials', action='store_const',
        const=True,
        default=True, 
        help='Print exon differentials (+1s and -1s)')
    parser.add_argument('--exon-intervals', action='store_const',
        const=True,
        default=False, 
        help='Print exon intervals')
    parser.add_argument('--end-to-end-sam', action='store_const', const=True,
        default=True, 
        help='Print read-by-read SAM output of end-to-end alignments for '
             'later consolidation after realignment.')
    parser.add_argument('--manifest', type=str, required=False,
        default='manifest',
        help='Path to manifest file')
    parser.add_argument('--test', action='store_const', const=True,
        default=False,
        help='Run unit tests; DOES NOT NEED INPUT FROM STDIN, AND DOES NOT '
             'WRITE EXONS AND INTRONS TO STDOUT')

    # Add command-line arguments for dependencies
    partition.addArgs(parser)
    bowtie.addArgs(parser)

    # Collect Bowtie arguments, supplied in command line after the -- token
    argv = sys.argv
    bowtie_args = ''
    in_args = False
    for i, argument in enumerate(sys.argv[1:]):
        if in_args:
            bowtie_args += argument + ' '
        if argument == '--':
            argv = sys.argv[:i + 1]
            in_args = True

    '''Now collect other arguments. While the variable args declared below is
    global, properties of args are also arguments of the go() function so
    different command-line arguments can be passed to it for unit tests.'''
    args = parser.parse_args(argv[1:])

if __name__ == '__main__' and not args.test:
    import time
    start_time = time.time()
    go(bowtie_exe=args.bowtie_exe,
        bowtie_index_base=args.bowtie_idx,
        bowtie_args=bowtie_args,
        manifest_file=args.manifest,
        verbose=args.verbose, 
        bin_size=args.partition_length,
        exon_differentials=args.exon_differentials,
        exon_intervals=args.exon_intervals,
        end_to_end_sam=args.end_to_end_sam,
        report_multiplier=args.report_multiplier)
    print >>sys.stderr, 'DONE with align_reads.py; in/out=%d/%d; ' \
        'time=%0.3f s' % (_input_line_count, _output_line_count,
                            time.time() - start_time)
elif __name__ == '__main__':
    # Test units
    del sys.argv[1:] # Don't choke on extra command-line parameters
    import unittest

    class TestWriteReads(unittest.TestCase):
        """ Tests write_reads(). """
        def setUp(self):
            '''input_reads has 5-, 6-, and 3-token examples, each on a
            different line. Every read has 100 bases.'''
            input_reads = \
                'r1;LB:holder' \
                    '\tTTACATACCATACAGTGCGCTAGCGGGTGACAGATATAATGCAGATCCAT' \
                      'ACAGACCAGATGGCAGACATGTGTTGCAGSCTGCAAGTGCAACGCGGTGA' \
                    '\tFFB<9889340///29==:766234466666340///29==:76623446' \
                      '744442<<9889<888@?FFFFFFFFFFFFDB<4444340///29==:76' \
                    '\tGACCAGAGGTGCACCAAATGACAGTGCGCCACAGATGGCGCATGCAGTGG' \
                      'CCAGAAACGTGCGACCGATGACAGGTGCACAATGCCGCTGACGACGTGAA' \
                    '\t444744442<<9889<8880///29==:766230///29==:766230//' \
                      '<9889340///29==:766234466666340///<9889340///29==:\n' \
                'r2/1;LB:holder' \
                    '\tGCAGAGTGCCGCAATGACGTGCGCCAAAGCGGTGACAGGGTGACAGTGAA' \
                      'CCAAGTGACAAGTGAACAGGTGCCAGAGTGACCGAGTGACCAGTGGACCA' \
                    '\t442<<9889<8880///29==:766230//442<<9889<8880///29=' \
                      '44442<<9889<8880///29==:76623044442<<9889<8880///2' \
                '\tr2/2;LB:holder' \
                    '\tCAGAGTGCCGCAATGACGTGCGCCAAAGCGGACAAAGCACCATGACAAGT' \
                      'ACACAGGTGACAGTGACAAGACAGAGGTGACACAGAGAAAGtGGGTGTGA' \
                    '\t<<9889<8880///29==:766230//442<<<<9889<8880///29==' \
                      '44442<<9889<8880///29==:766230///2944442<<9889<888\n' \
                'r3;LB:holder' \
                    '\tATCGATTAAGCTATAACAGATAACATAGACATTGCGCCCATAATAGATAA' \
                      'CTGACACCTGACCAGTGCCAGATGACCAGTGCCAGATGGACGACAGTAGC' \
                    '\tFFFFFFFFFFFFDB<4444340///29==:766234466666777689<3' \
                      '44=<<;444744442<<9889<888@?FFFFFFFFFFFFDB<4444340/\n'
            self.temp_dir_path = tempfile.mkdtemp()
            self.input_file = os.path.join(self.temp_dir_path,
                                'sample_input.tsv')
            # Store reads in temporary file
            with open(self.input_file, 'w') as reads_stream:
                reads_stream.write(input_reads)
            # For storing output of write_reads()
            self.output_file = os.path.join(self.temp_dir_path,
                                'sample_output.tsv')

        def test_output(self):
            """ Fails if output of write_reads() is not in the right form. """
            with open(self.input_file) as input_stream:
                with open(self.output_file, 'w') as output_stream:
                    write_reads(output_stream, input_stream=input_stream)

            # Verify output
            with open(self.output_file) as processed_stream:
                first_read_from_pair = \
                    processed_stream.readline().rstrip().split('\t')
                second_read_from_pair = \
                    processed_stream.readline().rstrip().split('\t')
                self.assertEquals(
                    [
                        'r1;LB:holder',
                        'TTACATACCATACAGTGCGCTAGCGGGTGACAGATATAATGCAGATCCAT'
                        'ACAGACCAGATGGCAGACATGTGTTGCAGSCTGCAAGTGCAACGCGGTGA',
                        'FFB<9889340///29==:766234466666340///29==:76623446'
                        '744442<<9889<888@?FFFFFFFFFFFFDB<4444340///29==:76'
                    ],
                    first_read_from_pair
                )
                self.assertEquals(
                    [
                        'r1;LB:holder',
                        'GACCAGAGGTGCACCAAATGACAGTGCGCCACAGATGGCGCATGCAGTGG'
                        'CCAGAAACGTGCGACCGATGACAGGTGCACAATGCCGCTGACGACGTGAA',
                        '444744442<<9889<8880///29==:766230///29==:766230//'
                        '<9889340///29==:766234466666340///<9889340///29==:'
                    ],
                    second_read_from_pair
                )
                first_read_from_pair = \
                    processed_stream.readline().rstrip().split('\t')
                second_read_from_pair = \
                    processed_stream.readline().rstrip().split('\t')
                self.assertEquals(
                    [
                        'r2/1;LB:holder',
                        'GCAGAGTGCCGCAATGACGTGCGCCAAAGCGGTGACAGGGTGACAGTGAA'
                        'CCAAGTGACAAGTGAACAGGTGCCAGAGTGACCGAGTGACCAGTGGACCA',
                        '442<<9889<8880///29==:766230//442<<9889<8880///29='
                        '44442<<9889<8880///29==:76623044442<<9889<8880///2'
                    ],
                    first_read_from_pair
                )
                self.assertEquals(
                    [
                        'r2/2;LB:holder',
                        'CAGAGTGCCGCAATGACGTGCGCCAAAGCGGACAAAGCACCATGACAAGT'
                        'ACACAGGTGACAGTGACAAGACAGAGGTGACACAGAGAAAGtGGGTGTGA',
                        '<<9889<8880///29==:766230//442<<<<9889<8880///29=='
                        '44442<<9889<8880///29==:766230///2944442<<9889<888'
                    ],
                    second_read_from_pair
                )
                self.assertEquals(
                    [
                        'r3;LB:holder',
                        'ATCGATTAAGCTATAACAGATAACATAGACATTGCGCCCATAATAGATAA'
                        'CTGACACCTGACCAGTGCCAGATGACCAGTGCCAGATGGACGACAGTAGC',
                        'FFFFFFFFFFFFDB<4444340///29==:766234466666777689<3'
                        '44=<<;444744442<<9889<888@?FFFFFFFFFFFFDB<4444340/'
                    ],
                    processed_stream.readline().rstrip().split('\t')
                )

        def tearDown(self):
            # Kill temporary directory
            shutil.rmtree(self.temp_dir_path)
   
    unittest.main()