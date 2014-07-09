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

Introns/insertions/deletions

(bed); tab-delimited output tuple columns (bed):
1. 'I' or 'D' insertion or deletion line
2. Sample label
3. Number string representing RNAME
4. Start position (Last base before insertion or first base of deletion)
5. End position (Last base before insertion or last base of deletion 
                    (exclusive))
6. Inserted sequence for insertions or deleted sequence for deletions
----Next fields are for introns only; they are '\x1c' for indels----
7. '\x1c'
8. '\x1c'
--------------------------------------------------------------------
9. Number of instances of insertion or deletion in sample; this is
    always +1 before bed_pre combiner/reducer

Read whose primary alignment is not end-to-end

Tab-delimited output tuple columns (unmapped):
1. SEQ
2. QNAME
3. QUAL

Tab-delimited output tuple columns (readletize)
1. SEQ or its reversed complement, whichever is first in alphabetical order
2. The sample label if field 1 is the read sequence; else '\x1c'
3. The sample label if field 1 is the read sequence's reversed complement;
    else '\x1c'

Tab-delimited tuple columns (fasta):
1. Read sequence
2. '\x1c' + FASTA reference name including '>'. The following format is used:
    original RNAME + '\x1d' + start position of sequence + '\x1d\x1d\x1d' + 'p'
    if derived from primary alignment to genome; 's' if derived from secondary
    alignment to genome; 'i' if derived from cointron search
3. FASTA sequence

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
import time

base_path = os.path.abspath(
                    os.path.dirname(os.path.dirname(os.path.dirname(
                        os.path.realpath(__file__)))
                    )
                )
utils_path = os.path.join(base_path, 'rna', 'utils')
site.addsitedir(utils_path)
site.addsitedir(base_path)

import bowtie
import bowtie_index
import partition
import manifest
from cigar_parse import indels_introns_and_exons, reference_from_seq
from dooplicity.tools import xstream

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
        if len(tokens) == 3:
            to_write = '%s\t%s\t%s' \
                % (tokens[0], tokens[1], tokens[2])
        else:
            to_write = '%s\t%s\t%s\n%s\t%s\t%s' \
                % (tokens[0], tokens[1], tokens[2],
                    (tokens[0] if len(tokens) == 5 else tokens[3]),
                        tokens[-2], tokens[-1])
            assert len(tokens[1]) == len(tokens[2])
            assert len(tokens[-2]) == len(tokens[-1])
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
                that permits access to reference
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
        i = 0
        for (qname,), xpartition in xstream(self.input_stream, 1):
            # While labeled multiread, this list may end up simply a uniread
            multiread = []
            sample_label = self.manifest_object.label_to_index[
                        qname.rpartition('\x1d')[2]
                    ]
            # Handle primary alignment
            rest_of_line = xpartition.next()
            i += 1
            flag = int(rest_of_line[0])
            cigar = rest_of_line[4]
            rname = rest_of_line[1]
            pos = int(rest_of_line[2])
            seq, qual = rest_of_line[8], rest_of_line[9]
            if self.verbose and next_report_line == i:
                print >>sys.stderr, \
                    'SAM output record %d: rdname="%s", flag=%d' \
                    % (i, qname, flag)
                next_report_line = max(int(next_report_line
                    * self.report_multiplier), next_report_line + 1)
            multiread.append((qname,) + rest_of_line)
            for rest_of_line in xpartition:
                multiread.append((qname,) + rest_of_line)
            if flag & 4 or 'S' in cigar:
                '''Write unmapped/soft-clipped primary alignments for
                realignment in a reduce step.'''
                seq = seq.upper()
                if flag & 16:
                    '''If it's reverse-complemented, write seq the way it was
                    found.'''
                    seq = seq[::-1].translate(
                        _reversed_complement_translation_table
                    )
                    qual = qual[::-1]
                print >>self.output_stream, '%s\t%s\t%s\t%s' % ('unmapped',
                    seq, qname, qual)
                _output_line_count += 1
                '''To reduce alignment burden, find reversed complement and
                only write sequence that comes first in lexicographic order.'''
                reversed_complement_seq = seq[::-1].translate(
                        _reversed_complement_translation_table
                    )
                if seq < reversed_complement_seq:
                    print >>self.output_stream, \
                        'readletize\t%s\t\x1c\t%s' % (seq, sample_label)
                else:
                    print >>self.output_stream, \
                        'readletize\t%s\t\x1c\t%s' \
                        % (reversed_complement_seq, sample_label)
                _output_line_count += 1
                '''Write fasta lines for realignment; they will be combined
                with fasta lines for transcriptome elements obtained in later
                reduce step.'''
                try:
                    md = [field for field in multiread[0]
                          if field[:5] == 'MD:Z:'][0][5:]
                except IndexError:
                    # Unmapped read
                    continue
                print >>self.output_stream, \
                    'fasta\t%s\t\x1c>%s\x1d%s\x1d\x1d\x1dp\t%s' \
                    % (multiread[0][9], multiread[0][2], multiread[0][3],
                        reference_from_seq(multiread[0][5],
                                            md, multiread[0][9]))
                try:
                    for alignment in multiread[1:]:
                        md = [field for field in alignment
                              if field[:5] == 'MD:Z:'][0][5:]
                        print >>self.output_stream, \
                            'fasta\t%s\t\x1c>%s\x1d%s\x1d\x1d\x1ds\t%s' \
                            % (alignment[9], alignment[2], alignment[3],
                                reference_from_seq(alignment[5], md,
                                alignment[9]))
                except IndexError:
                    pass
            else:
                if self.end_to_end_sam:
                    '''End-to-end SAM is output for every line with at least
                    one possible alignment.'''
                    for alignment in multiread:
                        print >>self.output_stream, (
                            ('%s\t%s\t%s\t%012d\t%s\t%s\t'
                                % ('end_to_end_sam', sample_label, 
                                    self.reference_index.rname_to_string[
                                        alignment[2]
                                    ],
                                int(alignment[3]),
                                alignment[0].partition('\x1d')[0],
                                alignment[1])) 
                            + '\t'.join(alignment[4:])
                            + ('\tNH:i:%d' % len(multiread)))
                        _output_line_count += 1
                # Output exonic chunks and indels for primary alignment
                alignment = multiread[0]
                md = [field for field in alignment
                            if field[:5] == 'MD:Z:'][0][5:]
                insertions, deletions, _, exons \
                                            = indels_introns_and_exons(
                                                cigar, md, pos, seq
                                            )
                # Output indels
                for insert_pos, insert_seq in insertions:
                    print >>self.output_stream, (
                           ('bed\tI\t%s\t%s\t%012d\t%012d\t%s'
                            '\t\x1c\t\x1c\t1')
                            % (sample_label, self.reference_index.\
                                rname_to_string[rname],
                                insert_pos, insert_pos,
                                insert_seq)
                        )
                    _output_line_count += 1
                for del_pos, del_seq in deletions:
                    print >>self.output_stream, (
                           ('bed\tD\t%s\t%s\t%012d\t%012d\t%s'
                            '\t\x1c\t\x1c\t1')
                            % (sample_label, self.reference_index.\
                                rname_to_string[rname],
                                del_pos, del_pos + len(del_seq),
                                del_seq)
                        )
                    _output_line_count += 1
                # Output exonic chunks
                if self.exon_intervals:
                    for exon_pos, exon_end_pos in exons:
                        partitions = partition.partition(
                                rname, exon_pos, exon_end_pos,
                                self.bin_size)
                        for partition_id, _, _ in partitions:
                            print >>self.output_stream, \
                                'exon_ival\t%s\t%012d\t' \
                                '%012d\t%s' \
                                % (partition_id,
                                    exon_pos, exon_end_pos, 
                                    sample_label)
                            _output_line_count += 1
                if self.exon_differentials:
                    for exon_pos, exon_end_pos in exons:
                        partitions = partition.partition(
                            rname, exon_pos, exon_end_pos,
                            self.bin_size)
                        for (partition_id, partition_start, 
                                partition_end) in partitions:
                            assert exon_pos < partition_end
                            # Print increment at interval start
                            print >>self.output_stream, \
                                'exon_diff\t%s\t%s\t%012d\t1' \
                                % (partition_id,
                                    sample_label,
                                    max(partition_start, exon_pos))
                            _output_line_count += 1
                            assert exon_end_pos > partition_start
                            if exon_end_pos < partition_end:
                                '''Print decrement at interval end 
                                iff exon ends before partition
                                ends.'''
                                print >>self.output_stream, \
                                    'exon_diff\t%s\t%s\t' \
                                    '%012d\t-1' \
                                    % (partition_id, 
                                        sample_label,
                                        exon_end_pos)
                                _output_line_count += 1

def go(input_stream=sys.stdin, output_stream=sys.stdout, bowtie2_exe='bowtie2',
    bowtie_index_base='genome', bowtie2_index_base='genome2', 
    manifest_file='manifest', bowtie2_args=None, bin_size=10000, verbose=False,
    exon_differentials=True, exon_intervals=False, end_to_end_sam=True,
    report_multiplier=1.2, keep_alive=False):
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

        Read whose primary alignment is not end-to-end

        Tab-delimited output tuple columns (unmapped):
        1. SEQ
        2. QNAME
        3. QUAL

        Tab-delimited output tuple columns (readletize)
        1. SEQ or its reversed complement, whichever is first in alphabetical
            order
        2. The sample label if field 1 is the read sequence; else '\x1c'
        3. The sample label if field 1 is the read sequence's reversed
            complement; else '\x1c'

        Tab-delimited tuple columns (fasta):
        1. Read sequence
        2. '\x1c' + FASTA reference name including '>'. The following format is
            used:
            original RNAME + '\x1d' + start position of sequence
            + '\x1d\x1d\x1d' + 'p' if derived from primary alignment to genome;
            's' if derived from secondary alignment to genome; 'i' if derived
            from cointron search
        3. FASTA sequence

        ALL OUTPUT COORDINATES ARE 1-INDEXED.

        input_stream: where to find input reads.
        output_stream: where to emit exonic chunks and introns.
        bowtie2_exe: filename of Bowtie2 executable; include path if not in
            $PATH.
        bowtie_index_base: the basename of the Bowtie1 index files associated
            with the reference.
        bowtie2_index_base: the basename of the Bowtie2 index files associated
            with the reference.
        manifest_file: filename of manifest
        bowtie2_args: string containing precisely extra command-line arguments
            to pass to first-pass Bowtie2.
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
        keep_alive: True iff "reporter:status:alive" should be printed to
            stderr periodically to keep job alive while Bowtie is running

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
    output_file = os.path.join(temp_dir, 'out.sam')
    bowtie_command = ' '.join([bowtie2_exe,
        bowtie2_args if bowtie2_args is not None else '',
        ' --local -t --no-hd --mm -x', bowtie2_index_base, '--12',
        reads_file, '-S', output_file])
    print >>sys.stderr, 'Starting Bowtie2 with command: ' + bowtie_command
    # Because of problems with buffering, write output to file
    bowtie_process = subprocess.Popen(bowtie_command, bufsize=-1,
        stdout=subprocess.PIPE, stderr=sys.stderr, shell=True)
    if keep_alive:
        period_start = time.time()
        while bowtie_process.poll() is None:
            now = time.time()
            if now - period_start > 60:
                print >>sys.stderr, '\nreporter:status:alive'
                period_start = now
            time.sleep(.2)
    else:
        bowtie_process.wait()
    output_thread = BowtieOutputThread(
            open(output_file), reference_index, manifest_object,
            exon_differentials=exon_differentials, 
            exon_intervals=exon_intervals, 
            bin_size=bin_size,
            verbose=verbose, 
            output_stream=output_stream,
            end_to_end_sam=end_to_end_sam,
            report_multiplier=report_multiplier
        )
    output_thread.start()
    # Join thread to pause execution in main thread
    if verbose: print >>sys.stderr, 'Joining thread...'
    output_thread.join()

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
    parser.add_argument('--keep-alive', action='store_const', const=True,
        default=False,
        help='Periodically print Hadoop status messages to stderr to keep ' \
             'job alive')
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
    partition.add_args(parser)
    bowtie.add_args(parser)

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
    start_time = time.time()
    go(bowtie2_exe=args.bowtie2_exe,
        bowtie_index_base=args.bowtie_idx,
        bowtie2_index_base=args.bowtie2_idx,
        bowtie2_args=bowtie_args,
        manifest_file=args.manifest,
        verbose=args.verbose, 
        bin_size=args.partition_length,
        exon_differentials=args.exon_differentials,
        exon_intervals=args.exon_intervals,
        end_to_end_sam=args.end_to_end_sam,
        report_multiplier=args.report_multiplier,
        keep_alive=args.keep_alive)
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
