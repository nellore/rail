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
  1. Nucleotide sequence or its reversed complement, whichever is first in 
    alphabetical order
  2. 1 if sequence was reverse-complemented else 0
  3. Name
  4. Quality sequence or its reverse, whichever corresponds to field 1

Format 2 (paired, 2 lines, 3 columns each)
(so this is the same as single-end)
  1. Nucleotide sequence for mate 1 or its reversed complement, whichever is
    first in alphabetical order
  2. 1 if sequence was reverse-complemented else 0
  3. Name for mate 1
  4. Quality sequence for mate 1 or its reverse, whichever corresponds to
    field 1
    
    (new line)

  1. Nucleotide sequence for mate 2 or its reversed complement, whichever is
    first in alphabetical order
  2. 1 if sequence was reverse complemented else 0
  3. Name for mate 2
  4. Quality sequence for mate 2 or its reverse, whichever corresponds to
    field 1

Input is partitioned and sorted by field 1, the read sequence.

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
4. +1 or -1 * count, the number of instances of a read sequence for which to
    print exonic chunks

Note that only unique alignments are currently output as ivals and/or diffs.

Format 3 (sam); tab-delimited output tuple columns:
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

Insertions/deletions (indel_bed)

Tab-delimited output tuple columns:
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
    always +1 * count before bed_pre combiner/reducer

Read whose primary alignment is not end-to-end

Tab-delimited output tuple columns (unmapped):
1. SEQ
2. 1 if SEQ is reverse-complemented, else 0
3. QNAME
4. QUAL

Tab-delimited output tuple columns (readletize)
1. SEQ or its reversed complement, whichever is first in alphabetical order
2. '\x1d'-separated list of sample indexes for which field 1 is the read
    sequence else '\x1c'
3. '\x1d'-separated list of sample indexes for which field 1 is the read
    sequence's reversed complement else '\x1c'

Tab-delimited tuple columns (postponed_sam):
Standard 11+ -column raw SAM output

Single column (unique):
1. A unique read sequence

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
import copy

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
from alignment_handlers import AlignmentPrinter
from dooplicity.tools import xstream
import re

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

def handle_bowtie_output(input_stream, reference_index, manifest_object,
        k_value=1, align_stream=None, other_stream=None,
        output_stream=sys.stdout, exon_differentials=True,
        exon_intervals=False, verbose=False, bin_size=10000,
        report_multiplier=1.2, min_exon_size=8):
    """ Prints end-to-end alignments and selects reads to be realigned.

        input_stream: where to retrieve Bowtie's SAM output, typically a
            Bowtie process's stdout.
        reference_index: object of class bowtie_index.BowtieIndexReference
            that permits access to reference
        manifest_object: object of class LabelsAndIndices that maps indices
            to labels and back; used to shorten intermediate output.
        align_stream: where to write reads that are to be aligned on second
            pass
        other_stream: where reads with the same sequence as those from the
            input stream are stored; if the primary alignment of a read from
            the input stream is not a tie and the edit distance (NM:i) is 0,
            these reads are assigned the same alignments. If None, second-pass
            alignment is being performed because there were ties or a nonzero
            edit distance during first-pass alignment.
        output_stream: where to emit exon and intron tuples; typically,
            this is sys.stdout.
        k_value: argument of Bowtie 2's -k parameter
        exon_differentials: True iff EC differentials are to be emitted.
        exon_intervals: True iff EC intervals are to be emitted.
        verbose: True if alignments should occasionally be written 
            to stderr.
        bin_size: genome is partitioned in units of bin_size for later load
            balancing.
        report_multiplier: if verbose is True, the line number of an
            alignment written to stderr increases exponentially with base
            report_multiplier.
        min_exon_size: minimum exon size searched for in intron_search.py
            later in pipeline; used to determine how large a soft clip on
            one side of a read is necessary to pass it on to intron search
            pipeline

        No return value.
    """
    global _output_line_count
    next_report_line = 0
    i = 0
    alignment_printer = AlignmentPrinter(
            manifest_object,
            reference_index,
            bin_size=bin_size,
            output_stream=output_stream,
            exon_ivals=exon_intervals,
            exon_diffs=exon_differentials
        )
    if other_stream:
        # First-pass alignment
        if verbose:
            print >>sys.stderr, 'Processing first-pass alignments.'
        other_xstream = xstream(other_stream, 1)
        for (qname,), xpartition in xstream(input_stream, 1):
            is_reverse, _, qname = qname.partition('\x1d')
            is_reverse = int(is_reverse)
            _, other_xpartition = other_xstream.next()
            sample_label = manifest_object.label_to_index[
                        qname.rpartition('\x1d')[2]
                    ]
            # Handle primary alignment
            rest_of_line = xpartition.next()
            i += 1
            flag = int(rest_of_line[0])
            if is_reverse:
                true_flag = flag ^ 16
            else:
                true_flag = flag
            cigar = rest_of_line[4]
            rname = rest_of_line[1]
            pos = int(rest_of_line[2])
            seq, qual = rest_of_line[8], rest_of_line[9]
            if verbose and next_report_line == i:
                print >>sys.stderr, \
                    'SAM output record %d: rdname="%s", flag=%d' \
                    % (i, qname, true_flag)
                next_report_line = max(int(next_report_line
                    * report_multiplier), next_report_line + 1)
            if flag & 4:
                print >>output_stream, 'unique\t%s' % seq
                print >>output_stream, 'unmapped\t%s\t%d\t%s\t%s' % (
                                                            seq,
                                                            is_reverse,
                                                            qname,
                                                            qual
                                                        )
                try:
                    for is_reverse, qname, qual in other_xpartition:
                        print >>output_stream, 'unmapped\t%s\t%s\t%s\t%s' % (
                                                                seq,
                                                                is_reverse,
                                                                qname,
                                                                qual
                                                            )
                except ValueError:
                    pass
                continue
            multiread = [(qname,) + rest_of_line] + \
                [(qname,) + next_line for next_line in xpartition]
            scores = \
                [field[5:] for field in multiread if field[:5] == 'XS:i:'] + \
                [field[5:] for field in multiread if field[:5] == 'AS:i:']
            tie_present = (len(scores) == 2) and (scores[0] == scores[1])
            clip_present = 'S' in cigar
            exact_match = 'NM:i:0' in multiread[0]
            alignment_reversed = [
                                    int(alignment[1]) & 16
                                    for alignment in multiread
                                ]
            if exact_match and not clip_present:
                # All alignments for all read sequences can be written
                NH_field = 'NH:i:%d' % len(multiread)
                reversed_flags = None
                '''First set count=0 to avoid printing exon_diffs and indels;
                then set count = k + 1 to print presummed exon_diffs and indels
                so fewer lines are written to disk.'''
                k = 0
                if k_value == 1 and not tie_present:
                    try:
                        for current_is_reverse, current_qname, current_qual \
                            in other_xpartition:
                            k += 1
                            if current_is_reverse == '1':
                                if not reversed_flags:
                                    reversed_flags = [
                                            str(int(alignment[1]) ^ 16)
                                            for alignment in multiread
                                        ]
                                _output_line_count += \
                                    alignment_printer.print_alignment_data(
                                            ([(current_qname,
                                                reversed_flags[j])
                                               + alignment[2:10]
                                               + (current_qual[::-1]
                                                    if alignment_reversed[j]
                                                    else current_qual,)
                                               + tuple([
                                                    field for field in
                                                    alignment[11:] if 
                                                    field[:5] != 'XS:i:'
                                                ])
                                               + (NH_field,)
                                               for j, alignment
                                               in enumerate(multiread)],),
                                            count=0
                                        )
                            else:
                                _output_line_count += \
                                    alignment_printer.print_alignment_data(
                                            ([(current_qname,)
                                                + alignment[1:10]
                                                + (current_qual[::-1]
                                                  if alignment_reversed[j]
                                                  else current_qual,)
                                                + tuple([
                                                    field for field in
                                                    alignment[11:] if 
                                                    field[:5] != 'XS:i:'
                                                ])
                                               + (NH_field,)
                                               for j, alignment
                                               in enumerate(multiread)],),
                                            count=0
                                        )
                    except ValueError:
                        pass
                if is_reverse:
                    _output_line_count += \
                        alignment_printer.print_alignment_data(
                                ([(alignment[0], str(int(alignment[1]) ^ 16))
                                   + alignment[2:] + (NH_field,)
                                   for alignment in multiread],),
                                count=(k+1)
                            )
                else:
                    _output_line_count += \
                        alignment_printer.print_alignment_data(
                                ([alignment + (NH_field,)
                                   for alignment in multiread],),
                                count=(k+1)
                            )
            if k_value != 1 \
                or tie_present or not (exact_match and not clip_present):
                # Prepare for second-pass Bowtie 2
                reversed_complement_seq = seq[::-1].translate(
                    _reversed_complement_translation_table
                )
                if seq < reversed_complement_seq:
                    seq_to_print = seq
                    qual_to_print = qual
                else:
                    seq_to_print = reversed_complement_seq
                    qual_to_print = qual[::-1]
                try:
                    for current_is_reverse, current_qname, current_qual \
                        in other_xpartition:
                        print >>align_stream, '\t'.join([
                                '%s\x1d%s' % (
                                        current_is_reverse, current_qname
                                    ), seq_to_print, current_qual
                            ])
                except ValueError:
                    # No other reads
                    pass
            if not (exact_match and not clip_present):
                print >>output_stream, 'unique\t%s' % seq_to_print
                # Write "postponed" SAM/readletize/unmapped lines
                split_cigar = re.split(r'([MINDS])', cigar)[:-1]
                try:
                    if ((split_cigar[1] == 'S'
                            and int(split_cigar[0]) >= min_exon_size) or
                        (split_cigar[-1] == 'S'
                            and int(split_cigar[-2]) >= min_exon_size)):
                        search_for_introns = True
                    else:
                        search_for_introns = False
                except IndexError:
                    search_for_introns = False
                if is_reverse:
                    for alignment in multiread:
                        print >>output_stream, \
                            '\t'.join(('postponed_sam', alignment[0])
                                + (str(int(alignment[1]) ^ 16),)
                                + alignment[2:])
                        _output_line_count += 1
                else:
                    for alignment in multiread:
                        print >>output_stream, \
                            '\t'.join(('postponed_sam',) + alignment)
                        _output_line_count += 1
                if search_for_introns:
                    if true_flag & 16:
                        print >>output_stream, \
                                'readletize\t%s\t\x1c\t%s' % (
                                        seq_to_print,
                                        sample_label
                                    )
                    else:
                        print >>output_stream, \
                                'readletize\t%s\t%s\t\x1c' % (
                                        seq_to_print,
                                        sample_label
                                    )
                    _output_line_count += 1
                print >>output_stream, 'unmapped\t%s\t%d\t%s\t%s' % (
                                                            seq_to_print,
                                                            is_reverse,
                                                            qname,
                                                            qual_to_print
                                                        )
    else:
        # Second-pass alignment
        if verbose:
            print >>sys.stderr, 'Processing second-pass alignments.'
        for (qname,), xpartition in xstream(input_stream, 1):
            is_reverse, _, qname = qname.partition('\x1d')
            is_reverse = int(is_reverse)
            sample_label = manifest_object.label_to_index[
                        qname.rpartition('\x1d')[2]
                    ]
            # Handle primary alignment
            rest_of_line = xpartition.next()
            i += 1
            flag = int(rest_of_line[0])
            if is_reverse:
                true_flag = flag ^ 16
            else:
                true_flag = flag
            cigar = rest_of_line[4]
            rname = rest_of_line[1]
            pos = int(rest_of_line[2])
            seq, qual = rest_of_line[8], rest_of_line[9]
            if verbose and next_report_line == i:
                print >>sys.stderr, \
                    'SAM output record %d: rdname="%s", flag=%d' \
                    % (i, qname, true_flag)
                next_report_line = max(int(next_report_line
                    * report_multiplier), next_report_line + 1)
            if flag & 4:
                print >>output_stream, 'unmapped\t%s\t%d\t%s\t%s' % (
                                                            seq,
                                                            is_reverse,
                                                            qname,
                                                            qual
                                                        )
                continue
            multiread = [(qname,) + rest_of_line] + \
                [(qname,) + next_line for next_line in xpartition]
            clip_present = 'S' in cigar
            exact_match = 'NM:i:0' in multiread[0]
            if exact_match and not clip_present:
                # All alignments for all read sequences can be written
                NH_field = 'NH:i:%d' % len(multiread)
                if is_reverse:
                    _output_line_count += \
                        alignment_printer.print_alignment_data(
                                ([(alignment[0], str(int(alignment[1]) ^ 16))
                                   + alignment[2:] + (NH_field,)
                                   for alignment in multiread],),
                                count=1
                            )
                else:
                    _output_line_count += \
                        alignment_printer.print_alignment_data(
                                ([alignment + (NH_field,)
                                   for alignment in multiread],),
                                count=1
                            )
            else:
                # Write "postponed" SAM/readletize/unmapped lines
                reversed_complement_seq = seq[::-1].translate(
                    _reversed_complement_translation_table
                )
                if seq < reversed_complement_seq:
                    seq_to_print = seq
                    qual_to_print = qual
                else:
                    seq_to_print = reversed_complement_seq
                    qual_to_print = qual[::-1]
                split_cigar = re.split(r'([MINDS])', cigar)[:-1]
                try:
                    if ((split_cigar[1] == 'S'
                            and int(split_cigar[0]) >= min_exon_size) or
                        (split_cigar[-1] == 'S'
                            and int(split_cigar[-2]) >= min_exon_size)):
                        search_for_introns = True
                    else:
                        search_for_introns = False
                except IndexError:
                    search_for_introns = False
                if is_reverse:
                    for alignment in multiread:
                        print >>output_stream, \
                            '\t'.join(('postponed_sam', alignment[0])
                                + (str(int(alignment[1]) ^ 16),)
                                + alignment[2:])
                        _output_line_count += 1
                else:
                    for alignment in multiread:
                        print >>output_stream, \
                            '\t'.join(('postponed_sam',) + alignment)
                        _output_line_count += 1
                if search_for_introns:
                    if true_flag & 16:
                        print >>output_stream, \
                                'readletize\t%s\t\x1c\t%s' % (
                                        seq_to_print,
                                        sample_label
                                    )
                    else:
                        print >>output_stream, \
                                'readletize\t%s\t%s\t\x1c' % (
                                        seq_to_print,
                                        sample_label
                                    )
                    _output_line_count += 1
                print >>output_stream, 'unmapped\t%s\t%d\t%s\t%s' % (
                                                            seq_to_print,
                                                            is_reverse,
                                                            qname,
                                                            qual_to_print
                                                        )
    output_stream.flush()

def go(input_stream=sys.stdin, output_stream=sys.stdout, bowtie2_exe='bowtie2',
    bowtie_index_base='genome', bowtie2_index_base='genome2', 
    manifest_file='manifest', bowtie2_args=None, bin_size=10000, verbose=False,
    exon_differentials=True, exon_intervals=False, report_multiplier=1.2,
    min_exon_size=8):
    """ Runs Rail-RNA-align_reads.

        A single pass of Bowtie is run to find end-to-end alignments. Unmapped
        reads are saved for readletizing to determine introns in sucessive
        reduce steps as well as for realignment in a later map step.

        Input (read from stdin)
        ----------------------------
        Tab-delimited input tuple columns in a mix of any of the following
        three formats:
        Format 1 (single-end, 3-column):
          1. Nucleotide sequence or its reversed complement, whichever is first
            in alphabetical order
          2. 1 if sequence was reverse-complemented else 0
          3. Name
          4. Quality sequence or its reverse, whichever corresponds to field 1

        Format 2 (paired, 2 lines, 3 columns each)
        (so this is the same as single-end)
          1. Nucleotide sequence for mate 1 or its reversed complement,
            whichever is first in alphabetical order
          2. 1 if sequence was reverse-complemented else 0
          3. Name for mate 1
          4. Quality sequence for mate 1 or its reverse, whichever corresponds
            to field 1
            
            (new line)

          1. Nucleotide sequence for mate 2 or its reversed complement,
            whichever is first in alphabetical order
          2. 1 if sequence was reverse complemented else 0
          3. Name for mate 2
          4. Quality sequence for mate 2 or its reverse, whichever corresponds
            to field 1

        Input is partitioned and sorted by field 1, the read sequence.

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
        4. +1 or -1 * count, the number of instances of a read sequence for
            which to print exonic chunks

        Note that only unique alignments are currently output as ivals and/or
        diffs.

        Format 3 (sam); tab-delimited output tuple columns:
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

        Insertions/deletions (indel_bed)

        tab-delimited output tuple columns:
        1. 'I' or 'D' insertion or deletion line
        2. Sample label
        3. Number string representing RNAME
        4. Start position (Last base before insertion or 
            first base of deletion)
        5. End position (Last base before insertion or last base of deletion 
                            (exclusive))
        6. Inserted sequence for insertions or deleted sequence for deletions
        ----Next fields are for introns only; they are '\x1c' for indels----
        7. '\x1c'
        8. '\x1c'
        --------------------------------------------------------------------
        9. Number of instances of insertion or deletion in sample; this is
            always +1 * count before bed_pre combiner/reducer

        Read whose primary alignment is not end-to-end

        Tab-delimited output tuple columns (unmapped):
        1. SEQ
        2. 1 if SEQ is reverse-complemented, else 0
        3. QNAME
        4. QUAL

        Tab-delimited output tuple columns (readletize)
        1. SEQ or its reversed complement, whichever is first in alphabetical
            order
        2. '\x1d'-separated list of sample indexes for which field 1 is the
            read sequence else '\x1c'
        3. '\x1d'-separated list of sample indexes for which field 1 is the
            read sequence's reversed complement else '\x1c'

        Tab-delimited tuple columns (postponed_sam):
        Standard 11+ -column raw SAM output

        Single column (unique):
        1. A unique read sequence

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
        report_multiplier: if verbose is True, the line number of an alignment
            or read written to stderr increases exponentially with base
            report_multiplier.
        min_exon_size: minimum exon size searched for in intron_search.py later
            in pipeline; used to determine how large a soft clip on one side of
            a read is necessary to pass it on to intron search pipeline

        No return value.
    """
    temp_dir = tempfile.mkdtemp()
    atexit.register(handle_temporary_directory, temp_dir)
    align_file = os.path.join(temp_dir, 'first_pass_reads.temp')
    other_reads_file = os.path.join(temp_dir, 'other_reads.temp')
    k_value, _, _ = bowtie.parsed_bowtie_args(bowtie2_args)
    with open(align_file, 'w') as align_stream, \
        open(other_reads_file, 'w') as other_stream:
        for seq_number, ((seq,), xpartition) in enumerate(
                                                        xstream(sys.stdin, 1)
                                                    ):
            is_reversed, name, qual = next(xpartition)
            print >>align_stream, '\t'.join([
                    '%s\x1d%s' % (is_reversed, name),
                    seq, qual
                ])
            others_printed = False
            for is_reversed, name, qual in xpartition:
                print >>other_stream, '\t'.join([
                        str(seq_number), is_reversed, name, qual
                    ])
                others_printed = True
            if not others_printed:
                print >>other_stream, str(seq_number)
    output_file = os.path.join(temp_dir, 'first_pass_out.sam')
    bowtie_command = ' '.join([bowtie2_exe,
        bowtie2_args if bowtie2_args is not None else '',
        ' --local -t --no-hd --mm -x', bowtie2_index_base, '--12',
        align_file, '-S', output_file])
    print >>sys.stderr, 'Starting first-pass Bowtie2 with command: ' \
        + bowtie_command
    # Because of problems with buffering, write output to file
    bowtie_process = subprocess.Popen(bowtie_command, bufsize=-1,
        stdout=subprocess.PIPE, stderr=sys.stderr, shell=True)
    bowtie_process.wait()
    reference_index = bowtie_index.BowtieIndexReference(bowtie_index_base)
    manifest_object = manifest.LabelsAndIndices(manifest_file)
    os.remove(align_file)
    align_file = os.path.join(temp_dir, 'second_pass_reads.temp')
    if os.path.exists(output_file):
        with open(output_file) as bowtie_stream, \
            open(other_reads_file) as other_stream, \
            open(align_file, 'w') as align_stream:
            handle_bowtie_output(
                    bowtie_stream,
                    reference_index,
                    manifest_object,
                    k_value=k_value,
                    align_stream=align_stream,
                    other_stream=other_stream,
                    exon_differentials=exon_differentials, 
                    exon_intervals=exon_intervals, 
                    bin_size=bin_size,
                    verbose=verbose, 
                    output_stream=output_stream,
                    report_multiplier=report_multiplier,
                    min_exon_size=min_exon_size
                )
        os.remove(output_file)
        output_file = os.path.join(temp_dir, 'second_pass_out.sam')
        bowtie_command = ' '.join([bowtie2_exe,
            bowtie2_args if bowtie2_args is not None else '',
            ' --local -t --no-hd --mm -x', bowtie2_index_base, '--12',
            align_file, '-S', output_file])
        print >>sys.stderr, 'Starting second-pass Bowtie2 with command: ' \
            + bowtie_command
        bowtie_process = subprocess.Popen(bowtie_command, bufsize=-1,
            stdout=subprocess.PIPE, stderr=sys.stderr, shell=True)
        bowtie_process.wait()
        if os.path.exists(output_file):
            with open(output_file) as bowtie_stream:
                handle_bowtie_output(
                        bowtie_stream,
                        reference_index,
                        manifest_object,
                        k_value=k_value,
                        exon_differentials=exon_differentials, 
                        exon_intervals=exon_intervals, 
                        bin_size=bin_size,
                        verbose=verbose, 
                        output_stream=output_stream,
                        report_multiplier=report_multiplier,
                        min_exon_size=min_exon_size
                    )

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
    parser.add_argument('--min-exon-size', type=int, required=False,
        default=8,
        help='Minimum size of exons searched for in intron_search.py')
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

    # Start keep_alive thread immediately
    if args.keep_alive:
        from dooplicity.tools import KeepAlive
        keep_alive_thread = KeepAlive(sys.stderr)
        keep_alive_thread.start()

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
        report_multiplier=args.report_multiplier,
        min_exon_size=args.min_exon_size)
    print >>sys.stderr, 'DONE with align_reads.py; in/out=%d/%d; ' \
        'time=%0.3f s' % (_input_line_count, _output_line_count,
                            time.time() - start_time)
elif __name__ == '__main__':
    # Test units
    del sys.argv[1:] # Don't choke on extra command-line parameters
    import unittest
    unittest.main()
