#!/usr/bin/env python
"""
Rail-RNA-align_reads

Follows Rail-RNA-preprocess
Precedes Rail-RNA-align_readlets, Rail-RNA-collapse, Rail-RNA-realign, and
    Rail-RNA-bam

Alignment script for MapReduce pipelines that wraps Bowtie. Obtains exonic 
chunks from end-to-end alignments. Outputs sequences that do not align for
readletizing and, later, junction inference.

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
2. Sample index
3. EC start (inclusive) on forward strand
4. EC end (exclusive) on forward strand

Format 2 (exon_diff); tab-delimited output tuple columns:
1. Reference name (RNAME in SAM format) + ';' + bin number
2. max(EC start, bin start) (inclusive) on forward strand IFF diff is
    positive and EC end (exclusive) on forward strand IFF diff is negative
3. Sample index
4. '1' if alignment from which diff originates is "unique" according to
    --tie-margin criterion; else '0'
5. +1 or -1 * count, the number of instances of a read sequence for which to
    print exonic chunks

Note that only unique alignments are currently output as ivals and/or diffs.

Format 3 (sam); tab-delimited output tuple columns:
Standard SAM output except fields are in different order, and the first field
corresponds to sample label. (Fields are reordered to facilitate partitioning
by sample name/RNAME and sorting by POS.) Each line corresponds to a
spliced alignment. The order of the fields is as follows.
1. Sample index if outputting BAMs by sample OR sample-rname index if
    outputting BAMs by chr
2. (Number string representing RNAME; see BowtieIndexReference class in
    bowtie_index for conversion information) OR '0' if outputting BAMs by chr
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
2. Number string representing RNAME
3. Start position (Last base before insertion or first base of deletion)
4. End position (Last base before insertion or last base of deletion 
                    (exclusive))
5. Inserted sequence for insertions or deleted sequence for deletions'
6. Sample index
----Next fields are for junctions only; they are '\x1c' for indels----
7. '\x1c'
8. '\x1c'
--------------------------------------------------------------------
9. Number of instances of insertion or deletion in sample; this is
    always +1 * count before bed_pre combiner/reducer

Read whose primary alignment is not end-to-end

Tab-delimited output tuple columns (unmapped):
1. Transcriptome Bowtie 2 index group number
2. SEQ
3. 2 if SEQ is reverse-complemented, else 1
4. QNAME
5. QUAL

Tab-delimited output tuple columns (readletized)
1. Readlet sequence or its reversed complement, whichever is first in
    alphabetical order
2. read sequence ID + ('-' if readlet sequence is
    reverse-complemented; else '+') + '\x1e' + displacement of readlet's 5' end
    from read's 5' end + '\x1e' + displacement of readlet's 3' end from read's
    3' end (+, for EXACTLY one readlet of a read sequence, '\x1e' +
    read sequence + '\x1e' + (an '\x1f'-separated list A of unique sample
    labels with read sequences that match the original read sequence) + '\x1e'
    + (an '\x1f'-separated list  of unique sample labels B with read sequences
    that match the reversed complement of the original read sequence))] + 
    '\x1e' + (an '\x1f'-separated list of the number of instances of the read
    sequence for each respective sample in list A) + '\x1e' + (an
    '\x1f'-separated list of the number of instances of the read
    sequence's reversed complement for each respective sample in list B). Here,
    a read sequence ID takes the form X:Y, where X is the
    "mapred_task_partition" environment variable -- a unique index for a task
    within a job -- and Y is the index of the read sequence relative to the
    beginning of the input stream.

Tab-delimited tuple columns (postponed_sam):
Standard 11+ -column raw SAM output

Single column (unique):
1. A unique read sequence

Two columns, exactly one line (dummy); ensures creation of junction index:
1. character "-"
2. the word "dummy"

ALL OUTPUT COORDINATES ARE 1-INDEXED.
"""
import sys
import os
import site
import subprocess
import time
import string

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
import tempdel
import group_reads
from dooplicity.tools import xstream, dlist, register_cleanup, xopen, \
    make_temp_dir
from alignment_handlers import AlignmentPrinter

# Initialize global variables for tracking number of input lines
_input_line_count = 0

_reversed_complement_translation_table = string.maketrans('ATCG', 'TAGC')

def go(input_stream=sys.stdin, output_stream=sys.stdout, bowtie2_exe='bowtie2',
    bowtie_index_base='genome', bowtie2_index_base='genome2', 
    manifest_file='manifest', bowtie2_args=None, bin_size=10000, verbose=False,
    exon_differentials=True, exon_intervals=False, report_multiplier=1.2,
    min_exon_size=8, search_filter=1, min_readlet_size=15, max_readlet_size=25,
    readlet_interval=12, capping_multiplier=1.5, drop_deletions=False,
    gzip_level=3, scratch=None, index_count=1, output_bam_by_chr=False,
    tie_margin=0, no_realign=False, no_polyA=False):
    """ Runs Rail-RNA-align_reads.

        A single pass of Bowtie is run to find end-to-end alignments. Unmapped
        reads are saved for readletizing to determine junctions in sucessive
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
        2. Sample index
        3. EC start (inclusive) on forward strand
        4. EC end (exclusive) on forward strand

        Format 2 (exon_diff); tab-delimited output tuple columns:
        1. Reference name (RNAME in SAM format) + ';' + bin number
        2. max(EC start, bin start) (inclusive) on forward strand IFF diff is
            positive and EC end (exclusive) on forward strand IFF diff is
            negative
        3. Sample index
        4. '1' if alignment from which diff originates is "unique" according to
            --tie-margin criterion; else '0'
        5. +1 or -1 * count, the number of instances of a read sequence for
            which to print exonic chunks

        Note that only unique alignments are currently output as ivals and/or
        diffs.

        Format 3 (sam); tab-delimited output tuple columns:
        Standard SAM output except fields are in different order, and the first
        field corresponds to sample label. (Fields are reordered to facilitate
        partitioning by sample name/RNAME and sorting by POS.) Each line
        corresponds to a spliced alignment. The order of the fields is as
        follows.
        1. Sample index if outputting BAMs by sample OR
                sample-rname index if outputting BAMs by chr
        2. (Number string representing RNAME; see BowtieIndexReference
            class in bowtie_index for conversion information) OR
            '0' if outputting BAMs by chr
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
        2. Number string representing RNAME
        3. Start position (Last base before insertion or 
            first base of deletion)
        4. End position (Last base before insertion or last base of deletion 
                            (exclusive))
        5. Inserted sequence for insertions or deleted sequence for deletions
        6. Sample index
        ----Next fields are for junctions only; they are '\x1c' for indels----
        7. '\x1c'
        8. '\x1c'
        --------------------------------------------------------------------
        9. Number of instances of insertion or deletion in sample; this is
            always +1 * count before bed_pre combiner/reducer

        Read whose primary alignment is not end-to-end

        Tab-delimited output tuple columns (unmapped):
        1. Transcriptome Bowtie 2 index group number
        2. SEQ
        3. 1 if SEQ is reverse-complemented, else 0
        4. QNAME
        5. QUAL

        Tab-delimited output tuple columns (readletized):
        1. Readlet sequence or its reversed complement, whichever is first in
            alphabetical order
        2. read sequence ID + ('-' if readlet
            sequence is reverse-complemented; else '+') + '\x1e' + displacement
            of readlet's 5' end from read's 5' end + '\x1e' + displacement of
            readlet's 3' end from read's 3' end (+, for EXACTLY one readlet of
            a read sequence, '\x1e' + read sequence + '\x1e' +
            (an '\x1f'-separated list A of unique sample labels with read
            sequences that match the original read sequence) + '\x1e' +
            (an '\x1f'-separated list  of unique sample labels B with read
            sequences that match the reversed complement of the original read
            sequence)) + '\x1e' + (an '\x1f'-separated list of the number of
            instances of the read sequence for each respective sample in list
            A) + '\x1e' + (an '\x1f'-separated list of the number of instances
            of the read sequence's reversed complement for each respective
            sample in list B). Here, a read sequence ID takes the form X:Y,
            where X is the "mapred_task_partition" environment variable -- a
            unique index for a task within a job -- and Y is the index of the
            read sequence relative to the beginning of the input stream.

        Tab-delimited tuple columns (postponed_sam):
        Standard 11+ -column raw SAM output

        Single column (unique):
        1. A unique read sequence

        Two columns, exactly one line (dummy); ensures creation of junction
            index:
        1. character "-"
        2. the word "dummy"

        ALL OUTPUT COORDINATES ARE 1-INDEXED.

        input_stream: where to find input reads.
        output_stream: where to emit exonic chunks and junctions.
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
        min_exon_size: minimum exon size searched for in junction_search.py
            later in pipeline; used to determine how large a soft clip on one
            side of a read is necessary to pass it on to junction search
            pipeline
        search_filter: how large a soft clip on one side of a read is necessary
            to pass it on to junction search pipeline
        min_readlet_size: "capping" readlets (that is, readlets that terminate
            at a given end of the read) are never smaller than this value
        max_readlet_size: size of every noncapping readlet
        readlet_interval: number of bases separating successive readlets along
            the read
        capping_multiplier: successive capping readlets on a given end of a
            read are increased in size exponentially with base
            capping_multiplier
        drop_deletions: True iff deletions should be dropped from coverage
            vector
        gzip_level: compression level to use for temporary files
        scratch: scratch directory for storing temporary files or None if 
            securely created temporary directory
        index_count: number of transcriptome Bowtie 2 indexes to which to
            assign unmapped reads for later realignment
        output_bam_by_chr: True iff final output BAMs will be by chromosome
        tie_margin: allowed score difference per 100 bases among ties in
            max score. For example, 150 and 144 are tied alignment scores
            for a 100-bp read when --tie-margin is 6.
        no_realign: True iff job flow does not need more than readlets: this
            usually means only a transcript index is being constructed
        no_polyA: kill noncapping readlets that are all As and write as
            unmapped all reads with polyA prefixes whose suffixes are <
            min_exon_size

        No return value.
    """
    global _input_line_count
    reference_index = bowtie_index.BowtieIndexReference(bowtie_index_base)
    manifest_object = manifest.LabelsAndIndices(manifest_file)
    alignment_printer = AlignmentPrinter(
            manifest_object,
            reference_index,
            bin_size=bin_size,
            output_stream=output_stream,
            exon_ivals=exon_intervals,
            exon_diffs=exon_differentials,
            drop_deletions=drop_deletions,
            output_bam_by_chr=output_bam_by_chr,
            tie_margin=tie_margin
        )
    # Get task partition to pass to align_reads_delegate.py
    try:
        task_partition = os.environ['mapred_task_partition']
    except KeyError:
        # Hadoop 2.x?
        try:
            task_partition = os.environ['mapreduce_task_partition']
        except KeyError:
            # A unit test is probably being run
            task_partition = '0'
    temp_dir = make_temp_dir(scratch)
    register_cleanup(tempdel.remove_temporary_directories, [temp_dir])
    align_file = os.path.join(temp_dir, 'first_pass_reads.temp.gz')
    other_reads_file = os.path.join(temp_dir, 'other_reads.temp.gz')
    second_pass_file = os.path.join(temp_dir, 'second_pass_reads.temp.gz')
    k_value, _, _ = bowtie.parsed_bowtie_args(bowtie2_args)
    nothing_doing = True
    # Required length of prefix after poly(A) is trimmed
    remaining_seq_size = max(min_exon_size - 1, 1)
    with xopen(True, align_file, 'w', gzip_level) as align_stream, \
        xopen(True, other_reads_file, 'w', gzip_level) as other_stream:
        for seq_number, ((seq,), xpartition) in enumerate(
                                                        xstream(sys.stdin, 1)
                                                    ):
            seq_length = len(seq)
            if no_polyA and (
                    all(seq[i] == 'A' 
                         for i in xrange(seq_length - remaining_seq_size))
                    or all(seq[i] == 'T' 
                         for i in xrange(remaining_seq_size, seq_length))
                    or all(seq[i] == 'A' 
                         for i in xrange(remaining_seq_size, seq_length))
                    or all(seq[i] == 'T' 
                         for i in xrange(seq_length - remaining_seq_size))
                ):
                if not no_realign:
                    '''If a sequence is too short without its poly(A) tail,
                    make all reads with that sequence unmapped. Technically,
                    this also kills poly(A)s at 5' ends, but we probably
                    couldn't align those sequences anyway.'''
                    reversed_complement_seq = seq[::-1].translate(
                                        _reversed_complement_translation_table
                                    )
                    for is_reversed, name, qual in xpartition:
                        if is_reversed == '0':
                            alignment_printer.print_unmapped_read(
                                                    name,
                                                    seq,
                                                    qual
                                                )
                        else:
                            alignment_printer.print_unmapped_read(
                                                    name,
                                                    reversed_complement_seq,
                                                    qual[::-1]
                                                )
                continue
            nothing_doing = False
            '''Select highest-quality read with alphabetically last qname
            for first-pass alignment.'''
            best_name, best_mean_qual, best_qual_index, i = None, None, 0, 0
            others_to_print = dlist()
            for is_reversed, name, qual in xpartition:
                _input_line_count += 1
                others_to_print.append(
                        '\t'.join([
                            str(seq_number), is_reversed, name, qual
                        ])
                    )
                mean_qual = (
                        float(sum([ord(score) for score in qual])) / len(qual)
                    )
                if (mean_qual > best_mean_qual
                        or mean_qual == best_mean_qual and name > best_name):
                    best_qual_index = i
                    best_mean_qual = mean_qual
                    best_name = name
                    to_align = '\t'.join([
                                        '%s\x1d%s' % (is_reversed, name),
                                        seq, qual
                                    ])
                i += 1
            assert i >= 1
            if i == 1:
                print >>other_stream, str(seq_number)
            else:
                for j, other_to_print in enumerate(others_to_print):
                    if j != best_qual_index:
                        print >>other_stream, other_to_print
            print >>align_stream, to_align
    # Print dummy line
    print 'dummy\t-\tdummy'
    sys.stdout.flush() # this is REALLY important b/c called script will stdout
    if nothing_doing:
        # No input
        sys.exit(0)
    input_command = 'gzip -cd %s' % align_file
    bowtie_command = ' '.join([bowtie2_exe,
        bowtie2_args if bowtie2_args is not None else '',
        ' --sam-no-qname-trunc --local -t --no-hd --mm -x',
        bowtie2_index_base, '--12 -'])
    delegate_command = ''.join(
                [sys.executable, ' ', os.path.realpath(__file__)[:-3],
                    ('_delegate.py --task-partition {task_partition} '
                     '--other-reads {other_reads} --second-pass-reads '
                     '{second_pass_reads} --min-readlet-size '
                     '{min_readlet_size} {drop_deletions} '
                     '--max-readlet-size {max_readlet_size} '
                     '--readlet-interval {readlet_interval} '
                     '--capping-multiplier {capping_multiplier:1.12f} '
                     '{verbose} --report-multiplier {report_multiplier:1.12f} '
                     '--k-value {k_value} '
                     '--bowtie-idx {bowtie_index_base} '
                     '--partition-length {bin_size} '
                     '--manifest {manifest_file} '
                     '{exon_differentials} {exon_intervals} '
                     '--gzip-level {gzip_level} '
                     '--search-filter {search_filter} '
                     '--index-count {index_count} '
                     '--tie-margin {tie_margin} '
                     '{no_realign} '
                     '{no_polyA} '
                     '{output_bam_by_chr}').format(
                        task_partition=task_partition,
                        other_reads=other_reads_file,
                        second_pass_reads=second_pass_file,
                        min_readlet_size=min_readlet_size,
                        drop_deletions=('--drop-deletions' if drop_deletions
                                            else ''),
                        max_readlet_size=max_readlet_size,
                        readlet_interval=readlet_interval,
                        capping_multiplier=capping_multiplier,
                        verbose=('--verbose' if verbose else ''),
                        report_multiplier=report_multiplier,
                        k_value=k_value,
                        bowtie_index_base=bowtie_index_base,
                        bin_size=bin_size,
                        manifest_file=manifest_file,
                        exon_differentials=('--exon-differentials'
                                            if exon_differentials else ''),
                        exon_intervals=('--exon-intervals'
                                        if exon_intervals else ''),
                        gzip_level=gzip_level,
                        search_filter=search_filter,
                        index_count=index_count,
                        tie_margin=tie_margin,
                        no_realign=('--no-realign' if no_realign else ''),
                        no_polyA=('--no-polyA' if no_polyA else ''),
                        output_bam_by_chr=('--output-bam-by-chr'
                                            if output_bam_by_chr
                                            else '')
                     )]
            )
    full_command = ' | '.join([input_command, 
                                bowtie_command, delegate_command])
    print >>sys.stderr, \
        'Starting first-pass Bowtie 2 with command: ' + full_command
    bowtie_process = subprocess.Popen(' '.join(
                    ['set -exo pipefail;', full_command]
                ),
            bufsize=-1, stdout=sys.stdout, stderr=sys.stderr, shell=True,
            executable='/bin/bash')
    return_code = bowtie_process.wait()
    if return_code:
        raise RuntimeError('Error occurred while reading first-pass Bowtie 2 '
                           'output; exitlevel was %d.' % return_code)
    os.remove(align_file)
    os.remove(other_reads_file)
    if not no_realign:
        input_command = 'gzip -cd %s' % second_pass_file
        bowtie_command = ' '.join([bowtie2_exe,
            bowtie2_args if bowtie2_args is not None else '',
            ' --sam-no-qname-trunc --local -t --no-hd --mm -x',
            bowtie2_index_base, '--12 -'])
        delegate_command = ''.join(
                    [sys.executable, ' ', os.path.realpath(__file__)[:-3],
                        ('_delegate.py --task-partition {task_partition} '
                         '--min-readlet-size {min_readlet_size} '
                         '{drop_deletions} '
                         '--max-readlet-size {max_readlet_size} '
                         '--readlet-interval {readlet_interval} '
                         '--capping-multiplier {capping_multiplier:012f} '
                         '{verbose} '
                         '--report-multiplier {report_multiplier:012f} '
                         '--k-value {k_value} '
                         '--bowtie-idx {bowtie_index_base} '
                         '--partition-length {bin_size} '
                         '--manifest {manifest_file} '
                         '{exon_differentials} {exon_intervals} '
                         '--gzip-level {gzip_level} '
                         '--search-filter {search_filter} ' 
                         '--index-count {index_count} '
                         '--tie-margin {tie_margin} '
                         '{output_bam_by_chr}').format(
                            task_partition=task_partition,
                            min_readlet_size=min_readlet_size,
                            drop_deletions=('--drop-deletions'
                                                if drop_deletions else ''),
                            readlet_interval=readlet_interval,
                            max_readlet_size=max_readlet_size,
                            capping_multiplier=capping_multiplier,
                            verbose=('--verbose' if verbose else ''),
                            report_multiplier=report_multiplier,
                            k_value=k_value,
                            bowtie_index_base=bowtie_index_base,
                            bin_size=bin_size,
                            manifest_file=manifest_file,
                            exon_differentials=('--exon-differentials'
                                                if exon_differentials else ''),
                            exon_intervals=('--exon-intervals'
                                            if exon_intervals else ''),
                            gzip_level=gzip_level,
                            search_filter=search_filter,
                            index_count=index_count,
                            tie_margin=tie_margin,
                            output_bam_by_chr=('--output-bam-by-chr'
                                                if output_bam_by_chr
                                                else '')
                         )]
                )
        full_command = ' | '.join([input_command, 
                                    bowtie_command, delegate_command])
        print >>sys.stderr, \
            'Starting second-pass Bowtie 2 with command: ' + full_command
        bowtie_process = subprocess.Popen(' '.join(
                        ['set -exo pipefail;', full_command]
                    ),
                bufsize=-1, stdout=sys.stdout, stderr=sys.stderr, shell=True,
                executable='/bin/bash')
        return_code = bowtie_process.wait()
        if return_code:
            raise RuntimeError('Error occurred while reading second-pass '
                               'Bowtie 2 output; exitlevel was %d.'
                                % return_code)
    sys.stdout.flush()

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
    parser.add_argument('--drop-deletions', action='store_const',
        const=True,
        default=False, 
        help='Drop deletions from coverage vectors')
    parser.add_argument('--output-bam-by-chr', action='store_const',
        const=True,
        default=False, 
        help='Final BAMs will be output by chromosome')
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
        help='Minimum size of exons searched for in junction_search.py')
    parser.add_argument('--search-filter', type=int, required=False,
        default=1,
        help=('Minimum size of soft-clipped end that dispatches a read for '
              'junction search'))
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
    parser.add_argument('--gzip-level', type=int, required=False,
        default=3,
        help='Level of gzip compression to use for temporary files')
    parser.add_argument('--no-realign', action='store_const',
        const=True,
        default=False, 
        help=('Suppresses some output and does not perform second-pass '
              'alignment if unnecessary for job flow'))
    parser.add_argument('--no-polyA', action='store_const',
        const=True,
        default=False, 
        help='Disallows any capping readlet that is a string of A nucleotides')
    
    # Add command-line arguments for dependencies
    partition.add_args(parser)
    bowtie.add_args(parser)
    manifest.add_args(parser)
    tempdel.add_args(parser)
    group_reads.add_args(parser)
    from alignment_handlers import add_args as alignment_handlers_add_args
    alignment_handlers_add_args(parser)

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

    start_time = time.time()
    go(bowtie2_exe=os.path.expandvars(args.bowtie2_exe),
        bowtie_index_base=os.path.expandvars(args.bowtie_idx),
        bowtie2_index_base=os.path.expandvars(args.bowtie2_idx),
        bowtie2_args=bowtie_args,
        manifest_file=os.path.expandvars(args.manifest),
        verbose=args.verbose, 
        bin_size=args.partition_length,
        exon_differentials=args.exon_differentials,
        exon_intervals=args.exon_intervals,
        report_multiplier=args.report_multiplier,
        min_exon_size=args.min_exon_size,
        search_filter=args.search_filter,
        min_readlet_size=args.min_readlet_size,
        max_readlet_size=args.max_readlet_size,
        readlet_interval=args.readlet_interval,
        capping_multiplier=args.capping_multiplier,
        drop_deletions=args.drop_deletions,
        gzip_level=args.gzip_level,
        scratch=tempdel.silentexpandvars(args.scratch),
        index_count=args.index_count,
        output_bam_by_chr=args.output_bam_by_chr,
        tie_margin=args.tie_margin,
        no_realign=args.no_realign,
        no_polyA=args.no_polyA)

    print >>sys.stderr, 'DONE with align_reads.py; in=%d; ' \
        'time=%0.3f s' % (_input_line_count, time.time() - start_time)