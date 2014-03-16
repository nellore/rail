#!/usr/bin/env python
"""
Rail-RNA-realign

Follows Rail-RNA-intron_post
Precedes Rail-RNA-collapse / Rail-RNA-bam / Rail-RNA-bed_pre

Realignment script for MapReduce pipelines that wraps Bowtie. Uses Bowtie
index including only sequences framing introns to align only those reads for 
which Bowtie did not report at least one alignment in Rail-RNA-align. Reference
names in this index encode intron sizes and locations in the (presmably) exonic
sequences it records. Infers exonic chunks and introns from alignments.

Input (read from stdin)
----------------------------
Tab-delimited input tuple columns:
  1. Name + '\x1f' + ('0' if not paired-end; '1' if first read from pair; '2'
        if second read from pair.)
  2. Nucleotide sequence
  3. Quality sequence

Hadoop output (written to stdout)
----------------------------
A given RNAME sequence is partitioned into intervals ("bins") of some 
user-specified length (see partition.py).

Exonic chunks (aka ECs; two formats, any or both of which may be emitted):

Format 1 (exon_ival); tab-delimited output tuple columns:
1. Reference name (RNAME in SAM format) + ';' + bin number
2. Sample label
3. EC start (inclusive) on forward strand
4. EC end (exclusive) on forward strand

Format 2 (exon_diff); tab-delimited output tuple columns:
1. Reference name (RNAME in SAM format) + ';' + 
    max(EC start, bin start) (inclusive) on forward strand IFF diff is
    positive and EC end (exclusive) on forward strand IFF diff is negative
2. Bin number
3. Sample label
4. +1 or -1.

Note that only unique alignments are currently output as ivals and/or diffs.

Exonic chunks / introns

Format 3 (splice_sam); tab-delimited output tuple columns:
Standard 11-column SAM output except fields are in different order, and the
first field corresponds to sample label. (Fields are reordered to facilitate
partitioning by sample name/RNAME and sorting by POS.) Each line corresponds to
read overlapping at least one intron in the reference. The CIGAR string
represents intronic bases with N's and exonic bases with M's.
The order of the fields is as follows.
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
... + optional fields, including -- for reads overlapping introns --:
XS:A:'+' or '-' depending on which strand is the sense strand

Introns

Tab-delimited output tuple columns (intron):
1. Reference name (RNAME in SAM format) 
    + ';'  + ('+' or '-' indicating which strand is the sense strand)
    + ';' + Sample label
    + ';' + Intron start (inclusive) on forward strand
    + ';' + Intron end (exclusive) on forward strand
3. Number of nucleotides between 5' end of intron and 5' end of read from which
it was inferred, ASSUMING THE SENSE STRAND IS THE FORWARD STRAND. That is, if
the sense strand is the reverse strand, this is the distance between the 3' end
of the read and the 3' end of the intron.
4. Number of nucleotides between 3' end of intron and 3' end of read from which
it was inferred, ASSUMING THE SENSE STRAND IS THE FORWARD STRAND.
5. Number of nucleotides spanned by EC on the left (that is, towards the 5'
end of the read) of the intron, ASSUMING THE SENSE STRAND IS THE FORWARD
STRAND.
6. Number of nucleotides spanned by EC on the right (that is, towards the 3'
end of the read) of the intron, ASSUMING THE SENSE STRAND IS THE FORWARD
STRAND.

ALL OUTPUT COORDINATES ARE 1-INDEXED.
"""
import sys
import os
import site
import threading
import string
import tempfile
import atexit
import subprocess
import itertools

base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for directory_name in ['bowtie', 'sample', 'interval']:
    site.addsitedir(os.path.join(base_path, directory_name))

import bowtie
import bowtie_index
import sample
import partition

# Initialize global variables for tracking number of input/output lines
_input_line_count = 0
_output_line_count = 0

def multiread_with_introns(multiread, stranded=False):
    """ Modifies read alignments to correct CIGARs and reference positions.

        An alignment takes the form of a line of SAM.
        Introns are encoded in a multiread's RNAME as follows:

        original RNAME + '+' or '-' indicating which strand is the sense
        strand + ';' + start position of sequence + ';' + comma-separated
        list of subsequence sizes framing introns + ';' + comma-separated
        list of intron sizes.

        More than one alignment output by Bowtie may correspond to the same
        alignment in reference space because at least two reference names in
        the intron Bowtie index contained overlapping bases. In this case, the
        alignments are collapsed into a single alignment. If an alignment is
        found to overlap introns, the XS:A:('+' or '-') field is appended to
        indicate which strand is the sense strand. An NH:i:(integer) field is
        also added to each alignment to indicate the number of alignments.
        These extra fields are used by Cufflinks.

        multiread: a list of lists, each of whose elements are the tokens
            from a line of SAM representing an alignment.
        stranded: if input reads are stranded, an alignment is returned only
            if the strand of any introns agrees with the strand of the 
            alignment.

        Return value: alignments modified according to the rules given above;
            a list of tuples, each of whose elements are the tokens from a line
            of SAM representing an alignment.
    """
    corrected_multiread = set()
    for i in xrange(len(multiread)):
        tokens = multiread[i][2].split(';')
        reverse_strand_string = tokens[-4][-1]
        assert reverse_strand_string == '+' or reverse_strand_string == '-'
        reverse_strand = (True if reverse_strand_string == '-' else False)
        flag = int(tokens[1])
        if stranded and (flag & 16 != 0) == reverse_strand:
            # Strand of alignment doesn't agree with strand of intron
            continue
        rname = tokens[-4][:-1]
        offset = int(multiread[i][3]) - 1
        seq_size = len(multiread[i][10])
        exon_sizes = [int(exon_size) for exon_size in tokens[-2].split(',')]
        intron_sizes = [int(intron_size) for intron_size
                            in tokens[-1].split(',')]
        assert len(intron_sizes) == len(exon_sizes) - 1
        cigar_sizes = [None]*(len(exon_sizes) + len(intron_sizes))
        cigar_sizes[::2] = exon_sizes
        cigar_sizes[1::2] = intron_sizes
        partial_sizes = [sum(exon_sizes[:j+1]) for j
                            in xrange(len(exon_sizes))]
        start_index, end_index = None, None
        for j, partial_size in enumerate(partial_sizes):
            if partial_size > offset:
                start_index = j
                new_offset = (offset - partial_sizes[j-1]
                                if j != 0 else offset)
                break
        end_offset = seq_size + offset
        for j, partial_size in enumerate(partial_sizes):
            if partial_size >= end_offset:
                end_index = j
                break
        pos = int(tokens[-3]) + new_offset + \
                sum(cigar_sizes[:start_index*2])
        if start_index is None or end_index is None:
            RuntimeError('Invalid SAM line; sum of exon sizes doesn\'t agree '
                         'with size of reference sequence.')
        if start_index == end_index:
            cigar = str(seq_size) + 'M'
            corrected_multiread.add(
                    (multiread[i][0], multiread[i][1], rname, str(pos),
                        multiread[i][4], cigar) + tuple(multiread[i][6:])
                )
        else:
            assert start_index < end_index
            cigar_sizes[start_index*2] = partial_sizes[start_index] - offset
            cigar_sizes[end_index*2] = end_offset - partial_sizes[end_index-1]
            cigar = ''.join([(str(size) + 'M') if j % 2 == 0
                                else (str(size) + 'N')
                                for j, size in enumerate(cigar_sizes)
                                if start_index*2 <= j <= end_index*2])
            corrected_multiread.add(
                    (multiread[i][0], multiread[i][1], rname, str(pos),
                        multiread[i][4], cigar) + tuple(multiread[i][6:]) 
                    + (('XS:A:' + reverse_strand_string),)
                )
    NH_field = 'NH:i:' + str(len(corrected_multiread))
    multiread_to_return = [alignment + (NH_field,) for alignment in
                            corrected_multiread]
    return multiread_to_return

class BowtieOutputThread(threading.Thread):
    """ Processes Bowtie alignments, emitting tuples for exons and introns. """
    
    def __init__(self, input_stream, reference_index, output_stream=sys.stdout,
        exon_differentials=True, exon_intervals=False, stranded=False,
        verbose=False, bin_size=10000, report_multiplier=1.2):
        """ Constructor for BowtieOutputThread.

            input_stream: where to retrieve Bowtie's SAM output, typically a
                Bowtie process's stdout.
            reference_index: object of class BowtieIndexReference; for
                outputing RNAME number strings to facilitate later sorting
            output_stream: where to emit exon and intron tuples; typically,
                this is sys.stdout.
            exon_differentials: True iff EC differentials are to be emitted.
            exon_intervals: True iff EC intervals are to be emitted.
            stranded: True iff input reads are strand-specific; this affects
                whether the algorithm imposes that the strand to which the read
                aligns agrees with the strand of the intron call.
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
        self.reference_index = reference_index
        self.output_stream = output_stream
        self.verbose = verbose
        self.bin_size = bin_size
        self.stranded = stranded
        self.report_multiplier = report_multiplier
        self.exon_differentials = exon_differentials
        self.exon_intervals = exon_intervals

    def run(self):
        """ Prints SAM, exon_ivals, exon_diffs, and introns.

            Overrides default method containing thread activity.

            No return value.
        """
        global _input_line_count, _output_line_count
        next_report_line = 0
        '''Next read must be known to tell if a read mapped to multiple
        locations, so always work with previous read.'''
        while True:
            line = self.input_stream.readline()
            if not line: return # Bowtie output nothing
            # Skip header line
            if line[0] == '@': continue
            last_tokens = line.rstrip().split('\t')
            last_qname = last_tokens[0]
            last_flag = int(last_tokens[1])
            break
        # Initialize counter
        i = 0
        # While labeled multiread, this list may end up simply a uniread
        multiread = []
        while True:
            line = self.input_stream.readline()
            if line:
                tokens = line.rstrip().split('\t')
                qname = tokens[0]
                flag = int(tokens[1])
                _input_line_count += 1
            if self.verbose and next_report_line == i:
                print >>sys.stderr, \
                    'SAM output record %d: rdname="%s", flag=%d' \
                    % (i, last_qname, last_flag)
                next_report_line = int((next_report_line + 1)
                    * self.report_multiplier + 1) - 1
            multiread.append(last_tokens)
            if not line or qname != last_qname:
                if (last_flag & 4):
                    # Write only the SAM output if the read was unmapped
                    print >>self.output_stream, 'splice_sam\t' \
                         + '\t'.join(
                                [sample.parseLab(last_tokens[0][:-2]),
                                    self.reference_index.rname_to_string[
                                            last_tokens[2]
                                        ], '%012d' % 
                                    int(last_tokens[3]),
                                    last_tokens[0][:-2],
                                    last_tokens[1]] + last_tokens[4:]
                            )
                else:
                    '''Correct positions to match original reference's, correct
                    CIGARs, and eliminate duplicates.'''
                    corrected_multiread = multiread_with_introns(
                                                multiread, self.stranded
                                            )
                    sample_label = sample.parseLab(multiread[0][0][:-2])
                    for alignment in corrected_multiread:
                        print >>self.output_stream, 'splice_sam\t' \
                            + '\t'.join(
                                (sample_label,
                                    self.reference_index.rname_to_string[
                                            alignment[2]
                                        ], '%012d' % int(alignment[3]),
                                    alignment[0][:-2],
                                    alignment[1]) + alignment[4:]
                            )
                        _output_line_count += 1
                    if len(corrected_multiread) == 1:
                        '''Output exonic chunks and introns only if the 
                        alignment was unique.'''
                        alignment = list(corrected_multiread)[0]
                        cigar = alignment[5]
                        assert 'M' == cigar[-1]
                        pos = int(alignment[3])
                        rname = alignment[2]
                        cigar_chars = [-1] + [i for i, char
                            in enumerate(cigar) if char == 'M' or char == 'N']
                        base_counts = \
                            [int(cigar[(cigar_chars[i]+1):cigar_chars[i+1]])
                                for i in xrange(len(cigar_chars)-1)]
                        if 'N' in cigar:
                            '''There's some chance the alignment was
                            entirely in a region called as exonic; output
                            introns only if they're present.'''
                            reverse_strand_string = alignment[-2][-1]
                            assert reverse_strand_string == '+' or \
                                reverse_strand_string == '-'
                            # Gather intron stats to be output
                            introns = [(pos + sum(base_counts[:i]),
                                            pos + sum(base_counts[:i+1]),
                                            base_counts[i-1],
                                            base_counts[i+1],
                                            sum(base_counts[:i:2]),
                                            sum(base_counts[i::2]))
                                            for i in xrange(1,
                                                len(base_counts), 2)]
                            for (intron_pos, intron_end_pos,
                                    left_anchor_size, right_anchor_size,
                                    left_displacement, right_displacement) \
                                in introns:
                                print >>self.output_stream, \
                                    'intron\t%s%s;%s;%012d;%012d\t' \
                                    '%d\t%d\t%d\t%d' % (rname,
                                            reverse_strand_string,
                                            sample_label,
                                            intron_pos,
                                            intron_end_pos,
                                            left_displacement,
                                            right_displacement,
                                            left_anchor_size,
                                            right_anchor_size
                                        )
                                _output_line_count += 1
                        # Output exonic chunks
                        exons = [(pos + sum(base_counts[:i]),
                                    pos + sum(base_counts[:i+1])) for i
                                    in xrange(0, len(base_counts), 2)]
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
                                    diff_rname, diff_bin \
                                        = partition_id.split(';')
                                    print >>self.output_stream, \
                                        'exon_diff\t%s;%d;%s\t%s\t1' \
                                        % (diff_rname,
                                            max(partition_start, exon_pos),
                                            sample_label,
                                            diff_bin)
                                    _output_line_count += 1
                                    assert exon_end_pos > partition_start
                                    if exon_end_pos < partition_end:
                                        '''Print decrement at interval end 
                                        iff exon ends before partition
                                        ends.'''
                                        print >>self.output_stream, \
                                            'exon_diff\t%s;%d;%s\t' \
                                            '%s\t-1' \
                                            % (diff_rname, 
                                                exon_end_pos,
                                                sample_label,
                                                diff_bin)
                                        _output_line_count += 1
                multiread = []
            if not line: break
            last_tokens, last_qname, last_flag = tokens, qname, flag
            i += 1

def handle_temporary_directory(archive, temp_dir_path):
    """ Archives or deletes temporary directory.

        archive: directory name, including path, to which temp_dir_path should
            be renamed when script is complete; or None if temp_dir_path
            should be trashed.
        temp_dir_path: path of temporary directory for storing intermediate
            alignments; archived if archive is not None.

        No return value.
    """
    if archive is not None:
        # Rename temporary directory for archiving
        os.rename(temp_dir_path, archive)
    else:
        # Kill temporary directory
        import shutil
        shutil.rmtree(temp_dir_path)

def go(input_stream=sys.stdin, output_stream=sys.stdout, bowtie_exe='bowtie',
    reference_index_base='genome', bowtie_index_base='intron',
    bowtie_args=None, temp_dir_path=tempfile.mkdtemp(), bin_size=10000,
    verbose=False, exon_differentials=True, exon_intervals=False,
    stranded=False, report_multiplier=1.2):
    """ Runs Rail-RNA-realign.

        Uses Bowtie index including only sequences framing introns to align
        only those reads for which Bowtie did not report at least one alignment
        in Rail-RNA-align. Referencenames in this index encode intron sizes and
        locations in the (presmably) exonic sequences it records. Infers exonic
        chunks and introns from alignments.

        Input (read from stdin)
        ----------------------------
        Tab-delimited input tuple columns:
          1. Name + '\x1f' + ('0' if not paired-end; '1' if first read from
            pair; '2' if second read from pair.)
          2. Nucleotide sequence
          3. Quality sequence

        Hadoop output (written to stdout)
        ----------------------------
        A given RNAME sequence is partitioned into intervals ("bins") of some 
        user-specified length (see partition.py).

        Exonic chunks (aka ECs; two formats, any or both of which may be
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

        Exonic chunks / introns

        Format 3 (splice_sam); tab-delimited output tuple columns:
        Standard 11-column SAM output except fields are in different order, and
        the first field corresponds to sample label. (Fields are reordered to
        facilitate partitioning by sample name/RNAME and sorting by POS.) Each
        line corresponds to read probably overlapping at least one intron in
        the reference. The CIGAR string represents intronic bases with N's and
        exonic bases with M's.
        The order of the fields is as follows.
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
        ... + optional fields, including -- for reads overlapping introns --:
        XS:A:'+' or '-' depending on which strand is the sense strand

        Introns

        Tab-delimited output tuple columns (intron):
        1. Reference name (RNAME in SAM format) + ';' + bin number +  
            ('+' or '-' indicating which strand is the sense strand)
        2. Sample label
        3. Intron start (inclusive) on forward strand
        4. Intron end (exclusive) on forward strand
        5. Number of nucleotides between 5' end of intron and 5' end of read
        from which it was inferred, ASSUMING THE SENSE STRAND IS THE FORWARD
        STRAND. That is, if the sense strand is the reverse strand, this is the
        distance between the 3' end of the read and the 3' end of the intron.
        6. Number of nucleotides between 3' end of intron and 3' end of read
        from which it was inferred, ASSUMING THE SENSE STRAND IS THE FORWARD
        STRAND.
        7. Number of nucleotides spanned by EC on the left (that is, towards
        the 5' end of the read) of the intron, ASSUMING THE SENSE STRAND IS THE
        FORWARD STRAND.
        8. Number of nucleotides spanned by EC on the right (that is, towards
        the 3' end of the read) of the intron, ASSUMING THE SENSE STRAND IS THE
        FORWARD STRAND.

        ALL OUTPUT COORDINATES ARE 1-INDEXED.

        input_stream: where to find input reads.
        output_stream: where to emit exonic chunks and introns.
        bowtie_exe: filename of Bowtie executable; include path if not in
            $PATH.
        reference_index_base: where to find original reference; for
            outputing RNAME number strings to facilitate later sorting
        bowtie_index_base: the basename of the Bowtie index files associated
            with only the introns.
        bowtie_args: string containing precisely extra command-line arguments
            to pass to Bowtie, e.g., "--tryhard --best"; or None.
        temp_dir_path: path of temporary directory for storing intermediate
            alignments
        bin_size: genome is partitioned in units of bin_size for later load
            balancing.
        verbose: True iff more informative messages should be written to
            stderr.
        exon_differentials: True iff EC differentials are to be emitted.
        exon_intervals: True iff EC intervals are to be emitted.
        stranded: True iff input reads are strand-specific; this affects
            whether an output partition has a terminal '+' or '-' indicating
            the sense strand. Further, if stranded is True, an alignment is
            returned only if its strand agrees with the intron's strand.
        report_multiplier: if verbose is True, the line number of an alignment,
            read, or first readlet of a read written to stderr increases
            exponentially with base report_multiplier.

        No return value.
    """
    import time
    start_time = time.time()
    reference_index = bowtie_index.BowtieIndexReference(reference_index_base)
    bowtie_process, bowtie_command, threads = bowtie.proc(
            bowtieExe=bowtie_exe, bowtieIdx=bowtie_index_base,
            readFn=None, bowtieArgs=bowtie_args, sam=True,
            stdoutPipe=True, stdinPipe=True
        )
    output_thread = BowtieOutputThread(
                        bowtie_process.stdout,
                        reference_index=reference_index,
                        exon_differentials=exon_differentials, 
                        exon_intervals=exon_intervals, 
                        bin_size=bin_size,
                        verbose=verbose, 
                        output_stream=output_stream,
                        stranded=stranded,
                        report_multiplier=report_multiplier
                    )
    threads.append(output_thread)
    output_thread.start()
    for line in input_stream:
        bowtie_process.stdin.write(line)
    bowtie_process.stdin.close()
    # Join threads to pause execution in main thread
    for thread in threads:
        if verbose: print >>sys.stderr, 'Joining thread...'
        thread.join()
    output_stream.flush()
    print >> sys.stderr, 'DONE with align.py; in/out=%d/%d; ' \
        'time=%0.3f s' % (_input_line_count, _output_line_count,
                                time.time() - start_time)

if __name__ == '__main__':
    import argparse
    # Print file's docstring if -h is invoked
    parser = argparse.ArgumentParser(description=__doc__, 
                formatter_class=argparse.RawDescriptionHelpFormatter)
    # Add command-line arguments
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
    parser.add_argument(\
        '--stranded', action='store_const', const=True, default=False,
        help='Assume input reads come from the sense strand; then partitions '
             'in output have terminal + and - indicating sense strand')
    parser.add_argument('--test', action='store_const', const=True,
        default=False,
        help='Run unit tests; DOES NOT NEED INPUT FROM STDIN, AND DOES NOT '
             'WRITE EXONS AND INTRONS TO STDOUT')
    parser.add_argument(\
        '--original-idx', metavar='INDEX', type=str, required=False,
        default='genome',
        help='Path to Bowtie index of original reference. Specify its '
             'basename.')
    parser.add_argument('--archive', metavar='PATH', type=str, 
        default=None,
        help='Save output and Bowtie command to a subdirectory (named using ' 
             'this process\'s PID) of PATH')

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
    temp_dir_path = tempfile.mkdtemp()
    archive = os.path.join(args.archive,
        str(os.getpid())) if args.archive is not None else None
    # Handle temporary directory if CTRL+C'd
    atexit.register(handle_temporary_directory, archive, temp_dir_path)
    if args.verbose:
        print >>sys.stderr, 'Creating temporary directory %s' \
            % temp_dir_path
    go(bowtie_exe=args.bowtie_exe,
        reference_index_base=args.original_idx,
        bowtie_index_base=args.bowtie_idx,
        bowtie_args=bowtie_args,
        temp_dir_path=temp_dir_path, 
        verbose=args.verbose, 
        bin_size=args.partition_length,
        exon_differentials=args.exon_differentials,
        exon_intervals=args.exon_intervals,
        stranded=args.stranded,
        report_multiplier=args.report_multiplier)
elif __name__ == '__main__':
    # Test units
    del sys.argv[1:] # Don't choke on extra command-line parameters
    import unittest

    class TestMultireadWithIntrons(unittest.TestCase):
        """ Tests composed_and_sorted_readlets(); needs no fixture. """
        def test_cigar_when_read_overlaps_multiple_introns(self):
            """ Fails if CIGAR is not corrected properly. """


    unittest.main()