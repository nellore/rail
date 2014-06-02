#!/usr/bin/env python
"""
Rail-RNA-realign_reads

Follows Rail-RNA-cointron_search
Precedes Rail-RNA-collapse / Rail-RNA-bam / Rail-RNA-bed_pre

Realignment script for MapReduce pipelines that wraps Bowtie2. Creates Bowtie2
index including only sequences framing introns to align only those reads for
which Bowtie2 did not report alignments in Rail-RNA-align. Reference names in
the index encode intron sizes and locations in the (presmably) exonic
sequences it records. Infers exonic chunks and introns from alignments.

Input (read from stdin)
----------------------------
Tab-delimited input tuple columns:
(two kinds)

Type 1:
1. Read sequence
2. '\x1c' + FASTA reference name including '>'. The following format is used:
    original RNAME + '+' or '-' indicating which strand is the sense strand
    + ';' + start position of sequence + ';' + comma-separated list of
    subsequence sizes framing introns + ';' + comma-separated list of intron
    sizes)
3. FASTA sequence

Type 2:
1. Read sequence
2. QNAME
3. Quality sequence

Type 1 corresponds to a FASTA line to index to which the read sequence is
predicted to align. Type 2 corresponds to a distinct read. Input is partitioned
by field 1.

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

Introns/insertions/deletions

Format 4 (bed); tab-delimited output tuple columns (bed):
1. 'I', 'D', or 'N' for insertion, deletion, or intron line
2. Sample label
3. Number string representing RNAME
4. Start position (Last base before insertion, first base of deletion,
                    or first base of intron)
5. End position (Last base before insertion, last base of deletion (exclusive),
                    or last base of intron (exclusive))
6. '+' or '-' indicating which strand is the sense strand for introns,
   inserted sequence for insertions, or deleted sequence for deletions
----Next fields are for introns only; they are '\x1c' for indels----
7. Number of nucleotides between 5' end of intron and 5' end of read from which
    it was inferred, ASSUMING THE SENSE STRAND IS THE FORWARD STRAND. That is,
    if the sense strand is the reverse strand, this is the distance between the
    3' end of the read and the 3' end of the intron.
8. Number of nucleotides between 3' end of intron and 3' end of read from which
    it was inferred, ASSUMING THE SENSE STRAND IS THE FORWARD STRAND.
--------------------------------------------------------------------
9. Number of instances of intron, insertion, or deletion in sample; this is
    always +1 before bed_pre combiner/reducer

ALL OUTPUT COORDINATES ARE 1-INDEXED.
"""
import sys
import os
# So intron index is accessible
sys.path.append(os.path.dirname(__file__))
import site
import threading
import string
import tempfile
import atexit
import subprocess
import itertools
import re
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
import bowtie
import bowtie_index
import partition
from cigar_parse import indels_introns_and_exons

# Initialize global variables for tracking number of input/output lines
_input_line_count = 0
_output_line_count = 0

def input_files_from_input_stream(input_stream,
                                    temp_dir_path=tempfile.mkdtemp(),
                                    verbose=False):
    """ Creates FASTA reference to index and separate file with reads.

        input_stream: where to find Hadoop input
        temp_dir_path: where to store files
        verbose: output extra debugging messages

        Return value: tuple (path to FASTA reference file, path to read file)
    """
    global _input_line_count
    prefasta_filename = os.path.join(temp_dir_path, 'temp.prefa')
    deduped_fasta_filename = os.path.join(temp_dir_path, 'temp.deduped.prefa')
    final_fasta_filename = os.path.join(temp_dir_path, 'temp.fa')
    reads_filename = os.path.join(temp_dir_path, 'reads.temp')
    if verbose:
        print >>sys.stderr, 'Writing prefasta and input reads...'
    with open(prefasta_filename, 'w') as fasta_stream:
        with open(reads_filename, 'w') as read_stream:
            for read_seq, xpartition in xstream(input_stream, 1):
                for value in xpartition:
                    _input_line_count += 1
                    if value[0][0] == '\x1c':
                        # Add to FASTA reference
                        print >>fasta_stream, '\t'.join([value[0][1:],
                                                         value[1]])
                    else:
                        # Add to reads
                        print >>read_stream, '\t'.join([value[0],
                                                        read_seq[0],
                                                        value[1]])
    if verbose:
        print >>sys.stderr, 'Done! Sorting and deduplicating prefasta...'
    # Sort prefasta and eliminate duplicate lines
    dedup_process_return = subprocess.call(
            r'''sort %s | uniq > %s'''
            % (prefasta_filename, deduped_fasta_filename), shell=True
        )
    if dedup_process_return != 0:
        raise RuntimeError('Problem encountered deduplicating FASTA reference')
    if verbose:
        print >>sys.stderr, 'Done! Writing final FASTA.'
    with open(final_fasta_filename, 'w') as output_stream:
        with open(deduped_fasta_filename) as fasta_stream:
            for line in fasta_stream:
                rname, seq = line.strip().split('\t')
                print >>output_stream, rname
                output_stream.write(
                    '\n'.join([seq[i:i+80] for i 
                                in xrange(0, len(seq), 80)])
                )
                output_stream.write('\n')
    return final_fasta_filename, reads_filename

def create_index_from_reference_fasta(bowtie2_build_exe, fasta_file,
        index_basename, keep_alive=False):
    """ Creates Bowtie2 index from reference fasta.

        bowtie2_build_exe: Path to Bowtie2
        fasta_file: Path to reference FASTA to index
        index_basename: Path to index basename to be created

        No return value.
    """
    if args.keep_alive:
        class BowtieBuildThread(threading.Thread):
            """ Wrapper class for bowtie-build to poll for completion.
            """
            def __init__(self, command_list):
                super(BowtieBuildThread, self).__init__()
                self.command_list = command_list
                self.bowtie_build_process = None
            def run(self):
                self.bowtie_build_process = subprocess.Popen(self.command_list,
                                                stdout=sys.stderr).wait()
        bowtie_build_thread = BowtieBuildThread([args.bowtie2_build_exe,
                                                    fasta_file,
                                                    index_basename])
        bowtie_build_thread.start()
        while bowtie_build_thread.is_alive():
            print >>sys.stderr, 'reporter:status:alive'
            sys.stderr.flush()
            time.sleep(5)
        if bowtie_build_thread.bowtie_build_process:
            raise RuntimeError('Bowtie index construction failed '
                               'w/ exitlevel %d.'
                                % bowtie_build_thread.bowtie_build_process)
    else:
        bowtie_build_process = subprocess.Popen(
                                    [args.bowtie2_build_exe,
                                        fasta_file,
                                        index_basename],
                                    stderr=sys.stderr,
                                    stdout=sys.stderr
                                )
        bowtie_build_process.wait()
        if bowtie_build_process.returncode:
            raise RuntimeError('Bowtie index construction failed w/ '
                               'exitlevel %d.'
                               % bowtie_build_process.returncode)

def running_sum(iterable):
    """ Generates a running sum of the numbers in an iterable

        iterable: some iterable with numbers

        Yield value: next value in running sum
    """
    total = 0
    for number in iterable:
        total += number
        yield total

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
    new_multiread = set()
    for alignment in multiread:
        tokens = alignment[2].split(';')
        offset = int(alignment[3]) - 1
        cigar = re.split(r'([MINDS])', alignment[5])[:-1]
        flag = int(alignment[1])
        reverse_strand_string = tokens[0][-1]
        assert reverse_strand_string in '+-'
        reverse_strand = (True if reverse_strand_string == '-' else False)
        if stranded and (flag & 16 != 0) == reverse_strand:
            # Strand of alignment doesn't agree with strand of intron
            continue
        rname = tokens[0][:-1]
        exon_sizes = map(int, tokens[-2].split(','))
        intron_sizes = map(int, tokens[-1].split(','))
        for i, exon_sum in enumerate(running_sum(exon_sizes)):
            if exon_sum > offset: break
        # Compute start position of alignment
        pos = offset + sum(intron_sizes[:i]) + int(tokens[1])
        # Adjust exon/intron lists so they start where alignment starts
        exon_sizes = exon_sizes[i:]
        exon_sizes[0] = exon_sum - offset
        intron_sizes = intron_sizes[i:]
        new_cigar = []
        for i in xrange(0, len(cigar), 2):
            char_type = cigar[i+1]
            base_count = int(cigar[i])
            if char_type in 'MD':
                for j, exon_sum in enumerate(running_sum(exon_sizes)):
                    if exon_sum >= base_count: break
                for k in xrange(j):
                    new_cigar.extend(
                            [(str(exon_sizes[k]) + char_type)
                                if exon_sizes[k] != 0 else '',
                             str(intron_sizes[k]), 'N']
                        )
                last_size = base_count - (exon_sum - exon_sizes[j])
                new_cigar.extend([str(last_size), char_type])
                new_size = exon_sum - base_count
                exon_sizes = exon_sizes[j:]
                exon_sizes[0] = new_size
                intron_sizes = intron_sizes[j:]
            elif char_type in 'IS':
                new_cigar.extend([cigar[i], char_type])
            else:
                raise RuntimeError('Bowtie2 CIGAR chars are expected to be '
                                   'in set (DIMS).')
        new_multiread.add(
                    (alignment[0], alignment[1], rname, str(pos),
                        alignment[4], ''.join(new_cigar))
                    + tuple(alignment[6:])
                    + ('XS:A:' + reverse_strand_string,)
                )
    NH_field = 'NH:i:' + str(len(new_multiread))
    multiread_to_return = [alignment + (NH_field,) for alignment in
                            new_multiread]
    return multiread_to_return

class BowtieOutputThread(threading.Thread):
    """ Processes Bowtie alignments, emitting tuples for exons and introns. """
    
    def __init__(self, input_stream, reference_index, manifest_object,
        output_stream=sys.stdout, exon_differentials=True, 
        exon_intervals=False, stranded=False, verbose=False, bin_size=10000,
        report_multiplier=1.2):
        """ Constructor for BowtieOutputThread.

            input_stream: where to retrieve Bowtie's SAM output, typically a
                Bowtie process's stdout.
            reference_index: object of class BowtieIndexReference; for
                outputing RNAME number strings to facilitate later sorting
            manifest_object: object of class LabelsAndIndices that maps indices
                to labels and back; used to shorten intermediate output.
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
        self.manifest_object = manifest_object

    def run(self):
        """ Prints SAM, exon_ivals, exon_diffs, and introns.

            Overrides default method containing thread activity.

            No return value.
        """
        global _output_line_count
        next_report_line = 0
        i = 0
        for (qname,), xpartition in xstream(self.input_stream, 1):
            # While labeled multiread, this list may end up simply a uniread
            multiread = []
            for rest_of_line in xpartition:
                i += 1
                flag = int(rest_of_line[0])
                if self.verbose and next_report_line == i:
                    print >>sys.stderr, \
                        'SAM output record %d: rdname="%s", flag=%d' \
                        % (i, qname, flag)
                    next_report_line = max(int(next_report_line
                        * self.report_multiplier), next_report_line + 1)
                multiread.append((qname,) + rest_of_line)
            if flag & 4:
                # Write only the SAM output if the read was unmapped
                print >>self.output_stream, 'splice_sam\t' \
                     + '\t'.join(
                            list((self.manifest_object.label_to_index[
                                    qname.rpartition('\x1d')[2]
                                ],
                                self.reference_index.rname_to_string[
                                        rest_of_line[1]
                                    ], '%012d' % 
                                int(rest_of_line[2]),
                                qname.partition('\x1d')[0],
                                rest_of_line[0]) + rest_of_line[3:])
                        )
            else:
                '''Correct positions to match original reference's, correct
                CIGARs, and eliminate duplicates.'''
                corrected_multiread = multiread_with_introns(
                                            multiread, self.stranded
                                        )
                sample_label = self.manifest_object.label_to_index[
                                    multiread[0][0].rpartition('\x1d')[2]
                                    ]
                for alignment in corrected_multiread:
                    print >>self.output_stream, 'splice_sam\t' \
                        + '\t'.join(
                            (sample_label,
                                self.reference_index.rname_to_string[
                                        alignment[2]
                                    ], '%012d' % int(alignment[3]),
                                alignment[0].partition('\x1d')[0],
                                alignment[1]) + alignment[4:]
                        )
                    _output_line_count += 1
                if len(corrected_multiread) == 1:
                    '''Output exonic chunks and introns only if there is
                    exactly one alignment.'''
                    alignment = list(corrected_multiread)[0]
                    cigar = alignment[5]
                    rname = alignment[2]
                    pos = int(alignment[3])
                    seq = alignment[9]
                    md = [field for field in alignment
                                if field[:5] == 'MD:Z:'][0][5:]
                    insertions, deletions, introns, exons \
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
                    try:
                        reverse_strand_string = [field for field in alignment
                                        if field[:5] == 'XS:A:'][0][5:]
                    except IndexError:
                        # No introns
                        continue
                    # Output introns
                    for (intron_pos, intron_end_pos,
                            left_displacement, right_displacement) \
                        in introns:
                        print >>self.output_stream, (
                                ('bed\tN\t%s\t%s\t%012d\t%012d\t%s\t'
                                 '%d\t%d\t1')
                                 % (sample_label, 
                                    self.reference_index.\
                                    rname_to_string[rname],
                                    intron_pos, intron_end_pos,
                                    reverse_strand_string,
                                    left_displacement,
                                    right_displacement)
                            )
                        _output_line_count += 1

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

def go(input_stream=sys.stdin, output_stream=sys.stdout, bowtie2_exe='bowtie2',
    bowtie2_build_exe='bowtie2-build', reference_index_base='genome',
    bowtie2_args=None, temp_dir_path=tempfile.mkdtemp(), bin_size=10000,
    manifest_file='manifest', verbose=False, exon_differentials=True,
    exon_intervals=False, stranded=False, report_multiplier=1.2,
    keep_alive=False):
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

        Introns/insertions/deletions

        Format 4 (bed); tab-delimited output tuple columns (bed):
        1. 'I', 'D', or 'N' for insertion, deletion, or intron line
        2. Sample label
        3. Number string representing RNAME
        4. Start position (Last base before insertion, first base of deletion,
                            or first base of intron)
        5. End position (Last base before insertion, last base of deletion
                            (exclusive), or last base of intron (exclusive))
        6. '+' or '-' indicating which strand is the sense strand for introns,
           inserted sequence for insertions, or deleted sequence for deletions
        ----Next fields are for introns only; they are '\x1c' for indels----
        7. Number of nucleotides between 5' end of intron and 5' end of read
            from which it was inferred, ASSUMING THE SENSE STRAND IS THE
            FORWARD STRAND. That is, if the sense strand is the reverse strand,
            this is the distance between the 3' end of the read and the 3' end
            of the intron.
        8. Number of nucleotides between 3' end of intron and 3' end of read
            from which it was inferred, ASSUMING THE SENSE STRAND IS THE
            FORWARD STRAND.
        ---------------------------------------------------------------------
        9. Number of instances of intron, insertion, or deletion in sample;
            this is always +1 before bed_pre combiner/reducer

        ALL OUTPUT COORDINATES ARE 1-INDEXED.

        input_stream: where to find input reads.
        output_stream: where to emit exonic chunks and introns.
        bowtie2_exe: filename of Bowtie executable; include path if not in
            $PATH.
        bowtie2_build_exe: path to bowtie2-build executable
        reference_index_base: where to find original reference; for
            outputing RNAME number strings to facilitate later sorting
        bowtie2_args: string containing precisely extra command-line arguments
            to pass to Bowtie2
        manifest_file: path to manifest file
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
        keep_alive: prints reporter:status:alive messages to stderr to keep EMR 
            task alive.

        No return value.
    """
    start_time = time.time()
    fasta_file, reads_file = input_files_from_input_stream(input_stream,
                                        verbose=verbose,
                                        temp_dir_path=temp_dir_path
                                    )
    bowtie2_index_base = os.path.join(temp_dir_path, 'tempidx')
    create_index_from_reference_fasta(bowtie2_build_exe, fasta_file,
        bowtie2_index_base)
    reference_index = bowtie_index.BowtieIndexReference(reference_index_base)
    manifest_object = manifest.LabelsAndIndices(manifest_file)
    bowtie_command = [bowtie2_exe,
        bowtie2_args if bowtie2_args is not None else '',
        '-t --no-hd --mm -x', bowtie2_index_base, '--12', reads_file]
    print >>sys.stderr, 'Starting Bowtie2 with command: ' \
         + ' '.join(bowtie_command)
    bowtie_process = subprocess.Popen(bowtie_command, bufsize=-1,
        stdout=subprocess.PIPE, stderr=sys.stderr)
    output_thread = BowtieOutputThread(
                        bowtie_process.stdout,
                        reference_index=reference_index,
                        manifest_object=manifest_object,
                        exon_differentials=exon_differentials, 
                        exon_intervals=exon_intervals, 
                        bin_size=bin_size,
                        verbose=verbose, 
                        output_stream=output_stream,
                        stranded=stranded,
                        report_multiplier=report_multiplier
                    )
    output_thread.start()
    # Join thread to pause execution in main thread
    if verbose: print >>sys.stderr, 'Joining thread...'
    output_thread.join()
    bowtie_process.wait()
    output_stream.flush()
    print >> sys.stderr, 'DONE with realign_reads.py; in/out=%d/%d; ' \
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
    parser.add_argument(\
        '--keep-alive', action='store_const', const=True, default=False,
        help='Prints reporter:status:alive messages to stderr to keep EMR '
             'task alive')
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
             'basename')
    parser.add_argument('--manifest', type=str, required=False,
        default='manifest',
        help='Path to manifest file')
    parser.add_argument('--archive', metavar='PATH', type=str, 
        default=None,
        help='Save output and Bowtie command to a subdirectory (named using ' 
             'this process\'s PID) of PATH')

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
    temp_dir_path = tempfile.mkdtemp()
    archive = os.path.join(args.archive,
        str(os.getpid())) if args.archive is not None else None
    # Handle temporary directory if CTRL+C'd
    atexit.register(handle_temporary_directory, archive, temp_dir_path)
    if args.verbose:
        print >>sys.stderr, 'Creating temporary directory %s' \
            % temp_dir_path
    go(bowtie2_exe=args.bowtie2_exe,
        reference_index_base=args.original_idx,
        bowtie2_build_exe=args.bowtie2_build_exe,
        bowtie2_args=bowtie_args,
        manifest_file=args.manifest,
        temp_dir_path=temp_dir_path,
        verbose=args.verbose, 
        bin_size=args.partition_length,
        exon_differentials=args.exon_differentials,
        exon_intervals=args.exon_intervals,
        stranded=args.stranded,
        report_multiplier=args.report_multiplier,
        keep_alive=args.keep_alive)
elif __name__ == '__main__':
    # Test units
    del sys.argv[1:] # Don't choke on extra command-line parameters
    import unittest

    class MultireadWithIntrons(unittest.TestCase):
        """ Tests multiread_with_introns(); needs no fixture 

            Some examples are ripped from:
            http://onetipperday.blogspot.com/2012/07/
            deeply-understanding-sam-tags.html; others are from actual
            SAM output of a dmel simulation
        """
        def test_multiread_1(self):
            """ Fails if SAM output is not adjusted properly. """
            pass

    unittest.main()
