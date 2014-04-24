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
9. Number of nucleotides spanned by EC on the left (that is, towards the 5'
    end of the read) of the intron, ASSUMING THE SENSE STRAND IS THE FORWARD
    STRAND.
10. Number of nucleotides spanned by EC on the right (that is, towards the 3'
    end of the read) of the intron, ASSUMING THE SENSE STRAND IS THE FORWARD
    STRAND.
--------------------------------------------------------------------
11. Number of instances of intron, insertion, or deletion in sample; this is
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

base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for directory_name in ['bowtie', 'dooplicity', 'sample',
                        'interval', 'manifest']:
    site.addsitedir(os.path.join(base_path, directory_name))

import dooplicity as dp
import manifest
import bowtie
import bowtie_index
import sample
import partition

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
    prefasta_filename = os.path.join(temp_dir_path, 'temp.prefa')
    deduped_fasta_filename = os.path.join(temp_dir_path, 'temp.deduped.prefa')
    final_fasta_filename = os.path.join(temp_dir_path, 'temp.fa')
    reads_filename = os.path.join(temp_dir_path, 'reads.temp')
    if verbose:
        print >>sys.stderr, 'Writing prefasta and input reads...'
    with open(prefasta_filename, 'w') as fasta_stream:
        with open(reads_filename, 'w') as read_stream:
            for read_seq, xpartition in dp.xstream(input_stream, 1):
                for value in xpartition:
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
            r'''sort %s | sed '$!N; /^\(.*\)\n\1$/!P; D' > %s'''
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
        raise RuntimeError('Bowtie index construction failed w/ exitlevel %d.'
                                % bowtie_build_process.returncode)

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
        assert reverse_strand_string == '+' or reverse_strand_string == '-', \
            'Strand is %s, which is invalid' % tokens[-4]
        reverse_strand = (True if reverse_strand_string == '-' else False)
        flag = int(multiread[i][1])
        if stranded and (flag & 16 != 0) == reverse_strand:
            # Strand of alignment doesn't agree with strand of intron
            continue
        rname = tokens[-4][:-1]
        offset = int(multiread[i][3]) - 1
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
        old_cigar_list = re.split(r'([MIDS])', multiread[i][5])[:-1]
        old_cigar_chars = old_cigar_list[1::2]
        old_cigar_sizes = map(int, old_cigar_list[::2])
        assert len(old_cigar_chars) == len(old_cigar_sizes)
        aligned_base_count = sum([old_cigar_sizes[j]
                                    for j in xrange(len(old_cigar_sizes))
                                    if old_cigar_chars[j] == 'M'])
        end_offset = aligned_base_count + offset
        for j, partial_size in enumerate(partial_sizes):
            if partial_size >= end_offset:
                end_index = j
                break
        cigar_sizes[start_index*2] = partial_sizes[start_index] - offset
        cigar_sizes[end_index*2] = end_offset - partial_sizes[end_index-1]
        pos = int(tokens[-3]) + new_offset + \
                sum(cigar_sizes[:start_index*2])
        if start_index is None or end_index is None:
            raise RuntimeError('Invalid SAM line; sum of exon sizes doesn\'t '
                               'agree with size of reference sequence.')
        if start_index == end_index:
            # Just use old cigar
            corrected_multiread.add(
                    (multiread[i][0], multiread[i][1], rname, str(pos))
                        + tuple(multiread[i][4:])
                )
            continue
        assert start_index < end_index
        matches = []
        skips = []
        for j in xrange(start_index*2, end_index*2+1):
            if j % 2 == 0:
                matches.append(cigar_sizes[j])
            else:
                skips.append(cigar_sizes[j])
        # Transform original cigar into final cigar
        match_index, exon_sum, new_exon_sum = 0, 0, 0
        new_cigar = []
        for j, old_cigar_size in enumerate(old_cigar_sizes):
            if old_cigar_chars[j] == 'M':
                exon_sum += old_cigar_size
                bases_to_align = exon_sum - new_exon_sum
                if bases_to_align > 0:
                    while True:
                        if matches[match_index] < bases_to_align:
                            if matches[match_index] > 0:
                                new_cigar.append(
                                        str(matches[match_index]) + 'M' +
                                        str(skips[match_index]) + 'N'
                                    )
                            else:
                                assert matches[match_index] == 0
                                new_cigar.append(
                                        str(skips[match_index]) + 'N'
                                    )
                            bases_to_align -= matches[match_index]
                            match_index += 1
                        else:
                            new_cigar.append(str(bases_to_align) + 'M')
                            matches[match_index] -= bases_to_align
                            break
                    new_exon_sum = exon_sum
                else:
                    new_cigar.append(
                            str(old_cigar_size) + 'M'
                        )
            elif old_cigar_chars[j] in 'DIS':
                new_cigar.append(str(old_cigar_size) + old_cigar_chars[j])
            else:
                raise RuntimeError('Bowtie2 CIGAR chars are expected to be '
                                   'in set (DIMS).')
        corrected_multiread.add(
                    (multiread[i][0], multiread[i][1], rname, str(pos),
                        multiread[i][4], ''.join(new_cigar))
                    + tuple(multiread[i][6:])
                    + (('XS:A:' + reverse_strand_string),)
                )
    NH_field = 'NH:i:' + str(len(corrected_multiread))
    multiread_to_return = [alignment + (NH_field,) for alignment in
                            corrected_multiread]
    return multiread_to_return

def filtered_cigar(cigar):
    """ Filters a cigar string so it contains only Ns and Ms.

        This facilitates outputting exon_diffs and introns.

        cigar: a cigar string with Ms, Ds, Is, Ns, and Ss.

        Return value: a cigar string with Ms and Ns.
    """
    cigar_to_filter = re.split(r'([MINDS])', cigar)[:-1]
    cigar_chars, cigar_sizes = [], []
    for i in xrange(0, len(cigar_to_filter), 2):
        if cigar_to_filter[i+1] in 'MN':
            cigar_chars.append(cigar_to_filter[i+1])
            cigar_sizes.append(int(cigar_to_filter[i]))
    last_cigar_char = cigar_chars[0]
    last_cigar_size = cigar_sizes[0]
    new_cigar = []
    for i, char in enumerate(cigar_chars[1:]):
        if char != last_cigar_char:
            new_cigar.append(str(last_cigar_size) + last_cigar_char)
            last_cigar_size = cigar_sizes[i+1]
            last_cigar_char = char
        else:
            last_cigar_size += cigar_sizes[i+1]
    new_cigar.append(str(last_cigar_size) + last_cigar_char)
    return ''.join(new_cigar)

def parsed_md(md):
    """ Divides an MD string up by boundaries between ^, letters, and numbers

        md: an MD string (example: 33A^CC).

        Return value: MD string split by boundaries described above.
    """
    md_to_parse = []
    md_group = [md[0]]
    for i, char in enumerate(md):
        if i == 0: continue
        if (re.match('[A-Za-z]', char) is not None) \
            != (re.match('[A-Za-z]', md[i-1]) is not None) or \
            (re.match('[0-9]', char) is not None) \
            != (re.match('[0-9]', md[i-1]) is not None):
            if md_group:
                md_to_parse.append(''.join(md_group))
            md_group = [char]
        else:
            md_group.append(char)
    if md_group:
        md_to_parse.append(''.join(md_group))
    return md_to_parse

def indels(cigar, md, pos, seq):
    """ Computes indels from CIGAR, MD string, and POS of a given alignment.

        cigar: CIGAR string
        md: MD:Z string
        pos: position of first aligned base
        seq: read sequence

        Return value: tuple (insertions, deletions), where insertions is a
            list of tuples (last genomic position before insertion, 
                                string of inserted bases), and deletions
            is a list of tuples (first genomic position of deletion,
                                    string of deleted bases)
    """
    insertions, deletions = [], []
    cigar = re.split(r'([MINDS])', cigar)[:-1]
    md = parsed_md(md)
    cigar_chars, cigar_sizes = [], []
    cigar_index, md_index, seq_index = 0, 0, 0
    max_cigar_index = len(cigar)
    while cigar_index != max_cigar_index:
        if cigar[cigar_index+1] == 'M':
            aligned_base_cap = int(cigar[cigar_index])
            aligned_bases = 0
            while True:
                try:
                    aligned_bases += int(md[md_index])
                    if aligned_bases <= aligned_base_cap:
                        md_index += 1
                except ValueError:
                    # Not an int, but should not have reached a deletion
                    assert md[md_index] != '^'
                    if aligned_bases + len(md[md_index]) > aligned_base_cap:
                        md[md_index] = md[md_index][
                                            :aligned_base_cap-aligned_bases
                                        ]
                        aligned_bases = aligned_base_cap
                    else:
                        aligned_bases += len(md[md_index])
                        md_index += 1
                if aligned_bases > aligned_base_cap:
                    md[md_index] = aligned_bases - aligned_base_cap
                    break
                elif aligned_bases == aligned_base_cap:
                    break
            pos += aligned_base_cap
            seq_index += aligned_base_cap
        elif cigar[cigar_index+1] == 'N':
            # Skip region of reference
            pos += int(cigar[cigar_index])
        elif cigar[cigar_index+1] == 'I':
            # Insertion
            insert_size = int(cigar[cigar_index])
            insertions.append(
                    (pos - 1, seq[seq_index:seq_index+insert_size])
                )
            seq_index += insert_size
        elif cigar[cigar_index+1] == 'D':
            assert md[md_index] == '^'
            # Deletion
            deletions.append((pos, md[md_index+1]))
            delete_size = int(cigar[cigar_index])
            assert len(md[md_index+1]) == delete_size
            # Skip deleted part of reference
            pos += delete_size
            md_index += 2
        else:
            # Soft clip; ignore
            assert cigar[cigar_index+1] == 'S'
        cigar_index += 2
    return insertions, deletions

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
                                [self.manifest_object.label_to_index[
                                        sample.parseLab(last_tokens[0][:-2])
                                    ],
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
                    sample_label = self.manifest_object.label_to_index[
                                        sample.parseLab(multiread[0][0][:-2])
                                        ]
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
                        '''Output exonic chunks and introns only if there is
                        exactly one alignment.'''
                        alignment = list(corrected_multiread)[0]
                        cigar = alignment[5]
                        cigar_for_introns = filtered_cigar(cigar)
                        pos = int(alignment[3])
                        rname = alignment[2]
                        cigar_chars = [-1] + [i for i, char
                            in enumerate(cigar_for_introns) 
                            if char == 'M' or char == 'N']
                        base_counts = \
                            [int(cigar_for_introns[
                                    (cigar_chars[i]+1):cigar_chars[i+1]
                                ]) for i in xrange(len(cigar_chars)-1)]
                        if 'I' in cigar or 'D' in cigar:
                            # Output indels
                            md = [field for field in alignment
                                    if field[:5] == 'MD:Z:'][0][5:]
                            insertions, deletions = indels(cigar, md, pos,
                                                            alignment[9])
                            for insert_pos, insert_seq in insertions:
                                print >>self.output_stream, (
                                       ('bed\tI\t%s\t%s\t%012d\t%012d\t%s'
                                        '\t\x1c\t\x1c\t\x1c\t\x1c\t1')
                                        % (sample_label, self.reference_index.\
                                            rname_to_string[rname],
                                            insert_pos, insert_pos,
                                            insert_seq)
                                    )
                                _output_line_count += 1
                            for del_pos, del_seq in deletions:
                                print >>self.output_stream, (
                                       ('bed\tD\t%s\t%s\t%012d\t%012d\t%s'
                                        '\t\x1c\t\x1c\t\x1c\t\x1c\t1')
                                        % (sample_label, self.reference_index.\
                                            rname_to_string[rname],
                                            del_pos, del_pos + len(del_seq),
                                            del_seq)
                                    )
                                _output_line_count += 1
                        if 'N' in cigar_for_introns:
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
                                            sum(base_counts[i+1::2]))
                                            for i in xrange(1,
                                                len(base_counts), 2)]
                            for (intron_pos, intron_end_pos,
                                    left_anchor_size, right_anchor_size,
                                    left_displacement, right_displacement) \
                                in introns:
                                print >>self.output_stream, (
                                        ('bed\tN\t%s\t%s\t%012d\t%012d\t%s\t'
                                         '%d\t%d\t%d\t%d\t1')
                                         % (sample_label, 
                                            self.reference_index.\
                                            rname_to_string[rname],
                                            intron_pos, intron_end_pos,
                                            reverse_strand_string,
                                            left_displacement,
                                            right_displacement,
                                            left_anchor_size,
                                            right_anchor_size)
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
        9. Number of nucleotides spanned by EC on the left (that is, towards
            the 5' end of the read) of the intron, ASSUMING THE SENSE STRAND IS
            THE FORWARD STRAND.
        10. Number of nucleotides spanned by EC on the right (that is, towards
            the 3' end of the read) of the intron, ASSUMING THE SENSE STRAND IS
            THE FORWARD STRAND.
        ---------------------------------------------------------------------
        11. Number of instances of intron, insertion, or deletion in sample;
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
    import time
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
    bowtie_command = ' '.join([bowtie2_exe,
        bowtie2_args if bowtie2_args is not None else '',
        '-t --no-hd --mm', bowtie2_index_base, '--12', reads_file])
    print >>sys.stderr, 'Starting Bowtie2 with command: ' + bowtie_command
    bowtie_process = subprocess.Popen(bowtie_command, bufsize=-1, shell=True,
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
    #bowtie_process.wait()
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
    class TestIndels(unittest.TestCase):
        """ Tests indels(); needs no fixture 

            Some examples are ripped from:
            http://onetipperday.blogspot.com/2012/07/
            deeply-understanding-sam-tags.html; others are from actual
            SAM output of a dmel simulation
        """
        def test_read_1(self):
            """ Fails if example doesn't give expected indels."""
            self.assertEquals(([], [(18909816, 'GG')]),
                         indels('20M2D9M', '20^GG7A1', 18909796,
                                'TAGCCTCTGTCAGCACTCCTGAGTTCAGA')
                    )

        def test_read_2(self):
            """ Fails if example doesn't give expected indels."""
            self.assertEquals(([], [(73888560, 'GG')]),
                         indels('20M2D9M', '20^GG8C0', 73888540,
                                'TAGCCTCTGTCAGCACTCCTGAGTTCAGA')
                    )

        def test_read_3(self):
            """ Fails if example doesn't give expected indels."""
            self.assertEquals(([(20620369, 'CA')], [(20620365, 'GT')]),
                         indels('20M151N47M2D3M2I4M', '67^GT3T2C0', 20620147,
                                'CCGCACCCGTACTGCTACAGATTTCCATCATCGCCACCCGCGGGC'
                                'ATTCTGAAAAAGAGCGACGAAGAAGCAACCT')
                    )

        def test_read_4(self):
            """ Fails if example doesn't give expected indels."""
            self.assertEquals(([(20620155, 'CT')], []),
                         indels('9M2I63M70N2M', '1A2C1A0G1G1C1C0C1G2A54',
                                 20620147,
                                'TTCTNCCTGCTTGTATGACCGTGTTGGGCGTGAGTGGCTTGTCCC'
                                'TCAAGTAGAGACCATAGCGAGATGGGTACCT')
                    )

    unittest.main()