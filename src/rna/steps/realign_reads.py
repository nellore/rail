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

THIS CODE ASSUMES THAT BOWTIE RUNS ON JUST ONE THREAD; it exploits that Bowtie
returns reads in the order in which they were sent, which is guaranteed only
when on a single thread.

Input (read from stdin)
----------------------------
Tab-delimited input tuple columns:
(two kinds)

Type 1:
1. Read sequence
2. '\x1c' + FASTA reference name including '>'. The following format is used:
    original RNAME + '+' or '-' indicating which strand is the sense strand
    + '\x1d' + start position of sequence + '\x1d' + comma-separated list of
    subsequence sizes framing introns + '\x1d' + comma-separated list of intron
    sizes) + '\x1d' + 'p' if derived from primary alignment to genome; 's' if
    derived from secondary alignment to genome; 'i' if derived from cointron
    search
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

Format 3 (sam_ties) output only for ties in alignment score (which are within
    some tie margin as decided by multiread_to_report);
tab-delimited output tuple columns:
Standard SAM output except fields are in different order -- and the first four
fields include sample/intron information. If an alignment overlaps k introns, k
lines are output.
1. The character 'N' so the line can be matched up
    with intron bed lines
2. Sample index
3. Number string representing RNAME; see BowtieIndexReference class
    in bowtie_index for conversion information
4. Intron start position or '\x1c' if no introns present
5. Intron end position or '\x1c' if no introns present
6. '+' or '-' indicating which strand is sense strand
7. POS
8. QNAME
9. FLAG
10. MAPQ
11. CIGAR
12. RNEXT
13. PNEXT
14. TLEN
15. SEQ
16. QUAL
... + optional fields

Format 4 (splice_sam); tab-delimited output tuple columns:
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

Introns (intron_bed), insertions/deletions (indel_bed)

Format 5; tab-delimited output tuple columns:
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
from alignment_handlers import multiread_with_introns, AlignmentPrinter
import argparse
import random

# Initialize global variables for tracking number of input/output lines
_input_line_count = 0
_output_line_count = 0

_reversed_complement_translation_table = string.maketrans('ATCG', 'TAGC')

def input_files_from_input_stream(input_stream,
                                    temp_dir_path=tempfile.mkdtemp(),
                                    verbose=False):
    """ Creates FASTA reference to index, rname file, and file with reads.

        Each line of the read file is in the following format:

        read number + TAB + SEQ + TAB + QUAL

        Each line of the rname file is in the following format:

        (number of qnames associated with rnames) + '\x1e'
        + ('\x1e'-separated list of valid rnames) .

        The RNAME file contains those RNAMEs to which a given sequence is
        "allowed" to align. When Bowtie2 is run in -a mode and only these
        alignments are permitted, this ensures the reproducibility of
        results regardless of partitioning (and hence # cores).

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
    rname_filename = os.path.join(temp_dir_path, 'rnames.temp')
    if verbose:
        print >>sys.stderr, 'Writing prefasta and input reads...'
    with open(prefasta_filename, 'w') as fasta_stream:
        with open(reads_filename, 'w') as read_stream:
            with open(rname_filename, 'w') as rname_stream:
                for read_seq, xpartition in xstream(input_stream, 1):
                    rnames = []
                    read_count = 0
                    for value in xpartition:
                        _input_line_count += 1
                        if value[0][0] == '\x1c':
                            # Print FASTA line
                            print >>fasta_stream, '\t'.join([value[0][1:-2],
                                                             value[1]])
                            rnames.append(value[0][2:-2])
                        else:
                            # Add to temporary seq stream
                            print >>read_stream, '\t'.join([value[0], read_seq,
                                                                value[1]])
                            read_count += 1
                    rname_stream.write(str(read_count) + '\x1e')
                    print >>rname_stream, '\x1e'.join(rnames)
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
    return final_fasta_filename, reads_filename, rname_filename

def create_index_from_reference_fasta(bowtie2_build_exe, fasta_file,
        index_basename, keep_alive=False):
    """ Creates Bowtie2 index from reference fasta.

        bowtie2_build_exe: Path to Bowtie2
        fasta_file: Path to reference FASTA to index
        index_basename: Path to index basename to be created

        Return value: return value of bowtie-build process
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
        return bowtie_build_thread.bowtie_build_process
    else:
        bowtie_build_process = subprocess.Popen(
                                    [args.bowtie2_build_exe,
                                        fasta_file,
                                        index_basename],
                                    stderr=sys.stderr,
                                    stdout=sys.stderr
                                )
        bowtie_build_process.wait()
        return bowtie_build_process.returncode

class BowtieOutputThread(threading.Thread):
    """ Processes Bowtie alignments, emitting tuples for exons and introns. """
    
    def __init__(self, input_stream, rname_stream, reference_index,
        manifest_object, bowtie2_args, return_set, tie_margin=6,
        output_stream=sys.stdout, exon_differentials=True,
        exon_intervals=False, stranded=False,
        verbose=False, bin_size=10000, report_multiplier=1.2):
        """ Constructor for BowtieOutputThread.

            input_stream: where to retrieve Bowtie's SAM output, typically a
                Bowtie process's stdout.
            rname_stream: where to retrieve RNAMES, each of which is in
                the form (number of qnames associated with rnames) + '\x1e'
                + ('\x1e'-separated list of valid rnames)
            reference_index: object of class BowtieIndexReference; for
                outputing RNAME number strings to facilitate later sorting
            manifest_object: object of class LabelsAndIndicesWithClusters that
                maps indices to labels and back; used to shorten intermediate
                output
            bowtie2_args: Bowtie 2 arguments to be parsed
            tie_margin: allowed score difference per 100 bases among ties in
                max score. For example, 150 and 144 are tied alignment scores
                for a 100-bp read when tie_margin is 6.
            return_set: 0 is added to this set if the thread finishes
                successfully
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
        self.rname_stream = rname_stream
        self.reference_index = reference_index
        self.output_stream = output_stream
        self.verbose = verbose
        self.bin_size = bin_size
        self.stranded = stranded
        self.report_multiplier = report_multiplier
        self.exon_differentials = exon_differentials
        self.exon_intervals = exon_intervals
        self.manifest_object = manifest_object
        self.return_set = return_set
        self.bowtie2_args = bowtie2_args
        self.tie_margin = tie_margin

    def run(self):
        """ Prints SAM, exon_ivals, exon_diffs, and introns.

            Overrides default method containing thread activity.

            No return value.
        """
        global _output_line_count
        next_report_line = 0
        i = 0
        tokens = rname_stream.readline().strip().split('\x1e')
        qname_total, rnames = int(tokens[0]), set(tokens[1:])
        qname_count = 0
        alignment_printer = AlignmentPrinter(
                self.manifest_object,
                self.reference_index,
                output_stream=sys.stdout,
                exon_ivals=exon_intervals,
                exon_diffs=exon_differentials
            )
        alignment_count_to_report, seed, non_deterministic \
            = bowtie.parsed_bowtie_args(bowtie2_args)
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
                rname = rest_of_line[1]
                if rname in rnames:
                    multiread.append(list((qname,) + rest_of_line))
            if flag & 4:
                # Write only the SAM output if the read was unmapped
                _output_line_count += alignment_printer.print_unmapped_read(
                                                        qname,
                                                        rest_of_line[8],
                                                        rest_of_line[9]
                                                    )
            else:
                '''Correct positions to match original reference's, correct
                CIGARs, eliminate duplicates, and decide primary alignment.'''
                sample_index = self.manifest_object.label_to_index[
                                        multiread[0][0].rpartition('\x1d')[2]
                                    ]
                corrected_multiread = multiread_with_introns(
                                            multiread,
                                            stranded=self.stranded
                                        )
                if not corrected_multiread:
                    '''This is effectively an unmapped read; write
                    corresponding SAM output.'''
                    if (int(multiread[0][1]) & 16):
                        multiread[0][9] = multiread[0][9][::-1].translate(
                                        _reversed_complement_translation_table
                                    )
                        multiread[0][10] = multiread[0][10][::-1]
                    _output_line_count \
                        += alignment_printer.print_unmapped_read(
                                                        qname,
                                                        multiread[0][9],
                                                        multiread[0][10]
                                                    )
                    continue
                _output_line_count += alignment_printer.print_alignment_data(
                    multiread_to_report(
                        corrected_multiread,
                        alignment_count_to_report=alignment_count_to_report,
                        seed=0,
                        non_deterministic=non_deterministic,
                        tie_margin=self.tie_margin
                    )
                )
            qname_count += 1
            if qname_count == qname_total:
                tokens = rname_stream.readline().strip().split('\x1e')
                qname_total, rnames = int(tokens[0]), set(tokens[1:])
                qname_count = 0
        self.return_set.add(0)

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
    bowtie2_args=None, tie_margin=6,
    temp_dir_path=tempfile.mkdtemp(), bin_size=10000,
    manifest_file='manifest', verbose=False, exon_differentials=True,
    exon_intervals=False, stranded=False, report_multiplier=1.2,
    keep_alive=False):
    """ Runs Rail-RNA-realign.

        Uses Bowtie index including only sequences framing introns to align
        only those reads for which Bowtie did not report at least one alignment
        in Rail-RNA-align. Referencenames in this index encode intron sizes and
        locations in the (presmably) exonic sequences it records. Infers exonic
        chunks and introns from alignments.

        THIS CODE ASSUMES THAT BOWTIE RUNS ON JUST ONE THREAD; it exploits that
        Bowtie returns reads in the order in which they were sent, which is
        guaranteed only when on a single thread.

        Input (read from stdin)
        ----------------------------
        Tab-delimited input tuple columns:
        (two kinds)

        Type 1:
        1. Read sequence
        2. '\x1c' + FASTA reference name including '>'. The following format is
            used: original RNAME + '+' or '-' indicating which strand is the
            sense strand + '\x1d' + start position of sequence + '\x1d'
            + comma-separated list of subsequence sizes framing introns
            + '\x1d' + comma-separated list of intron sizes) + '\x1d'
            + 'p' if derived from primary alignment to genome; 's' if derived
            from secondary alignment to genome; 'i' if derived from cointron
            search
        3. FASTA sequence

        Type 2:
        1. Read sequence
        2. QNAME
        3. Quality sequence

        Type 1 corresponds to a FASTA line to index to which the read sequence
        is predicted to align. Type 2 corresponds to a distinct read. Input is
        partitioned by field 1.

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

        Format 3 (sam_ties) output only for ties in alignment score (which are
            within some tie margin as decided by multiread_to_report);
        tab-delimited output tuple columns:
        Standard SAM output except fields are in different order -- and the
        first four fields include sample/intron information. If an alignment
        overlaps k introns, k lines are output.
        1. The character 'N' so the line can be matched up
            with intron bed lines
        2. Sample index
        3. Number string representing RNAME; see BowtieIndexReference class
            in bowtie_index for conversion information
        4. Intron start position or '\x1c' if no introns present
        5. Intron end position or '\x1c' if no introns present
        6. '+' or '-' indicating which strand is sense strand
        7. POS
        8. QNAME
        9. FLAG
        10. MAPQ
        11. CIGAR
        12. RNEXT
        13. PNEXT
        14. TLEN
        15. SEQ
        16. QUAL
        ... + optional fields

        Format 4 (sam); tab-delimited output tuple columns:
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

        Introns (intron_bed), insertions/deletions (indel_bed)

        Format 5; tab-delimited output tuple columns:
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
        tie_margin: allowed score difference per 100 bases among ties in
            max score. For example, 150 and 144 are tied alignment scores
            for a 100-bp read when tie_margin is 6.
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
    fasta_file, reads_file, rnames_file = input_files_from_input_stream(
                                                input_stream,
                                                verbose=verbose,
                                                temp_dir_path=temp_dir_path
                                            )
    bowtie2_index_base = os.path.join(temp_dir_path, 'tempidx')
    bowtie_build_return_code = create_index_from_reference_fasta(
                                    bowtie2_build_exe,
                                    fasta_file,
                                    bowtie2_index_base)
    if bowtie_build_return_code == 0:
        reference_index = bowtie_index.BowtieIndexReference(
                                                reference_index_base
                                        )
        manifest_object = manifest.LabelsAndIndices(manifest_file)
        output_file = os.path.join(temp_dir_path, 'out.sam')
        bowtie_command = ' ' .join([bowtie2_exe,
            bowtie2_args if bowtie2_args is not None else '',
            '--local -t -a --no-hd --mm -x', bowtie2_index_base, '--12',
            reads_file, '-S', output_file])
        print >>sys.stderr, 'Starting Bowtie2 with command: ' + bowtie_command
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
        if os.path.exists(output_file):
            return_set = set()
            output_thread = BowtieOutputThread(
                                open(output_file),
                                open(rnames_file),
                                reference_index=reference_index,
                                manifest_object=manifest_object,
                                bowtie2_args=bowtie2_args,
                                return_set=return_set,
                                exon_differentials=exon_differentials, 
                                exon_intervals=exon_intervals, 
                                bin_size=bin_size,
                                verbose=verbose, 
                                output_stream=output_stream,
                                stranded=stranded,
                                report_multiplier=report_multiplier,
                                tie_margin=tie_margin,
                            )
            output_thread.start()
            # Join thread to pause execution in main thread
            if verbose: print >>sys.stderr, 'Joining thread...'
            output_thread.join()
            if not return_set:
                raise RuntimeError('Error occurred in BowtieOutputThread.')
    elif bowtie_build_return_code == 1:
        print >>sys.stderr, ('Bowtie build failed, but probably because '
                             'FASTA file was empty. Continuing...')
    else:
        raise RuntimeError('Bowtie build process failed with exitlevel %d.'
                            % bowtie_build_return_code)
    print >> sys.stderr, 'DONE with realign_reads.py; in/out=%d/%d; ' \
        'time=%0.3f s' % (_input_line_count, _output_line_count,
                                time.time() - start_time)

if __name__ == '__main__':
    # Print file's docstring if -h is invoked
    parser = argparse.ArgumentParser(description=__doc__, 
                formatter_class=argparse.RawDescriptionHelpFormatter)
    # Add command-line arguments
    parser.add_argument('--report_multiplier', type=float, required=False,
        default=1.2,
        help='When --verbose is also invoked, the only lines of lengthy '
             'intermediate output written to stderr have line number that '
             'increases exponentially with this base')
    parser.add_argument('--tie-margin', type=int, required=False,
        default=6,
        help='Allowed score difference per 100 bases among ties in '
             'max score. For example, 150 and 144 are tied alignment scores '
             'for a 100-bp read when --tie-margin is 6.')
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
        tie_margin=args.tie_margin,
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