#!/usr/bin/env python
"""
Rail-RNA-compare_alignments
Follows Rail-RNA-align_reads / Rail-RNA-realign_reads
Precedes Rail-RNA-collapse / Rail-RNA-bam / Rail-RNA-bed_pre
    / Rail-RNA-intron_coverage

Combines alignments of same reads from previous steps (so far,
Rail-RNA-align_reads and Rail-RNA-realign_reads) so alignment scores can be 
compared. If there are no ties in the top alignment score, alignments are
printed. Otherwise, alignment printing is postponed until Rail-RNA-break_ties,
when coverage of alignments by introns is used to break ties.

Input (read from stdin)
----------------------------
Tab-delimited input tuple columns:
Standard 11+ - column SAM input.

Input is partitioned by QNAME, the first field of SAM.

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

Format 3 (sam_intron_ties) output only for ties in alignment score (which are
    within some tie margin as decided by multiread_to_report);
tab-delimited output tuple columns:
Standard SAM output except fields are in different order -- and the first four
fields include sample/intron information. If an alignment overlaps k introns, k
lines are output.
1. The character 'N' so the line can be matched up
    with intron bed lines
2. Sample index
3. Number string representing RNAME; see BowtieIndexReference class
    in bowtie_index for conversion information
4. Intron start position
5. Intron end position
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

Format 4 (sam_clip_ties) output only for ties in alignment score when
    no introns are overlapped -- these alignments are almost invariably
    soft-clipped
[SAME AS SAM FIELDS; see SAM format specification]

Format 5 (sam); tab-delimited output tuple columns:
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

Format 6; tab-delimited output tuple columns:
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
import site
import string
import time
import argparse
import random

base_path = os.path.abspath(
                    os.path.dirname(os.path.dirname(os.path.dirname(
                        os.path.realpath(__file__)))
                    )
                )
utils_path = os.path.join(base_path, 'rna', 'utils')
site.addsitedir(utils_path)
site.addsitedir(base_path)

from dooplicity.tools import xstream
from alignment_handlers \
    import multiread_with_introns, AlignmentPrinter, multiread_to_report
import partition
import manifest
import bowtie
import bowtie_index

if __name__ == '__main__':
    # Print file's docstring if -h is invoked
    parser = argparse.ArgumentParser(description=__doc__, 
                formatter_class=argparse.RawDescriptionHelpFormatter)
    # Add command-line arguments
    parser.add_argument('--tie-margin', type=int, required=False,
        default=6,
        help='Allowed score difference per 100 bases among ties in '
             'max score. For example, 150 and 144 are tied alignment scores '
             'for a 100-bp read when --tie-margin is 6.')
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

    # Add command-line arguments for dependencies
    partition.add_args(parser)
    bowtie.add_args(parser)
    manifest.add_args(parser)

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

    args = parser.parse_args(argv[1:])

    reversed_complement_translation_table = string.maketrans('ATCG', 'TAGC')
    manifest_object = manifest.LabelsAndIndices(args.manifest)
    reference_index = bowtie_index.BowtieIndexReference(
                                            args.bowtie_idx
                                        )
    alignment_printer = AlignmentPrinter(
                    manifest_object,
                    reference_index,
                    bin_size=args.partition_length,
                    output_stream=sys.stdout,
                    exon_ivals=args.exon_intervals,
                    exon_diffs=args.exon_differentials
                )
    alignment_count_to_report, seed, non_deterministic \
                = bowtie.parsed_bowtie_args(bowtie_args)

    input_line_count = 0
    output_line_count = 0

    start_time = time.time()
    for (qname,), xpartition in xstream(sys.stdin, 1):
        # While labeled multiread, this list may end up simply a uniread
        initial_multiread \
            = [(qname,) + rest_of_line for rest_of_line in xpartition]
        input_line_count += len(initial_multiread)
        multiread = [alignment for alignment in initial_multiread
                        if not (int(alignment[1]) & 4)]
        flag = int(initial_multiread[0][1])
        if not multiread:
            # Write only the SAM output if the read was unmapped
            output_line_count += alignment_printer.print_unmapped_read(
                                                    qname,
                                                    initial_multiread[0][9],
                                                    initial_multiread[0][10]
                                                )
        else:
            '''Correct positions to match original reference's, correct
            CIGARs, eliminate duplicates, and decide primary alignment.'''
            try:
                corrected_multiread = multiread_with_introns(
                                            multiread,
                                            stranded=args.stranded
                                        )
            except:
                print >>sys.stderr, ('Error encountered interpreting '
                                     'multiread %s' % (multiread,))
                raise
            if not corrected_multiread:
                '''This is effectively an unmapped read; write
                corresponding SAM output.'''
                if flag & 16:
                    seq_to_write = initial_multiread[0][9][::-1].translate(
                                    reversed_complement_translation_table
                                )
                    qual_to_write = initial_multiread[0][10][::-1]
                else:
                    seq_to_write = initial_multiread[0][9]
                    qual_to_write = initial_multiread[0][10]
                output_line_count \
                    += alignment_printer.print_unmapped_read(
                                                    qname,
                                                    seq_to_write,
                                                    qual_to_write
                                                )
                continue
            output_line_count += alignment_printer.print_alignment_data(
                multiread_to_report(
                    corrected_multiread,
                    alignment_count_to_report=alignment_count_to_report,
                    seed=seed,
                    non_deterministic=non_deterministic,
                    tie_margin=args.tie_margin
                )
            )

    print >>sys.stderr, 'DONE with compare_alignments.py; in/out=%d/%d; ' \
        'time=%0.3f s' % (input_line_count, output_line_count,
                            time.time() - start_time)