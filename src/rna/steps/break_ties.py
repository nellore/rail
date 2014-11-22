#!/usr/bin/env python
"""
Rail-RNA-break-ties
Follows Rail-RNA-intron_coverage, Rail-RNA-realign_reads
Precedes Rail-RNA-bed_pre, Rail-RNA-bam

Decides primary alignments from among ties according to the read coverage
(by "uniquely" aligning reads) of the minimally covered intron overlapped
by the alignment. The rules are as follows.
-Break ties first by selecting from among alignments that overlap the fewest
introns.
-If no introns are overlapped in this group of alignments, break ties
uniformly at random from among max scores
-Otherwise, each alignment is weighted by the read coverage of its minimally
covered intron, and the tie is broken at random according to these weights.
Examples: if all weights are 0, the tie is broken uniformly at random. If there
are two alignments, and one is weighted 0 while the other is weighted 2, the
second alignment is always chosen. If there are two alignments, and one is
weighted 1 while the other is weighted 2, the first is chosen with probability
1/3 and the second is chosen with probability 1/2.

The seed qname + seq + qual is used.

Input (read from stdin)
----------------------------
Two formats -- format 1's tab-delimited input columns (introns):
Format 1's tab-deliminted input columns (introns):
[ALL SAM FIELDS; see SAM format specification for details]
Last field (only for alignments overlapping introns): XC:i:(coverage of an
intron from the read); there are as many lines for a given alignment
(i.e., QNAMEs) as there are introns overlapped by the alignment

Input is partitioned by QNAME (field 1).

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

Format 3 (sam); tab-delimited output tuple columns:
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

Format 4; tab-delimited output tuple columns:
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
import time
from collections import defaultdict

base_path = os.path.abspath(
                    os.path.dirname(os.path.dirname(os.path.dirname(
                        os.path.realpath(__file__)))
                    )
                )
utils_path = os.path.join(base_path, 'rna', 'utils')
site.addsitedir(utils_path)
site.addsitedir(base_path)

from dooplicity.tools import xstream
from alignment_handlers import AlignmentPrinter, multiread_to_report
import bowtie
import bowtie_index
import manifest
import partition

import argparse
# Print file's docstring if -h is invoked
parser = argparse.ArgumentParser(description=__doc__, 
            formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('--exon-differentials', action='store_const',
    const=True,
    default=True, 
    help='Print exon differentials (+1s and -1s)')
parser.add_argument('--exon-intervals', action='store_const',
    const=True,
    default=False, 
    help='Print exon intervals')
parser.add_argument('--drop-deletions', action='store_const',
        const=True,
        default=False, 
        help='Drop deletions from coverage vectors')

bowtie.add_args(parser)
manifest.add_args(parser)
partition.add_args(parser)

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

reference_index = bowtie_index.BowtieIndexReference(args.bowtie_idx)
manifest_object = manifest.LabelsAndIndices(args.manifest)
alignment_count_to_report, seed, non_deterministic \
    = bowtie.parsed_bowtie_args(bowtie_args)

alignment_printer = AlignmentPrinter(manifest_object,
                                        reference_index,
                                        output_stream=sys.stdout,
                                        bin_size=args.partition_length,
                                        exon_ivals=args.exon_intervals,
                                        exon_diffs=args.exon_differentials,
                                        drop_deletions=args.drop_deletions)
input_line_count, output_line_count = 0, 0
start_time = time.time()

for (qname,), xpartition in xstream(sys.stdin, 1):
    alignments = [(qname,) + alignment for alignment in xpartition]
    input_line_count += len(alignments)
    '''Separate into alignments that overlap the fewest introns and
    alignments that don't.'''
    intron_counts = [alignment[5].count('N') for alignment in alignments]
    min_intron_count = min(intron_counts)
    if not min_intron_count:
        '''There is at least one alignment that overlaps no introns; report 
        an alignment with highest score at random.'''
        clipped_alignments = [alignments[i] for i in xrange(len(intron_counts))
                                if intron_counts[i] == 0]
        alignments_and_scores = [(alignment, [int(tokens[5:])
                                                for tokens in alignment
                                                if tokens[:5] == 'AS:i:'][0])
                                    for alignment in clipped_alignments]
        alignments_and_scores.sort(key=lambda alignment: alignment[1],
                                    reverse=True)
        max_score = alignments_and_scores[0][1]
        alignments_to_report = [alignment for alignment, score
                                    in alignments_and_scores
                                    if score == max_score]
        output_line_count += alignment_printer.print_alignment_data(
                    multiread_to_report(
                        alignments_to_report,
                        alignment_count_to_report=alignment_count_to_report,
                        seed=seed,
                        non_deterministic=non_deterministic,
                        weights=[1]*len(alignments_to_report)
                    )
                )
    else:
        # All alignments overlap introns
        intron_alignments = [alignments[i] for i in xrange(len(intron_counts))
                                if intron_counts[i] == min_intron_count]
        # Compute min coverage for alignments that are the same
        alignment_dict = defaultdict(int)
        for alignment in intron_alignments:
            if alignment[:-1] not in alignment_dict:
                alignment_dict[alignment[:-1]] = int(alignment[-1][5:])
            else:
                alignment_dict[alignment[:-1]] \
                    = min(alignment_dict[alignment[:-1]],
                            int(alignment[-1][5:]))
        weights = alignment_dict.values()
        if not any(weights): weights = [1] * len(weights)
        output_line_count += alignment_printer.print_alignment_data(
                    multiread_to_report(
                        [alignment + ('XC:i:%s' % alignment_dict[alignment],)
                            for alignment in alignment_dict.keys()],
                        alignment_count_to_report=alignment_count_to_report,
                        seed=seed,
                        non_deterministic=non_deterministic,
                        weights=weights
                    )
                )

print >>sys.stderr, 'DONE with break_ties.py; in/out=%d/%d; ' \
                    'time=%0.3f s' % (input_line_count, output_line_count,
                                        time.time() - start_time)