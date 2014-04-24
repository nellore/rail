#!/usr/bin/env python
"""
Rail-RNA-junctions_bed_pre
Follows Rail-RNA-realign_reads
Precedes Rail-RNA-junctions_bed

Reduce step in MapReduce pipelines that computes coverage of introns and other
statistics by sample from Rail-RNA-realign and outputs BED lines for
consolidation in Rail-RNA-bed.

Input (read from stdin)
----------------------------
Tab-delimited input tuple columns:
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

Input is partitioned by column 1.

Hadoop output (written to stdout)
----------------------------
Tab-delimited input tuple columns are in sample label + 12-column BED format.
(See https://genome.ucsc.edu/FAQ/FAQformat.html#format1 for a detailed
description.)

1. Sample label
2. Number string representing chrom (chromosome name); see BowtieIndexReference
    class in bowtie_index for conversion information
3. chromStart (start position of region; 0-BASED)
4. chromEnd (end position of region; 0-BASED)
5. name (includes maximin anchor size and unique displacement count)
6. score (number of reads supporting this intron)
7. strand (forward is sense strand=+; reverse is sense strand=-)
8. thickStart (same as chromStart)
9. thickEnd (same as chromEnd)
10. itemRgb (always 227,29,118; that is, hot pink)
11. blockCount (always 2; each block denotes the maximal left/right overhang of
    any read spanning a junction)
12. blockSizes (comma-separated lengths of the maximal left/right overhangs of
    any read spanning a junction)
13. blockStarts (comma-separated list of the overhang start positions)

OUTPUT COORDINATES ARE 0-INDEXED.
"""
import os
import sys
import argparse
import site

base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
site.addsitedir(os.path.join(base_path, 'bowtie'))
import bowtie
import bowtie_index

# Print file's docstring if -h is invoked
parser = argparse.ArgumentParser(description=__doc__, 
            formatter_class=argparse.RawDescriptionHelpFormatter)
# So Bowtie index can be specified
bowtie.addArgs(parser)
args = parser.parse_args()

import time
start_time = time.time()

reference_index = bowtie_index.BowtieIndexReference(args.bowtie_idx)
maximin_anchor_size, unique_displacements, coverage = None, set(), 0
max_left_overhang, max_right_overhang = None, None
write_line = False
input_line_count, output_line_count = 0, 0
last_rname = None

while True:
    line = sys.stdin.readline().rstrip()
    if line:
        input_line_count += 1
        tokens = line.split('\t')
        rname, sample_label, pos, end_pos = tokens[0].split(';')
        pos, end_pos = int(pos), int(end_pos)
        left_overhang, right_overhang, left_anchor, right_anchor = \
            (int(tokens[1]), int(tokens[2]), int(tokens[3]),
                int(tokens[4]))
        reverse_strand_string = rname[-1]
        rname = reference_index.rname_to_string[rname[:-1]]
        if last_rname is not None and not (rname == last_rname 
            and sample_label == last_sample_label and pos == last_pos
            and end_pos == last_end_pos
            and reverse_strand_string == last_reverse_strand_string):
            write_line = True
    if not line: write_line = True
    if write_line and last_rname is not None:
        start_position = last_pos - max_left_overhang - 1
        end_position = last_end_pos + max_right_overhang - 1
        print ('bed\t%s\t%s\t%012d\t%012d\tmaximin_anchor_size=%d;' 
               'unique_displacement_count=%d\t%d\t%s\t%d\t%d\t227,29,118' 
               '\t2\t%d,%d\t0,%d') % (last_sample_label,
                                        last_rname,
                                        start_position,
                                        end_position,
                                        maximin_anchor_size,
                                        len(unique_displacements),
                                        coverage,
                                        last_reverse_strand_string,
                                        start_position,
                                        end_position,
                                        max_left_overhang,
                                        max_right_overhang,
                                        max_left_overhang
                                            + last_end_pos - last_pos)
        output_line_count += 1
        maximin_anchor_size, unique_displacements, coverage = None, set(), 0
        max_left_overhang, max_right_overhang = None, None
        write_line = False
    if not line: break
    coverage += 1
    maximin_anchor_size = max(maximin_anchor_size,
                            min(left_anchor, right_anchor))
    unique_displacements.add(left_overhang)
    max_left_overhang = max(left_overhang, max_left_overhang)
    max_right_overhang = max(right_overhang, max_right_overhang)
    (last_rname, last_sample_label, last_pos, last_end_pos,
        last_reverse_strand_string) = (rname, sample_label, pos, end_pos,
        reverse_strand_string)

print >>sys.stderr, 'DONE with bed_pre.py; in/out =%d/%d; time=%0.3f s' \
                        % (input_line_count, output_line_count,
                            time.time() - start_time)