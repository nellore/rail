#!/usr/bin/env python
"""
Rail-RNA-bed_pre
Follows Rail-RNA-realign_reads
Precedes Rail-RNA-bed

Combiner/reducer for MapReduce pipelines that sums coverage of
indels/introns and computes other statistics by sample from Rail-RNA-realign.

Input (read from stdin)
----------------------------
Tab-delimited input tuple columns:
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

Input is partitioned by fields 1-6.

Hadoop output (written to stdout)
----------------------------
Tab-delimited output tuple columns:
1. 'I', 'D', or 'N' for insertion, deletion, or intron line
2. Sample label
3. Number string representing RNAME (+ '+ or -' if intron; same as field 6)
4. Start position (Last base before insertion, first base of deletion,
                    or first base of intron)
5. End position (Last base before insertion, last base of deletion (exclusive),
                    or last base of intron (exclusive))
6. '+' or '-' indicating which strand is the sense strand for introns,
   inserted sequence for insertions, or deleted sequence for deletions
----Next fields are for introns only; they are '\x1c' for indels----
7. COMMA-SEPARATED LIST OF UNIQUE number of nucleotides between 5' end of
    intron and 5' end of read from which it was inferred, ASSUMING THE SENSE
    STRAND IS THE FORWARD STRAND. That is, if the sense strand is the reverse
    strand, this is the distance between the 3' end of the read and the 3' end
    of the intron.
8. COMMA-SEPARATED LIST OF UNIQUE number of nucleotides between 3' end of
    intron and 3' end of read from which it was inferred, ASSUMING THE SENSE
    STRAND IS THE FORWARD STRAND.
9. MIN mumber of nucleotides spanned by EC on the left (that is, towards the 5'
    end of the read) of the intron, ASSUMING THE SENSE STRAND IS THE FORWARD
    STRAND.
10. MIN number of nucleotides spanned by EC on the right (that is, towards the
    3' end of the read) of the intron, ASSUMING THE SENSE STRAND IS THE FORWARD
    STRAND.
--------------------------------------------------------------------
11. SUMMED number of instances of intron, insertion, or deletion in sample

OUTPUT COORDINATES ARE 0-INDEXED.
"""
import os
import sys
import argparse
import site

base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
site.addsitedir(os.path.join(base_path, 'dooplicity'))
import dooplicity as dp

# Print file's docstring if -h is invoked
parser = argparse.ArgumentParser(description=__doc__, 
            formatter_class=argparse.RawDescriptionHelpFormatter)
args = parser.parse_args()

import time
start_time = time.time()
input_line_count, output_line_count = 0, 0
for (line_type, sample_label, rname,
        pos, end_pos, strand_or_seq), xpartition in dp.xstream(sys.stdin, 6):
    if line_type == 'N':
        left_displacement_set = set()
        right_displacement_set = set()
        coverage_sum = 0
        min_left_anchor_size, min_right_anchor_size = None, None
        for (left_displacements, right_displacements,
                left_anchor_size, right_anchor_size, coverage) in xpartition:
            input_line_count += 1
            for displacement in left_displacements.split(','):
                left_displacement_set.add(displacement)
            for displacement in right_displacements.split(','):
                right_displacement_set.add(displacement)
            min_left_anchor_size = int(left_anchor_size) \
                if min_left_anchor_size is None \
                else min(int(left_anchor_size), min_left_anchor_size)
            min_right_anchor_size = int(right_anchor_size) \
                if min_right_anchor_size is None \
                else min(int(right_anchor_size), min_right_anchor_size)
            coverage_sum += int(coverage)
        assert min_left_anchor_size is not None
        assert min_right_anchor_size is not None
        left_displacement_list = ','.join(list(left_displacement_set))
        right_displacement_list = ','.join(list(right_displacement_set))
        print >>sys.stdout, 'N\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d' \
                            % (sample_label, rname, pos, end_pos,
                                strand_or_seq, left_displacement_list,
                                right_displacement_list,
                                min_left_anchor_size,
                                min_right_anchor_size,
                                coverage_sum)
        output_line_count += 1
    else:
        assert line_type in 'ID'
        coverage_sum = 0
        for _, _, _, _, coverage in xpartition:
            input_line_count += 1
            coverage_sum += int(coverage)
        print >>sys.stdout, '%s\t%s\t%s\t%s\t%s\t%s\t' \
                            '\x1c\t\x1c\t\x1c\t\x1c\t%d' \
                            % (line_type, sample_label, rname, pos, end_pos,
                                strand_or_seq, coverage_sum)
        output_line_count += 1

print >>sys.stderr, 'DONE with bed_pre.py; in/out =%d/%d; time=%0.3f s' \
                        % (input_line_count, output_line_count,
                            time.time() - start_time)