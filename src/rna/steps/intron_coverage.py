#!/usr/bin/env python
"""
Rail-RNA-intron_coverage
Follows Rail-RNA-realign_reads
Precedes Rail-RNA-break_ties

Associates each alignment that overlaps introns with their coverage by reads
aligning to those introns uniquely. Here, "uniquely" is defined by the
--tie-margin parameter in rail-rna_config: if more than one alignment of a read
falls within --tie-margin of the maximum alignment score, the read does not 
align uniquely.

Input (read from stdin)
----------------------------
Two formats -- format 1's tab-delimited input columns (introns):
Format 1's tab-deliminted input columns (introns):
1. The character 'N'
2. Sample index
3. Number string representing RNAME
4. Start position (Last base before insertion, first base of
                    deletion, or first base of intron)
5. End position (Last base before insertion, last base of deletion
    (exclusive), or last base of intron (exclusive))
6. '+' or '-' indicating which strand is the sense strand
7. Number of nucleotides between 5' end of intron and 5' end of
    read from which it was inferred, ASSUMING THE SENSE STRAND IS
    THE FORWARD STRAND. That is, if the sense strand is the reverse
    strand, this is the distance between the 3' end of the read and
    the 3' end of the intron.
8. Number of nucleotides between 3' end of intron and 3' end of
    read from which it was inferred, ASSUMING THE SENSE STRAND IS
    THE FORWARD STRAND.
9. Number of instances of intron, insertion, or deletion in sample;
    this is always +1 before bed_pre combiner/reducer

Format 2's tab-delimited input columns (alignments)
1. The character 'N' so the line can be matched up with intron bed lines
2. Sample index
3. Number string representing RNAME; see BowtieIndexReference class
    in bowtie_index for conversion information
4. Intron start position
5. Intron end position
6. '+' or '-' indicating which strand is sense strand
7. '-' TO ENSURE THAT THE LINE FOLLOWS ALL INTRON LINES
8. POS
9. QNAME
10. FLAG
11. MAPQ
12. CIGAR
13. RNEXT
14. PNEXT
15. TLEN
16. SEQ
17. QUAL
... + optional fields

Input is partitioned by fields 1-6 and sorted by field 7 to ensure that
all alignment lines follow coverage lines

Hadoop output (written to stdout)
----------------------------
Tab-delimited tuple columns:
[ALL SAM FIELDS; see SAM format specification for details]
Last field: XC:i:(coverage of an intron from the read); there are as many lines
for a given alignment (i.e., QNAMEs) as there are introns overlapped by the
alignment

ALL COORDINATES ARE 1-BASED.
"""
import sys
import os
import site
import time
import argparse

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
from dooplicity.tools import xstream

# Print file's docstring if -h is invoked
parser = argparse.ArgumentParser(description=__doc__, 
            formatter_class=argparse.RawDescriptionHelpFormatter)
bowtie.add_args(parser)
args = parser.parse_args()
input_line_count, output_line_count = 0, 0

start_time = time.time()

reference_index = bowtie_index.BowtieIndexReference(args.bowtie_idx)
for (_, sample_index, rname_string, intron_pos, intron_end_pos,
        sense), xpartition in xstream(sys.stdin, 6):
    coverage = 0
    for value in xpartition:
        input_line_count += 1
        try:
            # Assume intron line
            _, _, instance_count = value
        except ValueError:
            # Alignment line
            print ('\t'.join((value[2], value[3],
                reference_index.string_to_rname[rname_string],
                str(int(value[1])))
                + value[4:]) + ('\tXC:i:%d' % coverage))
            output_line_count += 1
        else:
            coverage += int(instance_count)

print >>sys.stderr, 'DONE with intron_coverage.py; in/out=%d/%d; ' \
                    'time=%0.3f s' % (input_line_count, output_line_count,
                                        time.time() - start_time)