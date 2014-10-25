#!/usr/bin/env python
"""
Rail-RNA-coverage_pre

Follows Rail-RNA-collapse
Precedes Rail-RNA-coverage

Reduce step in MapReduce pipelines that takes count differentials (exon_diff)
from Hadoop output of Rail-RNA-align and compiles per-sample coverage 
information.

Input (read from stdin)
----------------------------
Tab-delimited input tuple columns:
1. Reference name (RNAME in SAM format) + ';' + bin number
2. Sample label
3. Position at which diff should be subtracted or added to coverage
4. A diff -- that is, by how much coverage increases or decreases at the 
    given genomic position: +n or -n for some natural number n.
Input is partitioned by sample genome partition (fields 1-2) and sorted by
position (field 3).

Hadoop output (written to stdout)
----------------------------
Coverage (coverage)

Tab-delimited output tuple columns:
1. Sample label
2. Number string representing reference name (RNAME in SAM format; see 
    BowtieIndexReference class in bowtie_index for conversion information)
3. Position
4. Coverage (that is, the number of called ECs in the sample
    overlapping the position)

Partition statistics (partition_stats)
1. Number of exon_diff lines contributing to partition
2. Time it took to compute coverage for partition

Reducer statistics (reducer_stats); only one line per reducer
1. Number of partitions processed by reducer
2. Number of input lines processed by reducer
3. Number of output lines processed by reducer
4. Time taken by reducer

IF INPUT COORDINATES ARE (0-) 1-INDEXED, OUTPUT COORDINATES ARE (0-) 1-INDEXED.
"""
import os
import sys
import site
import argparse
import time
import itertools

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

parser = argparse.ArgumentParser(description=__doc__, 
            formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument(\
    '--partition-stats', action='store_const', const=True, 
    help='Output statistics about bin sizes, time taken per bin, number of '
         'bins per reducer')
bowtie.add_args(parser)

args = parser.parse_args()

start_time = time.time()
input_line_count, output_line_count = 0, 0
bin_count = 0
# For converting RNAMEs to number strings
reference_index = bowtie_index.BowtieIndexReference(args.bowtie_idx)
for (partition_id, sample_label), xpartition in xstream(sys.stdin, 2):
    bin_count += 1
    bin_start_time, bin_diff_count = time.time(), 0
    rname = reference_index.rname_to_string[
                    partition_id.rpartition(';')[0]
                ]
    coverage = 0
    for pos, diffs in itertools.groupby(xpartition, lambda val: val[0]):
        input_line_count += 1
        pos = int(pos)
        for _, diff in diffs:
            coverage += int(diff)
            bin_diff_count += 1
        print 'coverage\t%s\t%s\t%012d\t%d' % (sample_label, 
                rname, pos, coverage)
        output_line_count += 1
    if args.partition_stats:
        print 'partition_stats\t%d\t%s' % (bin_diff_count, 
                                            time.time() - bin_start_time)

end_time = time.time()
if args.partition_stats:
    print 'reducer_stats\t%d\t%d\t%d\t%d' % (bin_count, input_line_count,
        output_line_count, end_time - start_time)

print >>sys.stderr, "DONE with coverage_pre.py; in/out = %d/%d; time=%0.3f s" \
    % (input_line_count, output_line_count, start_time - end_time)
