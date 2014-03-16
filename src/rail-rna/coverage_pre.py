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

base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for directory_name in ['interval', 'bowtie']:
    site.addsitedir(os.path.join(base_path, directory_name))

import partition
import bowtie
import bowtie_index

parser = argparse.ArgumentParser(description=__doc__, 
            formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument(\
    '--partition-stats', action='store_const', const=True, 
    help='Output statistics about bin sizes, time taken per bin, number of '
         'bins per reducer')
bowtie.addArgs(parser)

args = parser.parse_args()

start_time = time.time()
input_line_count, output_line_count = 0, 0
last_partition_id, last_sample_label, last_pos = [None]*3
output_coverage = False
bin_count, bin_diff_count, bin_start_time = 0, 0, start_time
# For converting RNAMEs to number strings
reference_index = bowtie_index.BowtieIndexReference(args.bowtie_idx)
while True:
    line = sys.stdin.readline()
    if line:
        input_line_count += 1
        tokens = line.rstrip().split('\t')
        assert len(tokens) == 4, 'Bad input line:\n' + line
        partition_id, sample_label, pos, differential = (tokens[0],
            tokens[1], int(tokens[2]), int(tokens[3]))
        # Convert RNAME to number string
        rname = reference_index.rname_to_string[partition_id.split(';')[0]]
        if partition_id != last_partition_id \
            or sample_label != last_sample_label:
            '''Reset coverage count for new partitions/sample labels; the next
            partition gets a +1 if an EC continues into it.'''
            coverage_count = 0
            output_coverage = True
        if pos != last_pos:
            output_coverage = True
    if not line or (output_coverage and last_pos is not None):
        print 'coverage\t%s\t%s\t%012d\t%d' \
            % (last_sample_label, last_rname, last_pos, coverage_count)
        output_line_count += 1
        output_coverage = False
    if not line or (partition_id != last_partition_id 
        and last_partition_id is not None):
        if args.partition_stats:
            time_difference = time.time() - bin_start_time
            print 'partition_stats\t%d\t%s' % (bin_diff_count, 
                                                str(time_difference))
            bin_start_time, bin_diff_count = time.time(), 0
        bin_count += 1
    if not line: break
    coverage_count += differential
    bin_diff_count += 1
    (last_partition_id, last_sample_label, last_pos, last_differential,
        last_rname) = (partition_id, sample_label, pos, differential, rname)

end_time = time.time()
if args.partition_stats:
    print 'reducer_stats\t%d\t%d\t%d\t%s' % (bin_count, input_line_count,
        output_line_count, str(end_time - start_time))

print >>sys.stderr, "DONE with coverage_pre.py; in/out = %d/%d; time=%0.3f s" \
    % (input_line_count, output_line_count, start_time - end_time)
