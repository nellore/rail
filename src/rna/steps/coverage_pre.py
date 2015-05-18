#!/usr/bin/env python
"""
Rail-RNA-coverage_pre

Follows Rail-RNA-collapse
Precedes Rail-RNA-coverage

Reduce step in MapReduce pipelines that takes count differentials (exon_diff)
from Hadoop output of align steps and compiles per-sample coverage 
information, both for "uniquely mapping" reads only and for all primary
alignments.

Input (read from stdin)
----------------------------
Tab-delimited input tuple columns:
1. Reference name (RNAME in SAM format) + ';' + bin number
2. Position at which diff should be subtracted or added to coverage
3. Sample index
4. '1' if alignment from which diff originates is "unique" according to
    --tie-margin criterion; else '0'
5. A diff -- that is, by how much coverage increases or decreases at the 
    given genomic position: +n or -n for some natural number n.
Input is partitioned by genome partition (field 1) and sorted by position, then
sample index (fields 2-3).

Hadoop output (written to stdout)
----------------------------
Coverage (coverage)

Tab-delimited output tuple columns:
1. Sample index OR mean[.RNAME] OR median[.RNAME]
2. Number string representing reference name (RNAME in SAM format; see 
    BowtieIndexReference class in bowtie_index for conversion information)
3. Position
4. Coverage counting all primary alignments (that is, the number of called ECs
    in the sample overlapping the position) (mean or median of sample coverages
    normalized by --library-size if field 1 specifies)
5. Coverage counting only "uniquely mapping" reads; here, unique mappings are
    defined according to the criteria implied in --tie-margin (mean or median
    of sample coverages normalized by --library-size if field 1 specifies)

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
import manifest
from dooplicity.tools import xstream, xopen
from collections import defaultdict

parser = argparse.ArgumentParser(description=__doc__, 
            formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument(
    '--partition-stats', action='store_const', const=True, 
    help=('Output statistics about bin sizes, time taken per bin, number of '
          'bins per reducer'))
parser.add_argument(
    '--library-size', type=int, required=False,
    default=40,
    help=('Library size (in millions of reads) to which every sample\'s '
          'coverage should be normalized when computing average coverage'))
parser.add_argument(
        '--read-counts', type=str, required=True,
        help=('File with read counts by sample output by '
              'collect_read_stats.py')
    )
parser.add_argument(
        '--output-ave-bigwig-by-chr', action='store_const', const=True,
        default=False,
        help='Divides bigwigs storing average coverages up by chromosome'
    )
bowtie.add_args(parser)
manifest.add_args(parser)
args = parser.parse_args()

def median(a_list):
    """ Finds canonically defined median of list.

        a_list: a list

        Return value: median
    """
    if not a_list:
        # Median's nothing for empty input in this case
        return 0
    sorted_list = sorted(a_list)
    list_size = len(a_list)
    index = (list_size - 1) // 2
    if (list_size % 2):
        return sorted_list[index]
    return (sorted_list[index] + sorted_list[index + 1]) / 2.0

library_size = args.library_size * 1000000
start_time = time.time()
input_line_count, output_line_count = 0, 0
bin_count = 0
# For converting RNAMEs to number strings
reference_index = bowtie_index.BowtieIndexReference(args.bowtie_idx)
manifest_object = manifest.LabelsAndIndices(args.manifest)
# Grab read counts
mapped_read_counts, unique_mapped_read_counts = {}, {}
with xopen(None, args.read_counts) as read_count_stream:
    read_count_stream.readline()
    for line in read_count_stream:
        tokens = line.strip().split('\t')
        sample_index = manifest_object.label_to_index[tokens[0]]
        (mapped_read_counts[sample_index],
            unique_mapped_read_counts[sample_index]) = [
                                        int(token) for token
                                        in tokens[-2].split(',')
                                    ]
try:
    mean_weight = 1. / len([_ for _ in mapped_read_counts.values() if _])
except ZeroDivisionError:
    mean_weight = 0.0
try:
    unique_mean_weight = 1. / len(
                        [_ for _ in unique_mapped_read_counts.values() if _]
                    )
except ZeroDivisionError:
    unique_mean_weight = 0.0

for (partition_id,), xpartition in xstream(sys.stdin, 1):
    bin_count += 1
    bin_start_time, bin_diff_count = time.time(), 0
    rname_index = reference_index.rname_to_string[
                    partition_id.rpartition(';')[0]
                ]
    for (pos, sample_indexes_and_diffs) in itertools.groupby(
                                            xpartition, lambda val: val[0]
                                        ):
        input_line_count += 1
        pos = int(pos)
        coverages, unique_coverages = defaultdict(int), defaultdict(int)
        for sample_index, diffs in itertools.groupby(
                                sample_indexes_and_diffs, lambda val: val[1]
                            ):
            for _, _, uniqueness, diff in diffs:
                coverages[sample_index] += int(diff)
                if uniqueness == '1':
                    unique_coverages[sample_index] += int(diff)
                bin_diff_count += 1
            print 'coverage\t%s\t%s\t%012d\t%d\t%d' % (sample_index, 
                    rname_index, pos, coverages[sample_index],
                        unique_coverages[sample_index])
            output_line_count += 1
        # Now output measures of center
        coverage_row = [float(coverages[sample_index])
                         / mapped_read_counts[sample_index]
                         * library_size
                        for sample_index in manifest_object.index_to_label
                        if mapped_read_counts[sample_index]]
        unique_coverage_row = [float(unique_coverages[sample_index])
                                / unique_mapped_read_counts[sample_index]
                                * library_size
                               for sample_index
                               in manifest_object.index_to_label
                               if unique_mapped_read_counts[sample_index]]
        print 'coverage\t%s\t%s\t%012d\t%08f\t%08f' % (
                'mean' + (('.' + rname)
                            if args.output_ave_bigwig_by_chr else ''), 
                rname_index, pos,
                sum([cov * mean_weight for cov in coverage_row]),
                sum([cov * unique_mean_weight for cov in unique_coverage_row])
            )
        print 'coverage\t%s\t%s\t%012d\t%08f\t%08f' % (
                'median' + (('.' + rname)
                            if args.output_ave_bigwig_by_chr else ''), 
                rname_index, pos,
                median(coverage_row),
                median(unique_coverage_row)
            )

    if args.partition_stats:
        print 'partition_stats\t%d\t%s' % (bin_diff_count, 
                                            time.time() - bin_start_time)

end_time = time.time()
if args.partition_stats:
    print 'reducer_stats\t%d\t%d\t%d\t%d' % (bin_count, input_line_count,
        output_line_count, end_time - start_time)

print >>sys.stderr, "DONE with coverage_pre.py; in/out = %d/%d; time=%0.3f s" \
    % (input_line_count, output_line_count, start_time - end_time)
