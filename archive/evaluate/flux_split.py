"""
flux_split.py

This file divides the BED and FASTQ output of a Flux simulation into BED and
FASTQs corresponding to a specified number of technical replicates. Under
num-replicates reads may be dropped to arrange that the number of reads per
replicate is equal.
"""

import argparse
import sys
import random
import itertools

# Print file's docstring if -h is invoked
parser = argparse.ArgumentParser(description=__doc__, 
            formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument(\
        '-b', '--bed', type=str, required=True,
        help='Path of Flux BED file to split'
    )
parser.add_argument(\
        '-f', '--fastq', type=str, required=True,
        help='Path of Flux FASTQ file to split'
    )
parser.add_argument(\
        '-n', '--num-replicates', type=int, required=True,
        help='Number of replicates'
    )
parser.add_argument(\
        '-s', '--random-seed', type=int, required=False, default=0,
        help='Random seed'
    )
parser.add_argument(\
        '-p', '--paired-end', action='store_const', const=True,
        default=False,
        help='Assume input reads are paired-end'
    )

args = parser.parse_args()

if args.num_replicates < 2:
    print >>sys.stderr, 'Number of replicates should be an integer >= 2; ' \
                        'you entered %d.' % args.num_replicates

random.seed(args.random_seed)

fastq_output_streams = [open(args.fastq + '.%d' % i, 'w')
                            for i in xrange(args.num_replicates)]
bed_output_streams = [open(args.bed + '.%d' % i, 'w')
                            for i in xrange(args.num_replicates)]

fastq_read_line_count = 4
bed_read_line_count = 1
if args.paired_end:
    fastq_read_line_count *= 2
    bed_read_line_count *= 2

i = 0
with open(args.bed) as bed_stream:
    with open(args.fastq) as fastq_stream:
        bed_lines = []
        fastq_lines = []
        while True:
            bed_lines.append(
                    list(itertools.islice(bed_stream, bed_read_line_count))
                )
            fastq_lines.append(
                    list(itertools.islice(fastq_stream, fastq_read_line_count))
                )
            i += 1
            if i == 1: continue
            if i % args.num_replicates == 0:
                line_order = range(args.num_replicates)
                random.shuffle(line_order)
                for k, l in enumerate(line_order):
                    for m in range(bed_read_line_count):
                        bed_output_streams[k].write(
                                bed_lines[l][m]
                            )
                    for m in range(fastq_read_line_count):
                        fastq_output_streams[k].write(
                                fastq_lines[l][m]
                            )
                bed_lines = []
                fastq_lines = []