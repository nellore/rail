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

args = parser.parse_args()

if args.num_replicates < 2:
	print >>sys.stderr, 'Number of replicates should be an integer >= 2; ' \
						'you entered %d.' % args.num_replicates

random.seed(args.random_seed)

fastq_output_streams = [open(args.fastq + '.%d' % i, 'w')
							for i in xrange(args.num_replicates)]
bed_output_streams = [open(args.bed + '.%d' % i, 'w')
							for i in xrange(args.num_replicates)]

with open(args.bed) as bed_stream:
	with open(args.fastq) as fastq_stream:
		bed_lines = []
		fastq_grouped_lines = []
		for i, line in enumerate(bed_stream):
			bed_lines.append(line)
			fastq_lines = []
			for j in range(4):
				fastq_lines.append(fastq_stream.readline())
			fastq_grouped_lines.append(fastq_lines)
			if i == 0: continue
			if (i + 1) % args.num_replicates == 0:
				line_order = range(args.num_replicates)
				random.shuffle(line_order)
				for k, l in enumerate(line_order):
					bed_output_streams[k].write(bed_lines[l])
					for m in range(4):
						fastq_output_streams[k].write(
								fastq_grouped_lines[l][m]
							)