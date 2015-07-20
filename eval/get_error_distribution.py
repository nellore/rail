#!/usr/bin/env python
"""
get_error_distribution.py

Gets error distribution of Flux BED/FASTQ. Considers only reads from the FASTQ
with the length specified at the command line.
"""

if __name__ == '__main__':
	import argparse
	# Print file's docstring if -h is invoked
    parser = argparse.ArgumentParser(description=__doc__, 
                formatter_class=argparse.RawDescriptionHelpFormatter)
    # Add command-line arguments
    parser.add_argument('-l', '--length', required=False, type=int,
        default=76,
        help='Length of reads to consider from input FASTQ'
    )
    parser.add_argument('--bed', type=str, required=True,
        help='Path to Flux bed'
    )
    parser.add_argument('--fastq', type=str, required=True,
    	help='Path to Flux FASTQ'
    )
    parser.add_argument('-x', '--bowtie-idx', required=True,
    	help='Path to Bowtie 1 index. Specify its basename.'
    )
    with open(args.bed) as bed_stream:
    	with open(args.fastq) as fastq_stream:
    		while True