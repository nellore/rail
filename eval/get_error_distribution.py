#!/usr/bin/env python
"""
get_error_distribution.py

Gets error distribution of Flux BED/FASTQ. Considers only reads from the FASTQ
with the length specified at the command line.

We executed

pypy get_error_distribution.py --bed <(cat *.bed) --fastq <(cat *_sim.fastq)
     -l 76 -x /dcl01/leek/data/railsims/indexes_for_paper/genome
    >error_distribution

from the directory containing the 112 Flux BEDs/FASTQs to obtain the overall
error distribution and

pypy get_error_distribution.py --bed <(cat *.bed) --fastq <(cat *_sim.fastq)
     -l 76 -x /dcl01/leek/data/railsims/indexes_for_paper/genome
     --ignore-b-tails >error_distribution_no_b_tails
     
"""

from count_introns import BowtieIndexReference
import string
_reversed_complement_translation_table = string.maketrans('ATCG', 'TAGC')

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
    parser.add_argument('--ignore-b-tails', action='store_const', const=True,
        default=False,
        help='Don\'t count reads for which the last 10 bases have Q <= 15'
    )
    parser.add_argument('-x', '--bowtie-idx', required=True,
        help='Path to Bowtie 1 index. Specify its basename.'
    )
    args = parser.parse_args()
    reference_index = BowtieIndexReference(args.bowtie_idx)
    truths = [0 for i in xrange(args.length)]
    recalls = [0 for i in xrange(args.length)]
    k = 0
    b_tail_count, read_count = 0, 0
    with open(args.bed) as bed_stream:
        with open(args.fastq) as fastq_stream:
            while True:
                bed_line = bed_stream.readline().strip()
                if not bed_line: break
                read_name = fastq_stream.readline().strip()[1:]
                seq = fastq_stream.readline().strip().upper()
                fastq_stream.readline()
                qual = fastq_stream.readline().strip()
                read_count += 1
                if args.ignore_b_tails and all(
                            [ord(score) - 33 <= 15 for score in qual[-10:]]
                        ):
                    b_tail_count += 1
                    continue
                if len(seq) != args.length: continue
                tokens = bed_line.rstrip().split('\t')
                assert len(tokens) >= 12
                assert read_name == tokens[3], k
                chrom = tokens[0]
                chrom_start = int(tokens[1])
                chrom_end = int(tokens[2])
                block_sizes = tokens[10].split(',')
                block_starts = tokens[11].split(',')
                # Handle trailing commas
                try:
                    int(block_sizes[-1])
                except ValueError:
                    block_sizes = block_sizes[:-1]
                try:
                    int(block_starts[-1])
                except ValueError:
                    block_starts = block_starts[:-1]
                block_count = len(block_sizes)
                assert block_count == len(block_starts)
                true_seq = ''.join([reference_index.get_stretch(
                                        chrom,
                                        chrom_start + int(block_starts[i]),
                                        int(block_sizes[i])
                                    ) for i in xrange(block_count)]).upper()
                if tokens[5] == '-':
                    true_seq = true_seq[::-1].translate(
                            _reversed_complement_translation_table
                        )
                for i in xrange(args.length):
                    truths[i] += 1
                    if true_seq[i] == seq[i]:
                        recalls[i] += 1
                k += 1
    overall_truths = sum(truths)
    overall_recalls = sum(recalls)
    print 'Read count: %d' % read_count
    print 'Overall error rate: %.12f' % (
                float(overall_truths - overall_recalls) / overall_truths
            )
    distribution = [float(truths[i] - recalls[i]) / truths[i]
                        for i in xrange(args.length)]
    print 'Error rate distribution:'
    print '\t'.join([str(val) for val in distribution])
    if args.ignore_b_tails:
        print 'B-tail count: %d' % b_tail_count
        print 'B-tail proportion: %.12f' % (float(b_tail_count) / read_count)