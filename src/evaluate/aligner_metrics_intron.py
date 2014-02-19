#!/usr/bin/env python
"""
aligner_metrics_intron.py

Outputs precision, recall, and related performance measurements of a spliced
aligner given a BED file T with true introns (typically from a simulator like
Flux) and a BED file Y with introns retrieved by the aligner.

A line in an input BED file characterizes an intron by decking it with blocks.
A single intron defines two junctions, one on either side of it. A true
positive is an intron recovered by the aligner, a false positive is an intron
that appears in Y but not in T, and a false negative is an intron that appears
in T but not in Y.

Output (written to stdout)
----------------------------
Two columns delimited by tabs, where the first column characterizes the numbers
in the second column. The first column is given below.

true intron count
retrieved intron count
true positive count
false positive count
false negative count
precision
recall
"""

import sys

def introns_from_bed_stream(bed_stream):
    """ Converts BED to dictionary that maps RNAMES to sets of introns.

        bed_stream: input stream containing lines of a BED file characterizing
            splice junctions.

        Return value: a dictionary. Each key is an RNAME, typically a
            chromosome, and its corresponding value is a set of tuples, each
            denoting an intron on RNAME. Each tuple is of the form
            (start position, end position).
    """
    introns = {}
    for line in bed_stream:
        tokens = line.rstrip().split('\t')
        if len(tokens) != 12:
            continue
        chrom = tokens[0]
        chrom_start = int(tokens[1])
        chrom_end = int(tokens[2])
        if chrom not in introns:
            introns[chrom] = set()
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
        if block_count < 2:
            # No introns
            continue
        assert block_count == len(block_starts)
        junctions = []
        # First block characterizes junction on left side of intron
        junctions.append(chrom_start + int(block_starts[0]) 
                                + int(block_sizes[0]))
        for i in xrange(1, block_count - 1):
            # Any intervening blocks characterize two junctions
            intron_start = chrom_start + int(block_starts[i])
            junctions.append(intron_start)
            junctions.append(intron_start + int(block_sizes[i]))
        # Final block characterizes junction on right side of intron
        junctions.append(chrom_start + int(block_starts[-1]))
        for i in xrange(len(junctions) - 1):
            introns[chrom].add((junctions[i], junctions[i+1]))
    return introns

if __name__ == '__main__':
    import argparse
    # Print file's docstring if -h is invoked
    parser = argparse.ArgumentParser(description=__doc__, 
            formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-t', '--true-introns-bed', type=str, required=True, 
        help='Full path of BED file containing true introns')
    parser.add_argument('-r', '--retrieved-introns-bed', type=str,
        required=True, 
        help='Full path of BED file containing introns retrieved by aligner')
    args = parser.parse_args(sys.argv[1:])
    with open(args.true_introns_bed) as true_introns_bed_stream:
        true_introns = introns_from_bed_stream(true_introns_bed_stream)
    with open(args.retrieved_introns_bed) as retrieved_introns_bed_stream:
        retrieved_introns \
            = introns_from_bed_stream(retrieved_introns_bed_stream)
    (false_positive_count, false_negative_count, true_intron_count, 
        retrieved_intron_count) = [0]*4
    # Use set differences to compute false negative and false positive counts
    for chrom in true_introns:
        true_intron_count += len(true_introns[chrom])
        false_negative_count += \
            len(true_introns[chrom] - retrieved_introns.get(chrom, set()))
    for chrom in retrieved_introns:
        retrieved_intron_count += len(retrieved_introns[chrom])
        false_positive_count += \
            len(retrieved_introns[chrom] - true_introns.get(chrom, set()))
    true_positive_count = retrieved_intron_count - false_positive_count
    precision = float(true_positive_count) / retrieved_intron_count
    recall = float(true_positive_count) / true_intron_count
    print 'true intron count\t%d' % true_intron_count
    print 'retrieved intron count\t%d' % retrieved_intron_count
    print 'true positive count\t%d' % true_positive_count
    print 'false positive count\t%d' % false_positive_count
    print 'false negative count\t%d' % false_negative_count
    print 'precision\t%.9f' % precision
    print 'recall\t%.9f' % recall