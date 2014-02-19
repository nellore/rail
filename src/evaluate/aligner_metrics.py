#!/usr/bin/env python
"""
aligner_metrics.py

Outputs precision, recall, and related performance measurements of a spliced
aligner given a BED file T with true junctions (typically from a simulator like
Flux) and a BED file Y with junctions retrieved by the aligner.

A line in an input BED file characterizes an intron by decking it with blocks.
A single intron defines two junctions, one on either side of it. A true
positive is a junction recovered by the aligner, a false positive is a junction
that appears in Y but not in T, and a false negative is a junction that appears
in T but not in Y.

Output (written to stdout)
----------------------------
Two columns delimited by tabs, where the first column characterizes the numbers
in the second column. The first column is given below.

true junction count
retrieved junction count
true positive count
false positive count
false negative count
precision
recall
"""

import sys

def junctions_from_bed_stream(bed_stream):
    """ Converts BED to dictionary that maps RNAMES to sets of junction POSes.

        bed_stream: input stream containing lines of a BED file characterizing
            splice junctions.

        Return value: a dictionary. Each key is an RNAME, typically a
            chromosome, and its corresponding value is a set of junction
            positions on the RNAME.
    """
    junctions = {}
    for line in bed_stream:
        tokens = line.rstrip().split('\t')
        if len(tokens) < 12:
            continue
        chrom = tokens[0]
        chrom_start = int(tokens[1])
        chrom_end = int(tokens[2])
        if chrom not in junctions:
            junctions[chrom] = set()
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
        # First block characterizes junction on left side of intron
        junctions[chrom].add(chrom_start + int(block_starts[0]) 
                                + int(block_sizes[0]))
        for i in xrange(1, block_count - 1):
            # Any intervening blocks characterize two junctions
            intron_start = chrom_start + int(block_starts[i])
            junctions[chrom].add(intron_start)
            junctions[chrom].add(intron_start + int(block_sizes[i]))
        # Final block characterizes junction on right side of intron
        junctions[chrom].add(chrom_start + int(block_starts[-1]))
    return junctions

if __name__ == '__main__':
    import argparse
    # Print file's docstring if -h is invoked
    parser = argparse.ArgumentParser(description=__doc__, 
            formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-t', '--true-junctions-bed', type=str, required=True, 
        help='Full path of BED file containing true junctions')
    parser.add_argument('-r', '--retrieved-junctions-bed', type=str,
        required=True, 
        help='Full path of BED file containing junctions retrieved by aligner')
    args = parser.parse_args(sys.argv[1:])
    with open(args.true_junctions_bed) as true_junctions_bed_stream:
        true_junctions = junctions_from_bed_stream(true_junctions_bed_stream)
    with open(args.retrieved_junctions_bed) as retrieved_junctions_bed_stream:
        retrieved_junctions \
            = junctions_from_bed_stream(retrieved_junctions_bed_stream)
    (false_positive_count, false_negative_count, true_junction_count, 
        retrieved_junction_count) = [0]*4
    # Use set differences to compute false negative and false positive counts
    for chrom in true_junctions:
        true_junction_count += len(true_junctions[chrom])
        false_negative_count += \
            len(true_junctions[chrom] - retrieved_junctions.get(chrom, set()))
    for chrom in retrieved_junctions:
        retrieved_junction_count += len(retrieved_junctions[chrom])
        false_positive_count += \
            len(retrieved_junctions[chrom] - true_junctions.get(chrom, set()))
    true_positive_count = retrieved_junction_count - false_positive_count
    precision = float(true_positive_count) / retrieved_junction_count
    recall = float(true_positive_count) / true_junction_count
    print 'true junction count\t%d' % true_junction_count
    print 'retrieved junction count\t%d' % retrieved_junction_count
    print 'true positive count\t%d' % true_positive_count
    print 'false positive count\t%d' % false_positive_count
    print 'false negative count\t%d' % false_negative_count
    print 'precision\t%.9f' % precision
    print 'recall\t%.9f' % recall