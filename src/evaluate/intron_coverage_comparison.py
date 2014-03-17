#!/usr/bin/env python
"""
intron_coverage_comparison.py

Facilitates comparison of true coverage of introns from a Flux BED with
coverages from spliced aligners that output TopHat-like junctions BEDs.

A line in an input BED file characterizes an intron by decking it with blocks.
A single intron defines two junctions, one on either side of it. A true
positive is a junction recovered by the aligner, a false positive is a junction
that appears in Y but not in T, and a false negative is a junction that appears
in T but not in Y. The same intron occurs on multiple lines of the Flux BED,
and coverage of the intron is obtained by summing the number of intron calls.
Coverages are taken from the score column of TopHat-like junctions BEDs.

Coverage is 0 if an intron was not recovered by the aligner.

Output (written to stdout)
----------------------------
Columns delimited by tabs:
1. True intron start
2. True intron end
3. True coverage
4. Coverage of aligner #1
5. Coverage of aligner #2
           .
           .
           .
... (for as many aligner BEDs specified at command line)
"""

import sys

def introns_from_bed_stream(bed_stream, flux=True):
    """ Converts BED to dictionary that maps introns to coverages.

        The same intron occurs on multiple lines of the Flux BED, and coverage
        of the intron is obtained by summing the number of intron calls when
        flux=True. Coverages are taken from the score column of TopHat-like
        junctions BEDs, for which flux should be False.

        bed_stream: input stream containing lines of a BED file characterizing
            splice junctions.

        Return value: a dictionary. Each key is an RNAME whose value is itself
            a dictionary with key the tuple (intron start position, intron end
            position) and corresponding value coverage.
    """
    coverages = {}
    for line in bed_stream:
        tokens = line.rstrip().split('\t')
        if len(tokens) < 12:
            continue
        chrom = tokens[0]
        chrom_start = int(tokens[1])
        chrom_end = int(tokens[2])
        if chrom not in coverages:
            coverages[chrom] = {}
        block_sizes = tokens[10].split(',')
        block_starts = tokens[11].split(',')
        coverage = tokens[4]
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
        for i in xrange(len(junctions) / 2):
            if flux:
                coverages[chrom][(junctions[2*i], junctions[2*i+1])] = \
                    coverages[chrom].get(
                                     (junctions[2*i], junctions[2*i+1]
                                    ), 0) + 1
            else:
                coverages[chrom][(junctions[2*i], junctions[2*i+1])] = \
                    int(coverage)
    return coverages

if __name__ == '__main__':
    import argparse
    # Print file's docstring if -h is invoked
    parser = argparse.ArgumentParser(description=__doc__, 
            formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-t', '--true-introns-bed', type=str, required=True, 
        help='Full path of Flux BED file containing true introns')
    parser.add_argument('-r', '--retrieved-intron-beds', type=str,
        required=True, 
        help='Comma-separated list of TopHat-like junctions BED files '
             'containing introns retrieved by aligners')
    args = parser.parse_args(sys.argv[1:])
    with open(args.true_introns_bed) as true_introns_bed_stream:
        print >>sys.stderr, \
            'Reading true introns BED %s....' % args.true_introns_bed
        true_introns = introns_from_bed_stream(true_introns_bed_stream)
    retrieved_bed_files = [filename.strip() 
                            for filename in 
                            args.retrieved_intron_beds.split(',')]
    retrieved_introns = []
    for filename in retrieved_bed_files:
        with open(filename) as retrieved_introns_bed_stream:
            print >>sys.stderr, \
                'Reading retrieved introns BED %s....' % filename
            retrieved_introns.append(
                    introns_from_bed_stream(retrieved_introns_bed_stream,
                                                flux=False)
                )
    for chrom in true_introns:
        for true_intron in true_introns[chrom]:
            sys.stdout.write('%s\t%d\t%d\t%d' % ((chrom,) + true_intron + 
                                (true_introns[chrom][true_intron],)))
            for retrieved_intron_set in retrieved_introns:
                if chrom not in retrieved_intron_set \
                    or true_intron not in retrieved_intron_set[chrom]:
                    sys.stdout.write('\t0')
                else:
                    sys.stdout.write('\t%d' 
                            % retrieved_intron_set[chrom][true_intron]
                        )
            sys.stdout.write('\n')
    sys.stdout.flush()