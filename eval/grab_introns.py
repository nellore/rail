#!/user/bin/env python
"""
grab_introns.py

Grabs introns from bed. Output format is

RNAME <tab> 1-based start position (inclusive) <tab>
    1-based end position (exclusive)
"""

def introns_from_bed(bed):
    """ Converts BED to dictionary that maps RNAMES to sets of introns.

        bed: input BED file characterizing splice junctions

        Return value: a dictionary. Each key is an RNAME, typically a
            chromosome, and its corresponding value is a set of tuples, each
            denoting an intron on RNAME. Each tuple is of the form
            (start position, end position).
    """
    with open(bed) as bed_stream:
        for line in bed_stream:
            tokens = line.rstrip().split('\t')
            if len(tokens) < 12:
                continue
            chrom = tokens[0]
            chrom_start = int(tokens[1])
            chrom_end = int(tokens[2])
            coverage = int(tokens[4])
            strand = tokens[5]
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
            for i in xrange(len(junctions)/2):
                print '%s\t%d\t%d\t%s\t%d' % (
                            chrom, junctions[2*i]+1, junctions[2*i+1]+1,
                            strand, coverage
                        )

import sys
introns_from_bed(sys.argv[1])