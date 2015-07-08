#!/usr/bin/env python
"""
unique_introns.py

Outputs introns unique to each bed file from a set on bed files supplied as
a space-separated list at the command line.

Introns are output as new files with the same prefixes as the old ones,
but with the suffix ".unique" rather than ".bed". ".unique" files are placed
in the same directory as the original bed files.

Output tuple columns per file
chrom <tab> intron start (1-based inclusive)
    <tab> intron end (1-based exclusive)
"""
import os
from collections import defaultdict

if __name__ == '__main__':
    # Print file's docstring if -h is invoked
    parser = argparse.ArgumentParser(description=__doc__, 
            formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-b', '--bed', type=str, required=True, nargs='+'
        help='Path to BED files whose exonic blocks deck introns')
    args = parser.parse_args()
    introns_to_sample_indexes = defaultdict(set)
    sample_indexes_to_filenames = {}
    for k, bed_file in enumerate(args.bed):
        bed_dir = os.path.dirname(bed_file)
        sample_indexes_to_filenames[i] = os.path.join(
                                            bed_dir,
                                            bed_file[:-4] + '.unique'
                                        )
        with open(bed_file) as bed_stream:
            for line in bed_stream:
                tokens = line.rstrip().split('\t')
                if len(tokens) < 12:
                    continue
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
                    introns_to_sample_indexes[
                            (chrom, junctions[2*i]+1, junctions[2*i+1]+1)
                        ].add(k)
    open_files = {}
    for intron in introns_to_sample_indexes:
        if len(introns_to_sample_indexes[intron]) == 1:
            sample_index = list(introns_to_sample_indexes[intron])[0]
            try:
                print >>open_files[sample_index], '\t'.join(intron)
            except KeyError:
                open_files[sample_index] = open(
                            sample_indexes_to_filenames[sample_index], 'w'
                        )
                print >>open_files[sample_index], '\t'.join(intron)
