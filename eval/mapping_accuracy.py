#!/usr/bin/env python
"""
mapping_accuracy.py

Outputs precision, recall, and related performance measurements of a spliced
aligner given a Flux-like BED file T with true spliced alignments and a SAM
file Y with the aligner's spliced alignments (read from stdin). Considers all
reads in the measurement. Two kinds of accuracy are studied: 1) basewise: the
proportion of aligned (non-soft-clipped) read bases placed correctly
(precision) and the proportion of all read bases placed correctly (recall);
2) read-level: the proportion of aligned reads for which >= K% of bases are
aligned correctly (precision) and the proportion of all reads for which >= K%
of bases are aligned correctly (recall). K is the argument of the
--base-threshold command-line parameter.

THIS FILE DEPENDS ON DOOPLICITY AND RAIL; don't move it in the Rail repo.

Output (written to stdout)
----------------------------
1. relevant instances <tab> value <tab> value
2. retrieved instances <tab> value <tab> value
3. intersection [of relevant and retrieved instances] <tab> value <tab> value
4. precision <tab> value <tab> value
5. recall <tab> value <tab> value

In each row, the first value is relevant to basewise accuracy, and the second
value is relevant to read-level accuracy.
"""

import sys
import site
import os
import re
from collections import defaultdict

base_path = os.path.dirname(os.path.dirname(
                        os.path.realpath(__file__))
                    )
utils_path = os.path.join(base_path, 'src', 'rna', 'utils')
src_path = os.path.join(base_path, 'src')

site.addsitedir(utils_path)
site.addsitedir(src_path)

from alignment_handlers import indels_introns_and_exons
from dooplicity.tools import xstream

def dummy_md_and_mapped_offsets(cigar):
    """ Creates dummy MD string from CIGAR in case of missing MD.

        cigar: cigar string
        mapped_bases: returns list of offsets of mapped bases from end of read
            rather than MD string if True

        Return value: tuple (dummy MD string, list of offsets of of mapped
                                bases from beginning of read)
    """
    cigar = re.split(r'([MINDS])', cigar)[:-1]
    cigar_index, offset = 0, 0
    max_cigar_index = len(cigar)
    md, mapped = [], []
    while cigar_index != max_cigar_index:
        if cigar[cigar_index] == 0:
            cigar_index += 2
            continue
        if cigar[cigar_index+1] == 'M':
            try:
                if type(md[-1]) is int:
                    md[-1] += base_count
                else:
                    md.append(base_count)
            except IndexError:
                md.append(base_count)
            cigar_index += 2
            mapped.extend(range(offset, base_count))
            offset += base_count
        elif cigar[cigar_index+1] in 'SI':
            offset += int(cigar[cigar_index])
            cigar_index += 2
        elif cigar[cigar_index+1] == 'N':
            cigar_index += 2
        elif cigar[cigar_index+1] == 'D':
            md.extend(['^', 'A'*int(cigar[cigar_index])])
            cigar_index += 2
        else:
            raise RuntimeError(
                        'Accepted CIGAR characters are only in [MINDS].'
                    )
    return (''.join(str(el) for el in md), mapped)

if __name__ == '__main__':
    import argparse
    # Print file's docstring if -h is invoked
    parser = argparse.ArgumentParser(description=__doc__, 
            formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-t', '--true-bed', type=str, required=True, 
        help='Full path of BED file containing true mappings')
    parser.add_argument('-g', '--generous', action='store_const', const=True,
        default=False,
        help='TopHat/STAR/HISAT cut off /1s and /2s from read names, even in '
             'unpaired mode. This loses information. In generous mode, '
             'this script provides extremely tight upper bounds on '
             'precision and recall for TopHat/STAR/HISAT')
    parser.add_argument('-b', '--base-threshold', type=float,
        required=False, default=0.5,
        help=('Proportion of a read\'s bases that must align correctly for '
              'the read to be considered a correct mapping'))
    args = parser.parse_args(sys.argv[1:])
    # Dict mapping read names to alignments (chrom, 1-based start, 1-based end)
    true_maps = {} # Warning: takes up a lot of memory for large samples!
    basewise_relevant, read_relevant = 0, 0
    with open(args.true_bed) as true_bed_stream:
        if args.generous:
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
                assert block_count == len(block_starts)
                exons = [(chrom,
                            chrom_start + block_starts[i],
                            chrom_start + block_starts[i] + block_sizes[i])
                            for i in xrange(block_count)]
                # Cut off /1 or /2
                true_maps[tokens[3][:-2]] = exons
        else:
            for line in true_bed_stream:
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
                assert block_count == len(block_starts)
                exons = [(chrom,
                            chrom_start + block_starts[i],
                            chrom_start + block_starts[i] + block_sizes[i])
                            for i in xrange(block_count)]
                # No cutting off read names
                true_maps[tokens[3]] = exons
    # Initialize counters for computing accuracy metrics
    basewise_retrieved, basewise_intersection = 0, 0, 0
    read_retrieved, read_intersection = 0, 0, 0
    # Read SAM for stdin
    for line in sys.stdin:
        tokens = line.strip().split('\t')
        flag = int(tokens[1])
        if flag & 256 or flag & 4:
            # Secondary alignment or unmapped and thus not retrieved; ignore
            continue 
        true_map = true_maps[tokens[0]]
        read_length = true_map[0][2] - true_map[0][1]
        basewise_retrieved += read_length
        read_retrieved += 1
        if tokens[2] != true_map[0][0]:
            # chr is wrong, but this is still counted as a retrieval above
            continue
        base_counter, base_truths = 0, set()
        # Each tuple in base_truths is (index of base in read, mapped location)
        for block in true_map:
            base_truths.update([(base_counter + i, j)
                                    for i, j in enumerate(xrange(
                                                        block[1], block[2]
                                                    ))])
            base_counter += block[2] - block[1]
        base_predictions = set()
        dummy_md, mapped = dummy_md_and_mapped_offsets(cigar)
        _, _, _, exons = indels_introns_and_exons(cigar, dummy_md, pos, seq)
        mapped_index = 0
        for exon in exons:
            base_predictions.update([(mapped[mapped_index + i], j)
                                      for i, j in enumerate(xrange(
                                                        exon[1], exon[2]
                                                    ))])
            mapped_index += exon[2] - exon[1]
        intersected_base_count = len(base_predictions & base_truths)
        basewise_intersection += intersected_base_count
        if intersected_base_count >= read_length * args.base_threshold:
            read_intersection += 1
    basewise_precision = float(basewise_intersection) / basewise_retrieved
    basewise_recall = float(basewise_intersection) / read_relevant
    read_precision = float(read_intersection) / read_retrieved
    read_recall = float(read_intersection) / basewise_relevant
    print 'relevant instances\t%d\t%d' % (basewise_relevant, read_relevant)
    print 'retrieved instances\t%d\t%d' % (basewise_retrieved, read_retrieved)
    print 'intersection\t%d\t%d' % (basewise_intersection, read_intersection)
    print 'precision\t%.12f\t%.12f' % (basewise_precision, read_precision)
    print 'recall\t%.12f\t%.12f' % (basewise_recall, read_recall)
