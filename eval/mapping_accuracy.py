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

Warning: this script is a memory hog. Eats ~15 GB of RAM on 40 million reads.
Be ready.

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

def dummy_md_and_mapped_offsets(cigar, clip_threshold=1.0):
    """ Creates dummy MD string from CIGAR in case of missing MD.

        cigar: cigar string
        mapped_bases: returns list of offsets of mapped bases from end of read
            rather than MD string if True
        clip_threshold: proportion of a read's bases that must be clipped
            for a read to be considered unmapped

        Return value: tuple (dummy MD string, list of offsets of of mapped
                                bases from beginning of read,
                                True iff read should be considered unmapped,
                                number of clipped bases)
    """
    cigar = re.split(r'([MINDS])', cigar)[:-1]
    cigar_index, offset, read_bases = 0, 0, 0
    max_cigar_index = len(cigar)
    md, mapped = [], []
    clip_count = 0
    while cigar_index != max_cigar_index:
        if cigar[cigar_index] == 0:
            cigar_index += 2
            continue
        if cigar[cigar_index+1] == 'M':
            base_count = int(cigar[cigar_index])
            read_bases += base_count
            try:
                if type(md[-1]) is int:
                    md[-1] += base_count
                else:
                    md.append(base_count)
            except IndexError:
                md.append(base_count)
            cigar_index += 2
            mapped.extend(range(offset, offset + base_count))
            offset += base_count
        elif cigar[cigar_index+1] == 'S':
            value = int(cigar[cigar_index])
            read_bases += value
            offset += value
            clip_count += value
            cigar_index += 2
        elif cigar[cigar_index+1] == 'I':
            value = int(cigar[cigar_index])
            offset += value
            read_bases += value
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
    return (''.join(str(el) for el in md), mapped,
                float(clip_count) / read_bases >= clip_threshold,
                clip_count)

def go(true_bed_stream, sam_stream=sys.stdin, generous=False,
        base_threshold=0.5, clip_threshold=1.0, dump_incorrect=False):
    """ Finds relevant and retrieved instance counts.

        true_bed_stream: file handle for BED output of Flux simulation
        sam_stream: where to read in aligner's mappings
        generous: True iff aligner cuts off /1 or /2 of a given read
        base_threshold: proportion of a read's bases that must align
            correctly for a read to be considered a correct mapping
        clip_threshold: proportion of a read's bases that must be clipped
            for a read to be considered unmapped
        dump_incorrect: write incorrect (read) alignments to stderr
    """
    from tempdel import remove_temporary_directories
    import tempfile
    import atexit
    temp_dir_path = tempfile.mkdtemp()
    #print >>sys.stderr, temp_dir_path
    atexit.register(remove_temporary_directories, [temp_dir_path])
    # Store everything in one file, then sort it on read name
    combined_file = os.path.join(temp_dir_path, 'combined.temp')
    with open(combined_file, 'w') as temp_stream:
        if generous:
            for line in true_bed_stream:
                tokens = line.strip().split('\t')
                print >>temp_stream, '\t'.join([tokens[3][:-2], '0']
                                                + tokens[:3] + tokens[4:])
        else:
            for line in true_bed_stream:
                tokens = line.strip().split('\t')
                print >>temp_stream, '\t'.join(
                                [tokens[3], '0'] + tokens[:3] + tokens[4:]
                            )
        for line in sam_stream:
            tokens = line.strip().split('\t')
            print >>temp_stream, '\t'.join([tokens[0], '1'] + tokens[1:])
    import subprocess
    sorted_combined_file = os.path.join(temp_dir_path, 'combined.sorted.temp')
    subprocess.check_call(' '.join(['sort -k1,1 -k2,2n', combined_file, 
                                    '>', sorted_combined_file]),
                            bufsize=-1, shell=True)
    basewise_relevant, read_relevant = 0, 0
    # Initialize counters for computing accuracy metrics
    basewise_retrieved, basewise_intersection = 0, 0
    read_retrieved, read_intersection = 0, 0
    with open(sorted_combined_file) as sorted_combined_stream:
        for (name,), xpartition in xstream(sorted_combined_stream, 1):
            '''Dict mapping read names to alignments
            (chrom, 1-based start, 1-based end)'''
            true_maps = []
            saved = []
            for tokens in xpartition:
                saved.append(tokens)
                if tokens[0] == '0':
                    if len(tokens) < 12:
                        continue
                    chrom = tokens[1]
                    chrom_start = int(tokens[2])
                    chrom_end = int(tokens[3])
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
                                chrom_start + int(block_starts[i]),
                                chrom_start + int(block_starts[i])
                                + int(block_sizes[i]))
                                for i in xrange(block_count)]
                    true_maps.append(exons)
                    basewise_relevant += sum([int(block_size) for block_size
                                                in block_sizes])
                    read_relevant += 1
                elif tokens[0] == '1':
                    flag = int(tokens[1])
                    if flag & 256 or flag & 4:
                        '''Secondary alignment or unmapped and thus not
                        retrieved; ignore'''
                        continue
                    cigar, pos, seq = tokens[5], int(tokens[3]), tokens[9]
                    (dummy_md, mapped,
                        unmapped, clip_count) = dummy_md_and_mapped_offsets(
                                            cigar,
                                            clip_threshold=clip_threshold
                                        )
                    if unmapped:
                        # Too much clipping
                        continue
                    # Assume same read length for read pair
                    read_length = sum(
                            true_maps[0][i][2]
                                - true_maps[0][i][1]
                                for i in xrange(len(true_maps[0])))
                    basewise_retrieved += read_length - clip_count
                    read_retrieved += 1
                    # Try both /1 and /2; choose the best basewise result
                    intersected_base_count = 0
                    for true_map in true_maps:
                        if tokens[2] != true_map[0][0]:
                            '''chr is wrong, but this is still counted as a
                            retrieval above'''
                            continue
                        base_counter, base_truths = 0, set()
                        '''Each tuple in base_truths is
                        (index of base in read, mapped location)'''
                        for block in true_map:
                            base_truths.update([(base_counter + i, j + 1)
                                                    for i, j in enumerate(
                                                        xrange(
                                                            block[1], block[2]
                                                        ))])
                            base_counter += block[2] - block[1]
                        base_predictions = set()
                        if unmapped:
                            # Too much clipping
                            continue
                        _, _, _, exons = indels_introns_and_exons(
                                                        cigar,
                                                        dummy_md, pos, seq,
                                                        drop_deletions=True
                                                    )
                        mapped_index = 0
                        for exon in exons:
                            base_predictions.update(
                                        [(mapped[mapped_index + i], j)
                                                  for i, j in enumerate(
                                                    xrange(
                                                        exon[0], exon[1]
                                                    ))])
                            mapped_index += exon[1] - exon[0]
                        intersected_base_count = max(intersected_base_count,
                                len(
                                    base_predictions.intersection(base_truths)
                                ))
                    basewise_intersection += intersected_base_count
                    if intersected_base_count >= read_length * base_threshold:
                        read_intersection += 1
                    elif dump_incorrect:
                        # Incorrect alignment; write to stderr
                        print >>sys.stderr, '\t'.join(
                                ['.'.join(line) for line in saved]
                            )
                else:
                    raise RuntimeError(
                                'Invalid intermediate line.'
                            )
    return (basewise_retrieved, basewise_relevant, basewise_intersection,
            read_retrieved, read_relevant, read_intersection)

if __name__ == '__main__':
    import argparse
    # Print file's docstring if -h is invoked
    parser = argparse.ArgumentParser(description=__doc__, 
            formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--test', action='store_const', const=True,
        default=False,
        help='Run unit tests')
    args = parser.parse_known_args(sys.argv[1:])[0]
    if not args.test:
        parser = argparse.ArgumentParser(description=__doc__, 
            formatter_class=argparse.RawDescriptionHelpFormatter)
        parser.add_argument('-t', '--true-bed', type=str, required=True, 
            help='Full path of BED file containing true mappings')
        parser.add_argument('-g', '--generous', action='store_const', const=True,
            default=False,
            help='TopHat/STAR/HISAT cut off /1s and /2s from read names, even '
                 'in unpaired mode. This loses information. In generous mode, '
                 'this script provides extremely tight upper bounds on '
                 'precision and recall for TopHat/STAR/HISAT')
        parser.add_argument('-b', '--base-threshold', type=float,
            required=False, default=0.5,
            help=('Proportion of a read\'s bases that must align correctly for '
                  'the read to be considered a correct mapping'))
        parser.add_argument('-c', '--clip-threshold', type=float,
            required=False, default=1.0,
            help=('Proportion of a read\'s bases that are clipped above '
                  'which the read is considered unmapped'))
        parser.add_argument('--dump-incorrect', action='store_const',
            const=True,
            default=False,
            help='Write lines that aligned incorrectly to stderr')
        args = parser.parse_known_args(sys.argv[1:])[0]
        with open(args.true_bed) as true_bed_stream:
            (basewise_retrieved, basewise_relevant, basewise_intersection,
                read_retrieved, read_relevant, read_intersection) = go(
                        true_bed_stream,
                        generous=args.generous,
                        base_threshold=args.base_threshold,
                        clip_threshold=args.clip_threshold,
                        dump_incorrect=args.dump_incorrect
                    )
        basewise_precision = float(basewise_intersection) / basewise_retrieved
        basewise_recall = float(basewise_intersection) / basewise_relevant
        read_precision = float(read_intersection) / read_retrieved
        read_recall = float(read_intersection) / read_relevant
        print 'relevant instances\t%d\t%d' % (basewise_relevant, read_relevant)
        print 'retrieved instances\t%d\t%d' % (basewise_retrieved,
                                                    read_retrieved)
        print 'intersection\t%d\t%d' % (basewise_intersection,
                                            read_intersection)
        print 'precision\t%.12f\t%.12f' % (basewise_precision, read_precision)
        print 'recall\t%.12f\t%.12f' % (basewise_recall, read_recall)
        print 'fscore\t%.12f\t%.12f' % (
                (2 * basewise_precision * basewise_recall) /
                    (basewise_precision + basewise_recall),
                (2 * read_precision * read_recall) /
                    (read_precision + read_recall)
            )
    else:
        # Test go()
        del sys.argv[1:] # Don't choke on extra command-line parameters
        import unittest
        import shutil
        import tempfile

        class TestGo(unittest.TestCase):
            """ Tests go(). """
            def setUp(self):
                # For storing output of print_readletized_output()
                self.temp_dir_path = tempfile.mkdtemp()
                self.true_bed = os.path.join(self.temp_dir_path,
                                                        'true.bed')
                self.sam = os.path.join(self.temp_dir_path,
                                            'recovered.sam')
                # These bed lines are selected from a Flux simulation
                with open(self.true_bed, 'w') as true_bed_stream:
                    print >>true_bed_stream, (
"""chr2L\t70326\t70402\tchr2L:66584-71390W:NM_001258881:2:3789:2629:2857:A\t0\t-\t.\t.\t0,0,0\t1\t76\t0
chr2L\t200609\t200751\tchr2L:159032-203250W:NM_001272876:177:18434:15833:15999:A\t0\t-\t.\t.\t0,0,0\t2\t66,10\t0,132
chrX\t5525473\t5527141\tchrX:5519709-5542520C:NM_001272325:24:5914:1276:1427:A\t0\t+\t.\t.\t0,0,0\t3\t55,18,3\t0,749,1665"""
                    )

            def test_perfection(self):
                """ Fails if basewise/read-level instance counts are wrong. """
                # SAM lines are based on Rail's output
                with open(self.sam, 'w') as sam_stream:
                    print >>sam_stream, (
"""chr2L:66584-71390W:NM_001258881:2:3789:2629:2857:A\t16\tchr2L\t70327\t44\t76M\t*\t0\t0\tTTAGAAGTGCGTTAAAGCGCTCTATAAAACAGGCCCAGGAGCAACAGCTTCAGCTGAAGCGAGGCAATCATGAAGAHHHHHHHHH5HHHH555HHHHhHHH555HHHHHHHHhHHH5HHHHHHHHHHhhHhhHHHHH5HHHHhhhhHHHHHH
chr2L:159032-203250W:NM_001272876:177:18434:15833:15999:A\t16\tchr2L\t200610\t255\t66M66N10M\t*\t0\t0\tATCCAGGCAAGGCAATGGATATGCAGTTGGATGAGATGGACCGCATGTCAATGATTGCAGCTGTCGTTCAACAACA\tH55HHHHHHhhhHHhH5HHHHHHHHHhhHhhHHhhHHHhhhhHHhHhhhhhHHHHHHH55HHHHHHHHHHHHHHHH
chrX:5519709-5542520C:NM_001272325:24:5914:1276:1427:A\t0\tchrX\t5525474\t255\t55M694N18M898N3M\t*\t0\t0\tCCACCACATGCATCAGGCCAAAGCGAGCGATCACCTTGAAGCGATGGATGTTCAACTGAAATGGCGGCCTAGTCCG\tHHHhhhhhhhhhhhhhHhhhhhhHHhhhhhhhhHHHHHHHHHHH5HHHhhHhHHHH55H55HHHHHHHHHHHH5HH"""
                    )
                with open(self.true_bed) as true_bed_stream, \
                        open(self.sam) as sam_stream:
                    (basewise_retrieved, basewise_relevant,
                        basewise_intersection,
                        read_retrieved, read_relevant,
                        read_intersection) = go(
                                true_bed_stream,
                                sam_stream=sam_stream,
                                generous=False,
                                base_threshold=0.5
                            )
                # This should be perfect
                self.assertEquals(basewise_retrieved, 228)
                self.assertEquals(basewise_relevant, 228)
                self.assertEquals(basewise_intersection, 228)
                self.assertEquals(read_retrieved, 3)
                self.assertEquals(read_relevant, 3)
                self.assertEquals(read_intersection, 3)

            def test_deletions(self):
                """ Fails if basewise/read-level instance counts are wrong. """
                # SAM lines are based on Rail's output
                with open(self.sam, 'w') as sam_stream:
                    print >>sam_stream, (
"""chr2L:66584-71390W:NM_001258881:2:3789:2629:2857:A\t16\tchr2L\t70327\t44\t35M1D41M\t*\t0\t0\tTTAGAAGTGCGTTAAAGCGCTCTATAAAACAGGCCCAGGAGCAACAGCTTCAGCTGAAGCGAGGCAATCATGAAGAHHHHHHHHH5HHHH555HHHHhHHH555HHHHHHHHhHHH5HHHHHHHHHHhhHhhHHHHH5HHHHhhhhHHHHHH
chr2L:159032-203250W:NM_001272876:177:18434:15833:15999:A\t16\tchr2L\t200610\t255\t66M66N10M\t*\t0\t0\tATCCAGGCAAGGCAATGGATATGCAGTTGGATGAGATGGACCGCATGTCAATGATTGCAGCTGTCGTTCAACAACA\tH55HHHHHHhhhHHhH5HHHHHHHHHhhHhhHHhhHHHhhhhHHhHhhhhhHHHHHHH55HHHHHHHHHHHHHHHH
chrX:5519709-5542520C:NM_001272325:24:5914:1276:1427:A\t0\tchrX\t5525474\t255\t55M694N18M898N1D3M\t*\t0\t0\tCCACCACATGCATCAGGCCAAAGCGAGCGATCACCTTGAAGCGATGGATGTTCAACTGAAATGGCGGCCTAGTCCG\tHHHhhhhhhhhhhhhhHhhhhhhHHhhhhhhhhHHHHHHHHHHH5HHHhhHhHHHH55H55HHHHHHHHHHHH5HH"""
                    )
                with open(self.true_bed) as true_bed_stream, \
                        open(self.sam) as sam_stream:
                    (basewise_retrieved, basewise_relevant,
                        basewise_intersection,
                        read_retrieved, read_relevant,
                        read_intersection) = go(
                                true_bed_stream,
                                sam_stream=sam_stream,
                                generous=False,
                                base_threshold=0.5
                            )
                self.assertEquals(basewise_retrieved, 228)
                self.assertEquals(basewise_relevant, 228)
                self.assertEquals(basewise_intersection, 184)
                self.assertEquals(read_retrieved, 3)
                self.assertEquals(read_relevant, 3)
                self.assertEquals(read_intersection, 2)

            def test_soft_clips(self):
                """ Fails if basewise/read-level instance counts are wrong. """
                # SAM lines are based on Rail's output
                with open(self.sam, 'w') as sam_stream:
                    print >>sam_stream, (
"""chr2L:66584-71390W:NM_001258881:2:3789:2629:2857:A\t16\tchr2L\t70332\t44\t5S71M\t*\t0\t0\tTTAGAAGTGCGTTAAAGCGCTCTATAAAACAGGCCCAGGAGCAACAGCTTCAGCTGAAGCGAGGCAATCATGAAGAHHHHHHHHH5HHHH555HHHHhHHH555HHHHHHHHhHHH5HHHHHHHHHHhhHhhHHHHH5HHHHhhhhHHHHHH
chr2L:159032-203250W:NM_001272876:177:18434:15833:15999:A\t16\tchr2L\t200610\t255\t66M66N1M9S\t*\t0\t0\tATCCAGGCAAGGCAATGGATATGCAGTTGGATGAGATGGACCGCATGTCAATGATTGCAGCTGTCGTTCAACAACA\tH55HHHHHHhhhHHhH5HHHHHHHHHhhHhhHHhhHHHhhhhHHhHhhhhhHHHHHHH55HHHHHHHHHHHHHHHH
chrX:5519709-5542520C:NM_001272325:24:5914:1276:1427:A\t0\tchrX\t5525474\t255\t55M694N18M898N1M2S\t*\t0\t0\tCCACCACATGCATCAGGCCAAAGCGAGCGATCACCTTGAAGCGATGGATGTTCAACTGAAATGGCGGCCTAGTCCG\tHHHhhhhhhhhhhhhhHhhhhhhHHhhhhhhhhHHHHHHHHHHH5HHHhhHhHHHH55H55HHHHHHHHHHHH5HH"""
                    )
                with open(self.true_bed) as true_bed_stream, \
                        open(self.sam) as sam_stream:
                    (basewise_retrieved, basewise_relevant,
                        basewise_intersection,
                        read_retrieved, read_relevant,
                        read_intersection) = go(
                                true_bed_stream,
                                sam_stream=sam_stream,
                                generous=False,
                                base_threshold=0.8
                            )
                # This should be perfect
                self.assertEquals(basewise_retrieved, 212)
                self.assertEquals(basewise_relevant, 228)
                self.assertEquals(basewise_intersection, 212)
                self.assertEquals(read_retrieved, 3)
                self.assertEquals(read_relevant, 3)
                self.assertEquals(read_intersection, 3)

            def test_unmapped(self):
                """ Fails if basewise/read-level instance counts are wrong. """
                # SAM lines are based on Rail's output
                with open(self.sam, 'w') as sam_stream:
                    print >>sam_stream, (
"""chr2L:66584-71390W:NM_001258881:2:3789:2629:2857:A\t16\tchr2L\t70327\t44\t76M\t*\t0\t0\tTTAGAAGTGCGTTAAAGCGCTCTATAAAACAGGCCCAGGAGCAACAGCTTCAGCTGAAGCGAGGCAATCATGAAGAHHHHHHHHH5HHHH555HHHHhHHH555HHHHHHHHhHHH5HHHHHHHHHHhhHhhHHHHH5HHHHhhhhHHHHHH
chr2L:159032-203250W:NM_001272876:177:18434:15833:15999:A\t16\tchr2L\t200610\t255\t66M66N10M\t*\t0\t0\tATCCAGGCAAGGCAATGGATATGCAGTTGGATGAGATGGACCGCATGTCAATGATTGCAGCTGTCGTTCAACAACA\tH55HHHHHHhhhHHhH5HHHHHHHHHhhHhhHHhhHHHhhhhHHhHhhhhhHHHHHHH55HHHHHHHHHHHHHHHH
chrX:5519709-5542520C:NM_001272325:24:5914:1276:1427:A\t4\t*\t0\t0\t*\t*\t0\t0\tCCACCACATGCATCAGGCCAAAGCGAGCGATCACCTTGAAGCGATGGATGTTCAACTGAAATGGCGGCCTAGTCCG\tHHHhhhhhhhhhhhhhHhhhhhhHHhhhhhhhhHHHHHHHHHHH5HHHhhHhHHHH55H55HHHHHHHHHHHH5HH"""
                    )
                with open(self.true_bed) as true_bed_stream, \
                        open(self.sam) as sam_stream:
                    (basewise_retrieved, basewise_relevant,
                        basewise_intersection,
                        read_retrieved, read_relevant,
                        read_intersection) = go(
                                true_bed_stream,
                                sam_stream=sam_stream,
                                generous=False,
                                base_threshold=0.8
                            )
                self.assertEquals(basewise_retrieved, 152)
                self.assertEquals(basewise_relevant, 228)
                self.assertEquals(basewise_intersection, 152)
                self.assertEquals(read_retrieved, 2)
                self.assertEquals(read_relevant, 3)
                self.assertEquals(read_intersection, 2)

            def test_clip_threshold(self):
                """ Fails if basewise/read-level instance counts are wrong. """
                # SAM lines are based on Rail's output
                with open(self.sam, 'w') as sam_stream:
                    print >>sam_stream, (
"""chr2L:66584-71390W:NM_001258881:2:3789:2629:2857:A\t16\tchr2L\t70327\t44\t72M4S\t*\t0\t0\tTTAGAAGTGCGTTAAAGCGCTCTATAAAACAGGCCCAGGAGCAACAGCTTCAGCTGAAGCGAGGCAATCATGAAGAHHHHHHHHH5HHHH555HHHHhHHH555HHHHHHHHhHHH5HHHHHHHHHHhhHhhHHHHH5HHHHhhhhHHHHHH
chr2L:159032-203250W:NM_001272876:177:18434:15833:15999:A\t16\tchr2L\t200610\t255\t66M66N10M\t*\t0\t0\tATCCAGGCAAGGCAATGGATATGCAGTTGGATGAGATGGACCGCATGTCAATGATTGCAGCTGTCGTTCAACAACA\tH55HHHHHHhhhHHhH5HHHHHHHHHhhHhhHHhhHHHhhhhHHhHhhhhhHHHHHHH55HHHHHHHHHHHHHHHH
chrX:5519709-5542520C:NM_001272325:24:5914:1276:1427:A\t0\tchrX\t5525474\t255\t55M694N18M898N3M\t*\t0\t0\tCCACCACATGCATCAGGCCAAAGCGAGCGATCACCTTGAAGCGATGGATGTTCAACTGAAATGGCGGCCTAGTCCG\tHHHhhhhhhhhhhhhhHhhhhhhHHhhhhhhhhHHHHHHHHHHH5HHHhhHhHHHH55H55HHHHHHHHHHHH5HH"""
                    )
                with open(self.true_bed) as true_bed_stream, \
                        open(self.sam) as sam_stream:
                    (basewise_retrieved, basewise_relevant,
                        basewise_intersection,
                        read_retrieved, read_relevant,
                        read_intersection) = go(
                                true_bed_stream,
                                sam_stream=sam_stream,
                                generous=False,
                                base_threshold=0.5,
                                clip_threshold=0.05
                            )
                # Clip should kill read 1 here
                self.assertEquals(basewise_retrieved, 152)
                self.assertEquals(basewise_relevant, 228)
                self.assertEquals(basewise_intersection, 152)
                self.assertEquals(read_retrieved, 2)
                self.assertEquals(read_relevant, 3)
                self.assertEquals(read_intersection, 2)
                with open(self.sam, 'w') as sam_stream:
                    print >>sam_stream, (
"""chr2L:66584-71390W:NM_001258881:2:3789:2629:2857:A\t16\tchr2L\t70327\t44\t73M3S\t*\t0\t0\tTTAGAAGTGCGTTAAAGCGCTCTATAAAACAGGCCCAGGAGCAACAGCTTCAGCTGAAGCGAGGCAATCATGAAGAHHHHHHHHH5HHHH555HHHHhHHH555HHHHHHHHhHHH5HHHHHHHHHHhhHhhHHHHH5HHHHhhhhHHHHHH
chr2L:159032-203250W:NM_001272876:177:18434:15833:15999:A\t16\tchr2L\t200610\t255\t66M66N10M\t*\t0\t0\tATCCAGGCAAGGCAATGGATATGCAGTTGGATGAGATGGACCGCATGTCAATGATTGCAGCTGTCGTTCAACAACA\tH55HHHHHHhhhHHhH5HHHHHHHHHhhHhhHHhhHHHhhhhHHhHhhhhhHHHHHHH55HHHHHHHHHHHHHHHH
chrX:5519709-5542520C:NM_001272325:24:5914:1276:1427:A\t0\tchrX\t5525474\t255\t55M694N18M898N3M\t*\t0\t0\tCCACCACATGCATCAGGCCAAAGCGAGCGATCACCTTGAAGCGATGGATGTTCAACTGAAATGGCGGCCTAGTCCG\tHHHhhhhhhhhhhhhhHhhhhhhHHhhhhhhhhHHHHHHHHHHH5HHHhhHhHHHH55H55HHHHHHHHHHHH5HH"""
                    )
                with open(self.true_bed) as true_bed_stream, \
                        open(self.sam) as sam_stream:
                    (basewise_retrieved, basewise_relevant,
                        basewise_intersection,
                        read_retrieved, read_relevant,
                        read_intersection) = go(
                                true_bed_stream,
                                sam_stream=sam_stream,
                                generous=False,
                                base_threshold=0.5,
                                clip_threshold=0.05
                            )
                # Clip should kill read 1 here
                self.assertEquals(basewise_retrieved, 225)
                self.assertEquals(basewise_relevant, 228)
                self.assertEquals(basewise_intersection, 225)
                self.assertEquals(read_retrieved, 3)
                self.assertEquals(read_relevant, 3)
                self.assertEquals(read_intersection, 3)

            def tearDown(self):
                # Kill temporary directory
                shutil.rmtree(self.temp_dir_path)

        unittest.main()
