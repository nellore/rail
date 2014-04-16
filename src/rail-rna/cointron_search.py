#!/usr/bin/env python
"""
Rail-RNA-cointron_search

Follows Rail-RNA-realign_readlets
Precedes Rail-RNA-cointron_fasta

Reduce step in MapReduce pipelines that builds a map of introns cooccurring
on reads.

Input (read from stdin)
----------------------------
Tab-delimited input tuple columns, where each line corresponds to a readlet
from a distinct read rather than a unique readlet sequence:
1. Read sequence ID
2. Displacement of readlet's 5' end from read's 5' end + '\x1e' +
    displacement of readlet's 3' end from read's 3' end (+, for EXACTLY one
    readlet of a read sequence, + '\x1e' + read sequence + '\x1e' + number of
    instances of read sequence + '\x1e' + number of instances of read
    sequence's reversed complement + '\x1e' + (an '\x1f'-separated set of
    unique sample labels with read sequences that match the original read
    sequence) + '\x1e' + (an '\x1f'-separated set of unique sample labels with
    read sequences that match the reversed complement of the original read
    sequence))
3. '\x1f'-separated list of alignment RNAMEs or '\x1c' if no alignments
4. '\x1f'-separated list of alignment FLAGs or '\x1c' if no alignments
5. '\x1f'-separated list of alignment POSes or '\x1c' if no alignments

Input is partitioned by field 1, the read sequence ID.

Hadoop output (written to stdout)
----------------------------
Tab-delimited tuple columns:
1. Reference name (RNAME in SAM format) + 
    '+' or '-' indicating which strand is the sense strand
2. Comma-separated list of intron start positions in configuration
3. Comma-separated list of intron end positions in configuration
4. left_extend_size: by how many bases on the left side of an intron the
    reference should extend
5. right_extend_size: by how many bases on the right side of an intron the
    reference should extend

ALL OUTPUT COORDINATES ARE 1-INDEXED.
"""
import sys
import os
import site
import re
import numpy as np

base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
site.addsitedir(os.path.join(base_path, 'dooplicity'))

import dooplicity as dp

# Initialize global variables for tracking number of input/output lines
_input_line_count = 0
_output_line_count = 0

def strand_and_introns(rname, offset, readlet_size):
    """ Extracts strand and introns overlapped by readlet from RNAME and POS.

        Introns are encoded in an alignment's RNAME as follows:

        original RNAME + '+' or '-' indicating which strand is the sense
        strand + ';' + start position of sequence + ';' + comma-separated
        list of subsequence sizes framing introns + ';' + comma-separated
        list of intron sizes.

        rname: RNAME from intron reference
        offset: offset from beginning of intron reference
        readlet_size: size of readlet

        Return value: tuple (strand, True if sense strand is forward strand
                                else False, alignment_start_position,
                                alignment_end_position, 
                                tuple of tuples (pos, end_pos) of start
                                and end positions of introns OR None if 
                                no introns are overlapped); here,
                            alignment_start_position and alignment_end_position
                            are the start and end positions of the 
                            alignment _along the original reference_
    """
    strand, original_pos, exon_sizes, intron_sizes = rname.split(';') 
    exon_sizes = [int(exon_size) for exon_size in exon_sizes.split(',')]
    intron_sizes = [int(intron_size) for intron_size
                        in intron_sizes.split(',')]
    assert len(intron_sizes) == len(exon_sizes) - 1
    cigar_sizes = [None]*(len(exon_sizes) + len(intron_sizes))
    cigar_sizes[::2] = exon_sizes
    cigar_sizes[1::2] = intron_sizes
    partial_sizes = [sum(exon_sizes[:j+1]) for j in xrange(len(exon_sizes))]
    start_index, end_index = None, None
    for j, partial_size in enumerate(partial_sizes):
        if partial_size > offset:
            start_index = j
            new_offset = (offset - partial_sizes[j-1] if j != 0 else offset)
            break
    end_offset = readlet_size + offset
    for j, partial_size in enumerate(partial_sizes):
        if partial_size >= end_offset:
            end_index = j
            break
    original_pos = int(original_pos) + new_offset \
                    + sum(cigar_sizes[:start_index*2])
    if start_index is None or end_index is None:
        raise RuntimeError('Invalid alignment; sum of exon sizes doesn\'t '
                           'agree with size of reference sequence.')
    if start_index == end_index:
        return None
    else:
        assert start_index < end_index
        cigar_sizes[start_index*2] = partial_sizes[start_index] - offset
        cigar_sizes[end_index*2] = end_offset - partial_sizes[end_index-1]
        last_pos = original_pos
        introns = []
        for j in xrange(start_index*2+1, end_index*2+1, 2):
            last_pos += cigar_sizes[j-1]
            introns.append((last_pos, last_pos + cigar_sizes[j]))
            last_pos += (cigar_sizes[j] + cigar_sizes[j+1])
        return (strand, True if strand[-1] == '-' else False, 
                original_pos, last_pos, tuple(introns), readlet_size)

def different_introns_overlap(intron_iterable_1, intron_iterable_2):
    """ Test whether different introns in two iterables overlap.

        intron_iterable_1: an iterable, each of whose items is an intron
            specified by the tuple (start position, end position)
        intron_iterable_2: takes the same form as intron_iterable_1

        Return value: True iff there is at least one pair of distinct introns,
            one from intron_iterable_1 and the other from intron_iterable_2,
            that overlap. Here, two introns overlap if there is not at least
            one exonic base between them.
    """
    for pos_1, end_pos_1 in intron_iterable_1:
        for pos_2, end_pos_2 in intron_iterable_2:
            if (pos_1, end_pos_1) == (pos_2, end_pos_2):
                continue
            if min(end_pos_1, end_pos_2) - max(pos_1, pos_2) >= 0:
                return True
    return False

def compatible(alignment_1, alignment_2, fudge=4):
    """ Determines whether alignments are compatible.

        Two alignments are compatible iff:
        1) They are on the same strand; that is, they have the same RNAME
        (chromosome), agree on which strand is the sense strand, and 
        both align to either the forward or the reverse strand.
        2) Their introns are consistent; that is, when an intron from
        alignment_1 and an intron from alignment_2 overlap, these introns
        are the same intron.
        3) The number of bases between alignment_1's start position and
        alignment_2's start position less any intervening intronic bases
        from either alignment is equal to the difference in corresponding
        readlet displacements from the 5' end of the read when the alignments
        overlap, and is within fudge of this difference when they don't.

        alignment_1: tuple (strand, True if sense strand is forward strand else
                            False, reference start position of alignment,
                            reference end position of alignment,
                            tuple of intron tuples overlapped by alignment,
                            readlet size,
                            True if alignment is to forward strand else False,
                            displacement of readlet from 5' end of read)

        alignment_2: takes the same form as alignment_1

        Return value: False if alignments are incompatible; True if alignments
            are compatible.
    """
    if alignment_1[0] != alignment_2[0] or alignment_1[6] != alignment_2[6]:
        # Strands don't agree or alignments aren't to same strand
        return False
    if different_introns_overlap(alignment_1[4], alignment_2[4]):
        # Two distinct introns overlap
        return False
    if alignment_1[-1] >= alignment_2[-1]:
        first_alignment = alignment_2
        second_alignment = alignment_1
    else:
        first_alignment = alignment_1
        second_alignment = alignment_2
    overlap_with_introns = second_alignment[2] - first_alignment[3]
    if overlap_with_introns >= 0:
        # No overlap; impose that readlet separation is within that tolerated
        if abs(second_alignment[-1] - first_alignment[-1] + first_alignment[-3]
                - second_alignment[2] + first_alignment[3]) <= fudge:
            return True
        else:
            return False
    else:
        '''Overlap; since there are no overlapping introns between alignments,
        compatibility necessarily holds.'''
        return True

def go(input_stream=sys.stdin, output_stream=sys.stdout, 
        verbose=False, stranded=False):
    """ Runs Rail-RNA-intron_search.

        Input (read from stdin)
        ----------------------------
        Tab-delimited input tuple columns, where each line corresponds to a
        readlet from a distinct read rather than a unique readlet sequence:
        1. Read sequence ID
        2. Displacement of readlet's 5' end from read's 5' end + '\x1e' +
            displacement of readlet's 3' end from read's 3' end (+, for EXACTLY
            one readlet of a read sequence, '\x1e' + read sequence + '\x1e' +
            number of instances of read sequence + '\x1e' + number of instances
            of read sequence's reversed complement + '\x1e' + (an
            '\x1f'-separated set of unique sample labels with read sequences
            that match the original read sequence) + '\x1e' + (an
            '\x1f'-separated set of unique sample labels with read sequences
            that match the reversed complement of the original read sequence))
        3. '\x1f'-separated list of alignment RNAMEs or '\x1c' if no alignments
        4. '\x1f'-separated list of alignment FLAGs or '\x1c' if no alignments
        5. '\x1f'-separated list of alignment POSes or '\x1c' if no alignments

        Input is partitioned by field 1, the read sequence ID.

        Hadoop output (written to stdout)
        ----------------------------
        Tab-delimited columns:
        1. Reference name (RNAME in SAM format) + ';' + partition number +  
            ('+' or '-' indicating which strand is the sense strand if input
                reads are strand-specific -- that is, --stranded in was
                invoked; otherwise, there is no terminal '+' or '-')
        2. Candidate intron start (inclusive) on forward strand (1-indexed)
        3. Candidate intron end (exclusive) on forward strand (1-indexed)
        4. '\x1f'-separated list of sample (label)s in which intron was
            detected
        5. Total number of reads supporting intron

        ALL OUTPUT COORDINATES ARE 1-INDEXED.

        input_stream: where to find input reads.
        output_stream: where to emit exonic chunks and introns.
        verbose: True iff more informative messages should be written to
            stderr.
        stranded: True iff input reads are strand-specific; this affects
            whether an output partition has a terminal '+' or '-' indicating
            the sense strand.

        No return value.
    """
    global _input_line_count, _output_line_count
    for (seq_id,), xpartition in dp.xstream(input_stream, 1):
        '''Input is readletized, and readlets must be composed.'''
        readlet_displacements, collected_readlets = [], []
        seq_info_captured = False
        for value in xpartition:
            assert len(value) == 4
            _input_line_count += 1
            seq_info = value[0].split('\x1e')
            seq_info_size = len(seq_info)
            assert seq_info_size >= 2
            if seq_info_size > 2:
                assert seq_info_size == 7
                seq = seq_info[2]
                seq_size = len(seq)
                seq_count = int(seq_info[3])
                reversed_complement_seq_count = int(seq_info[4])
                sample_labels = seq_info[5]
                reversed_complement_sample_labels = seq_info[6]
                seq_info_captured = True
            if value[-1] != '\x1c':
                # If there are alignments, add them to collected_readlets
                collected_readlets.append(
                    zip(
                            value[-3].split('\x1f'),
                            [(int(flag) & 16) != 0
                                for flag in value[-2].split('\x1f')],
                            [int(pos) for pos in value[-1].split('\x1f')]
                        )
                )
                readlet_displacements.append(
                    (int(seq_info[0]), int(seq_info[1]))
                )
        if len(collected_readlets):
            assert seq_info_captured, 'Sequence info is not in a collected ' \
                                      'readlet'
            '''Each multireadlet is a list of tuples representing readlet
            alignments. Each tuple takes the form
            (strand, True if sense strand is forward strand else False,
                reference start position of alignment, reference end position
                of alignment, tuple of intron tuples overlapped by alignment,
                readlet size,
                True if alignment is to forward strand else False, displacement
                of readlet from 5' end of read).'''
            multireadlets = [[strand_and_introns(rname, pos - 1, seq_size
                                - readlet_displacements[i][1]
                                - readlet_displacements[i][0]) + 
                               (reverse_strand,
                                readlet_displacements[i][1]
                                    if reverse_strand 
                                    else readlet_displacements[i][0]) 
                                for rname, reverse_strand, pos in multireadlet]
                                for i, multireadlet
                                in enumerate(collected_readlets)]
            # Eliminate alignments that do not overlap introns
            multireadlets = [[alignment for alignment in multireadlet if
                                alignment[0] is not None] for multireadlet in
                                multireadlets]
            # Eliminate alignments that do not overlap introns
            multireadlets = [set([alignment for alignment in multireadlet if
                                alignment[0] is not None])
                                for multireadlet in multireadlets]
            # Cull unique alignments and sort by strand/positions
            alignments = sorted(
                                        [list(multireadlet)[0]
                                            for multireadlet in multireadlets
                                            if len(multireadlet) == 1],
                                        key=lambda alignment: alignment[:4]
                                    )
            if alignments and False not in \
                                [compatible(alignments[i], alignments[i-1])
                                    for i in xrange(1, len(alignments))]:
                '''If all successive alignments are compatible (see that
                    function's docstring for a definition), output a
                    cointron.'''
                left_extend_size = (alignments[0][4][0][0] # first intron start
                                        - alignments[0][2] # alignment start
                                        + alignments[0][-1]) # displacement
                right_extend_size = (seq_size
                                        - alignments[-1][-1] # displacement
                                        - alignments[-1][5] # readlet size
                                        - alignments[-1][4][-1][1] # last
                                                                   # intron end
                                        + alignments[-1][3]) # alignment end
                introns_to_add = set()
                for alignment in alignments:
                    for intron in alignment[4]:
                        introns_to_add.add(intron)
                intron_starts_and_ends = zip(*sorted(list(introns_to_add)))
                print >>output_stream, 'intron\t%s\t%s\t%s\t%d\t%d' % (
                        alignments[0][0], # strand
                        ','.join(map(str, intron_starts_and_ends[0])),
                        ','.join(map(str, intron_starts_and_ends[1])),
                        left_extend_size,
                        right_extend_size
                    )
                _output_line_count += 1

if __name__ == '__main__':
    import argparse
    # Print file's docstring if -h is invoked
    parser = argparse.ArgumentParser(description=__doc__, 
                formatter_class=argparse.RawDescriptionHelpFormatter)
    # Add command-line arguments
    parser.add_argument('--verbose', action='store_const', const=True,
        default=False,
        help='Print out extra debugging statements')
    parser.add_argument(\
        '--stranded', action='store_const', const=True, default=False,
        help='Assume input reads come from the sense strand; then partitions '
             'in output have terminal + and - indicating sense strand')
    parser.add_argument('--test', action='store_const', const=True,
        default=False,
        help='Run unit tests; DOES NOT NEED INPUT FROM STDIN, AND DOES NOT '
             'WRITE EXONS AND INTRONS TO STDOUT')

    args = parser.parse_args(sys.argv[1:])

if __name__ == '__main__' and not args.test:
    import time
    start_time = time.time()
    go(verbose=args.verbose, stranded=args.stranded)
    print >> sys.stderr, 'DONE with intron_search.py; in/out=%d/%d; ' \
        'time=%0.3f s' % (_input_line_count, _output_line_count,
                                time.time() - start_time)
elif __name__ == '__main__':
    # Test units
    del sys.argv[1:] # Don't choke on extra command-line parameters
    # Insert unit tests here
    unittest.main()