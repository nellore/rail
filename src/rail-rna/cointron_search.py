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
2. First intron start position in configuration
3. Rest of intron start positions in configuration or '\x1c' if there are none
4. Comma-separated list of intron end positions in configuration
5. left_extend_size: by how many bases on the left side of an intron the
    reference should extend
6. right_extend_size: by how many bases on the right side of an intron the
    reference should extend
7. Read sequence

ALL OUTPUT COORDINATES ARE 1-INDEXED.
"""
import sys
import os
import site
import re
import numpy as np
from collections import defaultdict
import random
import argparse

base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
site.addsitedir(os.path.join(base_path, 'dooplicity'))

import dooplicity as dp

# Initialize global variables for tracking number of input/output lines
_input_line_count = 0
_output_line_count = 0

def rname_and_introns(rname, offset, readlet_size):
    """ Extracts strand and introns overlapped by readlet from RNAME and POS.

        Some alignments are purely exonic and have no semicolons in the RNAME;
        in this case, no introns are returned. Introns are encoded in an
        alignment's RNAME as follows:

        original RNAME + '+' or '-' indicating which strand is the sense
        strand + ';' + start position of sequence + ';' + comma-separated
        list of subsequence sizes framing introns + ';' + comma-separated
        list of intron sizes + ';' + distance to previous intron or 'NA' if
        beginning of strand + ';' distance to next intron or 'NA' if end of
        strand.

        rname: RNAME from intron reference
        offset: offset from beginning of intron reference
        readlet_size: size of readlet

        Return value: tuple (rname, True if sense strand is forward strand or
                                False if sense strand is reverse strand or
                                None if not known, alignment_start_position,
                                alignment_end_position, 
                                tuple of tuples (pos, end_pos) of start
                                and end positions of introns OR None if 
                                no introns are overlapped, readlet_size,
                                distance to previous intron or None if 
                                beginning of strand, distance to next
                                intron or None if end of strand); here,
                            alignment_start_position and alignment_end_position
                            are the start and end positions of the 
                            alignment _along the original reference_
    """
    try:
        (strand, original_pos, exon_sizes,
            intron_sizes, left_size, right_size) = rname.split(';')
    except ValueError:
        # Exonic alignment
        return (rname, None, offset + 1,
                offset + 1 + readlet_size, None, readlet_size)
    print >>sys.stderr, rname
    if left_size == 'NA':
        left_size = None
    if right_size == 'NA':
        right_size = None
    rname = strand[:-1]
    reverse_strand = True if strand[-1] == '+' else False
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
        return (rname, None, original_pos,
                original_pos + readlet_size, None, readlet_size)
    else:
        assert start_index < end_index
        if start_index != 0:
            left_size = cigar_sizes[start_index*2]
        cigar_sizes[start_index*2] = partial_sizes[start_index] - offset
        try:
            left_size -= cigar_sizes[start_index*2]
        except TypeError:
            left_size = None
        if end_index * 2 != len(cigar_sizes) - 1:
            right_size = cigar_sizes[end_index * 2]
        cigar_sizes[end_index*2] = end_offset - partial_sizes[end_index-1]
        try:
            right_size -= cigar_sizes[end_index*2]
        except TypeError:
            right_size = None
        last_pos = original_pos
        introns = []
        for j in xrange(start_index*2+1, end_index*2+1, 2):
            last_pos += cigar_sizes[j-1]
            introns.append((last_pos, last_pos + cigar_sizes[j]))
            last_pos += (cigar_sizes[j] + cigar_sizes[j+1])
        if start_index != 0:
            left_size = cigar_sizes[start_index*2 - 1]
        return (rname, reverse_strand, 
                original_pos, last_pos, tuple(introns), readlet_size,
                left_size, right_size)

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

def maximal_cliques(intron_alignments):
    """ Based on http://www.kuchaev.com/files/graph.py .
    """
    cliques = []
    alignment_count = len(intron_alignments)
    intron_alignment_list = list(intron_alignments)
    if alignment_count == 1:
        return [intron_alignment_list]
    if alignment_count == 2:
        if not different_introns_overlap(intron_alignment_list[0][4],
                                            intron_alignment_list[1][4]):
            return [sorted(intron_alignment_list,
                        key=lambda alignment: alignment[2])]
        else:
            return [[intron_alignment_list[0]], [intron_alignment_list[1]]]
    graph = defaultdict(set)
    for i in xrange(alignment_count):
        for j in xrange(i+1, alignment_count):
            if not different_introns_overlap(intron_alignment_list[i][4],
                                                intron_alignment_list[j][4]):
                graph[intron_alignment_list[i]].add(intron_alignment_list[j])
                graph[intron_alignment_list[j]].add(intron_alignment_list[i])
    nodes = set(graph.keys())
    stack = [(set(), intron_alignments,
                set(), None, len(intron_alignments))]
    while stack:
        (c_compsub, c_candidates, c_not, c_nd, c_disc_num) = stack.pop()
        if not len(c_candidates) and not len(c_not):
            if len(c_compsub) > 2:
                cliques.append(
                        sorted(list(c_compsub),
                                key=lambda alignment: alignment[2])
                    )
                continue
        for u in list(c_candidates):
            if c_nd is None or u not in graph[c_nd]:
                c_candidates.remove(u)
                Nu=graph[u]                               
                new_compsub=set(c_compsub)
                new_compsub.add(u)
                new_candidates=set(c_candidates.intersection(Nu))
                new_not=set(c_not.intersection(Nu))  
                if c_nd is not None:
                    if c_nd in new_not:
                        new_disc_num=c_disc_num-1
                        if new_disc_num>0:
                            new_search_node=(new_compsub,new_candidates,
                                                new_not,c_nd,new_disc_num)                        
                            stack.append(new_search_node)
                    else:
                        new_disc_num=len(nodes)
                        new_nd=c_nd
                        for cand_nd in new_not:
                            cand_disc_num=len(new_candidates)\
                                -len(new_candidates.intersection(
                                        graph[cand_nd])
                                    ) 
                            if cand_disc_num<new_disc_num:
                                new_disc_num=cand_disc_num
                                new_nd=cand_nd
                        new_search_node=(new_compsub,new_candidates,
                                            new_not,new_nd,new_disc_num)                        
                        stack.append(new_search_node)                
                else:
                    new_search_node=(new_compsub,new_candidates,
                                        new_not,c_nd,c_disc_num)
                    stack.append(new_search_node)
                c_not.add(u) 
                new_disc_num=0
                for x in c_candidates:
                    if u not in graph[x]:
                        new_disc_num+=1
                if new_disc_num<c_disc_num and new_disc_num>0:
                    new1_search_node=(c_compsub,c_candidates,
                                        c_not,u,new_disc_num)
                    stack.append(new1_search_node)
                else:
                    new1_search_node=(c_compsub,c_candidates,
                                        c_not,c_nd,c_disc_num)
                    stack.append(new1_search_node)    
    return cliques

def selected_introns_by_clustering(multireadlets, seed=0):
    '''(rname, True if sense strand is forward strand or
                                False if sense strand is reverse strand or
                                None if not known, alignment_start_position,
                                alignment_end_position, 
                                tuple of tuples (pos, end_pos) of start
                                and end positions of introns OR None if 
                                no introns are overlapped, readlet_size,
                                distance to previous intron or None if
        beginning of strand, distance to next intron or None if end of
        strand
                True if alignment is to forward strand else False, displacement
                of readlet from 5' end of read)'''
    random.seed(seed)
    alignments = [alignment + (i,)
                    for i, multireadlet in enumerate(multireadlets)
                    for alignment in multireadlet]
    unclustered_alignments = range(len(alignments))
    clustered_alignments = []
    while unclustered_alignments:
        pivot = i = random.choice(unclustered_alignments)
        new_unclustered_alignments = []
        alignment_cluster = defaultdict(list)
        for j in unclustered_alignments:
            if j == pivot: continue
            if not (alignments[i][-1] == alignments[j][-1] or
                alignments[i][0] != alignments[j][0] or
                (alignments[i][1] is not None and alignments[j][1] is not None
                 and alignments[i][1] != alignments[j][1]) or
                alignments[i][-3] != alignments[j][-3] or
                (alignments[i][-2] == alignments[j][-2] and
                   alignments[i][2] != alignments[j][2]) or
                ((alignments[i][-2] < alignments[j][-2]) !=
                    (alignments[i][2] < alignments[j][2]))):
                alignment_cluster[alignments[j][-1]].append(j)
            else:
                new_unclustered_alignments.append(j)
        # Choose alignments closest to pivot in each multireadlet group
        alignment_cluster_list = [pivot]
        for group in alignment_cluster:
            overlap_distances = [max(alignments[i][2], alignments[j][2])
                                    - min(alignments[i][3], alignments[j][3])
                                    for j in alignment_cluster[group]]
            min_overlap_distance = min(overlap_distances)
            for j in xrange(len(overlap_distances)):
                if overlap_distances[j] == min_overlap_distance:
                    alignment_cluster_list.append(alignment_cluster[group][j])
                else:
                    new_unclustered_alignments.append(
                            alignment_cluster[group][j]
                        )
        clustered_alignments.append(
                [alignments[j] for j in alignment_cluster_list]
            )
        unclustered_alignments = new_unclustered_alignments
    cluster_sizes = [len(set([alignment[-1] for alignment in cluster]))
                        for cluster in clustered_alignments]
    largest_cluster_size = max(cluster_sizes)
    largest_clusters = []
    for i, cluster_size in enumerate(cluster_sizes):
        if cluster_size == largest_cluster_size:
            largest_clusters.append(clustered_alignments[i])
    multimap_count = len(largest_clusters)
    if multimap_count:
        # Return ALL largest clusters
        return [set([alignment[:-1] for alignment in largest_clusters[i]
                    if alignment[1] is not None])
                for i in xrange(multimap_count)]
    else:
        return []

def go(input_stream=sys.stdin, output_stream=sys.stdout, 
        verbose=False, stranded=False, fudge=10):
    """ Runs Rail-RNA-intron_search.

        Reduce step in MapReduce pipelines that builds a map of introns
        cooccurring on reads.

        Input (read from stdin)
        ----------------------------
        Tab-delimited input tuple columns, where each line corresponds to a
        readlet from a distinct read rather than a unique readlet sequence:
        1. Read sequence ID
        2. Displacement of readlet's 5' end from read's 5' end + '\x1e' +
            displacement of readlet's 3' end from read's 3' end (+, for EXACTLY
            one readlet of a read sequence, + '\x1e' + read sequence + '\x1e'
            + number of instances of read sequence + '\x1e' + number of
            instances of read sequence's reversed complement + '\x1e'
            + (an '\x1f'-separated set of unique sample labels with read
            sequences that match the original read sequence) + '\x1e'
            + (an '\x1f'-separated set of unique sample labels with read
            sequences that match the reversed complement of the original read
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
        2. First intron start position in configuration
        3. Rest of intron start positions in configuration or '\x1c' if there
            are none
        4. Comma-separated list of intron end positions in configuration
        5. left_extend_size: by how many bases on the left side of an intron
            the reference should extend
        6. right_extend_size: by how many bases on the right side of an intron
            the reference should extend
        7. Read sequence

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
        collected_readlets = defaultdict(list)
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
                collected_readlets[
                        (int(seq_info[0]), int(seq_info[1]))
                    ].extend(
                            zip(
                                value[-3].split('\x1f'),
                                [(int(flag) & 16) != 0
                                    for flag in value[-2].split('\x1f')],
                                [int(pos) for pos in value[-1].split('\x1f')]
                            )
                        )
        if len(collected_readlets):
            assert seq_info_captured, 'Sequence info is not in a collected ' \
                                      'readlet'
            '''Each multireadlet is a list of tuples representing readlet
            alignments. Each tuple takes the form
            (rname, True if sense strand is forward strand else False,
                reference start position of alignment, reference end position
                of alignment, tuple of intron tuples overlapped by alignment,
                readlet size, distance to previous intron or None if
                beginning of strand, distance to next intron or None if end of
                strand, True if alignment is to forward strand else False,
                displacement of readlet from 5' end of read).'''
            multireadlets = [set([rname_and_introns(rname, pos - 1, seq_size
                                - displacement[1] - displacement[0])
                                + (reverse_strand, displacement[1]
                                    if reverse_strand else displacement[0]) 
                                for rname, reverse_strand, pos
                                in multireadlet])
                                for displacement, multireadlet
                                in collected_readlets.items()]
            if False not in [alignment[1] is None 
                                for multireadlet in multireadlets
                                for alignment in multireadlet]:
                # No introns to see here
                continue
            clusters = selected_introns_by_clustering(
                                multireadlets,
                                seed=seq
                            )
            if len(clusters):
                for selected_introns in clusters:
                    for alignments in maximal_cliques(selected_introns):
                        left_extend_size = (alignments[0][4][0][0] # first
                                                                   # intron 
                                                                   # start
                                            - alignments[0][2] # alignment
                                                               # start
                                            + alignments[0][-1]) # displacement
                        right_extend_size = (seq_size
                                            - alignments[-1][-1] # displacement
                                            - alignments[-1][5] # readlet size
                                            - alignments[-1][4][-1][1] # last
                                                                       # intron
                                                                       # end
                                            + alignments[-1][3]) # alignment
                                                                 # end
                        introns_to_add = set()
                        for alignment in alignments:
                            for intron in alignment[4]:
                                introns_to_add.add(intron)
                        intron_starts_and_ends = zip(*sorted(
                                                          list(introns_to_add))
                                                        )
                        '''Ensure number of exonic bases is within fudge of
                        read size.'''
                        if abs(seq_size 
                                - sum([intron_starts_and_ends[0][i+1] 
                                        - intron_starts_and_ends[1][i]
                                        for i 
                                        in xrange(
                                            len(intron_starts_and_ends[0])
                                                        - 1)])
                                - left_extend_size
                                - right_extend_size) > fudge:
                            if verbose:
                                print >>sys.stderr, 'Killing intron combo', \
                                    intron_starts_and_ends, 'because its ' \
                                    'exonic base sum was not within ' \
                                    'fudge=%d bases of the sequence ' \
                                    'size=%d' % (fudge, seq_size)
                            continue
                        '''Determine by how much reference should be extended
                        on either side of intron.'''
                        if left_extend_size > alignments[0][6]:
                            left_size = left_extend_size
                        else:
                            left_size = min(left_extend_size + fudge,
                                            alignments[0][6])
                        if right_extend_size > alignments[-1][7]:
                            right_size = right_extend_size
                        else:
                            right_size = min(right_extend_size + fudge,
                                             alignments[-1][7])
                        print >>output_stream, 'intron\t%s\t%d\t%s\t%s' \
                                '\t%d\t%d\t%s' % (alignments[0][0] + ('+' if 
                                                    alignments[0][1] else '-'),
                                intron_starts_and_ends[0][0],
                                ','.join(map(str,
                                                intron_starts_and_ends[0][1:]))
                                if len(intron_starts_and_ends[0][1:]) \
                                    else '\x1c',
                                ','.join(map(str, intron_starts_and_ends[1])),
                                left_size,
                                right_size,
                                seq
                            )
                        _output_line_count += 1

if __name__ == '__main__':
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
    parser.add_argument('--fudge', type=int, required=False,
        default=15,
        help='Permits a sum of exonic bases for an intron combo to be within '
             'the specified number of bases of a read sequence\'s size; '
             'this allows for indels with respect to the reference')

    args = parser.parse_args(sys.argv[1:])

if __name__ == '__main__' and not args.test:
    import time
    start_time = time.time()
    go(verbose=args.verbose, stranded=args.stranded, fudge=args.fudge)
    print >> sys.stderr, 'DONE with cointron_search.py; in/out=%d/%d; ' \
        'time=%0.3f s' % (_input_line_count, _output_line_count,
                                time.time() - start_time)
elif __name__ == '__main__':
    # Test units
    del sys.argv[1:] # Don't choke on extra command-line parameters
    # Insert unit tests here
    unittest.main()