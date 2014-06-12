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
from collections import defaultdict
import random
import argparse

base_path = os.path.abspath(
                    os.path.dirname(os.path.dirname(os.path.dirname(
                        os.path.realpath(__file__)))
                    )
                )
utils_path = os.path.join(base_path, 'rna', 'utils')
site.addsitedir(utils_path)
site.addsitedir(base_path)

from dooplicity.tools import xstream

# Initialize global variables for tracking number of input/output lines
_input_line_count = 0
_output_line_count = 0

def running_sum(iterable):
    """ Generates a running sum of the numbers in an iterable

        iterable: some iterable with numbers

        Yield value: next value in running sum
    """
    total = 0
    for number in iterable:
        total += number
        yield total

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
                                beginning of strand / not relevant (exonic
                                alignment), distance to next intron or None if
                                end of strand / not relevant
                                (exonic alignment)); here,
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
                offset + 1 + readlet_size, None, readlet_size,
                None, None)
    try:
        left_size = int(left_size)
    except ValueError:
        left_size = None
    try:
        right_size = int(right_size)
    except ValueError:
        right_size = None
    rname = strand[:-1]
    reverse_strand = True if strand[-1] == '-' else False
    exon_sizes = [int(exon_size) for exon_size in exon_sizes.split(',')]
    intron_sizes = [int(intron_size) for intron_size
                        in intron_sizes.split(',')]
    assert len(intron_sizes) == len(exon_sizes) - 1
    for i, exon_sum in enumerate(running_sum(exon_sizes)):
        if exon_sum > offset: break
    # Compute start position of alignment
    pos = offset + sum(intron_sizes[:i]) + int(original_pos)
    # Adjust exon/intron lists so they start where alignment starts
    exon_sizes = exon_sizes[i:]
    if i != 0:
        left_size = exon_sizes[0]
    exon_sizes[0] = exon_sum - offset
    intron_sizes = intron_sizes[i:]
    for i, exon_sum in enumerate(running_sum(exon_sizes)):
        if exon_sum >= readlet_size: break
    if i != len(exon_sizes) - 1:
        right_size = exon_sizes[i]
    intron_size_sum = 0
    introns = []
    for j in xrange(i):
        intron_pos = sum(exon_sizes[:j+1]) + sum(intron_sizes[:j]) + pos
        intron_end_pos = intron_pos + intron_sizes[j]
        introns.append((intron_pos, intron_end_pos))
        intron_size_sum += intron_sizes[j]
    end_pos = readlet_size + intron_size_sum + pos
    if introns:
        return (rname, reverse_strand, 
                    pos, end_pos, tuple(introns), readlet_size,
                    left_size, right_size)
    return (rname, None, pos, end_pos, None, readlet_size, None, None)

def different_introns_overlap(intron_iterable_1, intron_iterable_2):
    """ Test whether distinct introns in two iterables overlap.

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

def cointron_adjacencies(intron_alignments):
    """ Generates adjacency matrix for cointron graph described below.

        Consider an undirected graph where each node corresponds to
        an alignment overlapping at least one cointron, as specified by a
        tuple from the set intron_alignments. Place an edge between two nodes
        if no two _distinct_ introns overlap, where one intron associated with
        one node, and the other intron is associated with the other node.

        intron_alignments: A set of tuples (rname, True if sense strand is
            forward strand or False if sense strand is reverse strand,
            alignment_start_position, alignment_end_position, tuple of tuples
            (pos, end_pos) of start and end positions of introns, readlet_size,
            distance to previous intron or None if beginning of strand,
            distance to next intron or None if end of strand, True if alignment
            is to forward strand else False, displacement of readlet from 5'
            end of read)
        
        Yield value: tuple (node A, [list of nodes that connect to A])
    """
    for compared_alignment in intron_alignments:
        yield (compared_alignment, [alignment 
                                    for alignment in intron_alignments
                                    if not different_introns_overlap(
                                                    compared_alignment[4],
                                                    alignment[4]
                                                )])

def maximal_cliques(intron_alignments):
    """ Finds maximal cliques of graph of intron combinations.

        Consider an undirected graph where each node corresponds to
        an alignment overlapping at least one cointron, as specified by a
        tuple from the set intron_alignments. Place an edge between two nodes
        for which no two _distinct_ introns overlap, where one intron is
        associated with one node, and the other intron is associated with the
        other node. Now enumerate maximal cliques. This gives all possible
        valid clusters of _consistent_ alignments, and each cluster subsumes
        a set of introns the read could possibly overlap.

        This code is adapted from NetworkX's find_cliques(), which requires
        inclusion of the following copyright notice.

        --------
        Copyright (C) 2004-2012, NetworkX Developers
        Aric Hagberg <hagberg@lanl.gov>
        Dan Schult <dschult@colgate.edu>
        Pieter Swart <swart@lanl.gov>
        All rights reserved.

        Redistribution and use in source and binary forms, with or without
        modification, are permitted provided that the following conditions are
        met:

          * Redistributions of source code must retain the above copyright
            notice, this list of conditions and the following disclaimer.

          * Redistributions in binary form must reproduce the above
            copyright notice, this list of conditions and the following
            disclaimer in the documentation and/or other materials provided
            with the distribution.

          * Neither the name of the NetworkX Developers nor the names of its
            contributors may be used to endorse or promote products derived
            from this software without specific prior written permission.


        THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
        "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
        LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
        A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
        OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
        SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
        LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
        DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
        THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
        (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
        OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
        --------

        intron_alignments: A set of tuples (rname, True if sense strand is
            forward strand or False if sense strand is reverse strand,
            alignment_start_position, alignment_end_position, tuple of tuples
            (pos, end_pos) of start and end positions of introns, readlet_size,
            distance to previous intron or None if beginning of strand,
            distance to next intron or None if end of strand, True if alignment
            is to forward strand else False, displacement of readlet from 5'
            end of read)
        
        Yield value: A maximal clique -- a list of intron alignments.
    """
    # Cache nbrs and find first pivot (highest degree)
    maxconn=-1
    nnbrs={}
    pivotnbrs=set() # handle empty graph
    for n,nbrs in cointron_adjacencies(intron_alignments):
        nbrs=set(nbrs)
        nbrs.discard(n)
        conn = len(nbrs)
        if conn > maxconn:
            nnbrs[n] = pivotnbrs = nbrs
            maxconn = conn
        else:
            nnbrs[n] = nbrs
    # Initial setup
    cand=set(nnbrs)
    smallcand = set(cand - pivotnbrs)
    done=set()
    stack=[]
    clique_so_far=[]
    # Start main loop
    while smallcand or stack:
        try:
            # Any nodes left to check?
            n=smallcand.pop()
        except KeyError:
            # back out clique_so_far
            cand,done,smallcand = stack.pop()
            clique_so_far.pop()
            continue
        # Add next node to clique
        clique_so_far.append(n)
        cand.remove(n)
        done.add(n)
        nn=nnbrs[n]
        new_cand = cand & nn
        new_done = done & nn
        # check if we have more to search
        if not new_cand:
            if not new_done:
                # Found a clique!
                yield clique_so_far[:]
            clique_so_far.pop()
            continue
        # Shortcut--only one node left!
        if not new_done and len(new_cand)==1:
            yield clique_so_far + list(new_cand)
            clique_so_far.pop()
            continue
        # find pivot node (max connected in cand)
        # look in done nodes first
        numb_cand=len(new_cand)
        maxconndone=-1
        for n in new_done:
            cn = new_cand & nnbrs[n]
            conn=len(cn)
            if conn > maxconndone:
                pivotdonenbrs=cn
                maxconndone=conn
                if maxconndone==numb_cand:
                    break
        # Shortcut--this part of tree already searched
        if maxconndone == numb_cand:
            clique_so_far.pop()
            continue
        # still finding pivot node
        # look in cand nodes second
        maxconn=-1
        for n in new_cand:
            cn = new_cand & nnbrs[n]
            conn=len(cn)
            if conn > maxconn:
                pivotnbrs=cn
                maxconn=conn
                if maxconn == numb_cand-1:
                    break
        # pivot node is max connected in cand from done or cand
        if maxconndone > maxconn:
            pivotnbrs = pivotdonenbrs
        # save search status for later backout
        stack.append( (cand, done, smallcand) )
        cand=new_cand
        done=new_done
        smallcand = cand - pivotnbrs

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
    unclustered_intron_alignments = [j for j, alignment
                                        in enumerate(alignments)
                                        if alignment[4] is not None]
    if not unclustered_intron_alignments:
        return []
    unclustered_alignments = range(len(alignments))
    clustered_alignments = []
    while unclustered_intron_alignments:
        pivot = i = random.choice(unclustered_intron_alignments)
        new_unclustered_alignments = []
        new_unclustered_intron_alignments = []
        alignment_cluster = defaultdict(list)
        for j in unclustered_alignments:
            if j == pivot: continue
            pivot_rname, compared_rname = alignments[i][0], alignments[j][0]
            pivot_sense, compared_sense = alignments[i][1], alignments[j][1]
            pivot_start, compared_start = alignments[i][2], alignments[j][2]
            pivot_end, compared_end = alignments[i][3], alignments[j][3]
            pivot_introns, compared_introns \
                = alignments[i][4], alignments[j][4]
            pivot_sign, compared_sign = alignments[i][8], alignments[j][8]
            pivot_displacement, compared_displacement \
                = alignments[i][9], alignments[j][9]
            pivot_group, compared_group = alignments[i][10], alignments[j][10]
            if (pivot_group != compared_group and
                pivot_rname == compared_rname and
                (pivot_sense is None or compared_sense is None or
                    pivot_sense == compared_sense) and
                ((pivot_start == compared_start
                    and pivot_displacement == compared_displacement) or
                 (pivot_start < compared_start)
                    == (pivot_displacement < compared_displacement)) and
                pivot_sign == compared_sign):
                alignment_cluster[compared_group].append(j)
            else:
                new_unclustered_alignments.append(j)
                if alignments[j][4] is not None:
                    new_unclustered_intron_alignments.append(j)
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
        unclustered_intron_alignments = new_unclustered_intron_alignments
    cluster_sizes = [len(set([alignment[-1] for alignment in cluster]))
                        for cluster in clustered_alignments]
    largest_cluster_size = max(cluster_sizes)
    largest_clusters = []
    for i, cluster_size in enumerate(cluster_sizes):
        if cluster_size == largest_cluster_size:
            largest_clusters.append(clustered_alignments[i])
    multimap_count = len(largest_clusters)
    return [set([alignment[:-1] for alignment in largest_clusters[i]
                if alignment[1] is not None])
            for i in xrange(multimap_count)]

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
    for (seq_id,), xpartition in xstream(input_stream, 1):
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
        if collected_readlets:
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
            if clusters:
                for selected_introns in clusters:
                    for alignments in maximal_cliques(selected_introns):
                        # Get stats on alignment with smallest start position
                        (_, _, left_pos,
                            _, _, left_readlet_size,
                            left_intron_distance, _,
                            _, left_displacement) = min(alignments,
                                key=lambda alignment: alignment[2]
                            )
                        # Get stats on alignment with largest end position
                        (_, _, _,
                            right_end_pos, _, right_readlet_size,
                            right_intron_distance, _,
                            _, right_displacement) = max(alignments,
                                key=lambda alignment: alignment[3]
                            )
                        introns_to_add = set()
                        for alignment in alignments:
                            for intron in alignment[4]:
                                introns_to_add.add(intron)
                        intron_starts_and_ends = zip(*sorted(
                                                          list(introns_to_add))
                                                        )
                        # Find start position of first intron
                        left_intron_pos = intron_starts_and_ends[0][0]
                        # Find end position of last intron
                        right_intron_end_pos = intron_starts_and_ends[1][-1]
                        '''Compute distances to ends of read from first and
                        last introns.'''
                        left_extend_size = (left_intron_pos - left_pos
                                             + left_displacement)
                        right_extend_size = (seq_size - right_displacement
                                                - right_readlet_size
                                                + right_end_pos
                                                - right_intron_end_pos)
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
                        if left_extend_size > left_intron_distance:
                            left_size = left_extend_size
                        else:
                            left_size = min(left_extend_size + fudge,
                                            left_intron_distance)
                        if right_extend_size > right_intron_distance:
                            right_size = right_extend_size
                        else:
                            right_size = min(right_extend_size + fudge,
                                             right_intron_distance)
                        print >>output_stream, '%s\t%d\t%s\t%s' \
                                '\t%d\t%d\t%s' % (alignments[0][0] + ('-' if 
                                                    alignments[0][1] else '+'),
                                left_intron_pos,
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
    import unittest
    # Test units
    del sys.argv[1:] # Don't choke on extra command-line parameters
    class TestRnameAndIntrons(unittest.TestCase):
        """ Tests rname_and_introns(); needs no fixture. """
        
        def test_read_1(self):
            """ Fails if example doesn't give expected introns."""
            self.assertEquals(('chr10', True, 
                                101478239, 101480750,
                                ((101478258, 101480744),), 25,
                                155, 82),
                            rname_and_introns(
                                'chr10-;101478234;24,24;2486;155;82',
                                5,
                                25
                            )
                    )

        def test_read_2(self):
            """ Fails if example doesn't give expected introns."""
            self.assertEquals(('chr6', True, 
                                161054930, 161056165,
                                ((161054945, 161056155),), 25,
                                16, 11273),
                        rname_and_introns(
                            'chr6-;161049359;24,16,24;5546,1210;21874;11273',
                            25,
                            25
                        )
                    )

    unittest.main()
