#!/usr/bin/env python
"""
Rail-RNA-intron_search

Follows Rail-RNA-align_readlets
Precedes Rail-RNA-intron_call

Reduce step in MapReduce pipelines that -- very roughly -- infers splice
junctions between successive readlets that align noncontiguously.

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
Tab-delimited columns:
1. Reference name (RNAME in SAM format) + ';' + partition number +  
    ('+' or '-' indicating which strand is the sense strand if input reads are
        strand-specific -- that is, --stranded was invoked; otherwise, there is
        no terminal '+' or '-')
2. Candidate intron start (inclusive) on forward strand (1-indexed)
3. Candidate intron end (exclusive) on forward strand (1-indexed)
4. '\x1f'-separated list of sample (label)s in which intron was detected
5. Total number of reads supporting intron

ALL OUTPUT COORDINATES ARE 1-INDEXED.
"""
import sys
import os
import site
import threading
import string
import subprocess
import re
import random
import itertools
from collections import defaultdict

base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for directory_name in ['bowtie', 'sample', 'interval']:
    site.addsitedir(os.path.join(base_path, directory_name))

import bowtie
import bowtie_index
import sample
import partition
from collections import deque

_reversed_complement_translation_table = string.maketrans('ATCG', 'TAGC')

# Initialize global variables for tracking number of input/output lines
_input_line_count = 0
_output_line_count = 0

# Initialize forward- and reverse-strand motifs

_prioritized_reverse_strand_motifs = [('CT', 'AC'), ('CT', 'GC'), ('GT', 'AT')]
_reverse_strand_motifs = set(_prioritized_reverse_strand_motifs)
_prioritized_forward_strand_motifs = [('GT', 'AG'), ('GC', 'AG'), ('AT', 'AC')]
_forward_strand_motifs = set(_prioritized_forward_strand_motifs)
_all_motifs = _reverse_strand_motifs | _forward_strand_motifs

def maximal_suffix_match(query_seq, search_window,
                            min_cap_size=8, max_cap_count=5):
    """ Finds maximum matching suffix of query_seq closest to start of window.

        query_seq: sequence to search for.
        search_window: sequence to search in.
        min_cap_size: minimum size of sequence to search for; this is used to
            enumerate initial possible matches for query_seq (that is, matches
            of the min_cap_size-length suffix of query_seq)
        max_cap_count: initial list of initial possible matches is barred from
            exceeding this value.

        Return value: tuple (offset of match from beginning of search_window, 
                                length of maximum matching suffix of query_seq)
                      for the longest suffix CLOSEST to the start of the window
                      or None if (no suffix found or max_cap_count exceeded
                        or last min_cap_size characters of suffix match
                        last min_cap_size characters of first len(query_seq)
                        characters of search_window)
    """
    query_seq_size = len(query_seq)
    offset = query_seq_size - min_cap_size
    suffix_seq = query_seq[offset:]
    suffixes = []
    first_search = True
    while len(suffixes) <= max_cap_count:
        suffix = re.search(suffix_seq, search_window[offset:])
        if suffix is None:
            break
        elif first_search and not suffix.start():
            '''Suffix was found on first search at very beginning of window;
            it's VERY likely just an extension.'''
            return None
        else:
            extra_base_count = 0
            offset += suffix.start()
            while min_cap_size + extra_base_count < query_seq_size:
                if query_seq[-min_cap_size-extra_base_count-1] \
                    == search_window[offset-extra_base_count-1]:
                    extra_base_count += 1
                else:
                    break
            suffixes.append((-(extra_base_count + min_cap_size),
                                offset - extra_base_count))
            offset += 1
        first_search = False
    try:
        if len(suffixes) <= max_cap_count:
            suffix = min(suffixes)
            return (suffix[1], -suffix[0])
        else:
            return None
    except ValueError:
        # No suffixes found
        return None

def selected_readlet_alignments_by_clustering(readlets):
    """ Selects multireadlet alignment via a correlation clustering algorithm.
    
        Consider a list "readlets" whose items {R_i} correspond to the aligned
        readlets from a given read. Each R_i is itself a list of the possible
        alignments {R_ij} of a readlet. Each R_ij is a tuple
        (rname, reverse_strand, pos, end_pos, displacement), where rname is the
        SAM-format rname (typically a chromosome), reverse_strand is True iff
        the readlet's reversed complement aligns to the reference, and 
        displacement is the number of bases between the 5' (3')end of the
        readlet, which aligns to the forward (reverse) strand, and the 5' (3')
        end of the read. Let K_i be the number of alignments {R_ij} of a given
        readlet R_i.

        The algorithm first constructs a complete signed graph where each node
        represents a different R_ij. A '+' edge is placed between each pair
        of nodes (R_ij, R_kl) iff R_ij and R_kl are alignments to the same
        strand, i != k (that is, if the alignments don't correspond to the same
        multireadlet), and R_ij and R_kl are "order-distance-consistent." Order
        -distance consistency demands the following:
            1) If the displacement of R_kl is greater than the displacement of
            R_ij along the read, R_kl must occur at a position greater than
            R_ij along the reference, and vice versa.
            2) If the displacements of R_kl and R_ij are the same, the
            positions of R_kl and R_ij along the reference must also be the
            same.
            3) There are no R_ip or R_kp for any p between the start positions
            of R_ij and R_kl along the reference.
        Otherwise, there is a '-' edge between R_ij and R_kl. The stipulation
        3) guarantees that an alignment never has more than one '+' edge
        between it and any set of alignments corresponding to the same
        multireadlet. For multireadlets R_i and R_k, order-distance consistency
        associates the closest pair (R_ij, R_kl) (in reference space) for which
        the orders of R_ij and R_kl in both read and reference space agree;
        if no such pair exists, no association is made.

        After the graph is constructed, correlation clustering is performed
        using QuickClust, a simple 3-approximation algorithm. (See Ailon,
        Nir, and Charikar. "Aggregating inconsistent information: ranking and
        clustering." STOC 2005: Proceedings of the thirty-seventh annual
        ACM symposium on Theory of Computing. pp. 684-693.) At each iteration,
        a pivot node is chosen from unclustered nodes at random. The nodes to
        which this pivot is connected by a '+' edge are placed in the same
        cluster as the pivot, while the remaining nodes are considered
        unclustered. Because a given pivot is associated with no more than one
        alignment from the same multireadlet, each cluster contains at most one
        alignment from the same multireadlet. Alignments from the largest
        cluster are returned. If more than one cluster is largest, no
        alignments are returned; this scenario is analogous to having full
        reads map to multiple locations in the genome.

        readlets: a list whose items {R_i} correspond to the aligned readlets
            from a given read. Each R_i is itself a list of the possible
            alignments {R_ij} of a readlet. Each R_ij is a tuple
            (rname, reverse_strand, pos, end_pos, displacement).
            See above for a detailed explanation.

        Return value: a list of selected alignment tuples
            (rname, reverse_strand, pos, end_pos, displacement).
    """
    alignments = [alignment + (i,) for i, multireadlet
                                        in enumerate(readlets)
                                        for alignment in multireadlet]
    '''Sort alignments so closest multireadlet from given group can be found
    easily.'''
    alignments.sort()
    multireadlet_groups = [alignment[-1] for alignment in alignments]
    unclustered_alignments = range(len(alignments))
    clustered_alignments = []
    while unclustered_alignments:
        pivot = i = random.choice(unclustered_alignments)
        new_unclustered_alignments = []
        alignment_cluster = [pivot]
        for j in unclustered_alignments:
            if j == pivot: continue
            intervening_multireadlets = (multireadlet_groups[i+1:j] if i < j 
                                            else multireadlet_groups[j+1:i])
            if not (alignments[i][-1] == alignments[j][-1] or \
                alignments[i][:2] != alignments[j][:2] or \
                (alignments[i][4] == alignments[j][4] and
                   alignments[i][2] != alignments[j][2]) or \
                ((alignments[i][4] < alignments[j][4]) != \
                    (alignments[i][2] < alignments[j][2])) or \
                (alignments[i][-1] in intervening_multireadlets or
                    alignments[j][-1] in intervening_multireadlets)):
                alignment_cluster.append(j)
            else:
                new_unclustered_alignments.append(j)
        clustered_alignments.append(
                [alignments[j][:-1] for j in alignment_cluster]
            )
        unclustered_alignments = new_unclustered_alignments
    largest_cluster_size = max(map(len, clustered_alignments))
    largest_clusters = [cluster for cluster in clustered_alignments
                            if len(cluster) == largest_cluster_size]
    if len(largest_clusters) == 1:
        return largest_clusters[0]
    else:
        return []

def alignment_adjacencies(alignments):
    """ Generates adjacency matrix for graph described below.

        Consider an undirected graph where each node corresponds to
        an alignment of a distinct multireadlet. Place an edge between two
        nodes that are "order-consistent"; that is:

        1) They have the same strand label
        2) If the displacements of the corresponding readlets from the 5' end
        of the read are equal, their genomic positions on the reference are
        also equal.
        3) If the displacement of node A is greater than the displacement of
        node B, the genomic position of node A is greater than the genomic
        position of node B.

        cluster: a list of alignment tuples
            (rname, reverse_strand, pos, end_pos, displacement), each
            corresponding to a distinct readlet
        
        Yield value: tuple (node A, [list of nodes that connect to A])
    """
    for compared_alignment in alignments:
        yield (compared_alignment, [alignment 
                                    for alignment in alignments
                                    if (compared_alignment[:2] == alignment[:2]
                                    and (compared_alignment[2]
                                            == alignment[2]
                                            if compared_alignment[4]
                                            == alignment[4] else
                                            ((compared_alignment[2]
                                                < alignment[2])
                                            == (compared_alignment[4]
                                                < alignment[4]))))])

def largest_maximal_clique(cluster):
    """ Finds maximal cliques of graph of multireadlet alignment cluster.

        Consider an undirected graph where each node corresponds to
        an alignment of a distinct multireadlet. Place an edge between two
        nodes that are "order-consistent"; that is:

        1) They have the same strand label
        2) If the displacements of the corresponding readlets from the 5' end
        of the read are equal, their genomic positions on the reference are
        also equal.
        3) If the displacement of node A is greater than the displacement of
        node B, the genomic position of node A is greater than the genomic
        position of node B.

        Now enumerate maximal cliques. This gives all possible maximal groups
        of mutually consistent alignments. Return the largest group. If there's
        a tie, break it at random.

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

        cluster: a list of alignment tuples
            (rname, reverse_strand, pos, end_pos, displacement), each
            corresponding to a distinct readlet
        
        Return value: largest maximal clique -- a list of alignments.
    """
    clique_time = time.time()
    cliques = []
    # Cache nbrs and find first pivot (highest degree)
    maxconn=-1
    nnbrs={}
    pivotnbrs=set() # handle empty graph
    for n,nbrs in alignment_adjacencies(cluster):
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
                cliques.append(clique_so_far[:])
            clique_so_far.pop()
            continue
        # Shortcut--only one node left!
        if not new_done and len(new_cand)==1:
            cliques.append(clique_so_far + list(new_cand))
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
    try:
        max_clique_size = max(map(len, cliques))
    except ValueError:
        assert not cluster
        return []
    largest_cliques = [clique for clique in cliques
                        if len(clique) == max_clique_size]
    largest_clique_count = len(largest_cliques)
    assert largest_clique_count >= 1
    if largest_clique_count == 1:
        return largest_cliques[0]
    return random.choice(largest_cliques)

def pairwise(iterable):
    """ Iterates through iterable in pairs.

        If iterable's items are [a1, a2, a3, a4, ...], yields tuples (a1, a2),
        (a2, a3), (a3, a4), .... See
        https://docs.python.org/2/library/itertools.html.

        Return value: generator for pairs
    """
    left, right = itertools.tee(iterable)
    next(right, None)
    return itertools.izip(left, right)

def introns_from_clique(clique, read_seq, reference_index,
        min_exon_size=8, min_intron_size=10, max_intron_size=500000,
        search_window_size=1000, stranded=False, motif_radius=1,
        reverse_reverse_strand=False):
    """ 
        NOTE THAT clique LIST IS SORTED ASSUMING THE ONLY READLETS WHOSE
        DISPLACEMENTS ARE THE SAME ARE CAPPING READLETS. IF THE READLETIZING
        SCHEME IS CHANGED, THIS SORTING SCHEME SHOULD ALSO BE CHANGED.

        clique: list of alignment tuples
            (rname, reverse_strand, pos, end_pos, displacement), all of which
            are order-consistent
        read_seq: sequence of the original read from which the readlets are
            derived.
        reference_index: object of class bowtie_index.BowtieIndexReference that
            permits access to reference; used for realignment of unmapped
            regions between exonic chunks.
        min_exon_size: minimum size of exons to search for
        min_intron_size: do not return introns < this size
        max_intron_size: do not return introns > this size
        search_window_size: size of search window flanking an alignment
            in which to search for an exon
        stranded: True iff input reads are strand-specific
        motif_radius: radius (in bp) around unmapped region in which to search
            for motif; keep this small!
        reverse_reverse_strand: if True, original read sequence was
            reverse-complemented before alignment of constituent readlets
    """
    read_seq = read_seq.upper()
    if not clique:
        return []
    rname, reverse_strand = clique[0][:2]
    if reverse_reverse_strand:
        reverse_strand = not reverse_strand
    candidate_introns = []
    # Arrange sort so that longest aligning capping readlet is used as flank
    clique.sort(key=lambda alignment: (alignment[-1], alignment[-2] if
                                                        alignment[-1] == 0
                                                        else alignment[-3]))
    '''Call introns overlapped by read using constraint that all are on the
    same strand. First call introns between mapped regions by iterating 
    through flanking alignments.'''
    for ((_, _, left_pos, left_end_pos, left_displacement),
            (_, _, right_pos, right_end_pos, right_displacement)) \
        in pairwise(clique):
        '''Assume that flanking alignments provide the exact _size_ of an
        intron; this happens if there are no intervening exons.'''
        read_span = right_displacement + right_end_pos - right_pos \
                        - left_displacement
        intron_size = right_end_pos - left_pos - read_span
        unmapped_start = left_displacement + left_end_pos - left_pos
        unmapped_end = right_displacement
        if unmapped_start > unmapped_end:
            '''Alignments overlap on read; turn the unmapped region into
            the overlap.'''
            unmapped_start, unmapped_end = unmapped_end, unmapped_start
            unmapped_size = unmapped_end - unmapped_start
            right_pos += unmapped_size
            left_end_pos -= unmapped_size
        else:
            unmapped_size = unmapped_end - unmapped_start
        # Expand unmapped region
        right_pos += motif_radius
        left_end_pos -= motif_radius
        unmapped_size += 2 * motif_radius
        unmapped_start -= motif_radius
        unmapped_end += motif_radius
        if min_intron_size <= intron_size <= max_intron_size:
            candidate_introns.append(
                    (left_end_pos, 
                        intron_size, 
                        unmapped_size,
                        unmapped_start,
                        unmapped_end)
                )
    # Motifs must agree on strand
    introns = []
    if stranded:
        if reverse_strand:
            search_motifs = _reverse_strand_motifs
        else:
            search_motifs = _forward_strand_motifs
    else:
        search_motifs = _all_motifs
    for i, (pos, intron_size, unmapped_size, unmapped_start, 
            unmapped_end) in enumerate(candidate_introns):    
        left_motif_window = reference_index.get_stretch(
                                    rname,
                                    pos - 1,
                                    unmapped_size
                                )
        right_motif_window = reference_index.get_stretch(
                                    rname,
                                    pos - 1 + intron_size - 2,
                                    unmapped_size
                                )
        found_motifs = []
        for j in xrange(unmapped_size - 1):
            candidate_motif = (left_motif_window[j:j+2],
                                right_motif_window[j:j+2])
            if candidate_motif in search_motifs:
                introns.append(
                        (i, pos + j,
                         pos + intron_size + j,
                         intron_size,
                         unmapped_start,
                         unmapped_end,
                         candidate_motif,
                         )
                    )
    # Enumerate groups of introns that agree on strand
    forward_strand_introns = [intron for intron in introns if intron[-1]
                                    in _forward_strand_motifs]
    reverse_strand_introns = [intron for intron in introns if intron[-1]
                                    in _reverse_strand_motifs]
    forward_strand_intron_set = set([intron[0] for intron
                                        in forward_strand_introns])
    reverse_strand_intron_set = set([intron[0] for intron
                                        in reverse_strand_introns])
    forward_strand_intron_count = len(forward_strand_intron_set)
    reverse_strand_intron_count = len(reverse_strand_intron_set)
    filtered_introns = []
    if forward_strand_intron_count != 0 \
        and forward_strand_intron_count == reverse_strand_intron_count:
        '''Try to break tie with priority class. If this can't be done,
        be conservative and don't call any introns.'''
        forward_strand_priority_sum = 0
        final_forward_strand_introns = []
        for i in forward_strand_intron_set:
            candidate_motifs = [intron for intron in forward_strand_introns
                                if i == intron[0]]
            if len(candidate_motifs) >= 1:
                priorities \
                    = [_prioritized_forward_strand_motifs.index(motif[-1])
                        for motif in candidate_motifs]
                min_priority = min(priorities)
                priority_indices = [j for j, priority in enumerate(priorities)
                                        if priority == min_priority]
                if len(priority_indices) == 1:
                    final_forward_strand_introns.append(candidate_motifs[
                                                          priority_indices[0]
                                                        ][1:] + (False,))
                forward_strand_priority_sum += min_priority
        reverse_strand_priority_sum = 0
        final_reverse_strand_introns = []
        for i in reverse_strand_intron_set:
            candidate_motifs = [intron for intron in reverse_strand_introns
                                if i == intron[0]]
            if len(candidate_motifs) >= 1:
                priorities \
                    = [_prioritized_reverse_strand_motifs.index(motif[-1])
                        for motif in candidate_motifs]
                min_priority = min(priorities)
                priority_indices = [j for j, priority in enumerate(priorities)
                                        if priority == min_priority]
                if len(priority_indices) == 1:
                    final_reverse_strand_introns.append(candidate_motifs[
                                                          priority_indices[0]
                                                        ][1:] + (True,))
                reverse_strand_priority_sum += min_priority
        if reverse_strand_priority_sum > forward_strand_priority_sum:
            filtered_introns = final_forward_strand_introns
        elif reverse_strand_priority_sum < forward_strand_priority_sum:
            filtered_introns = final_reverse_strand_introns
    elif forward_strand_intron_count > reverse_strand_intron_count:
        for i in forward_strand_intron_set:
            candidate_motifs = [intron for intron in forward_strand_introns
                                if i == intron[0]]
            if len(candidate_motifs) >= 1:
                priorities \
                    = [_prioritized_forward_strand_motifs.index(motif[-1])
                        for motif in candidate_motifs]
                min_priority = min(priorities)
                priority_indices = [j for j, priority in enumerate(priorities)
                                        if priority == min_priority]
                if len(priority_indices) == 1:
                    filtered_introns.append(candidate_motifs[
                                                priority_indices[0]
                                            ][1:] + (False,))
    elif reverse_strand_intron_count > forward_strand_intron_count:
        for i in reverse_strand_intron_set:
            candidate_motifs = [intron for intron in reverse_strand_introns
                                if i == intron[0]]
            if len(candidate_motifs) >= 1:
                priorities \
                    = [_prioritized_reverse_strand_motifs.index(motif[-1])
                        for motif in candidate_motifs]
                min_priority = min(priorities)
                priority_indices = [j for j, priority in enumerate(priorities)
                                        if priority == min_priority]
                if len(priority_indices) == 1:
                    filtered_introns.append(candidate_motifs[
                                                priority_indices[0]
                                            ][1:] + (True,))
    final_introns = deque()
    # Search for small exons that could split introns
    for (pos, end_pos, intron_size, unmapped_start, unmapped_end,
            intron_motif, intron_reverse_strand) in filtered_introns:
        final_introns.append((rname, intron_reverse_strand, pos, end_pos))
        '''if unmapped_end - unmapped_start >= min_exon_size:
            # Search for just _one_ intervening exon
            if intron_reverse_strand:
                left_caps = [motif[1] for motif in _reverse_strand_motifs
                                if motif[0] == intron_motif[0]]
                right_caps = [motif[0] for motif in _reverse_strand_motifs
                                if motif[1] == intron_motif[1]]
            else:
                left_caps = [motif[1] for motif in _forward_strand_motifs
                                if motif[0] == intron_motif[0]]
                right_caps = [motif[0] for motif in _forward_strand_motifs
                                if motif[1] == intron_motif[1]]
            caps = itertools.product(left_caps, right_caps)
            search_strings = [cap[0] + read_seq[unmapped_start:unmapped_end]
                                + cap[1] for cap in caps]
            if intron_size <= search_window_size:
                search_displacement = 0
                search_window = reference_index.get_stretch(
                        rname,
                        pos - 1,
                        intron_size
                    )
            else:
                # Randomly select start position of search window 
                search_displacement = random.randint(
                                            0, 
                                            intron_size - search_window_size
                                        )
                search_window = reference_index.get_stretch(
                        rname,
                        pos - 1 + search_displacement,
                        search_window_size
                    )
            hits = []
            for search_string in search_strings:
                hits.extend(list(re.findall('(?=(%s))' % search_string, 
                                    search_window)))
            hit_count = len(hits)
            if hit_count == 0:
                # String not found; call single intron
                final_introns.append((rname, intron_reverse_strand,
                                            pos, end_pos))
            elif hit_count == 1:
                # Call two introns
                final_introns.extend([(rname, intron_reverse_strand,
                                            pos, hits[0].start()),
                                           (rname, intron_reverse_strand,
                                            hits[0].end(), end_pos)])
            else:
                # At least two pairs of possible introns, so add nothing
                pass'''
    # Search for prefix and suffix exons
    return final_introns

def go(input_stream=sys.stdin, output_stream=sys.stdout,
    bowtie_index_base='genome', verbose=False, stranded=False,
    min_exon_size=8, min_intron_size=15, max_intron_size=500000,
    motif_radius=1, search_window_size=1000):
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
        bin_size: genome is partitioned in units of bin_size for later load
            balancing.
        bowtie_index_base: the basename of the Bowtie index files associated
            with the reference.
        verbose: True iff more informative messages should be written to
            stderr.
        stranded: True iff input reads are strand-specific; this affects
            whether an output partition has a terminal '+' or '-' indicating
            the sense strand.
        max_discrepancy: if the difference in length between an unmapped region
            framed by two ECs and its corresponding gap in the reference is <=
            this value, the unmapped region is considered a candidate for
            incorporation into a single EC spanning the two original ECs via
            DP filling.
        min_seq_similarity: if the difference in length between an unmapped
            region framed by two ECs and its corresponding gap in the reference
            is <= max_discrepancy AND the score of global alignment is >=
            min_seq_similarity * (length of unmapped region), the unmapped
            region is incorporated into a single EC spanning the two original
            ECs via DP filling. See the GlobalAlignment class for the
            substitution matrix used.
        min_intron_size: introns smaller than this number of bases are filtered
            out.
        max_intron_size: an intron of that spans more than this number of bases 
            is suppressed from output.
        intron_partition_overlap: number of bases to subtract from
            reference start position of intron when determining genome
            partition it is in.
        search_for_caps: True iff reference should be searched for the segment
            of a read (a cap) that precedes the first EC and the segment of a
            read that follows the last EC. Such segments are subsequently added
            as an ECs themselves, and introns may be called between them.
        min_cap_size: the reference is not searched for a cap smaller
            than this size.
        cap_search_window_size: the size (in bp) of the reference subsequence
            in which to search for a cap.
        max_cap_count: maximum number of possible caps of size
            min_cap_size to consider when searching for caps.
        global_alignment: instance of GlobalAlignment class used for fast
                realignment of exonic chunks via Weave.
        report_multiplier: if verbose is True, the line number of an alignment,
            read, or first readlet of a read written to stderr increases
            exponentially with base report_multiplier.

        No return value.
    """
    global _input_line_count, _output_line_count
    reference_index = bowtie_index.BowtieIndexReference(bowtie_index_base)
    '''Input is readletized, and readlets must be composed.'''
    readlet_displacements, collected_readlets = [], []
    seq_info_captured = False
    line = input_stream.readline()
    if not line:
        # Empty input
        return
    tokens = line.rstrip().split('\t')
    assert len(tokens) == 5
    seq_id = tokens[0]
    while True:
        seq_info = tokens[1].split('\x1e')
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
        if tokens[-1] != '\x1c':
            # If there are alignments, add them to collected_readlets
            collected_readlets.append(
                    zip(
                            tokens[-3].split('\x1f'),
                            [(int(flag) & 16) != 0
                                for flag in tokens[-2].split('\x1f')],
                            [int(pos) for pos in tokens[-1].split('\x1f')]
                        )
                )
            readlet_displacements.append(
                (int(seq_info[0]), int(seq_info[1]))
            )
        line = input_stream.readline()
        if line:
            tokens = line.rstrip().split('\t')
            assert len(tokens) == 5
            _input_line_count += 1
            next_seq_id = tokens[0]
        if (not line or seq_id != next_seq_id) \
            and len(collected_readlets):
            assert seq_info_captured, \
                'Sequence info was not in a collected readlet'
            multireadlets = [[(rname, reverse_strand, pos, pos + seq_size
                                - readlet_displacements[i][1]
                                - readlet_displacements[i][0],
                                readlet_displacements[i][1]
                                    if reverse_strand 
                                    else readlet_displacements[i][0]) 
                                for rname, reverse_strand, pos in multireadlet]
                                for i, multireadlet
                                in enumerate(collected_readlets)]
            # Set seed for each read so results for read are reproducible
            random.seed(seq)
            sample_labels = set((sample_labels.split('\x1f')
                                    if len(sample_labels) else []))
            reversed_complement_sample_labels \
                = set((reversed_complement_sample_labels.split('\x1f')
                        if len(reversed_complement_sample_labels) else []))
            if stranded:
                if sample_labels:
                    introns = introns_from_clique(
                            largest_maximal_clique(
                                    selected_readlet_alignments_by_clustering(
                                            multireadlets
                                        )
                                ),
                            seq,
                            reference_index,
                            min_exon_size=min_exon_size,
                            min_intron_size=min_intron_size,
                            max_intron_size=max_intron_size,
                            search_window_size=search_window_size,
                            stranded=stranded,
                            motif_radius=motif_radius,
                            reverse_reverse_strand=False
                        )
                    for sample_label in sample_labels:
                        for (intron_rname, intron_reverse_strand,
                             intron_pos, intron_end_pos) in introns:
                            print '%s%s\t%s\t%012d\t%012d' % (
                                    intron_rname,
                                    '-' if intron_reverse_strand else '+',
                                    sample_label,
                                    intron_pos,
                                    intron_end_pos
                                )
                if reversed_complement_sample_labels:
                    introns = introns_from_clique(
                            largest_maximal_clique(
                                    selected_readlet_alignments_by_clustering(
                                            multireadlets
                                        )
                                ),
                            seq,
                            reference_index,
                            min_exon_size=min_exon_size,
                            min_intron_size=min_intron_size,
                            max_intron_size=max_intron_size,
                            search_window_size=search_window_size,
                            stranded=stranded,
                            motif_radius=motif_radius,
                            reverse_reverse_strand=True
                        )
                    for sample_label in reversed_complement_sample_labels:
                        for (intron_rname, intron_reverse_strand,
                             intron_pos, intron_end_pos) in introns:
                            print '%s%s\t%s\t%012d\t%012d' % (
                                    intron_rname,
                                    '-' if intron_reverse_strand else '+',
                                    sample_label,
                                    intron_pos,
                                    intron_end_pos
                                )
            else:
                introns = introns_from_clique(
                        largest_maximal_clique(
                                selected_readlet_alignments_by_clustering(
                                        multireadlets
                                    )
                            ),
                        seq,
                        reference_index,
                        min_exon_size=min_exon_size,
                        min_intron_size=min_intron_size,
                        max_intron_size=max_intron_size,
                        search_window_size=search_window_size,
                        motif_radius=motif_radius,
                        stranded=stranded
                    )
                for sample_label in (sample_labels
                                        | reversed_complement_sample_labels):
                    for (intron_rname, intron_reverse_strand,
                             intron_pos, intron_end_pos) in introns:
                        print '%s%s\t%s\t%012d\t%012d' % (
                                intron_rname,
                                '-' if intron_reverse_strand else '+',
                                sample_label,
                                intron_pos,
                                intron_end_pos
                            )
            seq_info_captured = False
            collected_readlets = []
            readlet_displacements = []
        if not line: break
        seq_id = next_seq_id

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
    parser.add_argument('--min-intron-size', type=int, required=False,
        default=5,
        help='Filters introns of length smaller than this value')
    parser.add_argument('--max-intron-size', type=int, required=False,
        default=500000, 
        help='Filters introns of length greater than this value')
    parser.add_argument('--min-exon-size', type=int, required=False,
        default=8,
        help='Minimum size of exons to search for')
    parser.add_argument('--search-window-size', type=int, required=False,
        default=1000,
        help='Size of window (in bp) in which to search for exons between '
             'anchoring alignments')
    parser.add_argument('--motif-radius', type=int, required=False,
        default=1,
        help='Number of bases to tack on to either end of an unmapped region '
             'when searching for motifs')

    # Add command-line arguments for dependencies
    partition.addArgs(parser)
    bowtie.addArgs(parser)

    '''Now collect arguments. While the variable args declared below is
    global, properties of args are also arguments of the go() function so
    different command-line arguments can be passed to it for unit tests.'''
    args = parser.parse_args(sys.argv[1:])

if __name__ == '__main__' and not args.test:
    import time
    start_time = time.time()
    go(bowtie_index_base=args.bowtie_idx,
        verbose=args.verbose, 
        stranded=args.stranded,
        min_intron_size=args.min_intron_size,
        max_intron_size=args.max_intron_size,
        min_exon_size=args.min_exon_size,
        motif_radius=args.motif_radius,
        search_window_size=args.search_window_size)
    print >> sys.stderr, 'DONE with intron_search.py; in/out=%d/%d; ' \
        'time=%0.3f s' % (_input_line_count, _output_line_count,
                                time.time() - start_time)
elif __name__ == '__main__':
    # Test units
    del sys.argv[1:] # Don't choke on extra command-line parameters
    import random
    import unittest
    import shutil
    import tempfile

    def random_sequence(seq_size):
        """ Gets random sequence of nucleotides.

            seq_size: number of bases to return.

            Return value: string of random nucleotides.
        """
        return ''.join([random.choice('ATCG') for _ in xrange(seq_size)])

    class TestComposedAndSortedReadlets(unittest.TestCase):
        """ Tests composed_and_sorted_readlets(); needs no fixture. """
        def test_overlapping_readlets_1(self):
            """ Fails if readlets are not consolidated as expected. """
            self.assertEqual(
                    composed_and_sorted_readlets([
                            ('chr1', True, 3, 20, 5), 
                            ('chr2', False, 2, 10, 4),
                            ('chr2', False, 3, 12, 5)
                        ]
                    ),
                    {
                        ('chr1', True): [(3, 20, 5)],
                        ('chr2', False): [(2, 12, 4)]
                    }
                )
        def test_overlapping_readlets_2(self):
            """ Fails if readlets are not consolidated as expected. """
            self.assertEqual(
                    composed_and_sorted_readlets([
                            ('chr1', True, 3, 20, 5), 
                            ('chr2', True, 2, 10, 4),
                            ('chr2', False, 2, 10, 4),
                            ('chr2', False, 3, 12, 5),
                            ('chr2', False, 20, 23, 10),
                            ('chr2', False, 21, 24, 11)
                        ]
                    ),
                    {
                        ('chr1', True): [(3, 20, 5)],
                        ('chr2', True): [(2, 10, 4)],
                        ('chr2', False): [(2, 12, 4), (20, 24, 10)]
                    }
                )
        def test_overlapping_readlets_and_sorting_1(self):
            """ Fails if readlets are not consolidated as expected. """
            self.assertEqual(
                    composed_and_sorted_readlets([
                            ('chr2', False, 3, 12, 5),
                            ('chr1', True, 3, 20, 5), 
                            ('chr2', False, 20, 23, 10),
                            ('chr2', True, 2, 10, 4),
                            ('chr2', False, 2, 10, 4),
                            ('chr2', False, 21, 24, 11)
                        ]
                    ),
                    {
                        ('chr2', False): [(2, 12, 4), (20, 24, 10)],
                        ('chr1', True): [(3, 20, 5)],
                        ('chr2', True): [(2, 10, 4)]
                    }
                )
        def test_overlapping_readlets_and_sorting_2(self):
            """ Fails if readlets are not consolidated as expected. """
            self.assertEqual(
                    composed_and_sorted_readlets([
                            ('chr3', False, 2, 42, 7),
                            ('chr1', True, 4, 7, 9),
                            ('chr1', True, 3, 5, 8),
                            ('chr3', False, 1, 51, 6),
                            ('chr3', False, 47, 59, 21),
                            ('chr3', False, 121, 145, 10),
                            ('chr4', False, 3, 6, 9),
                            ('chr2', False, 20, 23, 10),
                        ]
                    ),
                    {
                        ('chr1', True): [(3, 7, 8)],
                        ('chr2', False): [(20, 23, 10)],
                        ('chr3', False): [(1, 59, 6), (121, 145, 10)],
                        ('chr4', False): [(3, 6, 9)]
                    }
                )

    class TestUnmappedRegionSplits(unittest.TestCase):
        """ Tests unmapped_region_splits(); needs no fixture. """
        def test_1000_random_exact_cases(self):
            """ Fails if proper split of unmapped region is not identified.

                Test is run on 1000 random instances.
            """
            for i in xrange(1000):
                seq_size = random.randint(50, 150)
                left_reference_seq = random_sequence(seq_size)
                right_reference_seq = random_sequence(seq_size)
                split = random.randint(0, seq_size)
                unmapped_seq = left_reference_seq[:split] \
                                + right_reference_seq[split:]
                # split must be among the highest-scoring splits
                self.assertTrue(split in unmapped_region_splits(
                                            unmapped_seq,
                                            left_reference_seq,
                                            right_reference_seq
                                        )
                                    )
    
    class TestMaximalSuffixMatch(unittest.TestCase):
        """ Tests maximal_suffix_match(); needs no fixture. """
        def test_one_instance_1(self):
            """ Fails if maximal suffix match is not identified.
            """
            self.assertEqual(
                    maximal_suffix_match(
                            'ATAGCATTA', 'CAGTCAGACCCATACCAATAGCATTA'
                        ),
                    (17, 9)
                )

        def test_one_instance_2(self):
            """ Fails if maximal suffix match is not identified.
            """
            self.assertEqual(
                    maximal_suffix_match(
                            'CGATACGTCAGACCATG',
                            'ATGGCATACGATACGTCAGACCATGCAGGACCTTTACCTACATACTG'
                        ),
                    (8, 17)
                )

        def test_one_instance_3(self):
            """ Fails if maximal suffix match is not identified.
            """
            self.assertEqual(
                    maximal_suffix_match(
                            'CGATACGTCAGACCATG',
                            'ATGGCATAATACGTCAGACCATGCAGGACCTTTACCTACATACTG'
                        ),
                    (8, 15)
                )

        def test_filtering_of_more_than_max_cap_count_instances(self):
            """ Fails if maximal suffix matches are not filtered out.
            """
            self.assertEqual(
                    maximal_suffix_match(
                            'ATAGCATTA',
                            'CAGTCAGACCCATACCAATAGCATTAATAGCATTA',
                            max_cap_count=1
                        ),
                    None
                )

    class TestIntronsFromRead(unittest.TestCase):
        """ Tests introns_from_read(). """
        def setUp(self):
            """ Creates temporary directory and Bowtie index. """
            reference_seq = 'ATGGCATACGATACGTCAGACCATGCAggACctTTacCTACATACTG' \
                            'GCTACATAGTACATCTAGGCATACTACGTgcCATACGgaCTACGTAA' \
                            'ATCCAGATTACGATACAAaTACGAAcTCccATAGCAaCATaCTAGac' \
                            'CAttAaaGACTAGACTAACAGACAaAACTAGCATacGATCATGACaA' \
                            'ACGAGATCCATATAtTTAGCAaGACTAaACGATACGATACAGTACaA' \
                            'ATACAGaaATCAGaGCAGAAaATACAGATCAaAGCTAGCAaAAtAtA'
            self.temp_dir_path = tempfile.mkdtemp()
            fasta_file = os.path.join(self.temp_dir_path, 'test.fa')
            self.bowtie_build_base = os.path.join(self.temp_dir_path, 'test')
            fasta_stream = open(fasta_file, 'w')
            print >>fasta_stream, '>chr1'
            print >>fasta_stream, \
                '\n'.join([reference_seq[i:i+50] for i in range(0, 251, 50)])
            fasta_stream.close()
            bowtie_build_process = subprocess.call(
                    [args.bowtie_build_exe,
                    fasta_file,
                    self.bowtie_build_base],
                    stdout=open(os.devnull, 'w'),
                    stderr=subprocess.STDOUT
                )
            self.reference_index = bowtie_index.BowtieIndexReference(
                                    self.bowtie_build_base
                                   )

        def test_DP_framing_1(self):
            """ Fails if splice junction is not accurate.

                While this test is passed, it _didn't have to be_; DP framing
                could have yielded more than one possible optimal jump, and the
                jump chosen by the algo may not have been right.

                Takes read to be concatenation of first and third lines of 
                reference_seq from setUp(); the intron should be identified as
                the second line of reference_seq.
            """
            '''Each line of read_seq below and reference_seq above spans
            47 bases.'''
            read_seq = 'ATGGCATACGATACGTCAGACCATGCAggACctTTacCTACATACTG' \
                       'ATCCAGATTACGATACAAaTACGAAcTCccATAGCAaCATaCTAGac'
            # Fake mapping each line of read to respective line of reference
            readlets = [('chr1', False, 1, 48, 0),
                          ('chr1', False, 95, 142, 47)]
            introns = introns_from_read(
                                self.reference_index, read_seq, readlets
                            )
            self.assertEquals({('chr1', False) : [(48, 95)]}, 
                                introns)
            '''Now try truncating readlets so there is an unmapped region of
            the read. This tests the DP framing code in context.'''
            readlets = [('chr1', False, 1, 42, 0),
                          ('chr1', False, 100, 142, 52)]
            introns = introns_from_read(
                                self.reference_index, read_seq, readlets
                            )
            self.assertEquals({('chr1', False) : [(48, 95)]},
                                introns)
            # Truncate readlets again to test DP framing.
            readlets = [('chr1', False, 1, 37, 0),
                          ('chr1', False, 105, 142, 57)]
            introns = introns_from_read(
                                self.reference_index, read_seq, readlets
                            )
            self.assertEquals({('chr1', False) : [(48, 95)]},
                                introns)

        def test_DP_framing_from_overlapping_ECs(self):
            """ Fails if splice junction is not accurate.

                While this test is passed, it _didn't have to be_; DP framing
                could have yielded more than one possible optimal jump, and the
                jump chosen by the algo may not have been right.

                Takes read to be concatenation of first and third lines of 
                reference_seq from setUp(); the intron should be identified as
                the second line of reference_seq. Exonic intervals are taken to
                overlap slightly on the read to see if
                exon_and_introns_from_read()'s remapping for overlapping ECs
                works.
            """
            '''Each line of read_seq below and reference_seq above spans
            47 bases.'''
            read_seq = 'ATGGCATACGATACGTCAGACCATGCAggACctTTacCTACATACTG' \
                       'ATCCAGATTACGATACAAaTACGAAcTCccATAGCAaCATaCTAGac'
            # Fake mapping each line of read to respective line of reference
            readlets = [('chr1', False, 1, 52, 0),
                          ('chr1', False, 91, 142, 43)]
            introns = introns_from_read(
                                self.reference_index, read_seq, readlets
                            )
            self.assertEquals({('chr1', False) : [(48, 95)]},
                                introns)
            '''Now try truncating readlets so there is an unmapped region of
            the read. This tests the DP framing code in context.'''
            readlets = [('chr1', False, 1, 42, 0),
                          ('chr1', False, 100, 142, 52)]
            introns = introns_from_read(
                                self.reference_index, read_seq, readlets
                            )
            self.assertEquals({('chr1', False) : [(48, 95)]},
                                introns)
            # Truncate readlets again to test DP framing.
            readlets = [('chr1', False, 1, 37, 0),
                          ('chr1', False, 105, 142, 57)]
            introns = introns_from_read(
                                self.reference_index, read_seq, readlets
                            )
            self.assertEquals({('chr1', False) : [(48, 95)]},
                                introns)

        def test_reverse_strand_DP_framing_1(self):
            """ Fails if splice junction is not accurate.

                While this test is passed, it _didn't have to be_; DP framing
                could have yielded more than one possible optimal jump, and the
                jump chosen by the algo may not have been right.

                Takes read to be REVERSE COMPLEMENT of concatenation of first
                and third lines of reference_seq from setUp(); the intron
                should be identified as the second line of reference_seq.

                This test is identical to test_DP_framing_1(), only alignments
                are to reverse strand and read_seq is reverse-complemented.
            """
            '''Each line of read_seq below and reference_seq above spans
            47 bases.'''
            read_seq = 'ATGGCATACGATACGTCAGACCATGCAggACctTTacCTACATACTG' \
                       'ATCCAGATTACGATACAAaTACGAAcTCccATAGCAaCATaCTAGac'
            '''Reverse-complement read_seq above and make reverse_strand True
            for all faked alignments below.'''
            read_seq = read_seq[::-1].translate(
                                        _reversed_complement_translation_table
                                      )
            # Fake mapping each line of read to respective line of reference
            readlets = [('chr1', True, 1, 48, 0),
                          ('chr1', True, 95, 142, 47)]
            introns = introns_from_read(
                                self.reference_index, read_seq, readlets
                            )
            self.assertEquals({('chr1', True) : [(48, 95)]},
                                introns)
            '''Now try truncating readlets so there is an unmapped region of
            the read. This tests the DP framing code in context.'''
            readlets = [('chr1', True, 1, 42, 0),
                          ('chr1', True, 100, 142, 52)]
            introns = introns_from_read(
                                self.reference_index, read_seq, readlets
                            )
            self.assertEquals({('chr1', True) : [(48, 95)]},
                                introns)
            # Truncate readlets again to test DP framing.
            readlets = [('chr1', True, 1, 37, 0),
                          ('chr1', True, 105, 142, 57)]
            introns = introns_from_read(
                                self.reference_index, read_seq, readlets
                            )
            self.assertEquals({('chr1', True) : [(48, 95)]},
                                introns)

        def test_DP_framing_2(self):
            """ Fails if splice junction is not accurate.

                While this test is passed, it _didn't have to be_; DP framing
                could have yielded more than one possible optimal jump, and the
                jump chosen by the algo may not have been right.

                Takes read to be concatenation of second and fourth lines of 
                reference_seq from setUp(); the intron should be identified as
                the second line of reference_seq.
            """
            '''Each line of read_seq below and reference_seq above spans
            47 bases.'''
            read_seq = 'GCTACATAGTACATCTAGGCATACTACGTgaCATACGgaCTACGTAA' \
                       'CAttAaaGACTAGACTAACAGACAaAACTAGCATacGATCATGACaA'
            # Fake mapping each line of read to respective line of reference
            readlets = [('chr1', False, 48, 95, 0),
                          ('chr1', False, 142, 189, 47)]
            introns = introns_from_read(
                                self.reference_index, read_seq, readlets
                            )
            self.assertEquals({('chr1', False) : [(95, 142)]},
                                introns)
            '''Now try truncating readlets so there is an unmapped region of
            the read. This tests the DP framing code in context.'''
            readlets = [('chr1', False, 48, 90, 0),
                          ('chr1', False, 147, 189, 52)]
            introns = introns_from_read(
                                self.reference_index, read_seq, readlets
                            )
            self.assertEquals({('chr1', False) : [(95, 142)]},
                                introns)
            # Truncate readlets again to test DP framing.
            readlets = [('chr1', False, 48, 87, 0),
                          ('chr1', False, 150, 189, 55)]
            introns = introns_from_read(
                                self.reference_index, read_seq, readlets
                            )
            self.assertEquals({('chr1', False) : [(95, 142)]},
                                introns)

        def test_DP_filling_between_ECs(self):
            """ Fails if unmapped region between two ECs is not filled.

                Takes read to be fifth line of reference_seq, with a few bases
                in the middle unmapped a priori.
            """
            # Second line of read_seq below was unmapped and should get filled
            read_seq = 'ACGAGATCCATATAtTTAGC' \
                       'AaGACTA' \
                       'aACGATACGATACAGTACaA'
            readlets = [('chr1', False, 189, 209, 0),
                          ('chr1', False, 216, 236, 27)]
            introns = introns_from_read(
                                self.reference_index, read_seq, readlets
                            )
            self.assertEquals({}, introns)

        def test_DP_filling_before_first_EC_and_after_last_EC(self):
            """ Fails if unaligned prefix and suffix don't become ECs.
                
                Reference is searched for unmapped regions before first EC
                and after last EC. Here, the read has exactly one EC, so it's
                the first one and the last one on the read. Read is taken to be
                the sixth line of reference_seq.
            """
            '''Second line of read_seq below is the only EC identified at
            first; first and third lines also appear in reference, and introns
            separate each line.'''
            read_seq = 'ATACAGaaAT' \
                       'CAGaGCAGAAaATACAGATCA' \
                       'aAGCTAGCAaAAtAtA'
            readlets = [('chr1', False, 246, 267, 10)]
            introns = introns_from_read(
                                self.reference_index, read_seq, readlets, 
                                min_seq_similarity=0
                            )
            self.assertEquals({}, introns)

        def test_that_prefix_and_suffix_ECs_are_found_1(self):
            """ Fails if unaligned prefix and suffix don't become ECs.
                
                Reference is searched for unmapped regions before first EC
                and after last EC. Here, the read has exactly one EC, so it's
                the first one and the last one on the read. Read is taken to be
                the first line of reference_seq.
            """
            '''Second line of read_seq below is the only EC identified at
            first; first and third lines also appear in reference, and introns
            separate each line.'''
            # ATGGCATACGATACGTCAGACCATGCAggACctTTacCTACATACTG
            read_seq = 'ATGGCATA' \
                       'GTCAGACCATGCAg' \
                       'CTACATAC'
            readlets = [('chr1', False, 15, 29, 8)]
            introns = introns_from_read(
                                self.reference_index, read_seq, readlets
                             )
            self.assertEquals(introns, 
                                {('chr1', False) : [(9, 15),
                                                    (29, 38)]})

        def test_that_prefix_and_suffix_ECs_are_found_2(self):
            """ Fails if unaligned prefix and suffix don't become ECs.
                
                Here, the read has exactly one EC, so it's the first one and
                the last one on the read. Read is taken to be the sixth line of
                reference_seq.
            """
            # ATACAGaaATCAGaGCAGAAaATACAGATCAaAGCTAGCAaAAtAtA
            # Second line of read_seq below is the only EC identified at first
            read_seq = 'ATACAGaa' \
                       'GCAGAAaATACAGATCA' \
                       'GCAaAAtAtA'
            readlets = [('chr1', False, 250, 267, 8)]
            introns = introns_from_read(
                                self.reference_index, read_seq, readlets
                             )
            self.assertEquals(introns, 
                                {('chr1', False) : [(244, 250),
                                                    (267, 273)]})

        def test_that_prefix_and_suffix_ECs_are_found_3(self):
            """ Fails if unaligned prefix and suffix don't become ECs.
                
                Reference is searched for unmapped regions before first EC
                and after last EC. Here, the read has exactly one EC, so it's
                the first one and the last one on the read. Read is taken to be
                the first line of reference_seq.
            """
            '''Second line of read_seq below is the only EC identified at
            first; first and third lines also appear in reference, and introns
            separate each line.'''
            # ATGGCATACGATACGTCAGACCATGCAggACctTTacCTACATACTG
            read_seq = 'ATGGCATACG' \
                       'GTCAGACCATGCAg' \
                       'CCTACATAC'
            readlets = [('chr1', False, 15, 29, 10)]
            introns = introns_from_read(
                                self.reference_index, read_seq, readlets
                             )
            self.assertEquals(introns, 
                                {('chr1', False) : [(11, 15),
                                                    (29, 37)]})

        def test_that_strange_mapping_is_thrown_out(self):
            """ Fails if any exons or introns are returned.

                Here, the read_seq is taken to be some random sequence, and
                EC #2 occurs BEFORE EC #1 on the read. (That is, displacements
                on read give different ordering of ECs than positions on
                reference.)
            """
            read_seq = random_sequence(60)
            readlets = [('chr1', False, 200, 220, 10),
                           ('chr1', False, 160, 180, 35)]
            introns = introns_from_read(
                                self.reference_index, read_seq, readlets
                            )
            self.assertEquals({}, introns)

        def test_that_small_exon_is_found_1(self):
            reference_seq = 'ATGGCATACGATACGTCAGACCATGCAggACctTTacCTACATACTG' \
                            'GCTACATAGTACATCTAGGCATACTACGTgcCATACGgaCTACGTAA' \
                            'ATCCAGATTACGATACAAaTACGAAcTCccATAGCAaCATaCTAGac' \
                            'CAttAaaGACTAGACTAACAGACAaAACTAGCATacGATCATGACaA' \
                            'ACGAGATCCATATAtTTAGCAaGACTAaACGATACGATACAGTACaA' \
                            'ATACAGaaATCAGaGCAGAAaATACAGATCAaAGCTAGCAaAAtAtA'
            """ Fails if two introns aren't returned rather than one.

                Here, the read_seq is taken from the second and third lines
                of reference_seq.
            """
            read_seq = 'GCTACATAGTACA' \
                       'GTGCCATACG' \
                       'ATCCAGATTACGATACAAATACGA'
            readlets = [('chr1', False, 48, 61, 0),
                        ('chr1', False, 95, 119, 23)]
            introns = introns_from_read(
                            self.reference_index, read_seq, readlets,
                            min_cap_size=6, readlet_interval=1
                        )
            self.assertEquals({('chr1', False) : [(61, 75),
                                                  (85, 95)]}, introns)

        def test_that_small_exon_is_found_2(self):
            reference_seq = 'ATGGCATACGATACGTCAGACCATGCAggACctTTacCTACATACTG' \
                            'GCTACATAGTACATCTAGGCATACTACGTgcCATACGgaCTACGTAA' \
                            'ATCCAGATTACGATACAAaTACGAAcTCccATAGCAaCATaCTAGac' \
                            'CAttAaaGACTAGACTAACAGACAaAACTAGCATacGATCATGACaA' \
                            'ACGAGATCCATATAtTTAGCAaGACTAaACGATACGATACAGTACaA' \
                            'ATACAGaaATCAGaGCAGAAaATACAGATCAaAGCTAGCAaAAtAtA'
            """ Fails if two introns aren't returned rather than one.

                Here, the read_seq is taken from the first and second lines
                of reference_seq.
            """
            read_seq = 'ATGGCATACGATACGTCAG' \
                       'CTACATACTG' \
                       'GTGCCATACGGACTACGTAA'
            readlets = [('chr1', False, 1, 20, 0),
                        ('chr1', False, 75, 95, 29)]
            introns = introns_from_read(
                            self.reference_index, read_seq, readlets,
                            min_cap_size=6, readlet_interval=1
                        )
            self.assertEquals({('chr1', False) : [(20, 38),
                                                  (48, 75)]}, introns)

        def tearDown(self):
            # Kill temporary directory
            shutil.rmtree(self.temp_dir_path)
   
    unittest.main()