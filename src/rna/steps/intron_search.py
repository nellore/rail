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

base_path = os.path.abspath(
                    os.path.dirname(os.path.dirname(os.path.dirname(
                        os.path.realpath(__file__)))
                    )
                )
utils_path = os.path.join(base_path, 'rna', 'utils')
site.addsitedir(utils_path)
site.addsitedir(base_path)

import bowtie
import bowtie_index
import partition
from dooplicity.tools import xstream

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
_left_reverse_elements = set(['CT', 'GT'])
_left_forward_elements = set(['GT', 'GC', 'AT'])
_right_reverse_elements = set(['AC', 'GC', 'AT'])
_right_forward_elements = set(['AG', 'AC'])
_left_elements = _left_reverse_elements | _left_forward_elements
_right_elements = _right_reverse_elements | _right_forward_elements

if 'pypy' not in sys.version.lower():
    # For fast global alignment without PyPy
    from scipy import weave
    class GlobalAlignment(object):
        """ Invokes Weave to obtain alignment score matrix with C. """

        def __init__(self, substitution_matrix=[[ 0,-1,-1,-1,-1,-1],
                                                [-1, 0,-1,-1,-1,-1],
                                                [-1,-1, 0,-1,-1,-1],
                                                [-1,-1,-1, 0,-1,-1],
                                                [-1,-1,-1,-1,-1,-1],
                                                [-1,-1,-1,-1,-1,-1]]):
            """ Constructor for GlobalAlignment.

                Places substitution_matrix directly in string containing C code
                that obtains score matrix for global alignment. Executes C with
                Weave on trivial case (two sequences, each one character) to 
                ensure code is precompiled.

                substitution_matrix: 6 x 6 substitution matrix (list of
                    lists); rows and columns correspond to ACGTN-, where N is
                    aNy and - is a gap. Default: +1 for match, -5 for gap,
                    -1 for everything else.
            """
            '''Set supporting code including libraries and the function 
            integer_sequence(), which converts a sequence of nucleotides to
            an array of indices that correspond to indices of the substitution
            matrix---also set in the constructor.'''
            self.support_code = """ 
                        int* integer_sequence(const char* seq, int length) {
                            int i;
                            int* int_seq = (int *)malloc(sizeof(int)*length);
                            for (i = 0; seq[i]; i++) {
                                switch (seq[i]) {
                                    case 'A':
                                    case 'a':
                                        int_seq[i] = 0;
                                        break;
                                    case 'C':
                                    case 'c':
                                        int_seq[i] = 1;
                                        break;
                                    case 'G':
                                    case 'g':
                                        int_seq[i] = 2;
                                        break;
                                    case 'T':
                                    case 't':
                                        int_seq[i] = 3;
                                        break;
                                    case 'N':
                                        int_seq[i] = 4;
                                        break;
                                    case '-':
                                        int_seq[i] = 5;
                                        break;
                                    default:
                                        int_seq[i] = 4;
                                }
                            }
                            return int_seq;
                        }
                        """
            '''Score matrix code is written in C/C++; see self.support_code
            above for dependencies. substitution_matrix is embedded in code
            below.'''
            self.score_matrix_code = """
                const char* c_first_seq = first_seq.c_str();
                const char* c_second_seq = second_seq.c_str();
                int row_count = first_seq.length() + 1;
                int column_count = second_seq.length() + 1;
                npy_intp dims[2] = {row_count, column_count};
                PyObject* to_return = PyArray_SimpleNew(2, dims, NPY_INT);
                int* score_matrix = (int *)((PyArrayObject*) to_return)->data;
                int* substitution_matrix = (int *)malloc(sizeof(int)*6*6);
                int* int_first_seq = integer_sequence(c_first_seq,
                                                        row_count - 1);
                int* int_second_seq = integer_sequence(c_second_seq,
                                                        column_count - 1);
            """ + """""".join(
                            [("""substitution_matrix[%d] = %d;""" % ((i*6 + j), 
                                substitution_matrix[i][j])) 
                                for i in xrange(6) for j in xrange(6)]
                        ) \
            + """
                int i, j;
                score_matrix[0] = 0;
                for (i = 1; i < row_count; i++) {
                    score_matrix[i*column_count] = i 
                        * substitution_matrix[int_first_seq[i-1]*6+5];
                }
                for (j = 1; j < column_count; j++) {
                    score_matrix[j] = j 
                        * substitution_matrix[5*6+int_second_seq[j-1]];
                }
                int diagonal, vertical, horizontal;
                for (i = 1; i < row_count; i++) {
                    for (j = 1; j < column_count; j++) {
                        diagonal = score_matrix[(i-1)*column_count + j-1]
                            + substitution_matrix[int_first_seq[i-1]*6
                                                    +int_second_seq[j-1]];
                        vertical = score_matrix[(i-1)*column_count + j]
                            + substitution_matrix[int_first_seq[i-1]*6+5];
                        horizontal = score_matrix[i*column_count + j-1]
                            + substitution_matrix[5*6+int_second_seq[j-1]];
                        score_matrix[i*column_count + j] = std::max(
                                horizontal, std::max(diagonal, vertical)
                            );
                    }
                }
                free(int_first_seq);
                free(int_second_seq);
                free(substitution_matrix);
                return_val = to_return;
            """
            stdout_holder = os.dup(sys.stdout.fileno())
            try:
                # Suppress compiler output by redirecting stdout temporarily
                devnull = open(os.devnull, 'w')
                os.dup2(devnull.fileno(), sys.stdout.fileno())
                self.score_matrix('N', 'N')
            except Exception as e:
                print >>sys.stderr, 'C code for computing score matrix ' \
                                    'failed to compile. Get Python\'s ' \
                                    'distutils to use g++, and try again.'
                raise
            finally:
                # Turn the tap back on
                os.dup2(stdout_holder, sys.stdout.fileno())

        def score_matrix(self, first_seq, second_seq):
            """ Computes score matrix for global alignment of two sequences.

                The substitution matrix is specified when the GlobalAlignment
                class is instantiated.

                first_seq: first sequence (string).
                second_seq: second sequence (string).

                Return value: score_matrix, a numpy array whose dimensions are
                    (len(first_seq) + 1) x (len(second_seq) + 1). It can be
                    used to trace back the best global alignment.
            """
            return weave.inline(
                                    self.score_matrix_code,
                                    ['first_seq', 'second_seq'], 
                                    support_code=self.support_code,
                                    verbose=0
                                )
else:
    # Use Python version of class
    class GlobalAlignment(object):
        """ Uses Python to obtain alignment score matrix. """

        def __init__(self, substitution_matrix=[[ 0,-1,-1,-1,-1,-1],
                                                [-1, 0,-1,-1,-1,-1],
                                                [-1,-1, 0,-1,-1,-1],
                                                [-1,-1,-1, 0,-1,-1],
                                                [-1,-1,-1,-1,-1,-1],
                                                [-1,-1,-1,-1,-1,-1]]):
            """ Constructor for GlobalAlignment.

                substitution_matrix: 6 x 6 substitution matrix (list of
                    lists); rows and columns correspond to ACGTN-, where N is
                    aNy and - is a gap. Default: +1 for match, -5 for gap,
                    -1 for everything else.
            """
            self.substitution_matrix = substitution_matrix

        def score_matrix(self, first_seq, second_seq):
            """ Computes score matrix for global alignment of two sequences.

                The substitution matrix is specified when the GlobalAlignment
                class is instantiated.

                first_seq: first sequence (string).
                second_seq: second sequence (string).

                Return value: score_matrix, a numpy array whose dimensions are
                    (len(first_seq) + 1) x (len(second_seq) + 1). It can be
                    used to trace back the best global alignment.
            """
            first_seq = ['ACGTN-'.index(char) for char in first_seq]
            second_seq = ['ACGTN-'.index(char) for char in second_seq]
            row_count = len(first_seq) + 1
            column_count = len(second_seq) + 1
            score_matrix = [[0 for i in xrange(column_count)]
                                for j in xrange(row_count)]
            for j in xrange(1, column_count):
                score_matrix[0][j] = j \
                    * self.substitution_matrix[5][second_seq[j-1]]
            for i in xrange(1, row_count):
                score_matrix[i][0] = i \
                    * self.substitution_matrix[first_seq[i-1]][5]
            for i in xrange(1, row_count):
                for j in xrange(1, column_count):
                    score_matrix[i][j] = max(score_matrix[i-1][j-1] 
                                                + self.substitution_matrix[
                                                        first_seq[i-1]]
                                                        [second_seq[j-1]
                                                    ], # diagonal
                                             score_matrix[i-1][j]
                                                + self.substitution_matrix[
                                                        first_seq[i-1]][5
                                                    ], # vertical
                                             score_matrix[i][j-1]
                                                + self.substitution_matrix[
                                                        5][second_seq[j-1]
                                                    ] # horizontal
                                            )
            return score_matrix

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

def maximum_clique(cluster):
    """ Finds maximum clique of graph of multireadlet alignment cluster.

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
        
        Return value: maximum clique -- a list of alignments.
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
        alignment from the same multireadlet. The cluster with the largest
        maximum clique (as described is returned. If there is a tie, no clique
        is returned.

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
    maximum_cliques = []
    for cluster in clustered_alignments:
        maximum_cliques.append(maximum_clique(cluster))
    largest_maximum_clique_size = max(map(len, maximum_cliques))
    largest_maximum_cliques = [clique for clique in maximum_cliques
                               if len(clique) == largest_maximum_clique_size]
    if len(largest_maximum_cliques) == 1:
        return largest_maximum_cliques[0]
    else:
        '''Alignment is in general too repetitive; postpone treatment until
        a nice systematic way to handle this case is found.'''
        pass

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
        reverse_reverse_strand=False, global_alignment=GlobalAlignment(),
        max_gaps_mismatches=5):
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
        global_alignment: object of class GlobalAlignment used for fast
            realignment to reference
        max_gaps_mismatches: maximum number of (gaps + mismatches) to permit
            in realignments to reference minus intron per 100 bp or None if
            unlimited
    """
    if not clique:
        return
    read_seq = read_seq.upper()
    read_seq_size = len(read_seq)
    if stranded:
        if reverse_strand:
            search_motifs = _reverse_strand_motifs
            left_search_motifs = _left_reverse_elements
            right_search_motifs = _right_reverse_elements
        else:
            search_motifs = _forward_strand_motifs
            left_search_motifs = _left_forward_elements
            right_search_motifs = _right_forward_elements
    else:
        search_motifs = _all_motifs
        left_search_motifs = _left_elements
        right_search_motifs = _right_elements
    rname, reverse_strand = clique[0][:2]
    if reverse_strand:
        read_seq = read_seq[::-1].translate(
            _reversed_complement_translation_table
        )
    if reverse_reverse_strand:
        reverse_strand = not reverse_strand
    # Arrange sort so that longest aligning capping readlet is used as flank
    clique.sort(key=lambda alignment: (alignment[-1], alignment[-2] if
                                                        alignment[-1] == 0
                                                        else alignment[-3]))
    _, _, prefix_pos, prefix_end_pos, prefix_displacement = clique[0]
    _, _, suffix_pos, suffix_end_pos, suffix_displacement = clique[-1]
    new_prefix, new_suffix = [], []
    if clique[0][-1] >= min_exon_size and search_window_size:
        # Find possible maximal matching prefixes
        search_pos = max(0, prefix_pos - 1 - search_window_size)
        search_window = reference_index.get_stretch(
                rname,
                search_pos,
                prefix_pos - 1 - search_pos
            )
        try:
            new_prefix_offset, new_prefix_size \
                = maximal_suffix_match(
                        read_seq[:prefix_displacement][::-1],
                        search_window[::-1],
                        min_cap_size=min_exon_size
                    )
            new_prefix = [(rname, reverse_strand, 
                            prefix_pos - new_prefix_offset
                            - new_prefix_size,
                            prefix_pos - new_prefix_offset,
                            0)]
        except TypeError:
            # maximal_suffix_match returned None
            pass
    unmapped_displacement = suffix_end_pos - suffix_pos + suffix_displacement
    unmapped_base_count = read_seq_size - unmapped_displacement
    if suffix_end_pos + unmapped_base_count - 1 \
        > reference_index.length[rname]:
        # Skip this part if too close to edge of reference
        unmapped_base_count = 0
    if unmapped_base_count >= min_exon_size and search_window_size:
        search_pos = suffix_end_pos - 1
        search_window = reference_index.get_stretch(
                rname,
                search_pos,
                min(search_window_size, 
                     reference_index.rname_lengths[rname] - search_pos)
            )
        try:
            new_suffix_offset, new_suffix_size \
                = maximal_suffix_match(
                        read_seq[-unmapped_base_count:],
                        search_window,
                        min_cap_size=min_exon_size
                    )
            new_suffix = [(rname, reverse_strand,
                            suffix_end_pos + new_suffix_offset,
                            suffix_end_pos + new_suffix_offset
                            + new_suffix_size,
                            read_seq_size - new_suffix_size)]
        except TypeError:
            # maximal_suffix_match returned None
            pass
    for ((_, _, left_pos, left_end_pos, left_displacement),
            (_, _, right_pos, right_end_pos, right_displacement)) \
        in pairwise(new_prefix + clique + new_suffix):
        candidate_introns = []
        unmapped_start = left_displacement + left_end_pos - left_pos
        unmapped_end = right_displacement
        read_span = right_displacement + right_end_pos - right_pos \
                        - left_displacement
        # Store number of intronic bases in intevening region
        intronic_base_count = right_end_pos - left_pos - read_span
        if intronic_base_count <= min_intron_size:
            # Filter out introns that are too small
            continue
        left_motif_search_pos_1 = left_end_pos - motif_radius
        right_motif_search_pos_1 = right_pos - motif_radius
        left_motif_search_pos_2 = right_motif_search_pos_1 \
                                   - intronic_base_count
        right_motif_search_pos_2 = left_motif_search_pos_1 \
                                   + intronic_base_count
        left_motif_search_bounds = (
                min(left_motif_search_pos_1, left_motif_search_pos_2),
                max(left_motif_search_pos_1 + 2 * motif_radius,
                    left_motif_search_pos_2 + 2 * motif_radius)
            )
        right_motif_search_bounds = (
                min(right_motif_search_pos_1, right_motif_search_pos_2) - 2,
                max(right_motif_search_pos_1 + 2 * motif_radius,
                    right_motif_search_pos_2 + 2 * motif_radius)
            )
        left_motif_search_size = left_motif_search_bounds[1] \
                                    - left_motif_search_bounds[0]
        right_motif_search_size = right_motif_search_bounds[1] \
                                    - right_motif_search_bounds[0]
        left_motif_window = reference_index.get_stretch(
                                        rname,
                                        left_motif_search_bounds[0] - 1,
                                        left_motif_search_size
                                    )
        right_motif_window = reference_index.get_stretch(
                                        rname,
                                        right_motif_search_bounds[0] - 1,
                                        right_motif_search_size
                                    )
        left_offsets = []
        left_motifs = []
        right_offsets = []
        right_motifs = []
        for i in xrange(left_motif_search_size):
            if left_motif_window[i:i+2] in left_search_motifs:
                left_motifs.append(left_motif_window[i:i+2])
                left_offsets.append(i)
        for i in xrange(right_motif_search_size):
            if right_motif_window[i:i+2] in right_search_motifs:
                right_motifs.append(right_motif_window[i:i+2])
                right_offsets.append(i)
        product_offsets = list(itertools.product(left_offsets, right_offsets))
        product_motifs = list(itertools.product(left_motifs, right_motifs))
        candidate_introns = []
        for i in xrange(len(product_motifs)):
            if product_motifs[i] in search_motifs:
                '''Appropriate motif combo found! Compute intron size
                ASSUMING ONE INTRON and number of intervening exonic bases.'''
                intron_pos \
                    = left_motif_search_bounds[0] + product_offsets[i][0]
                intron_end_pos \
                    = right_motif_search_bounds[0] + product_offsets[i][1] + 2
                intron_size = intron_end_pos - intron_pos
                if intron_size < intronic_base_count:
                    '''Intervening exonic bases should only split up the single
                    intron and be mistakenly counted as intronic bases,
                    increasing intron_size. If intron_size
                    < intronic_base_count, a deletion or a mismap may be
                    responsible, but admitting the option increases the
                    false positive rate. Ignore.'''
                    continue
                # Assume at most ONE intermediate small exon
                small_exon_size = intron_size - intronic_base_count
                if min_intron_size <= intron_size <= max_intron_size:
                    if small_exon_size == 0:
                        '''No need to search for intermediate exon. Realign
                        and add to list of candidate exons.'''
                        left_stretch = reference_index.get_stretch(
                                                        rname,
                                                        left_pos - 1,
                                                        intron_pos
                                                        - left_pos
                                                    )
                        right_stretch = reference_index.get_stretch(
                                                        rname,
                                                        intron_end_pos - 1,
                                                        right_end_pos
                                                        - intron_end_pos
                                                    )
                        reference_minus_intron = left_stretch + right_stretch
                        alignment_score = global_alignment.score_matrix(
                                                read_seq[
                                                    left_displacement:
                                                    left_displacement+read_span
                                                ],
                                                reference_minus_intron)[-1][-1]
                        candidate_introns.append(
                                (
                                    rname,
                                    True if product_motifs[i] in
                                    _reverse_strand_motifs else False,
                                    intron_pos,
                                    intron_end_pos,
                                    alignment_score
                                )
                            )
                    elif small_exon_size >= min_exon_size \
                        and search_window_size:
                        '''Search only within 1000 bases of either end of the
                        intron to find where to split it. Enumerate possible
                        motif ends and starts, and match small exon. Its bases
                        from the read are well-defined assuming no complicating
                        indels.'''
                        small_exon_start = intron_pos - left_end_pos \
                                             + unmapped_start
                        small_exon_end = unmapped_end - right_pos \
                                          + intron_end_pos
                        small_exon = read_seq[small_exon_start:small_exon_end]
                        if len(small_exon) != small_exon_size:
                            # Motif search provided combo that's out of bounds
                            continue
                        left_caps = [motif[1] for motif in search_motifs
                                        if motif[0] == product_motifs[i][0]]
                        right_caps = [motif[0] for motif in search_motifs
                                        if motif[1] == product_motifs[i][1]]
                        cap_combos = list(
                                itertools.product(left_caps, right_caps)
                            )
                        small_exon = read_seq[small_exon_start:
                                                small_exon_end]
                        '''Only need to compute alignment score once; it's the
                        same for all small exons found since only exact
                        matches are admitted.'''
                        left_stretch = reference_index.get_stretch(
                                                        rname,
                                                        left_pos - 1,
                                                        intron_pos
                                                        - left_pos
                                                    )
                        right_stretch = reference_index.get_stretch(
                                                        rname,
                                                        intron_end_pos - 1,
                                                        right_end_pos
                                                        - intron_end_pos
                                                    )
                        reference_minus_introns = left_stretch + small_exon \
                                                    + right_stretch
                        if intron_size <= 2 * search_window_size:
                            search_windows = [reference_index.get_stretch(
                                                        rname,
                                                        intron_pos - 1,
                                                        intron_size
                                                    )]
                        else:
                            search_windows = [reference_index.get_stretch(
                                                        rname,
                                                        intron_pos - 1,
                                                        search_window_size
                                                ),
                                              reference_index.get_stretch(
                                                        rname,
                                                        intron_end_pos
                                                        - search_window_size
                                                        - 1,
                                                        search_window_size
                                                    )]
                        found_alignment_score = False
                        for cap_combo in cap_combos:
                            capped_exon \
                                = cap_combo[0] + small_exon + cap_combo[1]
                            for j, search_window in enumerate(search_windows):
                                for small_exon_match \
                                    in re.finditer(capped_exon, search_window):
                                    if not found_alignment_score:
                                        alignment_score \
                                            = global_alignment.score_matrix(
                                                    read_seq[
                                                        left_displacement:
                                                        left_displacement
                                                        + read_span
                                                    ],
                                                    reference_minus_introns
                                                )[-1][-1]
                                        found_alignment_score = True
                                    # Split original intron
                                    if j == 0:
                                        # First search-window type
                                        first_intron_end_pos \
                                            = intron_pos \
                                                + small_exon_match.start() + 2
                                        second_intron_pos \
                                            = first_intron_end_pos \
                                                + small_exon_size
                                    else:
                                        # j == 1: second search-window type
                                        second_intron_pos \
                                            = intron_end_pos \
                                                - search_window_size \
                                                + small_exon_match.end() - 2
                                        first_intron_end_pos \
                                            = second_intron_pos \
                                                - small_exon_size
                                    if min_intron_size \
                                        <= first_intron_end_pos - intron_pos \
                                        <= max_intron_size:
                                        candidate_introns.append(
                                                (rname,
                                                  True if
                                                  (product_motifs[i][0],
                                                   cap_combo[0])
                                                  in _reverse_strand_motifs
                                                  else False,
                                                  intron_pos,
                                                  first_intron_end_pos,
                                                  alignment_score)
                                            )
                                    if min_intron_size \
                                        <= intron_end_pos - second_intron_pos \
                                        <= max_intron_size:
                                        candidate_introns.append(
                                                (rname,
                                                  True if
                                                  (cap_combo[1],
                                                   product_motifs[i][1])
                                                  in _reverse_strand_motifs
                                                  else False,
                                                  second_intron_pos,
                                                  intron_end_pos,
                                                  alignment_score)
                                            )
        try:
            max_score = max([intron[-1] for intron in candidate_introns])
            if max_gaps_mismatches is None \
                or max_score >= -max_gaps_mismatches * read_seq_size / 100.:
                '''Filter out alignments with more than max_gap_mismatches
                gaps or mismatches per 100 bp.'''
                for (rname, intron_reverse_strand,
                        pos, end_pos, alignment_score) in candidate_introns:
                    if alignment_score == max_score:
                        yield (rname, intron_reverse_strand, pos, end_pos)
        except ValueError:
            # No introns found
            pass

def go(input_stream=sys.stdin, output_stream=sys.stdout,
    bowtie_index_base='genome', verbose=False, stranded=False, min_exon_size=8,
    min_intron_size=15, max_intron_size=500000, motif_radius=1,
    search_window_size=1000, global_alignment=GlobalAlignment(),
    max_gaps_mismatches=5):
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
        min_exon_size: the reference is not searched for an exon smaller
            than this size.
        search_window_size: the size (in bp) of the reference subsequence
            in which to search for an exon.
        max_cap_count: maximum number of possible caps of size
            min_cap_size to consider when searching for caps.
        global_alignment: instance of GlobalAlignment class used for fast
                realignment via Weave or Pypy.
        max_gaps_mismatches: maximum number of gaps/mismatches to permit in
            a realignment to reference without intron per 100 bp
            or None if unlimited
        report_multiplier: if verbose is True, the line number of an alignment,
            read, or first readlet of a read written to stderr increases
            exponentially with base report_multiplier.

        No return value.
    """
    global _input_line_count, _output_line_count
    reference_index = bowtie_index.BowtieIndexReference(bowtie_index_base)
    '''Input is readletized, and readlets must be composed.'''
    for (seq_id,), xpartition in xstream(input_stream, 1):
        readlet_displacements, collected_readlets = [], []
        seq_info_captured = False
        for seq_info, rnames, flags, poses in xpartition:
            _input_line_count += 1
            seq_info = seq_info.split('\x1e')
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
            if poses != '\x1c':
                # If there are alignments, add them to collected_readlets
                collected_readlets.append(
                        zip(
                                rnames.split('\x1f'),
                                [(int(flag) & 16) != 0
                                    for flag in flags.split('\x1f')],
                                [int(pos) for pos in poses.split('\x1f')]
                            )
                    )
                readlet_displacements.append(
                    (int(seq_info[0]), int(seq_info[1]))
                )
        assert seq_info_captured, \
            'Sequence info was not in a collected readlet'
        if collected_readlets:
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
                    introns = list(introns_from_clique(
                                    selected_readlet_alignments_by_clustering(
                                            multireadlets
                                        ),
                            seq,
                            reference_index,
                            min_exon_size=min_exon_size,
                            min_intron_size=min_intron_size,
                            max_intron_size=max_intron_size,
                            search_window_size=search_window_size,
                            stranded=stranded,
                            motif_radius=motif_radius,
                            reverse_reverse_strand=False,
                            global_alignment=global_alignment,
                            max_gaps_mismatches=max_gaps_mismatches
                        ))
                    for (intron_rname, intron_reverse_strand,
                            intron_pos, intron_end_pos) in introns:
                        for sample_label in sample_labels:
                            print '%s%s\t%s\t%012d\t%012d' % (
                                    intron_rname,
                                    '-' if intron_reverse_strand else '+',
                                    sample_label,
                                    intron_pos,
                                    intron_end_pos
                                )
                            _output_line_count += 1
                if reversed_complement_sample_labels:
                    introns = introns_from_clique(
                                    selected_readlet_alignments_by_clustering(
                                            multireadlets
                                        ),
                            seq,
                            reference_index,
                            min_exon_size=min_exon_size,
                            min_intron_size=min_intron_size,
                            max_intron_size=max_intron_size,
                            search_window_size=search_window_size,
                            stranded=stranded,
                            motif_radius=motif_radius,
                            reverse_reverse_strand=True,
                            global_alignment=global_alignment,
                            max_gaps_mismatches=max_gaps_mismatches
                        )
                    for (intron_rname, intron_reverse_strand,
                             intron_pos, intron_end_pos) in introns:
                        for sample_label in reversed_complement_sample_labels:
                            print '%s%s\t%s\t%012d\t%012d' % (
                                    intron_rname,
                                    '-' if intron_reverse_strand else '+',
                                    sample_label,
                                    intron_pos,
                                    intron_end_pos
                                )
                            _output_line_count += 1
            else:
                introns = introns_from_clique(
                                selected_readlet_alignments_by_clustering(
                                        multireadlets
                                    ),
                        seq,
                        reference_index,
                        min_exon_size=min_exon_size,
                        min_intron_size=min_intron_size,
                        max_intron_size=max_intron_size,
                        search_window_size=search_window_size,
                        motif_radius=motif_radius,
                        stranded=stranded,
                        global_alignment=global_alignment,
                        max_gaps_mismatches=max_gaps_mismatches
                    )
                for (intron_rname, intron_reverse_strand,
                        intron_pos, intron_end_pos) in introns:
                    for sample_label in (sample_labels
                                          | reversed_complement_sample_labels):
                        print '%s%s\t%s\t%012d\t%012d' % (
                                intron_rname,
                                '-' if intron_reverse_strand else '+',
                                sample_label,
                                intron_pos,
                                intron_end_pos
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
    parser.add_argument('--max-gaps-mismatches', type=int, required=False,
        default=3,
        help='Maximum number of (gaps + mismatches) to permit in a '
             'realignment to reference without intron per 100 bp')
    parser.add_argument('--motif-radius', type=int, required=False,
        default=1,
        help='Number of bases to tack on to either end of an unmapped region '
             'when searching for motifs')

    # Add command-line arguments for dependencies
    partition.add_args(parser)
    bowtie.add_args(parser)

    '''Now collect arguments. While the variable args declared below is
    global, properties of args are also arguments of the go() function so
    different command-line arguments can be passed to it for unit tests.'''
    args = parser.parse_args(sys.argv[1:])

if __name__ == '__main__' and not args.test:
    import time
    start_time = time.time()
    global_alignment = GlobalAlignment()
    go(bowtie_index_base=args.bowtie_idx,
        verbose=args.verbose, 
        stranded=args.stranded,
        min_intron_size=args.min_intron_size,
        max_intron_size=args.max_intron_size,
        min_exon_size=args.min_exon_size,
        motif_radius=args.motif_radius,
        search_window_size=args.search_window_size,
        max_gaps_mismatches=args.max_gaps_mismatches,
        global_alignment=global_alignment)
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

    # Precomile global_alignment
    if 'pypy' not in sys.version.lower():
        global_alignment = GlobalAlignment()

    def random_sequence(seq_size):
        """ Gets random sequence of nucleotides.

            seq_size: number of bases to return.

            Return value: string of random nucleotides.
        """
        return ''.join([random.choice('ATCG') for _ in xrange(seq_size)])
    
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

    class TestIntronsFromClique(unittest.TestCase):
        """ Tests introns_from_clique(). """
        def setUp(self):
            """ Creates temporary directory and Bowtie index. """
            reference_seq = 'ATGGCATACGATACGTCAGACCATGCAggACctTTacCTACATACTG' \
                            'GTTACATAGTACATATAGGCATACTAGGTgcCATACGgaCTACGTAG' \
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

        def test_split_intron_1(self):
            """ Fails if single intron is not split in two.

                Single intron that's too large spans second line of
                reference_seq above. There is a small exon in the middle
                whose bases are initially mistaken for intronic bases.
                read_seq below contains the small exon and a few bases outside
                it.
            """
            '''Each line of read_seq below and reference_seq above spans
            47 bases.'''
            read_seq = 'AGGACCTTTACCTACATACTGGCATACTAGATCCAGATTACGATAC'
            clique = [('chr1', False, 27, 48, 0), ('chr1', False, 95, 111, 30)]
            introns = introns_from_clique(clique, read_seq,
                                            self.reference_index,
                                            min_exon_size=8,
                                            min_intron_size=5,
                                            search_window_size=1000,
                                            stranded=False,
                                            motif_radius=0,
                                            global_alignment=global_alignment,
                                            max_gaps_mismatches=5)
            self.assertEquals(
                    list(introns),
                    [('chr1', False, 48, 66), ('chr1', False, 75, 95)]
                )
            

        def tearDown(self):
            # Kill temporary directory
            shutil.rmtree(self.temp_dir_path)
    unittest.main()
