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
import numpy as np
import random
import itertools
# For fast global alignment
from scipy import weave
from collections import defaultdict

base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for directory_name in ['bowtie', 'sample', 'interval']:
    site.addsitedir(os.path.join(base_path, directory_name))

import bowtie
import bowtie_index
import sample
import partition

_reversed_complement_translation_table = string.maketrans('ATCG', 'TAGC')

class GlobalAlignment:
    """ Invokes Weave to obtain alignment score matrix with C. """

    def __init__(self, substitution_matrix=[[ 1,-1,-1,-1,-1,-5],
                                            [-1, 1,-1,-1,-1,-5],
                                            [-1,-1, 1,-1,-1,-5],
                                            [-1,-1,-1, 1,-1,-5],
                                            [-1,-1,-1,-1,-1,-5],
                                            [-5,-5,-5,-5,-5,-5]]):
        """ Constructor for GlobalAlignment.

            Places substitution_matrix directly in string containing C code
            that obtains score matrix for global alignment. Executes C with
            Weave on trivial case (two sequences, each one character) to 
            ensure code is precompiled.

            substitution_matrix: 6 x 6 substitution matrix (list of
                lists); rows and columns correspond to ACGTN-, where N is aNy
                and - is a gap. Default: +1 for match, -2 for gap, -1 for
                everything else.
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
        '''Score matrix code is written in C/C++; see self.support_code above
        for dependencies. substitution_matrix is embedded in code below.'''
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
            print >>sys.stderr, 'C code for computing score matrix failed ' \
                                'to compile. Get Python\'s distutils to use ' \
                                'g++, and try again.'
            raise
        finally:
            # Turn the tap back on
            os.dup2(stdout_holder, sys.stdout.fileno())

    def score_matrix(self, first_seq, second_seq):
        """ Computes score matrix for global alignment of two DNA sequences.

            The substitution matrix is specified when the GlobalAlignment
            class is instantiated.

            first_seq: first sequence (string).
            second_seq: second sequence (string).

            Return value: score_matrix, a numpy array whose dimensions are
                (len(first_seq) + 1) x (len(second_seq) + 1). It can be used to
                trace back the best global alignment.
        """
        return weave.inline(
                                self.score_matrix_code,
                                ['first_seq', 'second_seq'], 
                                support_code=self.support_code,
                                verbose=0
                            )

# Initialize global variables for tracking number of input/output lines
_input_line_count = 0
_output_line_count = 0

def composed_and_sorted_readlets(readlets):
    """ Composes overlapping readlets and sorts them by reference position.

        This function assumes that for overlapping readlets mapping to, for
        example, the forward strand and covering the reference like so:

            Reference  5'---------------------3'
            Readlet 1       =========
            Readlet 2           ============
            Readlet 3      ===============    ,

        the displacement of the readlet whose pos is smallest
        (here, Readlet 3), is the displacement of the interval overlapped ---
        EVEN IF THE OTHER READLETS' DISPLACEMENTS SEEM INCONSISTENT WITH THAT
        DISPLACEMENT given readlet pos's AND end_pos's. See below for
        definitions of displacement, pos, and end_pos.

        readlets: list of tuples (rname, reverse_strand, pos, end_pos,
            displacement), each of which corresponds to a readlet. rname is the
            SAM-format RNAME, typically the chromosome to which the readlet
            maps. reverse_strand is True iff the readlet's reversed complement
            aligns to the reference. The readlet should span the reference
            interval [pos, end_pos). displacement is the number of bases
            between the 5' (3') end of the readlet, which aligns to the forward
            (reverse) strand, and the 5' (3') end of the read.

        Return value: dictionary with each key a tuple (rname, reverse_strand)
            uniquely identifying a strand and each value a list of tuples
            (pos, end_pos, displacement) sorted by pos in ascending order;
            here, each tuple from the list denotes an exonic interval spanning
            [pos, end_pos) and displacement is the number of bases between the
            5' (3') end of the interval, composed of readlets aligning to the
            forward (reverse) strand, and the 5' (3') end of the read.
    """
    # Create dictionary separating readlets by strands to which they align
    uncomposed = {}
    for rname, reverse_strand, pos, end_pos, displacement in readlets:
        assert end_pos > pos
        if not uncomposed.has_key((rname, reverse_strand)):
            uncomposed[(rname, reverse_strand)] = []
        uncomposed[(rname, reverse_strand)].append((pos, end_pos, 
                                                    displacement))
    to_delete = []
    for strand in uncomposed:
        # Sort by start positions on reference
        uncomposed[strand].sort()
    for strand in to_delete: del uncomposed[strand]
    
    # Create dictionary for merging readlets whose alignments overlap
    composed = {}
    for strand in uncomposed:
        composed[strand] = [uncomposed[strand][0]]
        for pos, end_pos, displacement in uncomposed[strand][1:]:
            if pos <= composed[strand][-1][1]:
                '''If start position of readlet is <= end position of last
                interval, merge into one interval by changing only the
                interval's end position.'''
                composed[strand][-1] = (composed[strand][-1][0],
                                        max(end_pos, composed[strand][-1][1]),
                                        composed[strand][-1][2])
            else:
                # Start a new interval
                composed[strand].append((pos, end_pos, displacement))
    '''Now composed[strand] is list of candidate exonic chunks
    (pos, end_pos, displacement) ordered by pos, the reference position.'''
    return composed

def unmapped_region_splits(unmapped_seq, left_reference_seq,
        right_reference_seq, global_alignment=GlobalAlignment()):
    """ Distributes a read's unmapped region between two framing mapped
        regions.

        The algorithm used is best illustrated with an example:
        Read         5' mmmmmmmmmmmmmmmmmmmm??????|????MMMMMMMMMMMMMMMMMMMMM 3'

        Reference    5' ...mmmmmmmmmLLLLLL|LLLLBB...BBRRRRRR|RRRRMMMMMMMM... 3'

        Besides the pipe bars, each character above represents a base. Call the 
        unmapped region of the read denoted by the string of question marks
        above unmapped_seq. Read M's were previously mapped to reference M's by
        Bowtie. Similarly, read m's were mapped to reference m's. THE ALGORITHM
        ASSUMES THAT THE M's and m's ARE FORWARD-STRAND ALIGNMENTS. B's denote
        extra (presumably intronic) bases of the reference not considered by
        the algorithm. Let K be the number of bases of the read spanned by the
        unmapped region (i.e., the number of question marks). Call the K bases
        of the reference following the lefthand mapped region (represented by
        L's above) left_reference_seq. Similarly, call the K bases of the
        reference preceding the righthand mapped region (represented by R's
        above) right_reference_seq. Note that the lengths of unmapped_seq,
        left_reference_seq, and right_reference_seq are the same. Also note
        that in general, left_reference_seq and right_reference_seq may
        overlap.

        The algo first fills two (K + 1) x (K + 1) score matrices according to
        rules encoded in global_alignment's substitution matrix:
        left_score_matrix scores alignment of unmapped_seq to
        left_reference_seq, while right_score_matrix scores alignment of
        REVERSED unmapped_seq to REVERSED right_reference_seq. Let all matrices
        be 0-indexed, and consider only left_score_matrix. Recall that the
        lower righthand corner element of the matrix is the score of the best
        global alignment. More generally, the (S_L, S_L) element of
        left_score_matrix is the score of the best global alignment of the
        FIRST S_L bases of unmapped_seq to the FIRST S_L bases of
        left_reference_seq. Similar logic would apply to right_score_matrix,
        but first flip right_score_matrix upside-down and leftside-right. Then
        the (K - S_R, K - S_R) element of right_score_matrix is the score of
        the best global alignment of the FINAL S_R bases of unmapped_seq to the
        FINAL S_R bases of right_reference_seq.

        The algo picks out the best "jump" (represented by the pipe bars) from
        left_reference_seq to right_reference_seq for alignment of
        unmapped_seq; that is, it finds the value of S_L that maximizes the sum
        of the scores of the best global alignments of the first S_L bases of
        unmapped_seq to the first S_L bases of left_reference_seq and the final
        (K - S_L) bases of unmapped_seq to the final (K - S_L) bases of
        right_reference_seq. The maximum diagonal element of the matrix
        (left_score_matrix + right_score_matrix) is indexed by (S_L, S_L).
        The maximum may be degenerate --- that is, there may be a tie among
        more than one S_L.

        unmapped_seq: unmapped sequence (string; see paragraph above).
        left_reference_seq: left reference sequence (string; see paragraphs 
            above).
        right_reference_seq: right reference sequence (string; see paragraphs
            above). unmapped_seq, left_reference_seq, and right_reference_seq
            all have the same length.
        global_alignment: instance of GlobalAlignment class used for fast
                alignment. Encapsulates substitution matrix.

        Return value: List of possible S_L (see paragraphs above for
            definition).
    """
    unmapped_seq_size = len(unmapped_seq)
    assert unmapped_seq_size == len(left_reference_seq)
    assert unmapped_seq_size == len(right_reference_seq)

    unmapped_seq = unmapped_seq.upper()
    left_reference_seq = left_reference_seq.upper()
    right_reference_seq = right_reference_seq.upper()

    left_score_matrix = global_alignment.score_matrix(
            unmapped_seq, left_reference_seq
        )
    right_score_matrix = global_alignment.score_matrix(
            unmapped_seq[::-1], right_reference_seq[::-1], 
        )[::-1,::-1]
    total_score_matrix_diagonal = np.diagonal(left_score_matrix) \
                                  + np.diagonal(right_score_matrix)
    return np.argwhere(np.amax(total_score_matrix_diagonal) 
                                == total_score_matrix_diagonal).flatten()

def introns_from_read(reference_index, read_seq, readlets,
    search_for_caps=True, max_discrepancy=2, min_seq_similarity=0.85,
    min_cap_query_size=8, cap_search_window_size=1000, seed=0,
    global_alignment=GlobalAlignment()):
    """ Composes a given read's aligned readlets and returns introns.

        reference_index: object of class bowtie_index.BowtieIndexReference that
            permits access to reference; used for realignment of unmapped
            regions between exonic chunks.
        read_seq: sequence of the original read from which the readlets are
            derived.
        readlets: list of tuples (rname, reverse_strand, pos, end_pos,
            displacement), each of which corresponds to a readlet. rname is the
            SAM-format RNAME, typically the chromosome to which the readlet
            maps. reverse_strand is True iff the readlet's reversed complement
            aligns to the reference. The readlet should span the interval
            [pos, end_pos). displacement is the number of bases between the
            5' (3') end of the readlet, which aligns to the forward (reverse)
            strand, and the 5' (3') end of the read.
        search_for_caps: True iff reference should be searched for the segment
            of a read (a cap) that precedes the first EC and the segment of a
            read that follows the last EC. Such segments are subsequently added
            as an ECs themselves, and introns may be called between them.
        max_discrepancy: if the difference in length between an unmapped region
            framed by two ECs and its corresponding gap in the reference is
            <= this value, the unmapped region is considered a candidate for
            incorporation into a single EC spanning the two original ECs via
            DP filling.
        min_seq_similarity: if the difference in length between an unmapped
            region framed by two ECs and its corresponding gap in the reference
            is <= max_discrepancy AND the score of global alignment is >= 
            min_seq_similarity * (length of unmapped region), the unmapped
            region is incorporated into a single EC spanning the two original
            ECs via DP filling. See the class GlobalAlignment for the
            substitution matrix used.
        min_cap_query_size: the reference is not searched for a cap smaller
            than this size.
        cap_search_window_size: the size (in bp) of the reference subsequence
            in which to search for a cap.
        seed: seed for random number generator; used to break ties in median 
            index when performing DP filling.
        global_alignment: instance of GlobalAlignment class used for fast
                realignment of exonic chunks via Weave.

        Return value: A dictionary where each key is a strand (i.e., a tuple
            (rname, reverse_strand)), and its corresponding value is a list
            of tuples (pos, end_pos)), each of which denotes an intron. rname
            contains the SAM-format RNAME --- typically a chromosome. When
            input reads are strand-specific, reverse_strand is True iff the
            sense strand is the reverse strand; otherwise, it merely denotes
            the strand to which the intron's flanking ECs were presumed to
            align. The intron spans the interval [pos, end_pos).
    """
    random.seed(seed)
    composed = composed_and_sorted_readlets(readlets)
    read_seq = read_seq.upper()
    reversed_complement_read_seq = read_seq[::-1].translate(
            _reversed_complement_translation_table
        )
    read_seq_size = len(read_seq)
    composed_and_capped = {}
    # Add caps
    for strand in composed:
        rname, reverse_strand = strand
        if reverse_strand:
            '''Handle reverse-strand reads the same way forward strands are
            handled.'''
            current_read_seq = reversed_complement_read_seq
        else:
            current_read_seq = read_seq
        new_prefix, new_suffix = [], []
        prefix_pos, prefix_end_pos, prefix_displacement = composed[strand][0]
        suffix_pos, suffix_end_pos, suffix_displacement = composed[strand][-1]
        if prefix_pos - prefix_displacement < 1:
            # Skip this part if we're close to an edge of the reference
            prefix_displacement = 0
        if prefix_displacement:
            if global_alignment.score_matrix(
                        current_read_seq[:prefix_displacement],
                        reference_index.get_stretch(rname,
                            prefix_pos - prefix_displacement - 1,
                            prefix_displacement)
                    )[-1, -1] >= \
                    min_seq_similarity * prefix_displacement:
                new_prefix = [(prefix_pos - prefix_displacement,
                                    prefix_pos, 0)]
            elif search_for_caps \
                and prefix_displacement >= min_cap_query_size \
                and cap_search_window_size > 0:
                '''If region shouldn't be filled and region to the left isn't
                too small, search for prefix.'''
                search_pos = max(0, prefix_pos - 1 - cap_search_window_size)
                search_window = reference_index.get_stretch(
                        rname,
                        search_pos,
                        prefix_pos - 1 - search_pos
                    )
                prefix = re.search(current_read_seq[:min_cap_query_size][::-1],
                                    search_window[::-1])
                if prefix is not None:
                    new_prefix_pos = prefix_pos - prefix.end()
                    new_prefix = [(new_prefix_pos,
                                        new_prefix_pos + min_cap_query_size,
                                        0)]
        unmapped_displacement = (suffix_end_pos - suffix_pos
                                    + suffix_displacement)
        unmapped_base_count = read_seq_size - unmapped_displacement
        if suffix_end_pos + unmapped_base_count - 1 \
            > reference_index.length[rname]:
            # Skip this part if we're too close to the edge of the reference
            unmapped_base_count = 0
        if unmapped_base_count > 0:
            if global_alignment.score_matrix(
                                current_read_seq[unmapped_displacement:],
                                reference_index.get_stretch(rname,
                                    suffix_end_pos - 1, 
                                    unmapped_base_count)
                            )[-1, -1] >= (min_seq_similarity
                                             * unmapped_base_count):
                new_suffix = [(suffix_end_pos,
                                suffix_end_pos + unmapped_base_count,
                                unmapped_displacement)]
            elif search_for_caps \
                and unmapped_base_count >= min_cap_query_size \
                and cap_search_window_size > 0:
                '''If region shouldn't be filled and region to the right isn't
                too small, search for suffix.'''
                search_pos = suffix_end_pos - 1
                search_window = reference_index.get_stretch(
                        rname,
                        search_pos,
                        min(cap_search_window_size, 
                             reference_index.rname_lengths[rname] - search_pos)
                    )
                suffix = re.search(
                                current_read_seq[-min_cap_query_size:],
                                search_window
                            )
                if suffix is not None:
                    new_suffix_pos = suffix_end_pos + suffix.start()
                    new_suffix = [(new_suffix_pos,
                                    new_suffix_pos + min_cap_query_size,
                                    read_seq_size - min_cap_query_size)]
        composed_and_capped[strand] = (new_prefix + composed[strand]
                                        + new_suffix)
    introns = {}
    '''To make continuing outer loop from inner loop below possible.'''
    continue_strand_loop = False
    for strand in composed_and_capped:
        introns[strand] = []
        rname, reverse_strand = strand
        if reverse_strand:
            '''Handle reverse-strand reads the same way forward strands are
            handled.'''
            current_read_seq = reversed_complement_read_seq
        else:
            current_read_seq = read_seq
        last_pos, last_end_pos, last_displacement \
            = composed_and_capped[strand][0]
        unmapped_displacement = last_displacement + last_end_pos - last_pos
        for pos, end_pos, displacement in composed_and_capped[strand][1:]:
            next_unmapped_displacement = displacement + end_pos - pos
            if last_displacement >= displacement or \
                unmapped_displacement >= next_unmapped_displacement or \
                unmapped_displacement > read_seq_size or \
                next_unmapped_displacement > read_seq_size:
                '''Example cases handled:
                The 5' or 3' end of EC #1 isn't to the left of the 5'  or 3'
                end of EC #2 along the read:

                        5'                                       3'
                                             EC #1
                Read      |--------------=================
                                        ===================-----|
                                              EC #2

                                             EC #1
                Read      |----------=======================
                                        ===================-----|
                                              EC #2


                Read      |-----=============------=============|
                                    EC #2              EC #1

                The 3' end of either EC is computed to be displaced beyond the
                3' end of the read:
                                                              EC
                                                     ===================
                Read      |-------------------------------------|

                These situations are pathological. Throw out the strand by
                deleting it from the exon and intron list after continuing.'''
                continue_strand_loop = True
                introns[strand] = []
                break
            reference_distance = pos - last_end_pos
            read_distance = displacement - unmapped_displacement
            discrepancy = reference_distance - read_distance
            call_exon, call_intron = False, False
            if not read_distance:
                if pos - last_end_pos == 0:
                    '''Example case handled:
                                        EC #1                   EC #2
                    Read       |=====================|=====================|
                                                 
                    Reference ...====================|====================...    
                                      EC #1                   EC #2
                                      (EC #1 and #2 are unseparated on
                                        read and reference)
                    The ECs should be merged.'''
                    unmapped_displacement = end_pos - pos + displacement
                    last_end_pos = end_pos
                elif pos - last_end_pos < 0:
                    '''Example case handled:
                                        EC #1                   EC #2
                    Read       |=====================|=====================|
                                                 
                    Reference ...=====================
                                                    =======================...    
                                      EC #1                   EC #2
                                      (EC #1 and #2 may overlap.)
                    This case covers the outside possibility that 
                    EC #1 overlaps EC #2 by as much the constraint that EC #2
                    doesn't begin before EC #1 on the reference allows. This is
                    likely an insertion in the read with respect to the
                    reference. Keep the ECs distinct. CAN BE MODIFIED LATER TO
                    ACCOMMODATE CALLING INDELS.'''
                    call_exon = True
                else:
                    '''pos - last_end_pos > 0
                    Example case handled:
                                       EC #1                   EC #2
                    Read       |=====================|=====================|
                                                     /\
                    Reference ...====================--====================...    
                                     EC #1         intron         EC #2
                    If there are no bases between the ECs in the read sequence 
                    and the size of the candidate intron is > 0, call EC #1 and
                    intron.'''
                    call_exon, call_intron = True, True
            else:
                if read_distance > 0:
                    '''Example case handled:
                                        EC #1        UMR      EC #2
                    Read          |==============-----------========|        
                                               /             \
                    Reference ...==============---------------========...
                                     EC #1         intron       EC #2
                    UMR = unmapped region.
                    The displacement of EC #2's left end is to the right of the
                    displacement of EC #1's right end, and there is an unmapped
                    read region between them. Extend the unmapped region by one
                    nucleotide on either side; splice junctions will have to be
                    determined more precisely. See past if-else.'''
                    last_end_pos -= 1
                    pos += 1
                    unmapped_displacement -= 1
                    displacement += 1
                    read_distance += 2
                else:
                    '''read_distance < 0
                    Example case handled:
                                      EC #1                
                    Read          |=====================        EC #2
                                                ========================|
                    Reference    ...=============------------------=======...
                                        EC #1          intron       EC #2
                    The displacement of EC #2's left end is to the left of the 
                    displacement of EC #1's right end; at least a pair of
                    readlets, one from each EC, overlaps on the read. Call the
                    overlap + one nucleotide on either side the unmapped
                    region; splice junctions will have to be determined more
                    precisely. See past if-else.'''
                    unmapped_displacement, displacement = \
                        displacement - 1, unmapped_displacement + 1
                    read_distance = -read_distance + 2
                    last_end_pos -= read_distance - 1
                    pos += read_distance - 1
                # Now decide whether to DP fill or DP frame
                if abs(discrepancy) <= max_discrepancy:
                    # DP Fill
                    if global_alignment.score_matrix(
                            current_read_seq[
                                unmapped_displacement:displacement
                            ],
                            reference_index.get_stretch(rname,
                                last_end_pos - 1, pos - last_end_pos) 
                        )[-1, -1] >= min_seq_similarity * read_distance:
                        '''If the unmapped region should be filled,
                        last_displacement and last_pos should remain the same
                        on next iteration, but last_end_pos must change; two
                        ECs are merged.'''
                        unmapped_displacement = end_pos - pos + displacement
                        last_end_pos = end_pos
                    else:
                        '''If the unmapped region should not be filled, ignore
                        it, and don't merge ECs. Call only EC #1.'''
                        call_exon = True
                elif discrepancy > max_discrepancy:
                    '''Do DP framing: compute how to distribute unmapped bases
                    between ECs.'''
                    splits = \
                        unmapped_region_splits(
                            current_read_seq[unmapped_displacement:
                                displacement], 
                            reference_index.get_stretch(rname,
                                last_end_pos - 1,
                                read_distance), 
                            reference_index.get_stretch(rname,
                                pos - read_distance - 1, read_distance)
                        )
                    '''Decide which split to use: choose the median index. If 
                    there is a tie in median indices, break it at random.'''
                    split_count = len(splits)
                    split_end = split_count / 2
                    if split_count % 2:
                        split = splits[split_end]
                    else:
                        # If split count is even, break tie
                        split = random.choice(
                                        splits[(split_end-1):(split_end+1)]
                                    )
                    pos += split - read_distance
                    displacement += split - read_distance
                    last_end_pos += split
                    call_exon, call_intron = True, True
                else: 
                    '''discrepancy < -max_discrepancy; is likely a large
                    insertion with respect to reference, but not too large:
                    EC #2 never begins before EC #1 on the reference.
                    Again, call distinct ECs. CAN BE MODIFIED LATER TO
                    ACCOMMODATE CALLING INDELS.

                    Example case handled:
                                      EC #1       UMR          EC #2               
                    Read          |==========--------------==============|

                    Reference        ...==========---==============...
                                           EC #1  gap     EC #2
                    UMR = unmapped region.
                    '''
                    call_exon = True
            if call_intron:
                # Call the reference region between the two ECs an intron.
                introns[strand].append((last_end_pos, pos))
            if call_exon:
                unmapped_displacement = end_pos - pos + displacement
                (last_pos, last_end_pos, last_displacement) = \
                    (pos, end_pos, displacement)
        if continue_strand_loop:
            continue_strand_loop = False
            continue
    strands_to_remove = []
    for strand in introns:
        if introns[strand] == []:
            strands_to_remove.append(strand)
    for strand in strands_to_remove:
        del introns[strand]
    return introns

def selected_readlet_alignments_by_clustering(readlets, seed=0):
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
        seed: seed for randomized algorithm QuickClust to ensure that results
            for a given read are reproducible.

        Return value: a list of selected alignment tuples
            (rname, reverse_strand, pos, end_pos, displacement).
    """
    random.seed(seed)
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
        pivot = i = random.sample(unclustered_alignments, 1)[0]
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

def go(input_stream=sys.stdin, output_stream=sys.stdout, bin_size=10000,
    bowtie_index_base='genome', verbose=False, stranded=False,
    max_discrepancy=2, min_seq_similarity=0.85, min_intron_size=5,
    max_intron_size=500000, intron_partition_overlap=20,
    global_alignment=GlobalAlignment(), search_for_caps=True,
    min_cap_query_size=8, cap_search_window_size=1000):
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
        min_cap_query_size: the reference is not searched for a cap smaller
            than this size.
        cap_search_window_size: the size (in bp) of the reference subsequence
            in which to search for a cap.
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
            filtered_alignments \
                = selected_readlet_alignments_by_clustering(
                                multireadlets,
                                seed=seq[::-1]
                            )
            introns = introns_from_read(
                    reference_index, seq,
                    filtered_alignments,
                    max_discrepancy=max_discrepancy,
                    min_seq_similarity=min_seq_similarity,
                    global_alignment=global_alignment,
                    search_for_caps=search_for_caps,
                    min_cap_query_size=min_cap_query_size,
                    cap_search_window_size=cap_search_window_size,
                    seed=seq
                )
            # Print introns
            for intron_strand in introns:
                intron_rname, intron_reverse_strand = intron_strand
                for intron_pos, intron_end_pos in introns[intron_strand]:
                    if intron_end_pos - intron_pos > max_intron_size:
                        if verbose: 
                            print >>sys.stderr, 'Intron of size > ' \
                                'max-intron-size = %d filtered at ' \
                                '%s:%d-%d' % (max_intron_size, intron_rname,
                                                intron_pos, intron_end_pos)
                        continue
                    if intron_end_pos - intron_pos < min_intron_size:
                        if verbose:
                            print >>sys.stderr, 'Intron of size < ' \
                                'min-intron-size = %d filtered at ' \
                                '%s:%d-%d' % (min_intron_size, intron_rname,
                                                intron_pos, intron_end_pos)
                        continue
                    partitions = partition.partition(intron_rname, intron_pos,
                            intron_pos + 1, bin_size,
                            fudge=intron_partition_overlap
                        )
                    if stranded:
                        intron_reverse_strand_string = '-' if \
                            intron_reverse_strand else '+'
                        intron_strand_string = '+' if \
                            intron_reverse_strand else '-'
                        '''If input reads are stranded, whether the intron
                        is on the forward or reverse strand depends on
                        whether the original sequence was reverse-complemented
                        before alignment.'''
                        if sample_labels:
                            for (partition_id, partition_start, 
                                    partition_end) in partitions:
                                print >>output_stream, \
                                    'intron\t%s%s\t%012d\t%012d\t%s\t%d' \
                                    % (partition_id,
                                        intron_reverse_strand_string,
                                        intron_pos,
                                        intron_end_pos,
                                        sample_labels,
                                        seq_count
                                    )
                                _output_line_count += 1
                        if reversed_complement_sample_labels:
                            for (partition_id, partition_start, 
                                    partition_end) in partitions:
                                print >>output_stream, \
                                    'intron\t%s%s\t%012d\t%012d\t%s\t%d' \
                                    % (partition_id,
                                        intron_strand_string,
                                        intron_pos,
                                        intron_end_pos,
                                        reversed_complement_sample_labels,
                                        reversed_complement_seq_count
                                    )
                                _output_line_count += 1
                    else:
                        '''Strand to which original sequence aligned isn't
                        significant.'''
                        final_labels = '\x1f'.join(
                                set((sample_labels.split('\x1f')
                            if len(sample_labels) else [])
                            + (reversed_complement_sample_labels.split('\x1f')
                            if len(reversed_complement_sample_labels) 
                            else []))
                            )
                        for (partition_id, partition_start, 
                                    partition_end) in partitions:
                            print >>output_stream, \
                                'intron\t%s\t%012d\t%012d\t%s\t%d' \
                                % (partition_id, intron_pos, intron_end_pos,
                                    final_labels, seq_count
                                    + reversed_complement_seq_count
                                )
                            _output_line_count += 1
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
    parser.add_argument('--intron-partition-overlap', type=int, required=False,
        default=20, 
        help='Amount by which partitions overlap their left and right '
             'neighbors')
    parser.add_argument('--min-intron-size', type=int, required=False,
        default=5,
        help='Filters introns of length smaller than this value')
    parser.add_argument('--max-intron-size', type=int, required=False,
        default=500000, 
        help='Filters introns of length greater than this value')
    parser.add_argument('--min-seq-similarity', type=float, required=False,
        default=1., 
        help='If the difference in length between an unmapped region framed '
             'by two ECs and its corresponding gap in the reference is <= the '
             'command-line option --max-discrepancy AND the score of global '
             'alignment is >= --min-seq-similarity * '
             '(length of unmapped region), the unmapped region is '
             'incorporated into a single EC spanning the two original ECs '
             'via DP filling')
    parser.add_argument('--do-not-search_for_caps',
        action='store_const',
        const=True,
        default=False,
        help='Ordinarily, reference is searched for the segment of a read (a '
             'cap) that precedes the first EC and the cap that follows the '
             'last EC. Such caps are subsequently added as ECs themselves. '
             'Use this command-line parameter to turn the feature off')
    parser.add_argument('--min-cap-query-size', type=int, required=False,
        default=8,
        help='The reference is not searched for a segment of a read that '
             'precedes the first EC or follows the last EC smaller than this '
             'size')
    parser.add_argument('--cap-search-window-size', type=int, required=False,
        default=1000,
        help='The size (in bp) of the reference subsequence in which to '
             'search for a cap --- i.e., a segment of a read that follows '
             'the last EC or precedes the first EC.')

    # Add command-line arguments for dependencies
    partition.addArgs(parser)
    bowtie.addArgs(parser)

    # Collect Bowtie arguments, supplied in command line after the -- token
    argv = sys.argv
    bowtie_args = ''
    in_args = False
    for i, argument in enumerate(sys.argv[1:]):
        if in_args:
            bowtie_args += argument + ' '
        if argument == '--':
            argv = sys.argv[:i + 1]
            in_args = True

    '''Now collect other arguments. While the variable args declared below is
    global, properties of args are also arguments of the go() function so
    different command-line arguments can be passed to it for unit tests.'''
    args = parser.parse_args(argv[1:])

if __name__ == '__main__' and not args.test:
    import time
    start_time = time.time()
    '''Compile global alignment code; substitution matrix can be set to
    nondefault here if desired.'''
    global_alignment = GlobalAlignment()
    go(bowtie_index_base=args.bowtie_idx,
        verbose=args.verbose, 
        bin_size=args.partition_length,
        stranded=args.stranded,
        min_seq_similarity=args.min_seq_similarity,
        min_intron_size=args.min_intron_size,
        max_intron_size=args.max_intron_size,
        intron_partition_overlap=args.intron_partition_overlap,
        search_for_caps=(not args.do_not_search_for_caps),
        min_cap_query_size=args.min_cap_query_size,
        cap_search_window_size=args.cap_search_window_size,
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
    
    class TestIntronsFromRead(unittest.TestCase):
        """ Tests introns_from_read(). """
        def setUp(self):
            """ Creates temporary directory and Bowtie index. """
            reference_seq = 'ATGGCATACGATACGTCAGACCATGCAggACctTTacCTACATACTG' \
                            'GCTACATAGTACATCTAGGCATACTACGTgaCATACGgaCTACGTAA' \
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

        def tearDown(self):
            # Kill temporary directory
            shutil.rmtree(self.temp_dir_path)
   
    unittest.main()