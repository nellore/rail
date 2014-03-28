#!/usr/bin/env python
"""
Rail-RNA-align

Follows Rail-RNA-preprocess
Precedes Rail-RNA-intron / Rail-RNA-collapse / Rail-RNA-bam

Alignment script for MapReduce pipelines that wraps Bowtie. Obtains exonic 
chunks from end-to-end alignments. Then obtains candidate introns by dividing
each RNA-seq read into several overlapping segments, or readlets, aligning the
readlets, and -- very roughly -- inferring splice junctions between successive
readlets that align noncontiguously.

Input (read from stdin)
----------------------------
Tab-delimited input tuple columns in a mix of any of the following three
formats:
 Format 1 (single-end, 3-column):
  1. Name
  2. Nucleotide sequence
  3. Quality sequence
 Format 2 (paired-end, 5-column):
  1. Name
  2. Nucleotide sequence for mate 1
  3. Quality sequence for mate 1
  4. Nucleotide sequence for mate 2
  5. Quality sequence for mate 2
 Format 3 (paired, 6-column):
  1. Name for mate 1
  2. Nucleotide sequence for mate 1
  3. Quality sequence for mate 1
  4. Name for mate 2
  5. Nucleotide sequence for mate 2
  6. Quality sequence for mate 2

Hadoop output (written to stdout)
----------------------------
A given RNAME sequence is partitioned into intervals ("bins") of some 
user-specified length (see partition.py).

Exonic chunks (aka ECs; three formats, any or all of which may be emitted):

Format 1 (exon_ival); tab-delimited output tuple columns:
1. Reference name (RNAME in SAM format) + ';' + bin number
2. Sample label
3. EC start (inclusive) on forward strand
4. EC end (exclusive) on forward strand

Format 2 (exon_diff); tab-delimited output tuple columns:
1. Reference name (RNAME in SAM format) + ';' + 
    max(EC start, bin start) (inclusive) on forward strand IFF diff is
    positive and EC end (exclusive) on forward strand IFF diff is negative
2. Bin number
3. Sample label
4. +1 or -1.

Note that only unique alignments are currently output as ivals and/or diffs.

Format 3 (end_to_end_sam); tab-delimited output tuple columns:
Standard SAM output except fields are in different order, and the first field
corresponds to sample label. (Fields are reordered to facilitate partitioning
by sample name/RNAME and sorting by POS.) Each line corresponds to a
spliced alignment. The order of the fields is as follows.
1. Sample label
2. Number string representing RNAME; see BowtieIndexReference class in
    bowtie_index for conversion information
3. POS
4. QNAME
5. FLAG
6. MAPQ
7. CIGAR
8. RNEXT
9. PNEXT
10. TLEN
11. SEQ
12. QUAL
... + optional fields

Introns

Tab-delimited output tuple columns (intron):
1. Reference name (RNAME in SAM format) + ';' + bin number +  
    ('+' or '-' indicating which strand is the sense strand if input reads are
    strand-specific -- that is, --stranded is invoked; otherwise, there is no
    terminal '+' or '-')
2. Intron start (inclusive) on forward strand
3. Intron end (exclusive) on forward strand

Reads with no end-to-end alignments

Tab-delimited output tuple columns (unmapped):
1. QNAME
2. SEQ
3. QUAL

Maximum read lengths found

Tab-delimited output tuple columns (max_len):
1. The character '-', enforcing a single partition.
2. The character 'a', which places it before 'i' in 
    lexicograhic sort order for reading in Rail-RNA-intron_post
3. Maximum read length found
4. The character '0'.
5. The character '0'.

ALL OUTPUT COORDINATES ARE 1-INDEXED.
"""
import sys
import os
import site
import threading
import string
import tempfile
import atexit
import subprocess
import re
import numpy as np
import random
import itertools
# For fast global alignment
from scipy import weave

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

def write_reads(output_stream, input_stream=sys.stdin, readletize=False, 
         min_readlet_size=8, max_readlet_size=25, readlet_interval=5, 
         capping_fraction=.75, verbose=False, report_multiplier=1.2):
    """ Writes input reads/readlets in tab-separated format parsable by Bowtie.

        Mate information and, for a readlet, its position with respect to the 
        parent read's 5' and 3' ends as well as the parent sequence are stored 
        in a semicolon-separated read ID. (This will be QNAME in SAM output.) 
        Unpaired reads are marked 0 after original read ID, while paired-end
        reads are marked 1 and 2.

        Input formats:
        3-token: name TAB seq TAB qual
        5-token: name TAB seq1 TAB qual1 TAB seq2 TAB qual2
        6-token: name1 TAB seq1 TAB qual1 TAB name2 TAB seq2 TAB qual2

        ONLY 3-TOKEN INPUT FORMAT IS VALID WHEN readletize=True.

        output_stream: where to write reads, typically a file stream.
        input_stream: where to retrieve reads, typically sys.stdin on first
            pass of Bowtie and a file stream on second pass.
        readletize: True if reads are to be readletized, otherwise False.
        min_readlet_size: "capping" readlets (that is, readlets that terminate
            at a given end of the read) are never smaller than this value.
            Ignored if readletize=False.
        max_readlet_size: size of every noncapping readlet. Ignored if 
            readletize=False.
        readlet_interval: number of bases separating successive readlets along
            the read. Ignored if readletize=False.
        capping_fraction: successive capping readlets on a given end of a read
            are tapered in size exponentially with fractional base
            capping_fraction. Ignored if readletize=False.
        verbose: True if reads/readlets should occasionally be written 
            to stderr.
        report_multiplier: if verbose is True, the line number of a read or its 
            first readlet written to stderr increases exponentially with base
            report_multiplier.

        No return value.
    """
    next_report_line = 0
    if readletize:
        # Build list of capping readlet sizes
        cap_sizes = []
        cap_size = max_readlet_size
        while cap_size >= min_readlet_size:
            cap_sizes.append(cap_size)
            cap_size = int(cap_size*capping_fraction)
            if cap_size == cap_sizes[-1]:
              cap_size -= 1
        for i, line in enumerate(input_stream):
            tokens = line.rstrip().split('\t')
            '''len(tokens) must = 3; write_reads() is called only by 
            BowtieOutputThread instance after a first pass of Bowtie on 
            full reads has failed. Further, name should already have ;0, ;1, 
            or ;2 tacked on depending on a read's single- or paired-end
            status.'''
            assert len(tokens) == 3
            '''A readlet name is in the following (slightly redundant) format,
                where all semicolons are ASCII unit separators (\x1f):
                original name;0/1/2 (single- or paired-end label);displace-
                ment of readlet's 5' end from read's 5' end;displacement of
                readlet's 3' end from read's 3' end;read sequence;
                quality sequence;number of readlets in read.'''
            to_write = []
            seq_size = len(tokens[1])
            # Add capping readlets
            to_write += [['%s\x1f%d\x1f%d\x1f%s\x1f%s' \
                            % (tokens[0], 0, seq_size - cap_size, tokens[1],
                                tokens[2]),
                            '\t%s\t%s' % (tokens[1][:cap_size],
                                            tokens[2][:cap_size])]
                            for cap_size in cap_sizes]
            to_write += [['%s\x1f%d\x1f%d\x1f%s\x1f%s' \
                            % (tokens[0], seq_size - cap_size, 0, tokens[1],
                                tokens[2]),
                            '\t%s\t%s' % (tokens[1][-cap_size:],
                                            tokens[2][-cap_size:])]
                            for cap_size in cap_sizes]
            # Add noncapping readlets
            to_write += [[
                            '%s\x1f%d\x1f%d\x1f%s\x1f%s' % (tokens[0], j, 
                                            seq_size - j - max_readlet_size,
                                            tokens[1], tokens[2]),
                            '\t%s\t%s' % (tokens[1][j:j + max_readlet_size],
                                            tokens[2][j:j + max_readlet_size])
                            ] for j in range(readlet_interval, 
                                seq_size - max_readlet_size, readlet_interval)]
            # Add total number of readlets to each string
            readlet_count = len(to_write)
            readlet_count_string = '\x1f' + str(readlet_count)
            for j in range(readlet_count):
                to_write[j] = to_write[j][0] + readlet_count_string \
                                + ''.join(to_write[j][1:])
            if verbose and next_report_line == i:
                print >>sys.stderr, 'First readlet from read %d: %s' \
                    % (i + 1, to_write[0])
                next_report_line = int((next_report_line + 1)
                    * report_multiplier + 1) - 1
            for readlet in to_write:
                print >>output_stream, readlet
    else:
        # Don't readletize
        global _input_line_count
        for i, line in enumerate(input_stream):
            tokens = line.rstrip().split('\t')
            if len(tokens) not in (3, 5, 6):
                raise RuntimeError('The following line has an invalid number'
                    ' of tab-separated tokens:\n%sA valid line has 3, 5,'
                    ' or 6 such tokens.' % line)
            # Check that a properly formed label is embedded in the read name
            sample.hasLab(tokens[0], mustHave=True)
            if len(tokens) == 3:
                to_write = '%s\t%s\t%s' \
                    % (tokens[0] + '\x1f0', tokens[1], tokens[2])
                print >>output_stream, to_write
            else:
                to_write = '%s\t%s\t%s\n%s\t%s\t%s' \
                    % (tokens[0] + '\x1f1', tokens[1], tokens[2],
                        (tokens[0] if len(tokens) == 5 else tokens[3]) \
                            + '\x1f2',
                            tokens[-2], tokens[-1])
                print >>output_stream, to_write
            if verbose and next_report_line == i:
                print >>sys.stderr, 'Read(s) %d: %s' % (i, to_write)
                next_report_line = int((next_report_line + 1)
                    * report_multiplier + 1) - 1
            _input_line_count += 1
    output_stream.flush()

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
                        split = random.sample(
                                        splits[(split_end - 1):split_end], 1
                                    )[0]
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

class BowtieOutputThread(threading.Thread):
    """ Processes Bowtie alignments, emitting tuples for exons and introns. """
    
    def __init__(self, input_stream, reference_index, readletized=False,
        unmapped_stream=None, output_stream=sys.stdout,
        exon_differentials=True, exon_intervals=False, stranded=False,
        end_to_end_sam=True, verbose=False, bin_size=10000, min_intron_size=5,
        max_intron_size=100000, max_discrepancy=2, min_seq_similarity=0.85,
        intron_partition_overlap=20, search_for_caps=True,
        min_cap_query_size=8, cap_search_window_size=1000,
        global_alignment=GlobalAlignment(), report_multiplier=1.2):
        """ Constructor for BowtieOutputThread.

            input_stream: where to retrieve Bowtie's SAM output, typically a
                Bowtie process's stdout.
            reference_index: object of class bowtie_index.BowtieIndexReference
                that permits access to reference; used for realignment of
                unmapped regions between exonic chunks.
            readletized: True if input contains readlets; if False, input
                contains whole reads.
            unmapped_stream: where to write reads Bowtie can't map, possibly
                for later readletizing; typically, this is a file stream. No
                unmapped reads are written if unmapped_stream is None. Ignored
                if readletized=True.
            output_stream: where to emit exon and intron tuples; typically,
                this is sys.stdout.
            exon_differentials: True iff EC differentials are to be emitted.
            exon_intervals: True iff EC intervals are to be emitted.
            stranded: True iff input reads are strand-specific; this affects
                whether an output partition has a terminal '+' or '-'
                indicating the sense strand.
            end_to_end_sam: True iff SAM with end_to_end alignments should
                be output. See docstring for more information.
            verbose: True if alignments should occasionally be written 
                to stderr.
            bin_size: genome is partitioned in units of bin_size for later load
                balancing.
            min_intron_size: introns smaller than this number of bases are 
                filtered out.
            max_intron_size: introns larger than this number of bases are
                filtered out.
            max_discrepancy: if the difference in length between an unmapped 
                region framed by two ECs and its corresponding gap in the
                reference is <= this value, the unmapped region is considered a
                candidate for incorporation into a single EC spanning the two
                original ECs via DP filling.
            min_seq_similarity: if the difference in length between an unmapped
                region framed by two ECs and its corresponding gap in the
                reference is <= max_discrepancy AND the score of global
                alignment is >=
                min_seq_similarity * (length of unmapped region), the unmapped
                region is incorporated into a single EC spanning the two
                original ECs via DP filling. See the GlobalAlignment class for
                the substitution matrix used.
            intron_partition_overlap: number of bases to subtract from
                reference start position of intron when determining genome
                partition it is in.
            search_for_caps: True iff reference should be searched for the
                segment of a read (a cap) that precedes the first EC and the
                segment of a read that follows the last EC. Such segments are
                subsequently added as an ECs themselves, and introns may be
                called between them.
            min_cap_query_size: the reference is not searched for a cap smaller
                than this size.
            cap_search_window_size: the size (in bp) of the reference
                subsequence in which to search for a cap.
            global_alignment: instance of GlobalAlignment class used for fast
                realignment of exonic chunks via Weave.
            report_multiplier: if verbose is True, the line number of an
                alignment written to stderr increases exponentially with base
                report_multiplier.
        """
        super(BowtieOutputThread, self).__init__()
        self.daemon = True
        self.input_stream = input_stream
        self.readletized = readletized
        self.unmapped_stream = unmapped_stream
        self.output_stream = output_stream
        self.reference_index = reference_index
        self.verbose = verbose
        self.bin_size = bin_size
        self.stranded = stranded
        self.max_discrepancy = max_discrepancy
        self.min_seq_similarity = min_seq_similarity
        self.report_multiplier = report_multiplier
        self.exon_differentials = exon_differentials
        self.exon_intervals = exon_intervals
        self.end_to_end_sam = end_to_end_sam
        self.max_intron_size = max_intron_size
        self.min_intron_size = min_intron_size
        self.search_for_caps = search_for_caps
        self.min_cap_query_size = min_cap_query_size
        self.cap_search_window_size = cap_search_window_size
        self.intron_partition_overlap = intron_partition_overlap
        self.global_alignment = global_alignment

    def run(self):
        """ Prints exons for reads and exons/introns for readlets.

            Overrides default method containing thread activity.

            No return value.
        """
        global _output_line_count
        next_report_line = 0
        if not self.readletized:
            max_read_size = 0
            '''Next read must be known to tell if a read mapped to multiple
            locations, so always work with previous read.'''
            while True:
                line = self.input_stream.readline()
                if not line: return # Bowtie output nothing
                # Skip header line
                if line[0] == '@': continue
                last_tokens = line.rstrip().split('\t')
                (last_qname, last_flag, last_rname, last_pos, last_mapq,
                    last_cigar, last_rnext, last_pnext,
                    last_tlen, last_seq, last_qual) = last_tokens[:11]
                last_flag = int(last_flag)
                last_pos = int(last_pos)
                last_seq_size = len(last_seq)
                last_end_pos = last_pos + last_seq_size
                '''Find XM:i field, which is > 0 if read had several valid
                alignments, but all were suppressed because bowtie -m X was
                invoked for X an integer >= 1. See Bowtie documentation.'''
                last_multimapped = False
                for field in last_tokens[::-1]:
                    if field[:5] == 'XM:i:':
                        if int(field[5:]) > 0 and (last_flag & 4):
                            '''If read is multimapped and all alignments were
                            suppressed.'''
                            last_multimapped = True
                        break
                break
            # Initialize counter
            i = 0
            # While labeled multiread, this list may end up simply a uniread
            multiread = []
            while True:
                line = self.input_stream.readline()
                if line:
                    tokens = line.rstrip().split('\t')
                    (qname, flag, rname, pos, mapq, cigar, rnext,
                        pnext, tlen, seq, qual) = tokens[:11]
                    max_read_size = max(max_read_size, len(seq))
                    flag = int(flag)
                    pos = int(pos)
                    if self.stranded:
                        '''A reverse-strand string is needed if and only if
                        input reads are strand-specific.'''
                        reverse_strand_string = '-' \
                            if (flag & 16) != 0 else '+'
                    else:
                        reverse_strand_string = ''
                    seq_size = len(seq)
                    end_pos = pos + seq_size
                    '''Find XM:i field, which is > 0 if read had several valid
                    alignments, but all were suppressed because bowtie -m X was
                    invoked for X an integer >= 1. See Bowtie documentation.'''
                    multimapped = False
                    for field in tokens[::-1]:
                        if field[:5] == 'XM:i:':
                            if int(field[5:]) > 0 and (flag & 4):
                                '''If read is multimapped and all alignments
                                were suppressed.'''
                                multimapped = True
                            break
                if self.verbose and next_report_line == i:
                    print >>sys.stderr, \
                        'SAM output record %d: rdname="%s", flag=%d' \
                        % (i, last_qname, last_flag)
                    next_report_line = int((next_report_line + 1)
                        * self.report_multiplier + 1) - 1
                last_sample_label = sample.parseLab(last_qname[:-2])
                multiread.append(last_tokens)
                if not line or qname != last_qname:
                    '''If the next qname doesn't match the last qname or there
                    are no more lines, all of a multiread's alignments have
                    been collected.'''
                    if (last_flag & 4) \
                        and not (len(multiread) > 1 or last_multimapped):
                        if self.unmapped_stream is not None:
                            '''Write only reads with no possible alignments to
                            unmapped_stream.'''
                            print >>self.unmapped_stream, '%s\t%s\t%s' \
                                % (last_qname, last_seq, last_qual)
                            '''Also write unmapped reads for realignment in 
                            another map step.'''
                            print >>self.output_stream, '%s\t%s\t%s\t%s' \
                                % ('unmapped', last_qname, last_seq, last_qual)
                            _output_line_count += 1
                    elif self.end_to_end_sam:
                        '''End-to-end SAM is output for every line with at
                        least one possible alignment.'''
                        for alignment_tokens in multiread:
                            print >>self.output_stream, (
                                ('%s\t%s\t%s\t%012d\t%s\t%s\t'
                                    % ('end_to_end_sam', last_sample_label, 
                                    self.reference_index.rname_to_string[
                                            alignment_tokens[2]
                                        ],
                                    int(alignment_tokens[3]),
                                    alignment_tokens[0][:-2],
                                    alignment_tokens[1])) 
                                + '\t'.join(alignment_tokens[4:])
                                + ('\tNH:i:%d' % len(multiread)))
                            _output_line_count += 1
                    if not (last_flag & 4) and len(multiread) == 1 \
                        and not last_multimapped:
                        '''Read maps uniquely; the full alignment is to be
                        called as an exonic chunk (EC).'''
                        partitions = partition.partition(last_rname, last_pos, 
                            last_end_pos, self.bin_size)
                        if self.exon_differentials:
                            for (partition_id, 
                                partition_start, partition_end)in partitions:
                                # Print increment at interval start
                                assert last_pos < partition_end \
                                    + self.intron_partition_overlap
                                diff_rname, diff_bin = partition_id.split(';')
                                print >>self.output_stream, \
                                    'exon_diff\t%s;%d;%s\t%s\t1' \
                                    % (diff_rname,
                                        max(partition_start, last_pos),
                                        last_sample_label,
                                        diff_bin)
                                _output_line_count += 1
                                assert last_end_pos > partition_start
                                if last_end_pos < partition_end:
                                    '''Print decrement at interval end iff exon
                                    ends before partition ends.'''
                                    print >>self.output_stream, \
                                        'exon_diff\t%s;%d;%s\t%s\t-1' \
                                        % (diff_rname,
                                            last_end_pos,
                                            last_sample_label,
                                            diff_bin)
                                    _output_line_count += 1
                        if self.exon_intervals:
                            for partition_id, _, _ in partitions:
                                print >>self.output_stream, \
                                    'exon_ival\t%s\t%012d\t%012d\t%s' \
                                    % (partition_id, last_pos, 
                                        last_end_pos, last_sample_label)
                                _output_line_count += 1
                    multiread = []
                if not line: 
                    # Write max read size for each strand
                    for max_len_rname in self.reference_index.length:
                        for max_len_strand in ['+', '-']:
                            print >>self.output_stream, \
                                'max_len\t%s%s\ta\t%d\t-' % (max_len_rname,
                                                                max_len_strand,
                                                                max_read_size)
                        _output_line_count += 1
                    break
                last_tokens = tokens
                (last_qname, last_flag, last_rname, last_pos, last_mapq,
                    last_cigar, last_rnext, last_pnext, last_tlen, last_seq,
                    last_qual) = (qname, flag, rname, pos, mapq, cigar,
                    rnext, pnext, tlen, seq, qual)
                (last_seq_size, last_end_pos, last_multimapped) = (seq_size,
                    end_pos, multimapped)
                i += 1
        else:
            # Input is readletized, and readlets must be composed
            '''Create dictionary for gathering aligned readlets belonging
            to same read.'''
            collected_readlets = {}
            '''Create dictionary for counting readlets, aligned or not,
            belonging to same read.'''
            readlet_count = {}
            '''Next readlet must be known to tell if a readlet mapped to
            multiple locations, so always work with previous read.'''
            while True:
                line = self.input_stream.readline()
                if not line: return # Bowtie output nothing
                # Skip header line
                if line[0] == '@': continue
                last_tokens = line.rstrip().split('\t')
                (last_full_qname, last_flag, last_rname, last_pos, last_mapq,
                    last_cigar, last_rnext, last_pnext,
                    last_tlen, last_seq, last_qual) = last_tokens[:11]
                last_flag = int(last_flag)
                last_pos = int(last_pos)
                (last_qname, last_paired_label,
                    last_five_prime_displacement, 
                    last_three_prime_displacement,
                    last_read_seq, last_qual_seq,
                    last_total_readlets) = last_full_qname.split('\x1f')
                last_reverse_strand = (last_flag & 16) != 0
                if last_reverse_strand:
                    '''last_seq is reverse-complemented; displacement of
                    readlet's 3' end from read's 3' end is displacement of
                    readlet's 5' end from read's 5' end on complementary
                    (forward) strand.'''
                    last_displacement = int(last_three_prime_displacement)
                else:
                    last_displacement = int(last_five_prime_displacement)
                break
            # Initialize counter
            i = 0
            '''While it's labeled multireadlet, this list may end up being
            simply a unireadlet.'''
            multireadlet = []
            while True:
                line = self.input_stream.readline()
                if line:
                    tokens = line.rstrip().split('\t')
                    (full_qname, flag, rname, pos, mapq, cigar, rnext,
                        pnext, tlen, seq, qual) = tokens[:11]
                    flag = int(flag)
                    pos = int(pos)
                    '''Recall that a readlet name is in the following format,
                    where all semicolons are ASCII unit separators (\x1f):
                    original name;0/1/2 (single- or paired-end label);
                    displacement of readlet's 5' end from read's 5' end;
                    displacement of readlet's 3' end from read's 3' end;
                    read sequence;quality sequence;
                    number of readlets in read.'''
                    (qname, paired_label, five_prime_displacement, 
                        three_prime_displacement, read_seq, qual_seq,
                        total_readlets) = full_qname.split('\x1f')
                    reverse_strand = (flag & 16) != 0
                    if reverse_strand:
                        '''seq is reverse-complemented; displacement of
                        readlet's 3' end from read's 3' end is displacement of
                        readlet's 5' end from read's 5' end on complementary
                        (forward) strand.'''
                        displacement = int(three_prime_displacement)
                    else:
                        displacement = int(five_prime_displacement)
                if self.verbose and next_report_line == i:
                    print >>sys.stderr, \
                        'SAM output record %d: rdname="%s", flag=%d' \
                        % (i, last_qname, last_flag)
                    next_report_line = int((next_report_line + 1)
                        * self.report_multiplier + 1) - 1
                multireadlet.append((last_rname, last_reverse_strand, last_pos,
                                        last_pos + len(last_seq),
                                        last_displacement))
                if not line or full_qname != last_full_qname:
                    '''If the next qname doesn't match the last qname or there
                    are no more lines, all of a multireadlet's alignments have
                    been collected.'''
                    if (last_qname, last_paired_label) not in readlet_count:
                        readlet_count[(last_qname, last_paired_label)] = 1
                    else:
                        readlet_count[(last_qname, last_paired_label)] += 1
                    if not (last_flag & 4):
                        '''Readlet maps, but perhaps not uniquely; decide which
                        alignment is to be called as an exonic chunk (EC) later
                        using a selected_readlet_alignments function. All
                        output positions will be with respect to 5' end of
                        forward strand. Note that last_seq is
                        reverse-complemented readlet if sense strand is reverse
                        strand, but last_read_seq is NEVER
                        reverse-complemented.'''
                        if (last_qname,
                                last_paired_label) not in collected_readlets:
                            collected_readlets[
                                    (last_qname, last_paired_label)
                                ] = []
                        collected_readlets[
                                    (last_qname, last_paired_label)
                                ].append(multireadlet)
                    last_total_readlets = int(last_total_readlets)
                    if readlet_count[(last_qname, last_paired_label)] \
                        == last_total_readlets and (last_qname,
                        last_paired_label) in collected_readlets:
                        '''Set seed for each read so results for read are
                        reproducible.'''
                        filtered_alignments \
                                = selected_readlet_alignments_by_clustering(
                                        collected_readlets[(last_qname,
                                            last_paired_label)],
                                        seed=(last_qual_seq + last_read_seq)
                                    )
                        introns = introns_from_read(
                            self.reference_index, last_read_seq,
                            filtered_alignments,
                            max_discrepancy=self.max_discrepancy,
                            min_seq_similarity=self.min_seq_similarity,
                            global_alignment=self.global_alignment,
                            search_for_caps=self.search_for_caps,
                            min_cap_query_size=self.min_cap_query_size,
                            cap_search_window_size=self.cap_search_window_size,
                            seed=(last_read_seq + last_qual_seq)
                        )
                        '''Kill mate information in last_qname so rest of
                        pipeline associates mates.'''
                        # Print introns
                        for intron_strand in introns:
                            intron_rname, intron_reverse_strand = intron_strand
                            intron_reverse_strand_string = '-' if \
                                    intron_reverse_strand else '+'
                            for (intron_pos, intron_end_pos) \
                                    in introns[intron_strand]:
                                if intron_end_pos - intron_pos \
                                    > self.max_intron_size:
                                    if self.verbose: 
                                        print >>sys.stderr, \
                                            'Intron of size > ' \
                                            'max-intron-size = %d' \
                                            ' filtered at %s:%d-%d' \
                                            % (self.max_intron_size,
                                                intron_rname, intron_pos,
                                                intron_end_pos)
                                    continue
                                if intron_end_pos - intron_pos \
                                    < self.min_intron_size:
                                    if self.verbose:
                                        print >>sys.stderr, \
                                            'Intron of size < ' \
                                            'min-intron-size = %d' \
                                            ' filtered at %s:%d-%d' \
                                            % (self.min_intron_size,
                                                intron_rname, intron_pos,
                                                intron_end_pos)
                                    continue
                                partitions = partition.partition(intron_rname,
                                    intron_pos, intron_pos + 1, self.bin_size,
                                    fudge=self.intron_partition_overlap)
                                for (partition_id, partition_start, 
                                        partition_end) in partitions:
                                    print >>self.output_stream, \
                                        'intron\t%s%s\t%012d\t%012d' \
                                        % (partition_id,
                                            intron_reverse_strand_string if 
                                            self.stranded else '',
                                            intron_pos,
                                            intron_end_pos)
                                    _output_line_count += 1
                        del readlet_count[(last_qname, last_paired_label)]
                        del collected_readlets[(last_qname, last_paired_label)]
                    multireadlet = []
                if not line: break
                last_tokens = tokens
                (last_full_qname, last_flag, last_rname, last_pos, last_mapq,
                    last_cigar, last_rnext, last_pnext, last_tlen, last_seq,
                    last_qual) = (full_qname, flag, rname, pos, mapq, cigar,
                    rnext, pnext, tlen, seq, qual)
                (last_qname, last_paired_label, last_five_prime_displacement, 
                    last_three_prime_displacement, last_read_seq,
                    last_qual_seq, last_total_readlets) = (qname, paired_label,
                    five_prime_displacement, three_prime_displacement, 
                    read_seq, qual_seq, total_readlets)
                (last_reverse_strand, last_displacement) = (reverse_strand,
                    displacement)
                i += 1

def handle_temporary_directory(archive, temp_dir_path):
    """ Archives or deletes temporary directory.

        archive: directory name, including path, to which temp_dir_path should
            be renamed when script is complete; or None if temp_dir_path
            should be trashed.
        temp_dir_path: path of temporary directory for storing intermediate
            alignments; archived if archive is not None.

        No return value.
    """
    if archive is not None:
        # Rename temporary directory for archiving
        os.rename(temp_dir_path, archive)
    else:
        # Kill temporary directory
        import shutil
        shutil.rmtree(temp_dir_path)

def go(input_stream=sys.stdin, output_stream=sys.stdout, bowtie_exe='bowtie',
    bowtie_index_base='genome', bowtie_args=None,
    bowtie_junction_args="-t --sam-nohead --startverbose -v 0 -a -m 80",
    temp_dir_path=tempfile.mkdtemp(), bin_size=10000, verbose=False,
    read_filename='reads.tsv', readlet_filename='readlets.tsv',
    unmapped_filename='unmapped.tsv', exon_differentials=True,
    exon_intervals=False, end_to_end_sam=True, stranded=False,
    min_readlet_size=8, max_readlet_size=25, readlet_interval=5,
    capping_fraction=.75, max_discrepancy=2, min_seq_similarity=0.85,
    min_intron_size=5, max_intron_size=100000, intron_partition_overlap=20,
    global_alignment=GlobalAlignment(), report_multiplier=1.2,
    search_for_caps=True, min_cap_query_size=8, cap_search_window_size=1000):
    """ Runs Rail-RNA-align.

        Two passes of Bowtie are run. The first attempts to align full reads.
        Mapped reads are emitted as exonic chunks. Unmapped reads are then
        segmented into readlets, which are aligned on second pass of Bowtie.
        Readlets belonging to the same read are used to infer introns and other
        exonic chunks, which are also emitted.

        Input (read from input_stream)---
        Tab-delimited input tuple columns in a mix of any of the following
        three formats:
        Format 1 (single-end, 3-column):
        1. Name
        2. Nucleotide sequence
        3. Quality sequence
        Format 2 (paired-end, 5-column):
        1. Name
        2. Nucleotide sequence for mate 1
        3. Quality sequence for mate 1
        4. Nucleotide sequence for mate 2
        5. Quality sequence for mate 2
        Format 3 (paired, 6-column):
        1. Name for mate 1
        2. Nucleotide sequence for mate 1
        3. Quality sequence for mate 1
        4. Name for mate 2
        5. Nucleotide sequence for mate 2
        6. Quality sequence for mate 2

        Hadoop output (written to stdout)
        ----------------------------
        A given RNAME sequence is partitioned into intervals ("bins") of some 
        user-specified length (see partition.py).

        Exonic chunks (aka ECs; three formats, any or all of which may be
        emitted:

        Format 1 (exon_ival); tab-delimited output tuple columns:
        1. Reference name (RNAME in SAM format) + ';' + bin number
        2. Sample label
        3. EC start (inclusive) on forward strand
        4. EC end (exclusive) on forward strand

        Format 2 (exon_diff); tab-delimited output tuple columns:
        1. Reference name (RNAME in SAM format) + ';' + 
            max(EC start, bin start) (inclusive) on forward strand IFF diff is
            positive and EC end (exclusive) on forward strand IFF diff is
            negative
        2. Bin number
        3. Sample label
        4. +1 or -1.

        Note that only unique alignments are currently output as ivals and/or
        diffs.

        Format 3 (end_to_end_sam); tab-delimited output tuple columns:
        Standard 11-column SAM output except fields are in different order, and
        the first field corresponds to sample label. (Fields are reordered to
        facilitate partitioning by sample name/RNAME and sorting by POS.) Each
        line corresponds to a spliced alignment.
        The order of the fields is as follows.
        1. Sample label
        2. Number string representing RNAME; see BowtieIndexReference class in
            bowtie_index for conversion information
        3. POS
        4. QNAME
        5. FLAG
        6. MAPQ
        7. CIGAR
        8. RNEXT
        9. PNEXT
        10. TLEN
        11. SEQ
        12. QUAL
        ... + optional fields

        Introns

        Tab-delimited output tuple columns (intron):
        1. Reference name (RNAME in SAM format) + ';' + bin number +  
            ('+' or '-' indicating which strand is the sense strand if input
            reads are strand-specific -- that is, --stranded is invoked;
            otherwise, there is no terminal '+' or '-')
        2. Intron start (inclusive) on forward strand
        3. Intron end (exclusive) on forward strand

        Reads with no end-to-end alignments

        Tab-delimited output tuple columns (unmapped):
        1. QNAME
        2. SEQ
        3. QUAL

        Maximum read lengths found

        Tab-delimited output tuple columns (max_len):
        1. The character '-', enforcing a single partition.
        2. The character 'a', which places it before 'i' in 
            lexicograhic sort order for reading in Rail-RNA-intron_post
        3. Maximum read length found
        4. The character '0'.
        5. The character '0'.

        ALL OUTPUT COORDINATES ARE 1-INDEXED.

        input_stream: where to find input reads.
        output_stream: where to emit exonic chunks and introns.
        bowtie_exe: filename of Bowtie executable; include path if not in
            $PATH.
        bowtie_index_base: the basename of the Bowtie index files associated
            with the reference.
        bowtie_args: string containing precisely extra command-line arguments
            to pass to first-pass Bowtie, e.g., "--tryhard --best"; or None.
        bowtie_junction_args: string containing extra command-line arguments
            to pass to second-pass Bowtie.
        temp_dir_path: path of temporary directory for storing intermediate
            alignments
        bin_size: genome is partitioned in units of bin_size for later load
            balancing.
        verbose: True iff more informative messages should be written to
            stderr.
        read_filename: The file, excluding path, for storing reads before
            alignment with Bowtie.
        readlet_filename: The file, excluding path, for storing readlets before
            alignment with Bowtie.
        unmapped_filename: The file, excluding path, for storing reads unmapped
            after Bowtie's first pass.
        min_readlet_size: "capping" readlets (that is, readlets that terminate
            at a given end of the read) are never smaller than this value.
            Ignored if readletize=False.
        max_readlet_size: size of every noncapping readlet. Ignored if 
            readletize=False.
        readlet_interval: number of bases separating successive readlets along
            the read. Ignored if readletize=False.
        capping_fraction: successive capping readlets on a given end of a read
            are tapered in size exponentially with fractional base
            capping_fraction. Ignored if readletize=False.
        exon_differentials: True iff EC differentials are to be emitted.
        exon_intervals: True iff EC intervals are to be emitted.
        end_to_end_sam: True iff SAM with end_to_end alignments should be
            output. See docstring for more information.
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
    import time
    start_time = time.time()

    read_filename = os.path.join(temp_dir_path, read_filename)
    with open(read_filename, 'w') as read_stream:
       write_reads(
            read_stream, input_stream=input_stream, readletize=False, 
            verbose=verbose, report_multiplier=report_multiplier
        )
    # First-pass Bowtie
    bowtie_process, bowtie_command, threads = bowtie.proc(
            bowtieExe=bowtie_exe, bowtieIdx=bowtie_index_base,
            readFn=read_filename, bowtieArgs=bowtie_args, sam=True,
            stdoutPipe=True, stdinPipe=False
        )
    reference_index = bowtie_index.BowtieIndexReference(bowtie_index_base)
    unmapped_filename = os.path.join(temp_dir_path, unmapped_filename)
    with open(unmapped_filename, 'w') as unmapped_stream:
        output_thread = BowtieOutputThread(
                bowtie_process.stdout, reference_index,
                readletized=False, 
                unmapped_stream=unmapped_stream, 
                exon_differentials=exon_differentials, 
                exon_intervals=exon_intervals, 
                bin_size=bin_size,
                verbose=verbose, 
                output_stream=output_stream,
                end_to_end_sam=end_to_end_sam,
                stranded=stranded,
                global_alignment=global_alignment,
                report_multiplier=report_multiplier
            )
        threads.append(output_thread)
        output_thread.start()
        # Join threads to pause execution in main thread
        for thread in threads:
            if verbose: print >>sys.stderr, 'Joining thread...'
            thread.join()

    if verbose: print >>sys.stderr, 'Bowtie\'s first pass is finished.'

    readlet_filename = os.path.join(temp_dir_path, readlet_filename)
    with open(readlet_filename, 'w') as readlet_stream:
        with open(unmapped_filename) as unmapped_stream:
            write_reads(
                readlet_stream, input_stream=unmapped_stream, readletize=True, 
                min_readlet_size=min_readlet_size, 
                max_readlet_size=max_readlet_size,
                readlet_interval=readlet_interval,
                capping_fraction=capping_fraction, verbose=verbose, 
                report_multiplier=report_multiplier
            )
        # Second-pass Bowtie
        bowtie_process, bowtie_command, threads = bowtie.proc(
                bowtieExe=bowtie_exe, bowtieIdx=bowtie_index_base,
                readFn=readlet_filename, bowtieArgs=bowtie_junction_args,
                sam=True, stdoutPipe=True, stdinPipe=False
            )

    output_thread = BowtieOutputThread(
            bowtie_process.stdout, reference_index, 
            readletized=True, unmapped_stream=None,
            exon_differentials=exon_differentials, 
            exon_intervals=exon_intervals,
            bin_size=bin_size,
            verbose=verbose, 
            output_stream=output_stream,
            min_intron_size=min_intron_size,
            max_discrepancy=max_discrepancy, 
            min_seq_similarity=min_seq_similarity, 
            max_intron_size=max_intron_size,
            stranded=stranded,
            intron_partition_overlap=intron_partition_overlap,
            search_for_caps=search_for_caps,
            min_cap_query_size=min_cap_query_size,
            cap_search_window_size=cap_search_window_size,
            global_alignment=global_alignment,
            report_multiplier=report_multiplier
        )
    threads.append(output_thread)
    output_thread.start()
    # Join threads to pause execution in main thread
    for thread in threads:
        if verbose: print >>sys.stderr, 'Joining thread...'
        thread.join()

    if verbose: print >>sys.stderr, 'Bowtie\'s second pass is finished.'
    output_stream.flush()

    print >> sys.stderr, 'DONE with align.py; in/out=%d/%d; ' \
        'time=%0.3f s' % (_input_line_count, _output_line_count,
                                time.time() - start_time)

if __name__ == '__main__':
    import argparse
    # Print file's docstring if -h is invoked
    parser = argparse.ArgumentParser(description=__doc__, 
                formatter_class=argparse.RawDescriptionHelpFormatter)
    # Add command-line arguments
    parser.add_argument('--min-readlet-size', type=int, required=False,
        default=12, 
        help='Capping readlets (that is, readlets that terminate '
             'at a given end of the read) are never smaller than this value')
    parser.add_argument('--max-readlet-size', type=int, required=False,
        default=25, 
        help='Size of every noncapping readlet')
    parser.add_argument('--readlet-interval', type=int, required=False,
        default=12, 
        help='Number of bases separating successive noncapping readlets along '
             'the read')
    parser.add_argument('--capping-fraction', type=float, required=False,
        default=0.85, 
        help='Successive capping readlets on a given end of a read are '
             'tapered in size exponentially with this fractional base')
    parser.add_argument('--report_multiplier', type=float, required=False,
        default=1.2,
        help='When --verbose is also invoked, the only lines of lengthy '
             'intermediate output written to stderr have line number that '
             'increases exponentially with this base')
    parser.add_argument('--verbose', action='store_const', const=True,
        default=False,
        help='Print out extra debugging statements')
    parser.add_argument('--exon-differentials', action='store_const',
        const=True,
        default=True, 
        help='Print exon differentials (+1s and -1s)')
    parser.add_argument('--exon-intervals', action='store_const',
        const=True,
        default=False, 
        help='Print exon intervals')
    parser.add_argument('--end-to-end-sam', action='store_const', const=True,
        default=True, 
        help='Print read-by-read SAM output of end-to-end alignments for '
             'later consolidation after realignment.')
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
    parser.add_argument('--archive', metavar='PATH', type=str, 
        default=None,
        help='Save output and Bowtie command to a subdirectory (named using ' 
             'this process\'s PID) of PATH')

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
    '''Compile global alignment code; substitution matrix can be set to
    nondefault here if desired.'''
    global_alignment = GlobalAlignment()
    temp_dir_path = tempfile.mkdtemp()
    archive = os.path.join(args.archive,
        str(os.getpid())) if args.archive is not None else None
    # Handle temporary directory if CTRL+C'd
    atexit.register(handle_temporary_directory, archive, temp_dir_path)
    if args.verbose:
        print >>sys.stderr, 'Creating temporary directory %s' \
            % temp_dir_path
    go(bowtie_exe=args.bowtie_exe,
        bowtie_index_base=args.bowtie_idx,
        bowtie_args=bowtie_args, temp_dir_path=temp_dir_path, 
        verbose=args.verbose, 
        bin_size=args.partition_length,
        read_filename='reads.tsv', 
        readlet_filename='readlets.tsv',
        unmapped_filename='unmapped.tsv', 
        exon_differentials=args.exon_differentials,
        exon_intervals=args.exon_intervals,
        end_to_end_sam=args.end_to_end_sam,
        stranded=args.stranded,
        min_readlet_size=args.min_readlet_size, 
        max_readlet_size=args.max_readlet_size,
        readlet_interval=args.readlet_interval,
        capping_fraction=args.capping_fraction,
        min_seq_similarity=args.min_seq_similarity, 
        max_intron_size=args.max_intron_size,
        intron_partition_overlap=args.intron_partition_overlap,
        search_for_caps=(not args.do_not_search_for_caps),
        min_cap_query_size=args.min_cap_query_size,
        cap_search_window_size=args.cap_search_window_size,
        global_alignment=global_alignment,
        report_multiplier=args.report_multiplier)
elif __name__ == '__main__':
    # Test units
    del sys.argv[1:] # Don't choke on extra command-line parameters
    import random
    import unittest
    import shutil

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
                            'GTCCAGATTACGATACAAaTACGAAcTCccATAGCAaCATaCTAGac' \
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
                       'GTCCAGATTACGATACAAaTACGAAcTCccATAGCAaCATaCTAGac'
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
                       'GTCCAGATTACGATACAAaTACGAAcTCccATAGCAaCATaCTAGac'
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
                       'GTCCAGATTACGATACAAaTACGAAcTCccATAGCAaCATaCTAGac'
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
    
    class TestWriteReads(unittest.TestCase):
        """ Tests write_reads(). """
        def setUp(self):
            '''input_reads has 5-, 6-, and 3-token examples, each on a
            different line. Every read has 100 bases.'''
            input_reads = \
                'r1;LB:holder' \
                    '\tTTACATACCATACAGTGCGCTAGCGGGTGACAGATATAATGCAGATCCAT' \
                      'ACAGACCAGATGGCAGACATGTGTTGCAGSCTGCAAGTGCAACGCGGTGA' \
                    '\tFFB<9889340///29==:766234466666340///29==:76623446' \
                      '744442<<9889<888@?FFFFFFFFFFFFDB<4444340///29==:76' \
                    '\tGACCAGAGGTGCACCAAATGACAGTGCGCCACAGATGGCGCATGCAGTGG' \
                      'CCAGAAACGTGCGACCGATGACAGGTGCACAATGCCGCTGACGACGTGAA' \
                    '\t444744442<<9889<8880///29==:766230///29==:766230//' \
                      '<9889340///29==:766234466666340///<9889340///29==:\n' \
                'r2/1;LB:holder' \
                    '\tGCAGAGTGCCGCAATGACGTGCGCCAAAGCGGTGACAGGGTGACAGTGAA' \
                      'CCAAGTGACAAGTGAACAGGTGCCAGAGTGACCGAGTGACCAGTGGACCA' \
                    '\t442<<9889<8880///29==:766230//442<<9889<8880///29=' \
                      '44442<<9889<8880///29==:76623044442<<9889<8880///2' \
                '\tr2/2;LB:holder' \
                    '\tCAGAGTGCCGCAATGACGTGCGCCAAAGCGGACAAAGCACCATGACAAGT' \
                      'ACACAGGTGACAGTGACAAGACAGAGGTGACACAGAGAAAGtGGGTGTGA' \
                    '\t<<9889<8880///29==:766230//442<<<<9889<8880///29==' \
                      '44442<<9889<8880///29==:766230///2944442<<9889<888\n' \
                'r3;LB:holder' \
                    '\tATCGATTAAGCTATAACAGATAACATAGACATTGCGCCCATAATAGATAA' \
                      'CTGACACCTGACCAGTGCCAGATGACCAGTGCCAGATGGACGACAGTAGC' \
                    '\tFFFFFFFFFFFFDB<4444340///29==:766234466666777689<3' \
                      '44=<<;444744442<<9889<888@?FFFFFFFFFFFFDB<4444340/\n'
            self.temp_dir_path = tempfile.mkdtemp()
            self.input_file = os.path.join(self.temp_dir_path,
                                'sample_input.tsv')
            # Store reads in temporary file
            with open(self.input_file, 'w') as reads_stream:
                reads_stream.write(input_reads)
            # For storing output of write_reads()
            self.output_file = os.path.join(self.temp_dir_path,
                                'sample_output.tsv')

        def test_unreadletized_output(self):
            """ Fails if output of write_reads() is not in the right form. """
            with open(self.input_file) as input_stream:
                with open(self.output_file, 'w') as output_stream:
                    write_reads(output_stream, input_stream=input_stream,
                                    readletize=False)

            # Verify output
            with open(self.output_file) as processed_stream:
                first_read_from_pair = \
                    processed_stream.readline().rstrip().split('\t')
                second_read_from_pair = \
                    processed_stream.readline().rstrip().split('\t')
                self.assertEquals(
                    [
                        'r1;LB:holder\x1f1',
                        'TTACATACCATACAGTGCGCTAGCGGGTGACAGATATAATGCAGATCCAT'
                        'ACAGACCAGATGGCAGACATGTGTTGCAGSCTGCAAGTGCAACGCGGTGA',
                        'FFB<9889340///29==:766234466666340///29==:76623446'
                        '744442<<9889<888@?FFFFFFFFFFFFDB<4444340///29==:76'
                    ],
                    first_read_from_pair
                )
                self.assertEquals(
                    [
                        'r1;LB:holder\x1f2',
                        'GACCAGAGGTGCACCAAATGACAGTGCGCCACAGATGGCGCATGCAGTGG'
                        'CCAGAAACGTGCGACCGATGACAGGTGCACAATGCCGCTGACGACGTGAA',
                        '444744442<<9889<8880///29==:766230///29==:766230//'
                        '<9889340///29==:766234466666340///<9889340///29==:'
                    ],
                    second_read_from_pair
                )
                first_read_from_pair = \
                    processed_stream.readline().rstrip().split('\t')
                second_read_from_pair = \
                    processed_stream.readline().rstrip().split('\t')
                self.assertEquals(
                    [
                        'r2/1;LB:holder\x1f1',
                        'GCAGAGTGCCGCAATGACGTGCGCCAAAGCGGTGACAGGGTGACAGTGAA'
                        'CCAAGTGACAAGTGAACAGGTGCCAGAGTGACCGAGTGACCAGTGGACCA',
                        '442<<9889<8880///29==:766230//442<<9889<8880///29='
                        '44442<<9889<8880///29==:76623044442<<9889<8880///2'
                    ],
                    first_read_from_pair
                )
                self.assertEquals(
                    [
                        'r2/2;LB:holder\x1f2',
                        'CAGAGTGCCGCAATGACGTGCGCCAAAGCGGACAAAGCACCATGACAAGT'
                        'ACACAGGTGACAGTGACAAGACAGAGGTGACACAGAGAAAGtGGGTGTGA',
                        '<<9889<8880///29==:766230//442<<<<9889<8880///29=='
                        '44442<<9889<8880///29==:766230///2944442<<9889<888'
                    ],
                    second_read_from_pair
                )
                self.assertEquals(
                    [
                        'r3;LB:holder\x1f0',
                        'ATCGATTAAGCTATAACAGATAACATAGACATTGCGCCCATAATAGATAA'
                        'CTGACACCTGACCAGTGCCAGATGACCAGTGCCAGATGGACGACAGTAGC',
                        'FFFFFFFFFFFFDB<4444340///29==:766234466666777689<3'
                        '44=<<;444744442<<9889<888@?FFFFFFFFFFFFDB<4444340/'
                    ],
                    processed_stream.readline().rstrip().split('\t')
                )

        def test_readletized_output(self):
            """ Fails if output of write_reads() is not in the right form. """
            with open(self.input_file) as input_stream:
                '''Read only three-token input line because write_reads()
                takes only three-token input when readletizing.'''
                input_stream.readline() # Kill 5-token input
                input_stream.readline() # Kill 6-token input
                with open(self.output_file, 'w') as output_stream:
                    write_reads(output_stream, input_stream=input_stream,
                                    readletize=True, capping_fraction=0.5,
                                    min_readlet_size=25, readlet_interval=5,
                                    max_readlet_size=50)

            collected_readlets = []
            with open(self.output_file) as processed_stream:
                for readlet in processed_stream:
                    collected_readlets.append(readlet.rstrip().split('\t'))
            '''r1 of input_reads spans 100 bases, and from the arguments passed
            to write_reads() above, noncapping readlets should span 50. The
            capping fraction should arrange for one readlet spanning 25 bases
            on either end of a given read. There should be 13 readlets in
            total. Spot-check some readlets.'''
            # Capping readlets
            self.assertTrue([
                                    'r3;LB:holder\x1f0\x1f50\x1f'
                                    'ATCGATTAAGCTATAACAGATAACA'
                                    'TAGACATTGCGCCCATAATAGATAA'
                                    'CTGACACCTGACCAGTGCCAGATGA'
                                    'CCAGTGCCAGATGGACGACAGTAGC\x1f'
                                    'FFFFFFFFFFFFDB<4444340///'
                                    '29==:766234466666777689<3'
                                    '44=<<;444744442<<9889<888'
                                    '@?FFFFFFFFFFFFDB<4444340/\x1f'
                                    '13',
                                    'ATCGATTAAGCTATAACAGATAACA'
                                    'TAGACATTGCGCCCATAATAGATAA',
                                    'FFFFFFFFFFFFDB<4444340///'
                                    '29==:766234466666777689<3'
                                ] in collected_readlets
                            )
            self.assertTrue([
                                    'r3;LB:holder\x1f75\x1f0\x1f'
                                    'ATCGATTAAGCTATAACAGATAACA'
                                    'TAGACATTGCGCCCATAATAGATAA'
                                    'CTGACACCTGACCAGTGCCAGATGA'
                                    'CCAGTGCCAGATGGACGACAGTAGC\x1f'
                                    'FFFFFFFFFFFFDB<4444340///'
                                    '29==:766234466666777689<3'
                                    '44=<<;444744442<<9889<888'
                                    '@?FFFFFFFFFFFFDB<4444340/\x1f'
                                    '13',
                                    'CCAGTGCCAGATGGACGACAGTAGC',
                                    '@?FFFFFFFFFFFFDB<4444340/'
                                ] in collected_readlets
                            )
            # Noncapping readlets
            self.assertTrue([
                                    'r3;LB:holder\x1f5\x1f45\x1f'
                                    'ATCGATTAAGCTATAACAGATAACA'
                                    'TAGACATTGCGCCCATAATAGATAA'
                                    'CTGACACCTGACCAGTGCCAGATGA'
                                    'CCAGTGCCAGATGGACGACAGTAGC\x1f'
                                    'FFFFFFFFFFFFDB<4444340///'
                                    '29==:766234466666777689<3'
                                    '44=<<;444744442<<9889<888'
                                    '@?FFFFFFFFFFFFDB<4444340/\x1f'
                                    '13',
                                    'TTAAGCTATAACAGATAACA'
                                    'TAGACATTGCGCCCATAATAGATAA'
                                    'CTGAC',
                                    'FFFFFFFDB<4444340///'
                                    '29==:766234466666777689<3'
                                    '44=<<'
                                ] in collected_readlets
                            )
            self.assertTrue([
                                    'r3;LB:holder\x1f40\x1f10\x1f'
                                    'ATCGATTAAGCTATAACAGATAACA'
                                    'TAGACATTGCGCCCATAATAGATAA'
                                    'CTGACACCTGACCAGTGCCAGATGA'
                                    'CCAGTGCCAGATGGACGACAGTAGC\x1f'
                                    'FFFFFFFFFFFFDB<4444340///'
                                    '29==:766234466666777689<3'
                                    '44=<<;444744442<<9889<888'
                                    '@?FFFFFFFFFFFFDB<4444340/\x1f'
                                    '13',
                                    'TAATAGATAA'
                                    'CTGACACCTGACCAGTGCCAGATGA'
                                    'CCAGTGCCAGATGGA',
                                    '66777689<3'
                                    '44=<<;444744442<<9889<888'
                                    '@?FFFFFFFFFFFFD'
                                ] in collected_readlets
                            )

        def tearDown(self):
            # Kill temporary directory
            shutil.rmtree(self.temp_dir_path)
   
    unittest.main()