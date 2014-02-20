#!/usr/bin/env python
"""
Rail-RNA-align

Follows Rail-RNA-preprocess
Precedes Rail-RNA-intron / Rail-RNA-coverage_pre

Alignment script for MapReduce pipelines that wraps Bowtie. Obtains a set of
of exonic chunks by aligning entire reads with Bowtie. Then obtains more exonic 
chunks as well as introns by dividing each RNA-seq read into several
overlapping segments, or readlets, aligning the readlets, and -- very roughly
-- inferring splice junctions between successive readlets that align 
noncontiguously.

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
----------------------------

Hadoop output (written to stdout)
----------------------------
A given RNAME sequence is partitioned into intervals ("bins") of some 
user-specified length (see partition.py).

Exonic chunks (aka ECs; two formats --- either or both may be emitted):

Format 1 (exon_ival); tab-delimited output tuple columns:
1. Reference name (RNAME in SAM format) + ';' + bin number
2. Sample label
3. EC start (inclusive) on forward strand
4. EC end (exclusive) on forward strand
5. '+' or '-' indicating which strand is the sense strand if input reads are
        strand-specific -- that is, --stranded is invoked; otherwise, there is
        no terminal '+' or '-'

Format 2 (exon_diff); tab-delimited output tuple columns:
1. Reference name (RNAME in SAM format) + ';' + bin number
2. Sample label
3. max(EC start, bin start) (inclusive) on forward strand IFF next column is +1 
   and EC end (exclusive) on forward strand IFF next column is -1.
4. '+' or '-' indicating which strand is the sense strand if input reads are
        strand-specific -- that is, --stranded is invoked; otherwise, there is
        no terminal '+' or '-'
5. +1 or -1.

Introns

Tab-delimited output tuple columns (intron):
1. Reference name (RNAME in SAM format) + ';' + bin number +  
    ('+' or '-' indicating which strand is the sense strand if input reads are
        strand-specific -- that is, --stranded is invoked; otherwise, there is
        no terminal '+' or '-')
2. Sample label
3. Intron start (inclusive) on forward strand
4. Intron end (exclusive) on forward strand
5. Number of nucleotides between 5' end of intron and 5' end of read from which
it was inferred, ASSUMING THE SENSE STRAND IS THE FORWARD STRAND. That is, if
the sense strand is the reverse strand, this is the distance between the 3' end
of the read and the 3' end of the intron.
6. Number of nucleotides between 3' end of intron and 3' end of read from which
it was inferred, ASSUMING THE SENSE STRAND IS THE FORWARD STRAND.
7. Number of nucleotides spanned by EC on the left (that is, towards the 5'
end of the read) of the intron, ASSUMING THE SENSE STRAND IS THE FORWARD
STRAND.
8. Number of nucleotides spanned by EC on the right (that is, towards the 3'
end of the read) of the intron, ASSUMING THE SENSE STRAND IS THE FORWARD
STRAND.
9. Match rate: the number of bases that match the reference per base of the
aligned readlets comprising the read from which the intron was inferred.
10. '-' if reversed complements of readlets from which intron was inferred
aligned to forward strand; else '+'

SAM (splice_sam):
Standard 11-column SAM output except fields are in different order. (Fields are
    reordered to facilitate partitioning by RNAME and sorting by POS.) Each
line corresponding to a spliced alignment.
The order of the fields is as follows.
1. RNAME
2. POS
3. QNAME
4. FLAG
5. MAPQ
6. CIGAR
7. RNEXT
8. PNEXT
9. TLEN
10. SEQ
11. QUAL

ALL OUTPUT COORDINATES ARE 1-INDEXED.
"""
import sys
import os
import site
import threading
import copy
import string
import tempfile
import atexit
import subprocess
import re
import numpy as np
# For fast global alignment
from scipy import weave

base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for directory_name in ['bowtie', 'sample', 'alignment', 'interval']:
    site.addsitedir(os.path.join(base_path, directory_name))

import bowtie
import bowtie_index
import sample
import partition

_reversed_complement_translation_table = string.maketrans('ATCG', 'TAGC')

class GlobalAlignment:
    """ Invokes Weave to obtain alignment score matrix with C. """

    def __init__(self, substitution_matrix=[[ 1,-1,-1,-1,-1,-2],
                                            [-1, 1,-1,-1,-1,-2],
                                            [-1,-1, 1,-1,-1,-2],
                                            [-1,-1,-1, 1,-1,-2],
                                            [-1,-1,-1,-1,-1,-2],
                                            [-2,-2,-2,-2,-2,-2]]):
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
                                'gcc, and try again.'
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
input_line_count = 0
output_line_count = 0

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
        global input_line_count
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
            input_line_count += 1
    output_stream.flush()

def composed_and_sorted_readlets(readlets, min_strand_readlets=1):
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
        min_strand_readlets: readlets are composed on a strand iff the number
            of readlets that align to the strand is >= this number.

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
        if len(uncomposed[strand]) < min_strand_readlets:
            '''If a strand has one or two stray readlets, they may be
            bad alignments, so kill them.'''
            to_delete.append(strand)
            continue
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

def exons_and_introns_from_read(reference_index, read_seq, readlets, 
    min_intron_size=5, min_strand_readlets=1, max_discrepancy=2, 
    min_seq_similarity=0.85, search_for_caps=True, min_cap_query_size=8,
    cap_search_window_size=1000,
    global_alignment=GlobalAlignment()):
    """ Composes a given read's aligned readlets and returns ECs and introns.

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
        min_intron_size: introns smaller than this number of bases are filtered
            out.
        min_strand_readlets: exons and introns are called on a strand (i.e., a 
            tuple (rname, reverse_strand)) iff the number of readlets that
            align to the strand is >= this number.
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

        Return value: tuple (exons, introns).
            -exons is a dictionary; each key is a strand (i.e., a tuple
                (rname, reverse_strand)), and its corresponding value is a list
                of tuples (pos, end_pos), each of which denotes an exonic chunk
                (EC). rname is the SAM-format RNAME---typically a chromosome.
                When input reads are strand-specific, reverse_strand is True
                iff the sense strand is the reverse strand; otherwise, it
                merely denotes the strand to which the EC was presumed
                to align. The EC spans the interval [pos, end_pos).
            -introns is a dictionary; each key is a strand (i.e., a tuple
                (rname, reverse_strand)), and its corresponding value is a list
                of tuples (pos, end_pos, five_prime_displacement,
                three_prime_displacement, left_EC_size, right_EC_size), each of
                which denotes an intron. rname contains the SAM-format RNAME
                --- typically a chromosome. When input reads are
                strand-specific, reverse_strand is True iff the sense strand is
                the reverse strand; otherwise, it merely denotes the strand to
                which the intron's flanking ECs were presumed to align. The
                intron spans the interval [pos, end_pos). Assume the sense
                strand is the forward strand; that is, if the sense strand is
                the reverse strand, consider the reversed complements of both
                the readlet and its parent read. five_prime_displacement is the
                displacement of the 5' end of the intron from the 5' end of the
                read, while three_prime_displacement is the displacement of the
                3' end of the intron from the 3' end of the read. left_EC_size
                (right_EC_size) is the number of nucleotides spanned by the EC
                on the left (right) of the intron. Here, "left" ("right")
                means "towards the 5' (3') end of the read," again assuming
                the sense strand is the forward strand.
    """
    composed = composed_and_sorted_readlets(readlets, min_strand_readlets)
    read_seq = read_seq.upper()
    reversed_complement_read_seq = read_seq[::-1].translate(
            _reversed_complement_translation_table
        )
    exons_and_introns = {}
    '''To make continuing outer loop from inner loop below possible.'''
    continue_strand_loop = False
    for strand in composed:
        exons_and_introns[strand] = []
        rname, reverse_strand = strand
        if reverse_strand:
            '''Handle reverse-strand reads the same way forward strands are
            handled.'''
            current_read_seq = reversed_complement_read_seq
        else:
            current_read_seq = read_seq
        read_seq_size = len(read_seq)
        last_pos, last_end_pos, last_displacement = composed[strand][0]
        loop_start_index = 1
        if last_displacement:
            '''If last_displacement isn't 0, check if region to the left should
            be filled.'''
            if global_alignment.score_matrix(
                        current_read_seq[:last_displacement],
                        reference_index.get_stretch(rname,
                            last_pos - last_displacement - 1,
                            last_displacement)
                    )[-1, -1] >= min_seq_similarity * last_displacement:
                # Fill
                last_pos -= last_displacement
                last_displacement = 0
            elif search_for_caps and last_displacement >= min_cap_query_size \
                and cap_search_window_size > 0:
                '''If region shouldn't be filled and region to the left isn't
                too small, search for it.'''
                search_pos = max(0, last_pos - 1 - cap_search_window_size)
                search_window = reference_index.get_stretch(
                        rname,
                        search_pos,
                        last_pos - 1 - search_pos
                    )
                prefix = re.search(current_read_seq[:last_displacement][::-1],
                                    search_window[::-1])
                if prefix is not None:
                    if prefix.start() != 0:
                        # If the cap is found, tack on an extra EC.
                        last_pos -= prefix.end()
                        last_end_pos = last_pos + last_displacement
                        last_displacement = 0
                        loop_start_index = 0
                    else:
                        '''Merge ECs (this should have been taken care of by
                        DP filling, but just in case....)'''
                        last_pos -= last_displacement
                        last_displacement = 0
                        loop_start_index = 1
        unmapped_displacement = last_displacement + last_end_pos - last_pos
        for pos, end_pos, displacement in composed[strand][loop_start_index:]:
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
                exons_and_introns[strand] = []
                break
            reference_distance = pos - last_end_pos
            read_distance = displacement - unmapped_displacement
            discrepancy = reference_distance - read_distance
            call_exon, call_intron = False, False
            if not read_distance:
                if pos - last_end_pos < min_intron_size:
                    '''Example case handled:
                                        EC #1                   EC #2
                    Read       |=====================|=====================|
                                                 
                    Reference ...====================|====================...    
                                      EC #1                   EC #2
                                      (EC #1 and #2 are allowed to
                                        overlap or be separated by up 
                                        to min_intron_size in reference)
                    Note that this case subsumes the outside possibility that 
                    EC #1 overlaps EC #2 by as much the constraint that EC #2
                    doesn't begin before EC #1 on the reference allows. This is
                    likely an insertion in the read with respect to the
                    reference. Keep the ECs distinct. CAN BE MODIFIED LATER TO
                    ACCOMMODATE CALLING INDELS.'''
                    call_exon = True
                else:
                    '''pos - last_end_pos >= min_intron_size
                    Example case handled:
                                       EC #1                   EC #2
                    Read       |=====================|=====================|
                                                     /\
                    Reference ...====================--====================...    
                                     EC #1         intron         EC #2
                    If there are no bases between the ECs in the read sequence 
                    and the size of the candidate intron is >= min_intron_size,
                    call EC #1 and intron.'''
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
                    '''Decide which split to use: choose the smallest
                    medoid index.'''
                    nonintegral_median = np.median(splits)
                    split = splits[np.argmin(abs(splits - nonintegral_median))]
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
            if call_exon:
                '''Push ONLY the first EC to the list. (The next EC
                may get merged into another EC on another iteration.)'''
                exons_and_introns[strand].append((True, last_pos,
                                                    last_end_pos))
            if call_intron:
                # Call the reference region between the two ECs an intron.
                exons_and_introns[strand].append((False, last_end_pos, pos,
                                        displacement, 
                                        len(read_seq) - displacement))
            if call_exon:
                unmapped_displacement = end_pos - pos + displacement
                (last_pos, last_end_pos, last_displacement) = \
                    (pos, end_pos, displacement)
        if continue_strand_loop:
            continue_strand_loop = False
            continue
        '''Final exonic chunk must still be called since it won't get merged 
        with any other chunks. Check to see if region to right of chunk needs 
        filling first. The method for determining unmapped_base_count may
        give a negative result. This pathological case is thrown out by
        skipping to the next strand. If the fill criterion is not met, 
        search for the residual unmapped bases downstream.'''
        unmapped_base_count = read_seq_size - unmapped_displacement
        if unmapped_base_count < 0:
            exons_and_introns[strand] = []
            continue
        if unmapped_base_count > 0:
            if global_alignment.score_matrix(
                                current_read_seq[unmapped_displacement:],
                                reference_index.get_stretch(rname,
                                    last_end_pos - 1, 
                                    unmapped_base_count)
                            )[-1, -1] >= (min_seq_similarity
                                             * unmapped_base_count):
                # Fill
                last_end_pos += unmapped_base_count
            elif search_for_caps \
                and unmapped_base_count >= min_cap_query_size \
                and cap_search_window_size > 0:
                '''If region shouldn't be filled and region to the right isn't
                too small, search for it.'''
                search_pos = last_end_pos - 1
                search_window = reference_index.get_stretch(
                        rname,
                        search_pos,
                        min(cap_search_window_size, 
                             reference_index.rname_lengths[rname] - search_pos)
                    )
                suffix = re.search(
                                current_read_seq[unmapped_displacement:],
                                search_window
                            )
                if suffix is not None:
                    if suffix.start() != 0:
                        # Call last EC from list composed[strand].
                        exons_and_introns[strand].append((True, last_pos,
                                                            last_end_pos))
                        # Call intron
                        exons_and_introns[strand].append(
                                (False, last_end_pos,
                                    last_end_pos + suffix.start(),
                                    len(read_seq) - unmapped_base_count, 
                                    unmapped_base_count)
                            )
                        last_pos  = last_end_pos + suffix.start()
                    # If suffix.start() is 0, the next assignment merges ECs
                    last_end_pos = last_pos + unmapped_base_count
        exons_and_introns[strand].append((True, last_pos, last_end_pos))
    # Separate exon and intron lists
    exons = {}
    introns = {}
    for strand, strand_exons_and_introns in exons_and_introns.items():
        exons_and_introns_count = len(strand_exons_and_introns)
        if exons_and_introns_count == 0:
            continue
        for i, exon_or_intron in enumerate(strand_exons_and_introns):
            if exon_or_intron[0]:
                # Call exon
                if strand not in exons:
                    exons[strand] = []
                exons[strand].append(exon_or_intron[1:])
            else:
                # Call intron; it should always be sandwiched between two exons
                assert i != 0 and i != exons_and_introns_count - 1
                assert strand_exons_and_introns[i-1][0] \
                        and strand_exons_and_introns[i+1][0]
                if strand not in introns:
                    introns[strand] = []
                introns[strand].append(exon_or_intron[1:] 
                      + (strand_exons_and_introns[i-1][2] 
                            - strand_exons_and_introns[i-1][1],
                         strand_exons_and_introns[i+1][2] 
                             - strand_exons_and_introns[i+1][1])
                    )
    return exons, introns

def selected_readlet_alignments_by_distance(readlets):
    """ Selects multireadlet alignment closest to set of unireadlet alignments.

        Consider a list "readlets" whose items {R_i} correspond to the aligned
        readlets from a given read. Each R_i is itself a list of the possible
        alignments {R_ij} of a readlet. Each R_ij is a tuple
        (rname, reverse_strand, pos, end_pos, displacement, mismatch_count),
        where rname is the SAM-format rname (typically a chromosome),
        reverse_strand is True iff the readlet's reversed complement aligns to
        the reference, displacement is the number of bases between the 5' (3')
        end of the readlet, which aligns to the forward (reverse) strand, and
        the 5' (3') end of the read, and mismatch_count is the number of
        mismatched bases in the alignment. Let K_i be the number of alignments
        {R_ij} of a given readlet R_i.

        The algo first separates readlets into two sets: those that align
        uniquely (unireadlets U_i, a subset of K_i) and those that don't
        (multireadlets M_i, a subset of K_i). It then chooses the multireadlet
        alignment closest (in bp) to the set of unireadlet alignments; if more
        than one alignment of a multireadlet overlaps any of the unireadlets,
        the alignment whose overlap with any of the unireadlets is greatest is
        chosen. More specifically, let {M_ij} be the set of alignments of M_i,
        and slightly abusing notation, let U_i correspond to not just a
        unireadlet, but also its alignment. Denote the interval spanned by some
        alignment X as [S_X, E_X). For a given multireadlet M_i, the algorithm
        selects M_ij by computing

        argmax_j max_k ( min(E_(U_k), E_(M_ij)) - max(S_(U_k), S_(M_ij)) ) .

        After multireadlet alignments are chosen, the function returns only
        those alignments from the strand that has accrued the most alignments.
        If there are no unireadlets, the algorithm returns an empty list of
        alignments.

        readlets: a list whose items {R_i} correspond to the aligned readlets
            from a given read. Each R_i is itself a list of the possible
            alignments {R_ij} of a readlet. Each R_ij is a tuple
            (rname, reverse_strand, pos, end_pos, displacement,
                mismatch_count). See above for a detailed explanation.

        Return value: a list of selected alignment tuples
            (rname, reverse_strand, pos, end_pos, displacement,
                mismatch_count).
    """
    # Final readlets is first populated with unireadlets
    unireadlets = {}
    multireadlets = []
    for readlet in readlets:
        if len(readlet) == 1:
            # Add unireadlet, where the key is (rname, reverse_strand)
            if (readlet[0][0], readlet[0][1]) not in unireadlets:
                unireadlets[(readlet[0][0], readlet[0][1])] = []
            unireadlets[(readlet[0][0], readlet[0][1])].append(readlet[0])
        else:
            multireadlets.append(readlet)
    if len(unireadlets.keys()) == 0:
        '''If there are no unireadlets, the algo doesn't work, so throw out
        read by returning an empty list.'''
        return []
    if len(multireadlets) == 0:
        '''If there are no multireadlets, no selection has to be performed,
        so just return alignments immediately.'''
        return sorted(unireadlets.values(), key=len, reverse=True)[0]
    final_readlets = copy.deepcopy(unireadlets)
    # Find multireadlet alignment with "closest" unireadlet
    for multireadlet in multireadlets:
        last_overlap = None
        alignment = None
        for (rname, reverse_strand, pos, end_pos,
                displacement, mismatch_count) in multireadlet:
            if (rname, reverse_strand) not in unireadlets: continue
            overlap = max([min(end_pos, compared_end_pos)
                            - max(pos, compared_pos)
                            for _, _, compared_pos, compared_end_pos, _, _
                            in unireadlets[(rname, reverse_strand)]])
            if last_overlap is None or overlap > last_overlap:
                alignment = (rname, reverse_strand, pos, end_pos,
                                displacement, mismatch_count)
            last_overlap = overlap
        if alignment is not None:
            final_readlets[(alignment[0], alignment[1])].append(alignment)
    return sorted(final_readlets.values(), key=len, reverse=True)[0]

def selected_readlet_alignments_by_coverage(readlets):
    """ Selects multireadlet alignment that spans region with highest coverage.
    
        Consider a list "readlets" whose items {R_i} correspond to the aligned
        readlets from a given read. Each R_i is itself a list of the possible
        alignments {R_ij} of a readlet. Each R_ij is a tuple
        (rname, reverse_strand, pos, end_pos, displacement, mismatch_count),
        where rname is the SAM-format rname (typically a chromosome),
        reverse_strand is True iff the readlet's reversed complement aligns to
        the reference, displacement is the number of bases between the 5' (3')
        end of the readlet, which aligns to the forward (reverse) strand, and
        the 5' (3') end of the read, and mismatch_count is the number of
        mismatched bases in the alignment. Let K_i be the number of alignments
        {R_ij} of a given readlet R_i.

        The algo first constructs a coverage distribution from the {R_i}. Each
        base position B of the reference spanned by a given R_ij contributes
        1 / K_i to the coverage at B. The algo then selects the j for which
        R_ij spans the region with the highest coverage, where a region's
        coverage is computed by summing the coverages over all its bases.
        If there is more than one such j, the corresponding R_ij are
        recorded in a new, typically smaller subset of possible multireadlet
        alignments. selected_readlet_alignments_by_distance() is then applied
        to the narrowed set of alignments for all readlets to resolve ties.
        After multireadlet alignments are chosen, the function returns only
        those alignments from the strand that has accrued the most alignments.

        readlets: a list whose items {R_i} correspond to the aligned readlets
            from a given read. Each R_i is itself a list of the possible
            alignments {R_ij} of a readlet. Each R_ij is a tuple
            (rname, reverse_strand, pos, end_pos, displacement,
                mismatch_count). See above for a detailed explanation.

        Return value: a list of selected alignment tuples
            (rname, reverse_strand, pos, end_pos, displacement,
                mismatch_count).
    """
    # Construct coverage distribution
    coverage = {}
    for multireadlet in readlets:
        coverage_unit = 1. / len(multireadlet)
        for (rname, reverse_strand, pos, end_pos,
                displacement, mismatch_count) in multireadlet:
            if (rname, reverse_strand) not in coverage:
                coverage[(rname, reverse_strand)] = {}
            for covered_base_pos in xrange(pos, end_pos):
                coverage[(rname, reverse_strand)][covered_base_pos] \
                    = coverage[(rname, reverse_strand)].get(
                            covered_base_pos, 0
                        ) + coverage_unit
    filtered_readlets = []
    '''Choose alignment of multireadlet with highest total coverage. If there
    is a tie among top alignments, .'''
    for multireadlet in readlets:
        assert len(multireadlet) >= 1
        if len(multireadlet) == 1:
            filtered_readlets.append([multireadlet[0]])
            continue
        alignments = []
        for (rname, reverse_strand, pos, end_pos,
                displacement, mismatch_count) in multireadlet:
            readlet_coverage = 0
            for covered_base_pos in xrange(pos, end_pos):
                readlet_coverage += coverage[(rname, reverse_strand)].get(
                    covered_base_pos, 0)
            alignments.append((readlet_coverage, (rname, reverse_strand, pos,
                                                    end_pos, displacement, 
                                                    mismatch_count)))
        # Find highest-coverage alignments
        alignments.sort(reverse=True)
        highest_coverage_alignments = [alignments[0][1]]
        for alignment in alignments[1:]:
            if alignment[0] == alignments[0][0]:
                highest_coverage_alignments.append(alignment[1])
            else:
                break
        filtered_readlets.append(highest_coverage_alignments)
    return selected_readlet_alignments_by_distance(filtered_readlets)

class BowtieOutputThread(threading.Thread):
    """ Processes Bowtie alignments, emitting tuples for exons and introns. """
    
    def __init__(self, input_stream, reference_index, readletized=False,
        unmapped_stream=None, output_stream=sys.stdout,
        exon_differentials=True, exon_intervals=False, stranded=False,
        splice_sam=True, verbose=False, bin_size=10000, min_intron_size=5,
        min_strand_readlets=1, max_discrepancy=2, min_seq_similarity=0.85,
        max_intron_size=100000, intron_partition_overlap=20,
        assign_multireadlets_by_distance=False,
        search_for_caps=True, min_cap_query_size=8,
        cap_search_window_size=1000,
        global_alignment=GlobalAlignment(), 
        report_multiplier=1.2):
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
            splice_sam: True iff SAM with spliced alignments is to be emitted.
                Ignored if readletized is False.
            verbose: True if alignments should occasionally be written 
                to stderr.
            bin_size: genome is partitioned in units of bin_size for later load
                balancing.
            min_intron_size: introns smaller than this number of bases are 
                filtered out.
            min_strand_readlets: exons and introns are called on a strand iff 
                the number of readlets that align to the strand is >= this
                number.
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
            max_intron_size: an intron of that spans more than this number of
                bases is suppressed from output.
            intron_partition_overlap: number of bases to subtract from
                reference start position of intron when determining genome
                partition it is in.
            assign_multireadlets_by_distance: False iff multireadlet alignments
                should be selected primarily based on readlet coverage via
                selected_readlet_alignments_by_coverage(). True iff
                multireadlet alignments should be selected based on distance to
                unireadlets via selected_readlet_alignments_by_distance().
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
        self.min_strand_readlets = min_strand_readlets
        self.max_discrepancy = max_discrepancy
        self.min_seq_similarity = min_seq_similarity
        self.report_multiplier = report_multiplier
        self.exon_differentials = exon_differentials
        self.exon_intervals = exon_intervals
        self.splice_sam = splice_sam
        self.max_intron_size = max_intron_size
        self.search_for_caps = search_for_caps
        self.min_cap_query_size = min_cap_query_size
        self.cap_search_window_size = cap_search_window_size
        self.assign_multireadlets_by_distance \
            = assign_multireadlets_by_distance
        self.intron_partition_overlap = intron_partition_overlap
        self.global_alignment = global_alignment

    def run(self):
        """ Prints exons for reads and exons/introns for readlets.

            Overrides default method containing thread activity.

            No return value.
        """
        global output_line_count
        next_report_line = 0
        if not self.readletized:
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
                if self.stranded:
                    '''A reverse-strand string is needed if and only if
                    input reads are strand-specific.'''
                    last_reverse_strand_string = '-' \
                        if (last_flag & 16) != 0 else '+'
                else:
                    last_reverse_strand_string = ''
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
                            '''If read is multimapped and all alignments were
                            suppressed.'''
                            multimapped = True
                        break
                if self.verbose and next_report_line == i:
                    print >>sys.stderr, \
                        'SAM output record %d: rdname="%s", flag=%d' \
                        % (i, last_qname, last_flag)
                    next_report_line = int((next_report_line + 1)
                        * self.report_multiplier + 1) - 1
                multiread.append((last_qname, last_flag, last_rname, last_pos,
                    last_mapq, last_cigar, last_rnext, last_pnext, last_tlen,
                    last_seq, last_qual))
                if not line or qname != last_qname:
                    '''If the next qname doesn't match the last qname or there
                    are no more lines, all of a multiread's alignments have
                    been collected.'''
                    if not (last_flag & 4) and len(multiread) == 1 \
                        and not last_multimapped:
                        '''Read maps uniquely; the full alignment is to be
                        called as an exonic chunk (EC). Kill mate information
                        in qname so rest of pipeline associates mates.'''
                        last_sample_label = sample.parseLab(last_qname[:-2])
                        partitions = partition.partition(last_rname, last_pos, 
                            last_end_pos, self.bin_size)
                        if self.exon_differentials:
                            for (partition_id, 
                                partition_start, partition_end) in partitions:
                                # Print increment at interval start
                                assert last_pos < partition_end \
                                    + self.intron_partition_overlap
                                print >>self.output_stream, \
                                    'exon_diff\t%s\t%s\t%012d\t%s\t1' \
                                    % (partition_id, 
                                        last_sample_label,
                                        max(partition_start, last_pos),
                                        last_reverse_strand_string)
                                output_line_count += 1
                                assert last_end_pos > partition_start
                                if last_end_pos < partition_end:
                                    '''Print decrement at interval end iff exon
                                    ends before partition ends.'''
                                    print >>self.output_stream, \
                                        'exon_diff\t%s\t%s\t%012d\t%s\t-1' \
                                        % (partition_id, 
                                            last_sample_label, last_end_pos,
                                            last_reverse_strand_string)
                                    output_line_count += 1
                        if self.exon_intervals:
                            for partition_id, _, _ in partitions:
                                print >>self.output_stream, \
                                    'exon_ival\t%s\t%012d\t%012d\t%s\t%s' \
                                    % (partition_id, 
                                        last_reverse_strand_string, last_pos, 
                                        last_end_pos, last_sample_label)
                                output_line_count += 1
                        if self.splice_sam:
                            print >>self.output_stream, ('%s\t'*2 + '%012d\t'
                                    + '%s\t'*8 + '%s') \
                                    % ('splice_sam', last_rname,
                                        last_pos, last_qname[:-2],
                                        last_flag, last_mapq,
                                        str(last_seq_size) + 'M', last_rnext,
                                        last_pnext, last_tlen,
                                        last_seq, last_qual)
                    elif self.unmapped_stream is not None \
                        and not (len(multiread) > 1 or last_multimapped):
                        '''Write only reads with no possible alignments to
                        unmapped_stream.'''
                        print >>self.unmapped_stream, '%s\t%s\t%s' \
                            % (last_qname, last_seq, last_qual)
                    multiread = []
                if not line: break
                last_tokens = tokens
                (last_qname, last_flag, last_rname, last_pos, last_mapq,
                    last_cigar, last_rnext, last_pnext, last_tlen, last_seq,
                    last_qual) = (qname, flag, rname, pos, mapq, cigar,
                    rnext, pnext, tlen, seq, qual)
                (last_seq_size, last_reverse_strand_string, last_end_pos,
                    last_multimapped) = (seq_size, reverse_strand_string,
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
                '''Find MD:Z field to obtain number of mismatches in
                alignment.'''
                last_mismatch_count = 0
                for field in last_tokens[::-1]:
                    if field[:5] == 'MD:Z:':
                        last_mismatch_count \
                            = len(re.findall(r'[A-Za-z]', field[5:]))
                        break
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
                    '''Find MD:Z field to obtain number of mismatches in
                    alignment.'''
                    mismatch_count = 0
                    for field in last_tokens[::-1]:
                        if field[:5] == 'MD:Z:':
                            mismatch_count \
                                = len(re.findall(r'[A-Za-z]', field[5:]))
                            break
                if self.verbose and next_report_line == i:
                    print >>sys.stderr, \
                        'SAM output record %d: rdname="%s", flag=%d' \
                        % (i, last_qname, last_flag)
                    next_report_line = int((next_report_line + 1)
                        * self.report_multiplier + 1) - 1
                multireadlet.append((last_rname, last_reverse_strand, last_pos,
                                        last_pos + len(last_seq),
                                        last_displacement,
                                        last_mismatch_count))
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
                        '''Choose algorithm for selecting alignments from
                        multireadlets.'''
                        if self.assign_multireadlets_by_distance:
                            filtered_alignments \
                                = selected_readlet_alignments_by_distance(
                                        collected_readlets[(last_qname,
                                                            last_paired_label)]
                                    )
                        else:
                            filtered_alignments \
                                = selected_readlet_alignments_by_coverage(
                                        collected_readlets[(last_qname,
                                                            last_paired_label)]
                                    )
                        '''Compute mean number of matched bases per aligned
                        base of the readlets from the read.'''
                        aligned_base_count = sum([alignment[3] - alignment[2]
                            for alignment in filtered_alignments])
                        mismatched_base_count = sum([alignment[-1]
                            for alignment in filtered_alignments])
                        try:
                            match_rate = (aligned_base_count 
                                             - float(mismatched_base_count)) \
                                            / aligned_base_count
                        except ZeroDivisionError: match_rate = 0
                        exons, introns = exons_and_introns_from_read(
                            self.reference_index, last_read_seq,
                            [alignment[:-1] 
                                for alignment in filtered_alignments],
                            min_strand_readlets=self.min_strand_readlets,
                            max_discrepancy=self.max_discrepancy,
                            min_seq_similarity=self.min_seq_similarity,
                            global_alignment=self.global_alignment,
                            search_for_caps=self.search_for_caps,
                            min_cap_query_size=self.min_cap_query_size,
                            cap_search_window_size=self.cap_search_window_size
                        )
                        '''Kill mate information in last_qname so rest of
                        pipeline associates mates.'''
                        last_sample_label = sample.parseLab(last_qname[:-2])
                        if self.exon_differentials:
                            for exon_strand in exons:
                                exon_rname, exon_reverse_strand = exon_strand
                                if self.stranded:
                                    '''A reverse-strand string is needed iff
                                    input reads are strand-specific.'''
                                    exon_reverse_strand_string = '-' if \
                                        exon_reverse_strand else '+'
                                else:
                                    exon_reverse_strand_string = ''
                                for exon_pos, exon_end_pos \
                                    in exons[exon_strand]:
                                    partitions = partition.partition(
                                        exon_rname, exon_pos, exon_end_pos,
                                        self.bin_size)
                                    for (partition_id, partition_start, 
                                            partition_end) in partitions:
                                        assert exon_pos < partition_end
                                        # Print increment at interval start
                                        print >>self.output_stream, \
                                            'exon_diff\t%s%s\t%s\t%012d\t1' \
                                            % (partition_id, 
                                                exon_reverse_strand_string,
                                                last_sample_label, 
                                                max(partition_start, exon_pos))
                                        output_line_count += 1
                                        assert exon_end_pos > partition_start
                                        if exon_end_pos < partition_end:
                                            '''Print decrement at interval end 
                                            iff exon ends before partition
                                            ends.'''
                                            print >>self.output_stream, \
                                                'exon_diff\t%s%s\t%s\t' \
                                                '%012d\t-1' \
                                                % (partition_id, 
                                                    exon_reverse_strand_string,
                                                    last_sample_label,
                                                    exon_end_pos)
                                            output_line_count += 1
                        if self.exon_intervals:
                            for exon_strand in exons:
                                exon_rname, exon_reverse_strand = exon_strand
                                if self.stranded:
                                    '''A reverse-strand string is needed iff
                                    input reads are strand-specific.'''
                                    exon_reverse_strand_string = '-' if \
                                        exon_reverse_strand else '+'
                                else:
                                    exon_reverse_strand_string = ''
                                for exon_pos, exon_end_pos \
                                    in exons[exon_strand]:
                                    partitions = partition.partition(
                                        exon_rname, exon_pos, exon_end_pos,
                                        self.bin_size)
                                    for partition_id, _, _ in partitions:
                                        print >>self.output_stream, \
                                            'exon_ival\t%s%s\t%012d\t' \
                                            '%012d\t%s' \
                                            % (partition_id, 
                                                exon_reverse_strand_string,
                                                exon_pos, exon_end_pos, 
                                                last_sample_label)
                                        output_line_count += 1
                        # Print introns
                        for intron_strand in introns:
                            intron_rname, intron_reverse_strand = intron_strand
                            intron_reverse_strand_string = '-' if \
                                    intron_reverse_strand else '+'
                            for (intron_pos, intron_end_pos,
                                    intron_five_prime_displacement, 
                                    intron_three_prime_displacement,
                                    intron_left_EC_size,
                                    intron_right_EC_size) \
                                    in introns[intron_strand]:
                                if intron_end_pos - intron_pos \
                                    > self.max_intron_size:
                                    if self.verbose: 
                                        print >> sys.stderr, \
                                            'Intron of size > ' \
                                            'max-intron-size = %d' \
                                            ' filtered at %s:%d-%d' \
                                            % (self.max_intron_size,
                                                intron_rname, intron_pos,
                                                intron_end_pos)
                                        output_line_count += 1
                                        continue
                                partitions = partition.partition(intron_rname,
                                    intron_pos, intron_pos + 1, self.bin_size,
                                    fudge=self.intron_partition_overlap)
                                for (partition_id, partition_start, 
                                        partition_end) in partitions:
                                    print >>self.output_stream, \
                                        'intron\t%s%s\t%s\t%012d\t%012d' \
                                        '\t%d\t%d\t%d\t%d\t%.12f\t%s' \
                                        % (partition_id,
                                            intron_reverse_strand_string if 
                                            self.stranded else '',
                                            last_sample_label, intron_pos,
                                            intron_end_pos, 
                                            intron_five_prime_displacement,
                                            intron_three_prime_displacement,
                                            intron_left_EC_size,
                                            intron_right_EC_size,
                                            match_rate, 
                                            intron_reverse_strand_string)
                                    output_line_count += 1
                        if self.splice_sam:
                            '''Print SAM with each line encoding candidate
                            introns from a single read. NO INTRONS ARE
                            FILTERED.'''
                            for intron_strand in introns:
                                (intron_rname,
                                    intron_reverse_strand) = intron_strand
                                strand_introns = introns[intron_strand]
                                cigar_list = []
                                cigar_list.append('%dM' \
                                                    % strand_introns[0][2])
                                last_intron = strand_introns[0]
                                for intron in strand_introns[1:]:
                                    cigar_list.append('%dN' \
                                        % (last_intron[1] - last_intron[0]))
                                    cigar_list.append('%dM' \
                                        % (intron[2] - last_intron[2]))
                                    last_intron = intron
                                cigar_list.append('%dN' \
                                    % (strand_introns[-1][1] 
                                        - strand_introns[-1][0]))
                                cigar_list.append('%dM' % \
                                                    strand_introns[-1][3])
                                '''First column is "splice_sam"; other 11
                                columns are standard SAM. See SAM format
                                specification for details.'''
                                if intron_reverse_strand:
                                    # Need to reverse-complement sequence
                                    splice_sam_tuple = ('splice_sam',
                                        intron_rname, strand_introns[0][0]
                                            - strand_introns[0][2],
                                        last_qname, '16', '255', 
                                        ''.join(cigar_list), '*', '0', '0', 
                                        last_read_seq[::-1].translate(
                                         _reversed_complement_translation_table
                                        ), last_qual_seq[::-1])
                                else:
                                    splice_sam_tuple = ('splice_sam',
                                        intron_rname, strand_introns[0][0]
                                            - strand_introns[0][2],
                                        last_qname, '0', '255', 
                                        ''.join(cigar_list), '*', '0', '0', 
                                        last_read_seq, last_qual_seq)
                                print >>self.output_stream, \
                                    ('%s\t'*2 + '%012d\t' + '%s\t'*8 + '%s') \
                                    % splice_sam_tuple
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
                    last_qual_seq, last_total_readlets,
                    last_mismatch_count) = (qname, paired_label,
                    five_prime_displacement, three_prime_displacement, 
                    read_seq, qual_seq, total_readlets, mismatch_count)
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

def go(input_stream=sys.stdin, output_stream=sys.stdout, bowtie_exe="bowtie",
    bowtie_index_base="genome", bowtie_args=None,
    temp_dir_path=tempfile.mkdtemp(), bin_size=10000, verbose=False,
    read_filename='reads.tsv', readlet_filename='readlets.tsv',
    unmapped_filename='unmapped.tsv', exon_differentials=True,
    exon_intervals=False, splice_sam=True, stranded=False, min_readlet_size=8,
    max_readlet_size=25, readlet_interval=5, capping_fraction=.75,
    min_strand_readlets=1, max_discrepancy=2, min_seq_similarity=0.85,
    min_intron_size=5, max_intron_size=100000, intron_partition_overlap=20, 
    assign_multireadlets_by_distance=False,
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

        Output (written to output_stream)---
        A given RNAME sequence is partitioned into intervals ("bins") of some 
        user-specified length (see partition.py) for later load balancing.

        Exonic chunks (aka ECs; two formats --- either or both may be emitted):

        Format 1 (exon_ival); tab-delimited output tuple columns:
        1. Reference name (RNAME in SAM format) + ';' + bin number +  
        ('+' or '-' indicating which strand is the sense strand if input reads
        are strand-specific -- that is, --stranded is invoked; otherwise, there
        is no terminal '+' or '-')
        2. Sample label
        3. EC start (inclusive) on forward strand
        4. EC end (exclusive) on forward strand

        Format 2 (exon_differential); tab-delimited output tuple columns:
        1. Reference name (RNAME in SAM format) + ';' + bin number +  
        ('+' or '-' indicating which strand is the sense strand if input reads
        are strand-specific -- that is, --stranded is invoked; otherwise, there
        is no terminal '+' or '-')
        2. Sample label
        3. max(EC start, bin start) (inclusive) on forward strand IFF next 
        column is +1 and EC end (exclusive) on forward strand IFF next column
        is -1.
        4. +1 or -1.

        Introns

        Tab-delimited output tuple columns:
        1. Reference name (RNAME in SAM format) + ';' + bin number +  
        ('+' or '-' indicating which strand is the sense strand if input reads
        are strand-specific -- that is, --stranded is invoked; otherwise, there
        is no terminal '+' or '-')
        2. Sample label
        3. Intron start (inclusive) on forward strand
        4. Intron end (exclusive) on forward strand
        5. Number of nucleotides between 5' end of intron and 5' end of read
        from which it was inferred, ASSUMING THE SENSE STRAND IS THE FORWARD
        STRAND. That is, if the sense strand is the reverse strand, this is the
        distance between the 3' end of the read and the 3' end of the intron.
        6. Number of nucleotides between 3' end of intron and 3' end of read
        from which it was inferred, ASSUMING THE SENSE STRAND IS THE FORWARD
        STRAND.
        7. Number of nucleotides spanned by EC on the left (that is, towards
        the 5' end of the read) of the intron, ASSUMING THE SENSE STRAND IS THE
        FORWARD STRAND.
        8. Number of nucleotides spanned by EC on the right (that is, towards
        the 3' end of the read) of the intron, ASSUMING THE SENSE STRAND IS THE
        FORWARD STRAND.
        9. Match rate: the number of bases that match the reference per base of
        the aligned readlets comprising the read from which the intron was
        inferred.
        10. '-' if reversed complements of readlets from which intron was
        inferred aligned to forward strand; else '+'

        SAM (splice_sam):
        Standard 11-column SAM output except fields are in different order.
        (Fields are reordered to facilitate partitioning by RNAME and sorting
            by POS.) Each line corresponding to a spliced alignment.
        The order of the fields is as follows.
        1. RNAME
        2. POS
        3. QNAME
        4. FLAG
        5. MAPQ
        6. CIGAR
        7. RNEXT
        8. PNEXT
        9. TLEN
        10. SEQ
        11. QUAL

        ALL OUTPUT COORDINATES ARE 1-INDEXED.

        input_stream: where to find input reads.
        output_stream: where to emit exonic chunks and introns.
        bowtie_exe: filename of Bowtie executable; include path if not in
            $PATH.
        bowtie_index_base: the basename of the Bowtie index files associated
            with the reference.
        bowtie_args: string containing precisely extra command-line arguments
            to pass to Bowtie, e.g., "--tryhard --best"; or None.
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
        splice_sam: True iff SAM with spliced alignments is to be emitted.
        stranded: True iff input reads are strand-specific; this affects
            whether an output partition has a terminal '+' or '-' indicating
            the sense strand.
        min_strand_readlets: exons and introns are called on a strand iff 
            the number of readlets that align to the strand is >= this number.
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
        assign_multireadlets_by_distance: False iff multireadlet alignments
            should be selected based on readlet coverage via
            selected_readlet_alignments_by_coverage(). True iff multireadlet
            alignments should be selected based on distance to unireadlets
            via selected_readlet_alignments_by_distance().
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
                min_intron_size=args.min_intron_size,
                min_strand_readlets=min_strand_readlets, 
                max_discrepancy=max_discrepancy, 
                min_seq_similarity=min_seq_similarity, 
                max_intron_size=max_intron_size,
                stranded=stranded,
                intron_partition_overlap=intron_partition_overlap,
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
                readFn=readlet_filename, bowtieArgs=bowtie_args, sam=True,
                stdoutPipe=True, stdinPipe=False
            )

    output_thread = BowtieOutputThread(
            bowtie_process.stdout, reference_index, 
            readletized=True, unmapped_stream=None,
            exon_differentials=exon_differentials, 
            exon_intervals=exon_intervals,
            splice_sam=splice_sam, 
            bin_size=bin_size,
            verbose=verbose, 
            output_stream=output_stream,
            min_intron_size=args.min_intron_size,
            min_strand_readlets=min_strand_readlets, 
            max_discrepancy=max_discrepancy, 
            min_seq_similarity=min_seq_similarity, 
            max_intron_size=max_intron_size,
            stranded=stranded,
            intron_partition_overlap=intron_partition_overlap,
            assign_multireadlets_by_distance=assign_multireadlets_by_distance,
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
        'time=%0.3f s' % (input_line_count, output_line_count,
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
    parser.add_argument('--splice-sam', action='store_const', const=True,
        default=True, 
        help='Print read-by-read SAM output of candidate introns for later '
             'consolidation in splice_sam.py. --max-intron-size is ignored.')
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
        default=200000, 
        help='Filters introns of length greater than this value')
    parser.add_argument('--min-strand-readlets', type=int, required=False,
        default=1, 
        help='Exons and introns are called on a strand (i.e., a tuple ' 
        '(rname, reverse_strand)) iff the number of readlets that align to '
        'the strand is >= this number')
    parser.add_argument('--max-discrepancy', type=int, required=False,
        default=2, 
        help='If the difference in length between an unmapped region framed '
             'by two ECs and its corresponding gap in the reference is <= '
             'the value, the unmapped region is considered a candidate for '
             'incorporation into a single EC spanning the two original ECs '
             'via DP filling')
    parser.add_argument('--min-seq-similarity', type=float, required=False,
        default=0.85, 
        help='If the difference in length between an unmapped region framed '
             'by two ECs and its corresponding gap in the reference is <= the '
             'command-line option --max-discrepancy AND the score of global '
             'alignment is >= --min-seq-similarity * '
             '(length of unmapped region), the unmapped region is '
             'incorporated into a single EC spanning the two original ECs '
             'via DP filling')
    parser.add_argument('--assign-multireadlets-by-distance',
        action='store_const',
        const=True,
        default=False, 
        help='Don\'t use readlet coverage to determine which alignment of '
             'multireadlet should be selected. Instead, choose alignment of '
             'multireadlet closest to set of unireadlets. If False, alignment '
             'that spans region with highest coverage is selected.')
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
        splice_sam=args.splice_sam,
        stranded=args.stranded,
        min_readlet_size=args.min_readlet_size, 
        max_readlet_size=args.max_readlet_size,
        readlet_interval=args.readlet_interval,
        capping_fraction=args.capping_fraction, 
        min_strand_readlets=args.min_strand_readlets,
        max_discrepancy=args.max_discrepancy,
        min_seq_similarity=args.min_seq_similarity, 
        max_intron_size=args.max_intron_size,
        intron_partition_overlap=args.intron_partition_overlap,
        assign_multireadlets_by_distance=args.assign_multireadlets_by_distance,
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
        def test_stray_readlet_filtering(self):
            """ Fails if all readlets are not filtered. """
            self.assertEqual(
                    composed_and_sorted_readlets([
                            ('chr1', True, 3, 20, 5), 
                            ('chr2', False, 2, 10, 4),
                            ('chr2', False, 3, 12, 5)
                        ],
                        min_strand_readlets=3
                    ),
                    {}
                )
        def test_overlapping_readlets_1(self):
            """ Fails if readlets are not consolidated as expected. """
            self.assertEqual(
                    composed_and_sorted_readlets([
                            ('chr1', True, 3, 20, 5), 
                            ('chr2', False, 2, 10, 4),
                            ('chr2', False, 3, 12, 5)
                        ],
                        min_strand_readlets=2
                    ),
                    {
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
                        ],
                        min_strand_readlets=2
                    ),
                    {
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
                        ],
                        min_strand_readlets=1
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

                        ],
                        min_strand_readlets=2
                    ),
                    {
                        ('chr1', True): [(3, 7, 8)],
                        ('chr3', False): [(1, 59, 6), (121, 145, 10)]
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
    
    class TestExonsAndIntronsFromRead(unittest.TestCase):
        """ Tests exons_and_introns_from_read(). """
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
            exons, introns = exons_and_introns_from_read(
                                self.reference_index, read_seq, readlets
                            )
            self.assertEquals({('chr1', False) : [(1, 48), (95, 142)]},
                                exons)
            self.assertEquals({('chr1', False) : [(48, 95, 47, 47, 47, 47)]}, 
                                introns)
            '''Now try truncating readlets so there is an unmapped region of
            the read. This tests the DP framing code in context.'''
            readlets = [('chr1', False, 1, 42, 0),
                          ('chr1', False, 100, 142, 52)]
            exons, introns = exons_and_introns_from_read(
                                self.reference_index, read_seq, readlets
                            )
            self.assertEquals({('chr1', False) : [(1, 48), (95, 142)]},
                                exons)
            self.assertEquals({('chr1', False) : [(48, 95, 47, 47, 47, 47)]},
                                introns)
            # Truncate readlets again to test DP framing.
            readlets = [('chr1', False, 1, 37, 0),
                          ('chr1', False, 105, 142, 57)]
            exons, introns = exons_and_introns_from_read(
                                self.reference_index, read_seq, readlets
                            )
            self.assertEquals({('chr1', False) : [(1, 48), (95, 142)]},
                                exons)
            self.assertEquals({('chr1', False) : [(48, 95, 47, 47, 47, 47)]},
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
            exons, introns = exons_and_introns_from_read(
                                self.reference_index, read_seq, readlets
                            )
            self.assertEquals({('chr1', False) : [(1, 48), (95, 142)]},
                                exons)
            self.assertEquals({('chr1', False) : [(48, 95, 47, 47, 47, 47)]},
                                introns)
            '''Now try truncating readlets so there is an unmapped region of
            the read. This tests the DP framing code in context.'''
            readlets = [('chr1', False, 1, 42, 0),
                          ('chr1', False, 100, 142, 52)]
            exons, introns = exons_and_introns_from_read(
                                self.reference_index, read_seq, readlets
                            )
            self.assertEquals({('chr1', False) : [(1, 48), (95, 142)]},
                                exons)
            self.assertEquals({('chr1', False) : [(48, 95, 47, 47, 47, 47)]},
                                introns)
            # Truncate readlets again to test DP framing.
            readlets = [('chr1', False, 1, 37, 0),
                          ('chr1', False, 105, 142, 57)]
            exons, introns = exons_and_introns_from_read(
                                self.reference_index, read_seq, readlets
                            )
            self.assertEquals({('chr1', False) : [(1, 48), (95, 142)]},
                                exons)
            self.assertEquals({('chr1', False) : [(48, 95, 47, 47, 47, 47)]},
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
            exons, introns = exons_and_introns_from_read(
                                self.reference_index, read_seq, readlets
                            )
            self.assertEquals({('chr1', True) : [(1, 48), (95, 142)]},
                                exons)
            self.assertEquals({('chr1', True) : [(48, 95, 47, 47, 47, 47)]},
                                introns)
            '''Now try truncating readlets so there is an unmapped region of
            the read. This tests the DP framing code in context.'''
            readlets = [('chr1', True, 1, 42, 0),
                          ('chr1', True, 100, 142, 52)]
            exons, introns = exons_and_introns_from_read(
                                self.reference_index, read_seq, readlets
                            )
            self.assertEquals({('chr1', True) : [(1, 48), (95, 142)]},
                                exons)
            self.assertEquals({('chr1', True) : [(48, 95, 47, 47, 47, 47)]},
                                introns)
            # Truncate readlets again to test DP framing.
            readlets = [('chr1', True, 1, 37, 0),
                          ('chr1', True, 105, 142, 57)]
            exons, introns = exons_and_introns_from_read(
                                self.reference_index, read_seq, readlets
                            )
            self.assertEquals({('chr1', True) : [(1, 48), (95, 142)]},
                                exons)
            self.assertEquals({('chr1', True) : [(48, 95, 47, 47, 47, 47)]},
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
            exons, introns = exons_and_introns_from_read(
                                self.reference_index, read_seq, readlets
                            )
            self.assertEquals({('chr1', False) : [(48, 95), (142, 189)]},
                                exons)
            self.assertEquals({('chr1', False) : [(95, 142, 47, 47, 47, 47)]},
                                introns)
            '''Now try truncating readlets so there is an unmapped region of
            the read. This tests the DP framing code in context.'''
            readlets = [('chr1', False, 48, 90, 0),
                          ('chr1', False, 147, 189, 52)]
            exons, introns = exons_and_introns_from_read(
                                self.reference_index, read_seq, readlets
                            )
            self.assertEquals({('chr1', False) : [(48, 95), (142, 189)]},
                                exons)
            self.assertEquals({('chr1', False) : [(95, 142, 47, 47, 47, 47)]},
                                introns)
            # Truncate readlets again to test DP framing.
            readlets = [('chr1', False, 48, 87, 0),
                          ('chr1', False, 150, 189, 55)]
            exons, introns = exons_and_introns_from_read(
                                self.reference_index, read_seq, readlets
                            )
            self.assertEquals({('chr1', False) : [(48, 95), (142, 189)]},
                                exons)
            self.assertEquals({('chr1', False) : [(95, 142, 47, 47, 47, 47)]},
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
            exons, introns = exons_and_introns_from_read(
                                self.reference_index, read_seq, readlets
                            )
            self.assertEquals({('chr1', False) : [(189, 236)]}, exons)
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
            exons, introns = exons_and_introns_from_read(
                                self.reference_index, read_seq, readlets, 
                                min_seq_similarity=0
                            )
            self.assertEquals({('chr1', False) : [(236, 283)]}, exons)
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
            exons, introns = exons_and_introns_from_read(
                                self.reference_index, read_seq, readlets
                             )
            self.assertEquals(exons, {('chr1', False) : [(1, 9),
                                                         (15, 29),
                                                         (38, 46)]})
            self.assertEquals(introns, 
                                {('chr1', False) : [(9, 15, 8, 22, 8, 14),
                                                    (29, 38, 22, 8, 14, 8)]})

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
            exons, introns = exons_and_introns_from_read(
                                self.reference_index, read_seq, readlets
                             )
            self.assertEquals(exons, {('chr1', False) : [(236, 244),
                                                         (250, 267),
                                                         (273, 283)]})
            self.assertEquals(introns, 
                                {('chr1', False) : [(244, 250, 8, 27, 8, 17),
                                                    (267, 273, 25, 10,
                                                        17, 10)]})

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
            exons, introns = exons_and_introns_from_read(
                                self.reference_index, read_seq, readlets
                            )
            self.assertEquals({}, exons)
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

    class TestSelectedReadletAlignmentsByDistance(unittest.TestCase):
        """ Tests selected_readlet_alignments_by_distance(). """

        def test_that_best_overlapped_alignment_is_chosen(self):
            """ Fails if proper alignment of multireadlet is not chosen.

                This test in particular checks if an alignment whose overlap
                with any of the unireadlets is greatest is chosen.
            """
            multireadlets = [[('chr1', False, 35, 55, 1, 0), 
                              ('chr1', False, 45, 96, 11, 0)],
                             [('chr1', False, 46, 90, 12, 0)]]
            final_alignments = selected_readlet_alignments_by_distance(
                                    multireadlets
                                )
            self.assertTrue(
                sorted(final_alignments) == [('chr1', False, 45, 96, 11, 0),
                                             ('chr1', False, 46, 90, 12, 0)]
            )

    class TestSelectedReadletAlignmentsByCoverage(unittest.TestCase):
        """ Tests selected_readlet_alignments_by_coverage(). """

        def test_that_best_covered_alignment_is_chosen(self):
            """ Fails if proper alignment of multireadlet is not chosen.
            """
            multireadlets = [[('chr1', False, 35, 55, 1, 0), 
                              ('chr1', False, 45, 96, 11, 0)],
                             [('chr1', False, 46, 90, 12, 0)]]
            final_alignments = selected_readlet_alignments_by_coverage(
                                    multireadlets
                                )
            self.assertTrue(
                sorted(final_alignments) == [('chr1', False, 45, 96, 11, 0),
                                             ('chr1', False, 46, 90, 12, 0)]
            )

        def test_that_no_alignment_is_chosen_in_case_of_tie(self):
            """ Fails if alignments of multireadlet aren't thrown out. """
            '''Below, the multireadlet at the first position in multireadlets
            contains alignments that overlap the unireadlet at the second
            position in multireadlets equally.'''
            multireadlets = [[('chr1', False, 25, 55, 1, 0), 
                              ('chr1', False, 55, 85, 31, 0)],
                             [('chr1', False, 52, 58, 31, 0)]]
            final_alignments = selected_readlet_alignments_by_coverage(
                                    multireadlets
                                )
            self.assertTrue(final_alignments == 
                            [('chr1', False, 52, 58, 31, 0)])

    unittest.main()