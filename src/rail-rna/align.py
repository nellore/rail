#!/usr/bin/env python
"""
Rail-RNA-align

Follows Rail-RNA-preprocess
Precedes Rail-RNA-intron / Rail-RNA-normalize_pre

Alignment script for MapReduce pipelines that wraps Bowtie. Obtains exonic 
chunks and introns by dividing each RNA-seq read into several overlapping 
segments, or readlets, aligning the readlets, and -- very roughly
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

Output (written to stdout)
----------------------------
A given RNAME sequence is partitioned into intervals ("bins") of some 
user-specified length (see partition.py).

Exonic chunks (aka ECs; two formats --- either or both may be emitted):

Format 1 (exon_ival); tab-delimited output tuple columns:
1. Reference name (RNAME in SAM format) + ';' + bin number +  
    '+' or '-' indicating which strand is the sense strand
2. Sample label
3. EC start (inclusive) on forward strand
4. EC end (exclusive) on forward strand

Format 2 (exon_differential); tab-delimited output tuple columns:
1. Reference name (RNAME in SAM format) + ';' + bin number +  
    '+' or '-' indicating which strand is the sense strand
2. Sample label
3. max(EC start, bin start) (inclusive) on forward strand IFF next column is +1 
   and EC end (exclusive) on forward strand IFF next column is -1.
4. +1 or -1.

Introns

Tab-delimited output tuple columns:
1. Reference name (RNAME in SAM format) + ';' + bin number +  
    '+' or '-' indicating which strand is the sense strand
2. Sample label
3. Intron start (inclusive) on forward strand
4. Intron end (exclusive) on forward strand
5. Number of nucleotides between 5' end of intron and 5' end of read from which
it was inferred, ASSUMING THE SENSE STRAND IS THE FORWARD STRAND. That is, if
the sense strand is the reverse strand, this is the distance between the 5' end
of the reverse-complemented read and the 5' end of the reverse-complemented 
intron.
6. Number of nucleotides between 3' end of intron and 3' end of read from which
it was inferred, ASSUMING THE SENSE STRAND IS THE FORWARD STRAND.

ALL OUTPUT COORDINATES ARE 1-INDEXED.
"""
import sys
import os
import site
import threading
import string
import tempfile
import argparse
import numpy as np

base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
site.addsitedir(os.path.join(base_path, "bowtie"))
site.addsitedir(os.path.join(base_path, "sample"))
site.addsitedir(os.path.join(base_path, "alignment"))
site.addsitedir(os.path.join(base_path, "fasta"))
site.addsitedir(os.path.join(base_path, "interval"))

import bowtie
import sample
import needlemanWunsch
import fasta
import partition

# Print file's docstring if -h is invoked
parser = argparse.ArgumentParser(description=__doc__, 
            formatter_class=argparse.RawDescriptionHelpFormatter)
# Add command-line arguments
parser.add_argument('--refseq', type=str, required=False, 
    help='The fasta sequence of the reference genome. The fasta index of the '
         'reference genome is also required to be built via samtools.')
# To be implemented; for now, index is always fasta filename + .fai
parser.add_argument('--faidx', type=str, required=False, 
    help='Fasta index file')
parser.add_argument('--max-intron-length', type=int, required=False,
    default=1000000, 
    help='Filters out all potential introns longer than this length')
parser.add_argument('--min-readlet-size', type=int, required=False, default=5, 
    help='Capping readlets (that is, readlets that terminate '
         'at a given end of the read) are never smaller than this value')
parser.add_argument('--max-readlet-size', type=int, required=False, default=25, 
    help='Size of every noncapping readlet')
parser.add_argument('--readlet-interval', type=int, required=False, default=5, 
    help='Number of bases separating successive noncapping readlets along '
         'the read')
parser.add_argument('--capping-fraction', type=float, required=False,
    default=0.75, 
    help='Successive capping readlets on a given end of a read are tapered '
         'in size exponentially with this fractional base')
parser.add_argument('--report_multiplier', type=float, required=False,
    default=1.2,
    help='When --verbose is also invoked, the only lines of lengthy '
         'intermediate output written to stderr have line number that '
         'increases exponentially with this base.')
parser.add_argument('--verbose', action='store_const', const=True,
    default=False,
    help='Print out extra debugging statements')
parser.add_argument('--exon-differentials', action='store_const', const=True,
    default=False, 
    help='Print exon differentials (+1s and -1s)')
parser.add_argument('--exon-intervals', action='store_const', const=True,
    default=False, 
    help='Print exon intervals')
parser.add_argument('--test', action='store_const', const=True, default=False,
    help='Run unit tests; DOES NOT NEED INPUT FROM STDIN, AND DOES NOT '
         'OUTPUT EXONS AND INTRONS TO STDOUT')
parser.add_argument('--intron-partition-overlap', type=int, required=False,
    default=20, 
    help='Amount by which partitions overlap their left and right neighbors')
parser.add_argument('--min-intron-size', type=int, required=False,
    default=5,
    help='Filters introns of length smaller than this value')
parser.add_argument('--max-intron-size', type=int, required=False,
    default=100000, 
    help='Filters introns of length greater than this value')
parser.add_argument('--min-strand-readlets', type=int, required=False,
    default=1, 
    help='Exons and introns are called on a strand (i.e., a tuple ' 
    '(rname, reverse_strand)) iff the number of readlets that align to the '
    ' strand is >= this number')
parser.add_argument('--max-discrepancy', type=int, required=False,
    default=2, 
    help='If the difference in length between an unmapped region framed by '
         'two ECs and its corresponding gap in the reference is <= this '
         'value, the unmapped region is considered a candidate for '
         'incorporation into a single EC spanning the two original ECs via '
         'DP filling')
parser.add_argument('--min-seq-similarity', type=float, required=False,
    default=0.5, 
    help='If the difference in length between an unmapped region framed by'
         'two ECs and its corresponding reference is <= the command-line '
         'option --max-discrepancy AND the score of Needleman-Wunsch '
         'alignment is >= --min-seq-similarity * (length of unmapped region). '
         'the unmapped region is incorporated into a single EC spanning the '
         'two original ECs via DP filling')
parser.add_argument('--archive', metavar="PATH", type=str, 
    default=None,
    help='Save output and Bowtie command to a subdirectory (named using this ' 
         'process\'s PID) of PATH')

# Add command-line arguments for dependency
partition.addArgs(parser)
bowtie.addArgs(parser)

# Collect the Bowtie arguments
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
global, properties of args are also arguments of the go() function so different
command-line arguments can be passed to it for unit tests.'''
args = parser.parse_args(argv[1:])

# Initialize global variables for tracking number of input/output lines
input_line_count = 0
output_line_count = 0

def write_reads(output_stream, input_stream=sys.stdin, readletize=False, 
         min_readlet_size=5, max_readlet_size=25, readlet_interval=5, 
         capping_fraction=.75, verbose=False, report_multiplier=1.2):
    """Writes input reads/readlets in tab-separated format parsable by Bowtie.

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
                         at a given end of the read) are never smaller than 
                         this value. Ignored if readletize=False.
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
                          first readlet written to stderr increases 
                          exponentially with base report_multiplier.

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
        for i, line in enumerate(input_stream):
            tokens = line.rstrip().split('\t')
            '''len(tokens) must = 3; write_reads() is called only by 
            BowtieOutputThread instance after a first pass of Bowtie on 
            full reads has failed. Further, name should already have ;0, ;1, 
            or ;2 tacked on depending on a read's single- or paired-end
            status.'''
            assert len(tokens) == 3
            '''A readlet name is in the following (slightly redundant) format:
                original name;0/1/2 (single- or paired-end label);displace-
                ment of readlet's 5' end from read's 5' end;displacement of
                readlet's 3' end from read's 3' end;read sequence;number of 
                readlets in read.'''
            to_write = []
            seq_size = len(tokens[1])
            # Add capping readlets
            to_write += [['%s;%d;%d;%s' % (tokens[0], 0, seq_size - cap_size,
                                            tokens[1]),
                            '\t%s\t%s' % (tokens[1][:cap_size],
                                            tokens[2][:cap_size])]
                            for cap_size in cap_sizes]
            to_write += [['%s;%d;%d;%s' % (tokens[0], seq_size - cap_size, 0,
                                            tokens[1]),
                            '\t%s\t%s' % (tokens[1][-cap_size:],
                                            tokens[2][-cap_size:])]
                            for cap_size in cap_sizes]
            # Add noncapping readlets
            to_write += [[
                            '%s;%d;%d;%s' % (tokens[0], j, 
                                            seq_size - j - max_readlet_size,
                                            tokens[1]),
                            '\t%s\t%s' % (tokens[1][j:j + readlet_interval],
                                            tokens[2][j:j + readlet_interval])
                            ] for j in range(readlet_interval, seq_size, 
                                                seq_size - readlet_interval)]
            # Add total number of readlets to each string
            readlet_count = len(to_write)
            readlet_count_string = ';' + str(readlet_count)
            for j in range(readlet_count):
                to_write[j][0] += readlet_count_string
                to_write[j] = ''.join(to_write[j])
            if verbose and next_report_line == i:
                print >>sys.stderr, 'First readlet from read %d: %s' \
                    % (i + 1, to_write[0])
                next_report_line = int((next_report_line + 1)
                    * report_multiplier + 1) - 1
            for readlet in to_write:
                print >>output_stream, readlet
    else:
        global input_line_count
        # Don't readletize
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
                    % (tokens[0] + ';0', tokens[1], tokens[2])
                print >>output_stream, to_write
            else:
                to_write = '%s\t%s\t%s\n%s\t%s\t%s' \
                    % (tokens[0] + ';1', tokens[1], tokens[2],
                        (tokens[0] if len(tokens) == 5 else tokens[3]) + ';2',
                        tokens[-2], tokens[-1])
                print >>output_stream, to_write
            if verbose and next_report_line == i:
                print >>sys.stderr, 'Read(s) %d: %s' % (i, to_write)
                next_report_line = int((next_report_line + 1)
                    * report_multiplier + 1) - 1
            input_line_count += 1
    output_stream.flush()

def composed_and_sorted_readlets(readlets, min_strand_readlets=1):
    """Composes overlapping readlets and sorts them by reference position.

       This function assumes that for overlapping readlets mapping to, for
       example, the forward strand and covering the reference like so:

             Reference  5'---------------------3'
             Readlet 1       =========
             Readlet 2           ============
             Readlet 3      ===============    ,

       the displacement of the readlet whose pos is smallest (here, Readlet 3),
       is the displacement of the interval overlapped --- EVEN IF THE OTHER 
       READLETS' DISPLACEMENTS SEEM INCONSISTENT WITH THAT DISPLACEMENT given
       readlet pos's AND end_pos's. See below for definitions of displacement,
       pos, and end_pos.

       readlets: list of tuples (rname, reverse_strand, pos, end_pos,
                 displacement), each of which corresponds to a readlet. rname 
                 is the SAM-format RNAME, typically the chromosome to which the 
                 readlet maps. reverse_strand is True iff the readlet's 
                 reversed complement aligns to the reference. The readlet 
                 should span the reference interval [pos, end_pos).
                 displacement is the number of bases between the 5' (3') end of
                 the readlet, which aligns to the forward (reverse) strand, and
                 the 5' (3') end of the read.
       min_strand_readlets: readlets are composed on a strand iff the number of
                            readlets that align to the strand is >= this
                            number.

       Return value: dictionary with each key a tuple (rname, reverse_strand)
                     uniquely identifying a strand and each value a list of
                     tuples (pos, end_pos, displacement) sorted by pos in 
                     ascending order; here, each tuple from the list denotes
                     an exonic interval spanning [pos, end_pos) and
                     displacement is the number of bases between the 5' (3')
                     end of the interval, composed of readlets aligning to the
                     forward (reverse) strand, and the 5' (3') end of the
                     read.
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
        uncomposed[strand].sort(key=lambda readlet: readlet[0])
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
        right_reference_seq, substitution_matrix=needlemanWunsch.matchCost()):
    """Distributes a read's unmapped region between two framing mapped regions.

       The algorithm used is best illustrated with an example:
       Read         5' mmmmmmmmmmmmmmmmmmmm??????|????MMMMMMMMMMMMMMMMMMMMM 3'

       Reference    5' ...mmmmmmmmmLLLLLL|LLLLBB...BBRRRRRR|RRRRMMMMMMMM... 3'

       Besides the pipe bars, each character above represents a base. Call the 
       unmapped region of the read denoted by the string of question marks
       above unmapped_seq. Read M's were previously mapped to reference M's by
       Bowtie. Similarly, read m's were mapped to reference m's. THE ALGORITHM
       ASSUMES THAT THE M's and m's ARE FORWARD-STRAND ALIGNMENTS. B's denote
       extra (presumably intronic) bases of the reference not considered by the
       algorithm. Let K be the number of bases of the read spanned by the
       unmapped region (i.e., the number of question marks). Call the K bases
       of the reference following the lefthand mapped region (represented by
       L's above) left_reference_seq. Similarly, call the K bases of the
       reference preceding the righthand mapped region (represented by R's
       above) right_reference_seq. Note that the lengths of unmapped_seq,
       left_reference_seq, and right_reference_seq are the same. Also note that
       in general, left_reference_seq and right_reference_seq may overlap.

       The algo first fills two (K + 1) x (K + 1) score matrices according to
       rules encoded in substitution_matrix: left_score_matrix scores alignment
       of unmapped_seq to left_reference_seq, while right_score_matrix scores
       alignment of REVERSED unmapped_seq to REVERSED right_reference_seq.
       Let all matrices be 0-indexed, and consider only left_score_matrix. 
       Recall that the lower righthand corner element of the matrix is the
       score of the best global alignment. More generally, the (S_L, S_L)
       element of left_score_matrix is the score of the best global alignment
       of the FIRST S_L bases of unmapped_seq to the FIRST S_L bases of
       left_reference_seq. Similar logic would apply to right_score_matrix, but
       first flip right_score_matrix upside-down and leftside-right. Then the
       (K - S_R, K - S_R) element of right_score_matrix is the score of the
       best global alignment of the FINAL S_R bases of unmapped_seq to the
       FINAL S_R bases of right_reference_seq.

       The algo picks out the best "jump" (represented by the pipe bars) from
       left_reference_seq to right_reference_seq for alignment of unmapped_seq;
       that is, it finds the value of S_L that maximizes the sum of the scores
       of the best global alignments of the first S_L bases of unmapped_seq to
       the first S_L bases of left_reference_seq and the final (K - S_L) bases
       of unmapped_seq to the final (K - S_L) bases of right_reference_seq.
       The maximum diagonal element of the matrix
       (left_score_matrix + right_score_matrix) is indexed by (S_L, S_L).
       The maximum may be degenerate --- that is, there may be a tie among more
       than one S_L.

       unmapped_seq: unmapped sequence (string; see paragraph above).
       left_reference_seq: left reference sequence (string; see paragraphs 
                           above).
       right_reference_seq: right reference sequence (string; see paragraphs
                            above). unmapped_seq, left_reference_seq, and
                            right_reference_seq all have the same length.
       substitution_matrix: 6 x 6 substitution matrix (numpy object or list of
                            lists; see paragraphs above); rows and columns
                            correspond to ACGTN-, where N is aNy and - is
                            a gap. Default: +1 for match, -2 for gap, -1 for 
                            everything else.

       Return value: List of possible S_L (see paragraphs above for
                     definition).
    """
    k = len(unmapped_seq)
    assert k == len(left_reference_seq)
    assert k == len(right_reference_seq)

    left_score_matrix = needlemanWunsch.needlemanWunsch(
            unmapped_seq, left_reference_seq, 
            substitution_matrix
        )[1]
    right_score_matrix = needlemanWunsch.needlemanWunsch(
            unmapped_seq[::-1], right_reference_seq[::-1], 
            substitution_matrix
        )[1][::-1,::-1]
    total_score_matrix_diagonal = np.diagonal(left_score_matrix) \
                                  + np.diagonal(right_score_matrix)
    return np.argwhere(np.amax(total_score_matrix_diagonal) 
                                == total_score_matrix_diagonal).flatten()

def exons_and_introns_from_read(fasta_object, read_seq, readlets, 
    min_intron_size=5, min_strand_readlets=1, max_discrepancy=2, 
    min_seq_similarity=0.5, substitution_matrix=needlemanWunsch.matchCost()):
    """Composes a given read's aligned readlets and returns ECs and introns.

       fasta_object: object of class fasta.fasta corresponding to 
                  FASTA reference; used for realignment of unmapped
                  regions between exonic chunks.
       read_seq: sequence of the original read from which the readlets are
                 derived.
       readlets: list of tuples (rname, reverse_strand, pos, end_pos,
                 displacement), each of which corresponds to a readlet. rname 
                 is the SAM-format RNAME, typically the chromosome to which the 
                 readlet maps. reverse_strand is True iff the readlet's 
                 reversed complement aligns to the reference. The readlet 
                 should span the interval [pos, end_pos). displacement is the 
                 number of bases between the 5' (3') end of the readlet, which
                 aligns to the forward (reverse) strand, and the 5' (3') end of
                 the read.
       min_intron_size: introns smaller than this number of bases are 
                        filtered out.
       min_strand_readlets: exons and introns are called on a strand (i.e., a 
                            tuple (rname, reverse_strand)) iff the number of
                            readlets that align to the strand is >= this
                            number.
       max_discrepancy: if the difference in length between an unmapped 
                        region framed by two ECs and its corresponding gap in
                        the reference is <= this value, the unmapped region is 
                        considered a candidate for incorporation into a single
                        EC spanning the two original ECs via DP filling.
       min_seq_similarity: if the difference in length between an unmapped
                           region framed by two ECs and its corresponding 
                           reference is <= max_discrepancy AND the score of 
                           Needleman-Wunsch alignment is >= 
                           min_seq_similarity * (length of unmapped region),
                           the unmapped region is incorporated into a single EC
                           spanning the two original ECs via DP filling. See 
                           needlemanWunsch.matchCost() for the substitution 
                           matrix used.
       substitution_matrix: 6 x 6 substitution matrix (numpy object or list
                            of lists) for scoring Needleman-Wunsch filling
                            and framing alignments; rows and columns correspond
                            to ACGTN-, where N is aNy and - is a gap. Default:
                            +1 for match, -2 for gap, -1 for everything else.

       Return value: tuple (exons, introns).
                       -exons is a list of tuples (rname, reverse_strand, 
                        pos, end_pos), each of
                       which denotes an exonic chunk (EC). rname is the 
                       SAM-format RNAME---typically a chromosome.
                       reverse_strand is True iff the EC's sense strand is the
                       reverse strand. The EC spans the interval
                       [pos, end_pos).
                       -introns is a list of tuples (rname, reverse_strand,
                         pos, end_pos, five_prime_displacement, 
                         three_prime_displacement), each of which denotes
                       an intron. rname contains the SAM-format RNAME --- 
                       typically a chromosome. reverse_strand is True iff the 
                       intron's sense strand is the reverse strand. The intron
                       spans the interval [pos, end_pos). Assume the sense 
                       strand is the forward strand; that is, if the sense 
                       strand is the reverse strand, consider the reversed
                       complements of both the readlet and its parent read.
                       five_prime_displacement is the displacement of the 5'
                       end of the intron from the 5' end of the read, while
                       three_prime_displacement is the displacement of the 3'
                       end of the intron from the 3' end of the read.
    """
    composed = composed_and_sorted_readlets(readlets, min_strand_readlets)
    read_seq = read_seq.upper()
    reversed_complement_read_seq = read_seq[::-1].translate(
            string.maketrans('ATCG', 'TAGC')
        )
    introns, exons = [], []
    for strand in composed:
        rname, reverse_strand = strand
        if reverse_strand:
            '''Handle reverse strand reads the same way forward strands are
            handled.'''
            current_read_seq = reversed_complement_read_seq
        else:
            current_read_seq = read_seq
        read_seq_size = len(read_seq)
        last_pos, last_end_pos, last_displacement = composed[strand][0]
        if last_displacement:
            '''If last_displacement isn't 0, check if region to the left should
            be filled.'''
            if needlemanWunsch.needlemanWunsch(
                        current_read_seq[0:last_displacement],
                        fasta_object.fetch_sequence(rname, 
                            last_pos - last_displacement, last_pos - 1),
                        substitution_matrix
                    )[0] >= min_seq_similarity * last_displacement:
                # Fill
                last_pos -= last_displacement
                last_displacement = 0
        for pos, end_pos, displacement in composed[strand][1:]:
            if last_displacement >= displacement:
                '''Example cases handled:
                The 5' end of EC #1 isn't to the left of the 5' end of EC #2
                along the read:
                                             EC #1
                Read      |--------------===================
                                        ===================-----|
                                              EC #2
                Read      |-----=============------=============|
                                    EC #2              EC #1
                This situation is pathological. Throw out the entire read by
                returning empty exon and intron lists.'''
                return [], []
            reference_distance = pos - last_end_pos
            unmapped_displacement = last_displacement + last_end_pos \
                - last_pos
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
                    If there are no bases between the ECs in the read
                    sequence and almost no bases between the ECs in 
                    the reference sequence, merge the ECs: last_displacement 
                    and last_pos should remain the same on next iteration, 
                    but last_end_pos must change. Note that this case
                    subsumes the outside possibility that EC #1 overlaps
                    EC #2 by as much the constraint that EC #2 doesn't begin 
                    before EC #1 on the reference allows. This is likely an 
                    insertion in the read with respect to the reference, and
                    the ECs are still merged.'''
                    last_end_pos = end_pos
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
                    The displacement of EC #2 is to the right of the
                    displacement of EC #1, and there is an unmapped read region
                    between them. Extend the unmapped region by one nucleotide
                    on either side; splice junctions will have to be
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
                    The displacement of EC #2 is to the left of the 
                    displacement of EC #1; at least a pair of readlets, one 
                    from each EC, overlaps on the read. Call the overlap + 
                    one nucleotide on either side the unmapped region; 
                    splice junctions will have to be determined more precisely.
                    See past if-else.'''
                    last_end_pos -= -read_distance + 1
                    pos += -read_distance + 1
                    unmapped_displacement -= -read_distance + 1
                    displacement += -read_distance + 1
                    read_distance += -read_distance + 2
                # Now decide whether to DP fill or DP frame
                if abs(discrepancy) <= max_discrepancy:
                    # DP Fill
                    if needlemanWunsch.needlemanWunsch(
                            current_read_seq[
                                unmapped_displacement:displacement
                            ],
                            fasta_object.fetch_sequence(rname, last_end_pos, 
                                pos - 1), 
                            substitution_matrix
                        )[0] >= min_seq_similarity * read_distance:
                        '''If the unmapped region should be filled,
                        last_displacement and last_pos should remain the same
                        on next iteration, but last_end_pos must change; two
                        ECs are merged.'''
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
                            fasta_object.fetch_sequence(rname,
                                last_end_pos,
                                last_end_pos + read_distance - 1), 
                            fasta_object.fetch_sequence(rname,
                                pos - read_distance, pos - 1)
                        )
                    '''Decide which split to use: choose the smallest
                    medoid index.'''
                    nonintegral_median = np.median(splits)
                    split = splits[np.argmin(abs(splits - nonintegral_median))]
                    pos += split - read_distance
                    displacement += split - read_distance
                    last_end_pos += split
                    unmapped_displacement += split
                    call_exon, call_intron = True, True
                else: 
                    '''discrepancy < -max_discrepancy; again is likely large
                    insertion with respect to reference, but not too large:
                    EC #2 never begins before EC #1 on the reference.
                    Again, merge ECs.

                    Example case handled:
                                      EC #1       UMR          EC #2               
                    Read          |==========--------------==============|

                    Reference        ...==========---==============...
                                           EC #1  gap     EC #2
                    UMR = unmapped region.
                    '''
                    last_end_pos = end_pos
            if call_intron:
                # Call the reference region between the two ECs an intron.
                introns.append((rname, reverse_strand, last_end_pos, pos, 
                    displacement, len(read_seq) - displacement))
            if call_exon:
                '''Push ONLY the first EC to the exon list. (The next EC
                may get merged into another EC on another iteration.)'''
                exons.append((rname, reverse_strand, last_pos, last_end_pos))
                (last_pos, last_end_pos, last_displacement) = \
                    (pos, end_pos, displacement)
        '''Final exonic chunk must still be called since it won't get merged 
        with any other chunks. Check to see if region to right of chunk needs 
        filling first.'''
        unmapped_displacement = last_displacement + last_end_pos - last_pos
        unmapped_base_count = read_seq_size - unmapped_displacement
        if unmapped_displacement != read_seq_size:
            if needlemanWunsch.needlemanWunsch(
                                current_read_seq[unmapped_displacement:],
                                fasta_object.fetch_sequence(rname,
                                    last_end_pos, 
                                    last_end_pos + unmapped_base_count - 1), 
                                substitution_matrix
                            )[0] >= min_seq_similarity * unmapped_base_count:
                # Fill
                last_end_pos += unmapped_base_count
        exons.append((rname, reverse_strand, last_pos, last_end_pos))
    return exons, introns

class BowtieOutputThread(threading.Thread):
    """Processes Bowtie alignments, emitting tuples for exons and introns."""
    
    def __init__(self, input_stream, fasta_object, readletized=False, 
        write_unmapped=True, unmapped_stream=open(os.devnull, 'w'),
        output_stream=sys.stdout, exon_differentials=True,
        exon_intervals=False, verbose=False, bin_size=10000, min_intron_size=5,
        min_strand_readlets=1, max_discrepancy=2, min_seq_similarity=0.5,
        max_intron_size=100000, intron_partition_overlap=20,
        substitution_matrix=needlemanWunsch.matchCost(), 
        report_multiplier=1.2):
        """Constructor for BowtieOutputThread.

           input_stream: where to retrieve Bowtie's SAM output, typically
                         a Bowtie process's stdout.
           fasta_object: object of class fasta.fasta corresponding to 
                         FASTA reference; used for realignment of unmapped
                         regions between exonic chunks.
           readletized: True if input SAM contains readlets; if False, input
                        SAM contains whole reads.
           write_unmapped: True iff unmapped reads are to be written to
                           unmapped_stream. Ignored if readletized=True.
           unmapped_stream: where to write reads Bowtie can't map, possibly for 
                            later readletizing; typically, this is a file 
                            stream. Ignored if write_unmapped=False or
                            readletized=True.
           output_stream: where to emit exon and intron tuples; typically,
                          this is sys.stdout.
           exon_differentials: True iff EC differentials are to be emitted.
           exon_intervals: True iff EC intervals are to be emitted.
           verbose: True if alignments should occasionally be written 
                    to stderr.
           bin_size: genome is partitioned in units of bin_size for later
                     load balancing.
           min_intron_size: introns smaller than this number of bases are 
                            filtered out.
           min_strand_readlets: exons and introns are called on a strand iff 
                                the number of readlets that align to the strand
                                is >= this number.
           max_discrepancy: if the difference in length between an unmapped 
                            region framed by two ECs and its corresponding gap
                            in the reference is <= this value, the unmapped
                            region is considered a candidate for incorporation
                            into a single EC spanning the two original ECs via
                            DP filling.
           min_seq_similarity: if the difference in length between an unmapped
                               region framed by two ECs and its corresponding 
                               reference is <= max_discrepancy AND the score of 
                               Needleman-Wunsch alignment is >= 
                               min_seq_similarity *
                               (length of unmapped region), the unmapped region
                               is incorporated into a single EC spanning the
                               two original ECs via DP filling. See 
                               needlemanWunsch.matchCost() for the substitution 
                               matrix used.
           max_intron_size: an intron of that spans more than this number 
                            of bases is suppressed from output.
           intron_partition_overlap: number of bases to subtract from reference
                                     start position and add to reference end
                                     position of intron when computing genome
                                     partitions it is in.
           substitution_matrix: 6 x 6 substitution matrix (numpy object or list
                                of lists) for scoring Needleman-Wunsch filling
                                and framing alignments; rows and columns
                                correspond to ACGTN-, where N is aNy and - is
                                a gap. Default: +1 for match, -2 for gap, -1
                                for everything else.
           report_multiplier: if verbose is True, the line number of an
                              alignment written to stderr increases 
                              exponentially with base report_multiplier.
        """
        super(BowtieOutputThread, self).__init__()
        self.daemon = True
        self.input_stream = input_stream
        self.readletized = readletized
        self.write_unmapped = write_unmapped
        self.unmapped_stream = unmapped_stream
        self.output_stream = output_stream
        self.fasta_object = fasta_object
        self.verbose = verbose
        self.bin_size = bin_size
        self.min_strand_readlets = min_strand_readlets
        self.max_discrepancy = max_discrepancy
        self.min_seq_similarity = min_seq_similarity
        self.report_multiplier = report_multiplier
        self.exon_differentials = exon_differentials
        self.exon_intervals = exon_intervals
        self.max_intron_size = max_intron_size
        self.intron_partition_overlap = intron_partition_overlap
        self.substitution_matrix = substitution_matrix

    def run(self):
        """Prints exons for reads and exons/introns for readlets.

           Overrides default method containing thread activity.

           No return value.
        """
        global output_line_count
        next_report_line = 0
        if not self.readletized:
            for i, line in enumerate(self.input_stream):
                # Skip header line
                if line[0] == '@': continue
                # Variable abbreviations below mirror SAM spec
                (qname, flag, rname, pos, mapq, cigar, rnext,
                    pnext, tlen, seq, qual) = line.rstrip().split('\t')[:11]
                flag = int(flag)
                pos = int(pos)
                if self.verbose and next_report_line == i:
                    print >>sys.stderr, \
                        'SAM output record %d: rdname="%s", flag=%d' \
                        % (i, qname, flag)
                    next_report_line = int((next_report_line + 1)
                        * self.report_multiplier + 1) - 1
                if flag != 4:
                    '''Read is mapped; the full alignment is to be called as 
                    an exonic chunk (EC). Kill mate information in qname 
                    so rest of pipeline associates mates.'''
                    reverse_strand_string = '-' if (flag & 16) != 0 else '+'
                    end_pos = pos + len(seq)
                    sample_label = sample.parseLab(qname[:-2])
                    partitions = partition.partition(rname, pos, 
                        end_pos, self.bin_size)
                    if self.exon_differentials:
                        for (partition_id, 
                            partition_start, partition_end) in partitions:
                            # Print increment at interval start
                            assert pos < partition_end
                            print >>self.output_stream, \
                                'exon_diff\t%s%s\t%s\t%d\t1' \
                                % (partition_id, reverse_strand_string, 
                                    sample_label, max(partition_start, pos))
                            output_line_count += 1
                            assert end_pos > partition_start
                            if end_pos < partition_end:
                                '''Print decrement at interval end iff exon
                                ends before partition ends.'''
                                print >>self.output_stream, \
                                    'exon_diff\t%s%s\t%s\t%d\t-1' \
                                    % (partition_id, reverse_strand_string, 
                                        sample_label, end_pos)
                                output_line_count += 1
                    if self.exon_intervals:
                            for partition_id, _, _ in partitions:
                                print >>self.output_stream, \
                                    'exon_ival\t%s%s\t%d\t%d\t%s' \
                                    % (partition_id, reverse_strand_string,
                                        pos, end_pos, sample_label)
                                output_line_count += 1
                elif self.write_unmapped:
                    print >>self.unmapped_stream, '%s\t%s\t%s' % (qname, 
                        seq, qual)
        else:
            # Input is readletized, and readlets must be composed
            '''Create dictionary for gathering aligned readlets belonging
            to same read.'''
            collected_readlets = {}
            '''Create dictionary for counting readlets, aligned or not,
            belonging to same read.'''
            readlet_count = {}
            for i, line in enumerate(self.input_stream):
                # Skip header line
                if line[0] == '@': continue
                # Variable abbreviations below mirror SAM spec
                (qname, flag, rname, pos, mapq, cigar, rnext,
                    pnext, tlen, seq, qual) = line.rstrip().split('\t')[:11]
                flag = int(flag)
                pos = int(pos)
                if self.verbose and next_report_line == i:
                    print >>sys.stderr, \
                        'SAM output record %d: rdname="%s", flag=%d' \
                        % (i, qname, flag)
                    next_report_line = int((next_report_line + 1)
                        * self.report_multiplier + 1) - 1
                '''Recall that readlet name is in the following format:
                original name;0/1/2 (single- or paired-end label);
                displacement of readlet's 5' end from read's 5' end;
                displacement of readlet's 3' end from read's 3' end;
                read sequence;number of readlets in read.'''
                split_qname = qname.split(';')
                (paired_label, five_prime_displacement, 
                    three_prime_displacement, 
                    read_seq, total_readlets) = split_qname[-5:]
                total_readlets = int(total_readlets)
                qname = ';'.join(split_qname[:-5])
                if not readlet_count.has_key((qname, paired_label)):
                    readlet_count[(qname, paired_label)] = 1
                else:
                    readlet_count[(qname, paired_label)] += 1
                if flag != 4:
                    # Readlet is mapped
                    '''All output positions will be with respect to 5' end of
                    forward strand. Note that seq is reverse-complemented
                    readlet if sense strand is reverse strand, but
                    read_seq is NEVER reverse-complemented.'''
                    reverse_strand = (flag & 16) != 0
                    if reverse_strand:
                        '''seq is reverse-complemented; displacement of 
                        readlet's 3' end from read's 3' end is displacement of 
                        readlet's 5' end from read's 5' end on complementary
                        (forward) strand.'''
                        displacement = int(three_prime_displacement)
                    else:
                        displacement = int(five_prime_displacement)
                    if not collected_readlets.has_key((qname, paired_label)):
                        collected_readlets[(qname, paired_label)] = []
                    collected_readlets[(qname, paired_label)].append(
                                (rname, reverse_strand, pos, pos + len(seq), 
                                    displacement, read_seq)
                            )
                if readlet_count[(qname, paired_label)] == total_readlets \
                    and collected_readlets.has_key((qname, paired_label)):
                    exons, introns = exons_and_introns_from_read(
                        self.fasta_object, read_seq, [el[0:5] for el in \
                            collected_readlets[(qname, paired_label)]],
                        min_strand_readlets=self.min_strand_readlets,
                        max_discrepancy=self.max_discrepancy,
                        min_seq_similarity=self.min_seq_similarity,
                        substitution_matrix=self.substitution_matrix
                    )
                    '''Kill mate information in qname so rest of pipeline
                    associates mates.'''
                    sample_label = sample.parseLab(qname[:-2])
                    if self.exon_differentials:
                        for (exon_rname, exon_reverse_strand, exon_pos, 
                                exon_end_pos) in exons:
                            partitions = partition.partition(exon_rname,
                                exon_pos, exon_end_pos, self.bin_size)
                            exon_reverse_strand_string = '-' if \
                                exon_reverse_strand else '+'
                            for (partition_id, partition_start, 
                                    partition_end) in partitions:
                                assert exon_pos < partition_end
                                # Print increment at interval start
                                print >>self.output_stream, \
                                    'exon_diff\t%s%s\t%s\t%d\t1' \
                                    % (partition_id, 
                                        exon_reverse_strand_string,
                                        sample_label, 
                                        max(partition_start, exon_pos))
                                output_line_count += 1
                                assert exon_end_pos > partition_start
                                if exon_end_pos < partition_end:
                                    '''Print decrement at interval end iff 
                                    exon ends before partition ends.'''
                                    print >>self.output_stream, \
                                        'exon_diff\t%s%s\t%s\t%d\t-1' \
                                        % (partition_id, 
                                            exon_reverse_strand_string,
                                            sample_label, exon_end_pos)
                                    output_line_count += 1
                    if self.exon_intervals:
                        for (exon_rname, exon_reverse_strand, exon_pos, 
                                exon_end_pos) in exons:
                            partitions = partition.partition(exon_rname,
                                exon_pos, exon_end_pos, self.bin_size)
                            exon_reverse_strand_string = '-' if \
                                exon_reverse_strand else '+'
                            for partition_id, _, _ in partitions:
                                print >>self.output_stream, \
                                    'exon_ival\t%s%s\t%d\t%d\t%s' \
                                    % (partition_id, 
                                        exon_reverse_strand_string,
                                        exon_pos, exon_end_pos, 
                                        sample_label)
                                output_line_count += 1
                    # Print introns
                    for (intron_rname, intron_reverse_strand, intron_pos, 
                            intron_end_pos, intron_five_prime_displacement,
                            intron_three_prime_displacement) in introns:
                        intron_reverse_strand_string = '-' if \
                            intron_reverse_strand else '+'
                        if intron_end_pos - intron_pos > self.max_intron_size:
                            if self.verbose: 
                                print >> sys.stderr, 'Intron of size > ' \
                                    'max-intron-size = %d' \
                                    ' filtered at %s:%d-%d' \
                                    % (self.max_intron_size, intron_rname, 
                                        intron_pos, intron_end_pos)
                                output_line_count += 1
                                continue
                        partitions = partition.partition(intron_rname,
                            intron_pos, intron_end_pos, self.bin_size,
                            fudge=self.intron_partition_overlap)
                        for (partition_id, partition_start, 
                                partition_end) in partitions:
                            print >>self.output_stream, \
                                'intron\t%s%s\t%s\t%d\t%d\t%d\t%d' \
                                % (partition_id, intron_reverse_strand_string,
                                    sample_label, intron_pos, intron_end_pos, 
                                    intron_five_prime_displacement,
                                    intron_three_prime_displacement)
                            output_line_count += 1
                    del readlet_count[(qname, paired_label)]
                    del collected_readlets[(qname, paired_label)]

def go(reference_fasta, input_stream=sys.stdin, output_stream=sys.stdout,
    bowtie_exe="bowtie", bowtie_index_base="genome", bowtie_args=None,
    temp_dir_path=tempfile.mkdtemp(), archive=None, bin_size=10000, 
    verbose=False, read_filename='reads.tsv', readlet_filename='readlets.tsv',
    unmapped_filename='unmapped.tsv', exon_differentials=True,
    exon_intervals=False, min_readlet_size=5, max_readlet_size=25, 
    readlet_interval=5, capping_fraction=.75, min_strand_readlets=1, 
    max_discrepancy=2, min_seq_similarity=0.5, min_intron_size=5,
    max_intron_size=100000, intron_partition_overlap=20, 
    substitution_matrix=needlemanWunsch.matchCost(), report_multiplier=1.2):
    """Runs Rail-RNA-align.

       Two passes of Bowtie are run. The first attempts to align full reads.
       Mapped reads are emitted as exonic chunks. Unmapped reads are then
       segmented into readlets, which are aligned on second pass of Bowtie.
       Readlets belonging to the same read are used to infer introns and other
       exonic chunks, which are also emitted.

       Input (read from input_stream)---
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

       Output (written to output_stream)---
       A given RNAME sequence is partitioned into intervals ("bins") of some 
       user-specified length (see partition.py) for later load balancing.

       Exonic chunks (aka ECs; two formats --- either or both may be emitted):

       Format 1 (exon_ival); tab-delimited output tuple columns:
       1. Reference name (RNAME in SAM format) + ';' + bin number +  
       '+' or '-' indicating which strand is the sense strand
       2. Sample label
       3. EC start (inclusive) on forward strand
       4. EC end (exclusive) on forward strand

       Format 2 (exon_differential); tab-delimited output tuple columns:
       1. Reference name (RNAME in SAM format) + ';' + bin number +  
       '+' or '-' indicating which strand is the sense strand
       2. Sample label
       3. max(EC start, bin start) (inclusive) on forward strand IFF next 
          column is +1 and EC end (exclusive) on forward strand IFF next column
          is -1.
       4. +1 or -1.

       reference_fasta: filename, including path, of the reference fasta
       input_stream: where to find input reads
       output_stream: where to emit exonic chunks and introns
       bowtie_exe: filename of Bowtie executable; include path if not in
                   $PATH.
       bowtie_index_base: the basename of the Bowtie index files associated
                          with reference_fasta.
       bowtie_args: string containing precisely extra command-line arguments
                    to pass to Bowtie, e.g., "--tryhard --best"; or None.
       temp_dir_path: path of temporary directory for storing intermediate
                      alignments; archived if archive is not None.
       archive: directory name, including path, to which temp_dir_path should
                be renamed when script is complete; or None if temp_dir_path
                should be trashed.
       bin_size: genome is partitioned in units of bin_size for later
                 load balancing.
       verbose: True iff more informative messages should be written to stderr.
       read_filename: The file, excluding path, for storing reads before
                      alignment with Bowtie. Archived if archive is not None.
       readlet_filename: The file, excluding path, for storing readlets
                         before alignment with Bowtie.
       unmapped_filename: The file, excluding path, for storing reads unmapped
                          after Bowtie's first pass.
       min_readlet_size: "capping" readlets (that is, readlets that terminate
                         at a given end of the read) are never smaller than 
                         this value. Ignored if readletize=False.
       max_readlet_size: size of every noncapping readlet. Ignored if 
                         readletize=False.
       readlet_interval: number of bases separating successive readlets along
                         the read. Ignored if readletize=False.
       capping_fraction: successive capping readlets on a given end of a read
                         are tapered in size exponentially with fractional base
                         capping_fraction. Ignored if readletize=False.
       exon_differentials: True iff EC differentials are to be emitted.
       exon_intervals: True iff EC intervals are to be emitted.
       min_strand_readlets: exons and introns are called on a strand iff 
                            the number of readlets that align to the strand
                            is >= this number.
       max_discrepancy: if the difference in length between an unmapped 
                        region framed by two ECs and its corresponding gap
                        in the reference is <= this value, the unmapped
                        region is considered a candidate for incorporation
                        into a single EC spanning the two original ECs via
                        DP filling.
       min_seq_similarity: if the difference in length between an unmapped
                           region framed by two ECs and its corresponding 
                           reference is <= max_discrepancy AND the score of 
                           Needleman-Wunsch alignment is >= min_seq_similarity
                           * (length of unmapped region), the unmapped region
                           is incorporated into a single EC spanning the
                           two original ECs via DP filling. See 
                           needlemanWunsch.matchCost() for the substitution 
                           matrix used.
       min_intron_size: introns smaller than this number of bases are 
                        filtered out.
       max_intron_size: an intron of that spans more than this number of bases 
                        is suppressed from output.
       intron_partition_overlap: number of bases to subtract from reference
                                 start position and add to reference end
                                 position of intron when computing genome
                                 partitions it is in.
       substitution_matrix: 6 x 6 substitution matrix (numpy object or list
                            of lists) for scoring Needleman-Wunsch filling
                            and framing alignments; rows and columns correspond
                            to ACGTN-, where N is aNy and - is a gap. Default:
                            +1 for match, -2 for gap, -1 for everything else.
       report_multiplier: if verbose is True, the line number of an alignment,
                          read, or first readlet of a read written to
                          stderr increases exponentially with base
                          report_multiplier.

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
    fasta_object = fasta.fasta(reference_fasta)
    unmapped_filename = os.path.join(temp_dir_path, unmapped_filename)
    with open(unmapped_filename, 'w') as unmapped_stream:
        output_thread = BowtieOutputThread(
                bowtie_process.stdout, fasta_object,
                readletized=False, write_unmapped=True, 
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
                intron_partition_overlap=intron_partition_overlap,
                substitution_matrix=substitution_matrix,
                report_multiplier=report_multiplier
            )
        threads.append(output_thread)
        output_thread.start()
        # Join threads to pause execution in main thread
        for thread in threads:
            if verbose: print >>sys.stderr, 'Joining thread...'
            thread.join()
            if verbose: print >>sys.stderr, "    ...joined!"

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
            bowtie_process.stdout, fasta_object, 
            readletized=True, write_unmapped=False, 
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
            intron_partition_overlap=intron_partition_overlap,
            substitution_matrix=substitution_matrix,
            report_multiplier=report_multiplier
        )
    threads.append(output_thread)
    output_thread.start()
    # Join threads to pause execution in main thread
    for thread in threads:
        if verbose: print >>sys.stderr, 'Joining thread...'
        thread.join()
        if verbose: print >>sys.stderr, '    ...joined!'

    if verbose: print >>sys.stderr, 'Bowtie\'s second pass is finished.'
    output_stream.flush()

    if archive is not None:
        # Rename temporary directory for archiving
        os.rename(temp_dir_path, archive)
    else:
        # Kill temporary directory
        import shutil
        shutil.rmtree(temp_dir_path)

    print >> sys.stderr, 'DONE with align.py; in/out = %d/%d; ' \
        'time=%0.3f secs' % (input_line_count, output_line_count,
                                time.time()-start_time)

if not args.test:
    temp_dir_path=tempfile.mkdtemp()
    if args.verbose:
        print >>sys.stderr, 'Creating temporary directory %s' \
            % temp_dir_path
    go(args.refseq, bowtie_exe=args.bowtie_exe,
        bowtie_index_base=args.bowtie_idx,
        bowtie_args=bowtie_args, temp_dir_path=temp_dir_path,
        archive=os.path.join(
                args.archive, str(os.getpid())
            ) if args.archive is not None else None, 
        verbose=args.verbose, 
        bin_size=args.partition_length,
        read_filename='reads.tsv', 
        readlet_filename='readlets.tsv', unmapped_filename='unmapped.tsv', 
        exon_differentials=args.exon_differentials,
        exon_intervals=args.exon_intervals,
        min_readlet_size=args.min_readlet_size, 
        max_readlet_size=args.max_readlet_size,
        readlet_interval=args.readlet_interval,
        capping_fraction=args.capping_fraction, 
        min_strand_readlets=args.min_strand_readlets,
        max_discrepancy=args.max_discrepancy,
        min_seq_similarity=args.min_seq_similarity, 
        max_intron_size=args.max_intron_size,
        intron_partition_overlap=args.intron_partition_overlap,
        substitution_matrix=needlemanWunsch.matchCost(),
        report_multiplier=args.report_multiplier)
else:
    # Test units
    del sys.argv[1:] # Don't choke on extra command-line parameters
    import random
    import unittest

    def random_sequence(seq_size):
        """Gets random sequence of nucleotides.

           seq_size: number of bases to return.

           Return value: string of random nucleotides.
        """
        to_return = ""
        for el in xrange(seq_size):
            to_return += random.choice('ATCG')
        return to_return

    class TestComposedAndSortedReadlets(unittest.TestCase):
        """Tests composed_and_sorted_readlets(); needs no fixture."""
        def test_stray_readlet_filtering(self):
            """Fails if all readlets are not filtered."""
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
            """Fails if readlets are not consolidated as expected."""
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
            """Fails if readlets are not consolidated as expected."""
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
            """Fails if readlets are not consolidated as expected."""
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
            """Fails if readlets are not consolidated as expected."""
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
        """Tests unmapped_region_splits(); needs no fixture."""
        def test_random_exact_case_1000x(self):
            """Fails if proper split of unmapped region is not identified.

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
    
    unittest.main()