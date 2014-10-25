#!/usr/bin/env python
"""
Rail-RNA-cointron_fasta
Follows Rail-RNA-cointron_search
Precedes Rail-RNA-realign_reads

Reduce step in MapReduce pipelines that outputs FASTA line for a "reference"
obtained by concatenating exonic sequences framing each intron in an intron
configuration.

Input (read from stdin)
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

Input is partitioned by the first two fields.

Hadoop output (written to stdout)
----------------------------
Tab-delimited tuple columns, one for each read sequence:
1. Read sequence
2. '\x1c' + FASTA reference name including '>'. The following format is used:
    original RNAME + '+' or '-' indicating which strand is the sense strand
    + '\x1d' + start position of sequence + '\x1d' + comma-separated list of
    subsequence sizes framing introns + '\x1d' + comma-separated list of intron
    sizes + '\x1d' + 'i' to indicate base string overlaps introns
3. FASTA sequence
"""
import sys
import time
import os
import site
import argparse
import copy
import string

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
from dooplicity.tools import xstream

parser = argparse.ArgumentParser(description=__doc__, 
            formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument(\
    '--verbose', action='store_const', const=True, default=False,
    help='Print out extra debugging statements')
bowtie.add_args(parser)
args = parser.parse_args()

start_time = time.time()
input_line_count = 0
reference_index = bowtie_index.BowtieIndexReference(args.bowtie_idx)
reversed_complement_translation_table = string.maketrans('ATCG', 'TAGC')
for key, xpartition in xstream(sys.stdin, 2, skip_duplicates=True):
    # Store max extend sizes and sequences for identical intron combos
    combos = {}
    first_start_position = int(key[1])
    seq_count = 0
    for value in xpartition:
        assert len(value) == 5
        input_line_count += 1
        seq_count += 1
        try:
            start_positions = (first_start_position,) + \
                                tuple(map(int, value[0].split(',')))
        except ValueError:
            # No other start positions
            start_positions = tuple([first_start_position])
        end_positions = tuple(map(int, value[1].split(',')))
        if (start_positions, end_positions) not in combos:
            combos[(start_positions, end_positions)] = \
                [int(value[2]), int(value[3]), set([value[4]])]
        else:
            combos[(start_positions, end_positions)][0] = \
                max(int(value[2]), combos[(start_positions,
                                                end_positions)][0])
            combos[(start_positions, end_positions)][1] = \
                max(int(value[3]), combos[(start_positions,
                                                end_positions)][1])
            combos[(start_positions, end_positions)][2].add(value[4])
    rname = key[0]
    reverse_strand_string = rname[-1]
    rname = rname[:-1]
    # The code below greedily assigns subsumed combos to larger combos
    final_combos = {}
    for combo in sorted(combos.keys(), key=len, reverse=True):
        subsumed = False
        intron_count = len(combo[0])
        assert intron_count == len(combo[1])
        for compared_combo in final_combos.keys():
            if (compared_combo[0][:intron_count],
                compared_combo[1][:intron_count]) == combo \
                and (compared_combo[0][intron_count]
                        - compared_combo[1][intron_count-1]
                        >= combos[combo][1]):
                # Subsume by adding sequences
                for seq_to_add in combos[combo][2]:
                    final_combos[compared_combo][2].add(seq_to_add)
                subsumed = True
                break
        if not subsumed:
            # Add as final combo
            final_combos[combo] = copy.deepcopy(combos[combo])
    for combo in final_combos:
        intron_combo = zip(combo[0], combo[1])
        left_extend_size = final_combos[combo][0]
        right_extend_size = final_combos[combo][1]
        reference_length = reference_index.length[rname]
        subseqs = []
        left_start = max(intron_combo[0][0] - left_extend_size, 1)
        # Add sequence before first intron
        subseqs.append(
                reference_index.get_stretch(rname, left_start - 1, 
                    intron_combo[0][0] - left_start)
            )
        # Add sequences between introns
        for i in xrange(1, len(intron_combo)):
            subseqs.append(
                    reference_index.get_stretch(rname, 
                        intron_combo[i-1][1] - 1,
                        intron_combo[i][0]
                        - intron_combo[i-1][1]
                    )
                )
        # Add final sequence
        subseqs.append(
                reference_index.get_stretch(rname,
                    intron_combo[-1][1] - 1,
                    min(right_extend_size, reference_length - 
                                        intron_combo[-1][1] + 1))
            )
        '''A given reference name in the index will be in the following format:
        original RNAME + '+' or '-' indicating which strand is the sense strand
        + ';' + start position of sequence + ';' + comma-separated list of
        subsequence sizes framing introns + ';' + comma-separated list of
        intron sizes'''
        fasta_info = ('\x1c>' + rname + reverse_strand_string 
                     + '\x1d' + str(left_start) + '\x1d'
                     + ','.join([str(len(subseq)) for subseq in subseqs])
                     + '\x1d' + ','.join([str(intron_end_pos - intron_pos)
                                  for intron_pos, intron_end_pos
                                  in intron_combo])
                     + '\x1di\t' + ''.join(subseqs))
        for read_seq in final_combos[combo][2]:
            print read_seq + '\t' + fasta_info
            reversed_complement_read_seq = read_seq[::-1].translate(
                    reversed_complement_translation_table
                )
            print reversed_complement_read_seq + '\t' + fasta_info
    if args.verbose:
        print >>sys.stderr, '%d potential FASTA reference sequences ' \
                            'condensed to %d' % (seq_count, len(final_combos))
print >>sys.stderr, 'DONE with cointron_fasta.py; in=%d; ' \
                    'time=%0.3f s' % (input_line_count,
                                        time.time() - start_time)
