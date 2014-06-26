#!/usr/bin/env python
"""
Rail-RNA-intron_fasta

Follows Rail-RNA-intron_config
Precedes Rail-RNA-intron_index

Reduce step in MapReduce pipelines that outputs FASTA line for a "reference"
obtained by concatenating exonic sequences framing each intron in an intron
configuration. 

Input (read from stdin)
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
6. By how many bases on the left side of an intron the reference COULD extend,
    or NA if beginning of strand
7. By how many bases on the right side of an intron the reference COULD extend,
    or NA if end of strand
8. Sample index

Input is partitioned by the first three fields.

Hadoop output (written to stdout)
----------------------------
Tab-delimited tuple columns:
1. '-' to enforce that all lines end up in the same partition
2. FASTA reference name including '>'. The following format is used:
    original RNAME + '+' or '-' indicating which strand is the sense strand
    + '\x1d' + start position of sequence + '\x1d' + comma-separated list of
    subsequence sizes framing introns + '\x1d' + comma-separated list of intron
    sizes + '\x1d' + base-36-encoded integer A such that A & 2^sample index
    != 0 iff sample contains intron combo
3. Sequence
"""
import sys
import time
import os
import site
import argparse

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
from manifest import string_from_int

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
for key, xpartition in xstream(sys.stdin, 3, skip_duplicates=True):
    '''For computing maximum left and right extend sizes for every key --
    that is, every intron combo (fields 1-3 of input).'''
    left_extend_size, right_extend_size = None, None
    left_size, right_size = None, None
    sample_indexes = set()
    for value in xpartition:
        assert len(value) == 5
        input_line_count += 1
        left_extend_size = max(left_extend_size, int(value[-5]))
        right_extend_size = max(right_extend_size, int(value[-4]))
        try:
            left_size = max(left_size, int(value[-3]))
        except ValueError:
            left_size = 'NA'
        try:
            right_size = max(right_size, int(value[-2]))
        except ValueError:
            right_size = 'NA'
        sample_indexes.add(int(value[-1]))
    rname = key[0]
    reverse_strand_string = rname[-1]
    rname = rname[:-1]
    intron_combo = \
            zip([int(pos) for pos in key[1].split(',')],
                    [int(end_pos) for end_pos in key[2].split(',')])
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
    + '\x1d' + start position of sequence + '\x1d' + comma-separated list of
    subsequence sizes framing introns + '\x1d' + comma-separated list of
    intron sizes + '\x1d' + distance to previous intron or 'NA' if beginning of
    strand + '\x1d' + distance to next intron or 'NA' if end of strand
    + '\x1d' + base-36-encoded integer A such that A & 2^sample index != 0
    iff sample contains intron combo.'''
    encoded = ['0'] * (max(sample_indexes) + 1)
    for sample_index in sample_indexes:
        encoded[-(sample_index + 1)] = '1'
    encoded = string_from_int(int(''.join(encoded), base=2))
    print ('-\t>' + rname + reverse_strand_string 
            + '\x1d' + str(left_start) + '\x1d'
            + ','.join([str(len(subseq)) for subseq in subseqs]) + '\x1d'
            + ','.join([str(intron_end_pos - intron_pos)
                        for intron_pos, intron_end_pos
                        in intron_combo])
            + '\x1d' + str(left_size) + '\x1d' + str(right_size)
            + '\x1d' + encoded
            + '\t' + ''.join(subseqs)
        )

print >>sys.stderr, 'DONE with intron_fasta.py; in=%d; ' \
                    'time=%0.3f s' % (input_line_count,
                                        time.time() - start_time)