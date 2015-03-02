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
2. Comma-separated list of intron start positions in configuration
3. Comma-separated list of intron end positions in configuration
4. left_extend_size: by how many bases on the left side of an intron the
    reference should extend
5. right_extend_size: by how many bases on the right side of an intron the
    reference should extend
6. Read sequence

Input is partitioned by the first two fields.

Hadoop output (written to stdout)
----------------------------
Tab-delimited tuple columns, one for each read sequence:
1. Transcriptome Bowtie 2 index group number
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
from dooplicity.tools import xstream, dlist
import group_reads

parser = argparse.ArgumentParser(description=__doc__, 
            formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument(\
    '--verbose', action='store_const', const=True, default=False,
    help='Print out extra debugging statements')
bowtie.add_args(parser)
group_reads.add_args(parser)
args = parser.parse_args()

start_time = time.time()
input_line_count = 0
reference_index = bowtie_index.BowtieIndexReference(args.bowtie_idx)
group_reads_object = group_reads.IndexGroup(args.index_count)
for (rname, poses, end_poses), xpartition in xstream(sys.stdin, 3,
                                                        skip_duplicates=True):
    reverse_strand_string = rname[-1]
    rname = rname[:-1]
    read_seqs = dlist()
    poses = [int(pos) for pos in poses.split(',')]
    end_poses = [int(end_pos) for end_pos in end_poses.split(',')]
    max_left_extend_size, max_right_extend_size = None, None
    for left_extend_size, right_extend_size, read_seq in xpartition:
        input_line_count += 1
        max_left_extend_size = max(max_left_extend_size, int(left_extend_size))
        max_right_extend_size \
            = max(max_right_extend_size, int(right_extend_size))
        read_seqs.append(read_seq)
    intron_combo = zip(poses, end_poses)
    assert max_left_extend_size is not None
    assert max_right_extend_size is not None
    reference_length = reference_index.length[rname]
    subseqs = []
    left_start = max(intron_combo[0][0] - max_left_extend_size, 1)
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
                min(max_right_extend_size, reference_length - 
                                    intron_combo[-1][1] + 1))
        )
    '''A given reference name in the index will be in the following format:
    original RNAME + '+' or '-' indicating which strand is the sense strand
    + ';' + start position of sequence + ';' + comma-separated list of
    subsequence sizes framing introns + ';' + comma-separated list of intron
    sizes'''
    fasta_info = ('\x1c>' + rname + reverse_strand_string 
                 + '\x1d' + str(left_start) + '\x1d'
                 + ','.join([str(len(subseq)) for subseq in subseqs])
                 + '\x1d' + ','.join([str(intron_end_pos - intron_pos)
                              for intron_pos, intron_end_pos
                              in intron_combo])
                 + '\x1di\t' + ''.join(subseqs))
    if args.verbose:
        read_seq_count = 0
        for read_seq in read_seqs:
            read_seq_count += 1
            print '\t'.join([group_reads_object.index_group(read_seq),
                                read_seq, fasta_info])
        print >>sys.stderr, ('Printed %d read seqs for transcript fragment '
                             'on %s with introns %s.') % (read_seq_count,
                                                            rname,
                                                            str(intron_combo))
    else:
        for read_seq in read_seqs:
            print '\t'.join([group_reads_object.index_group(read_seq),
                                read_seq, fasta_info])

print >>sys.stderr, 'DONE with cointron_fasta.py; in=%d; ' \
                    'time=%0.3f s' % (input_line_count,
                                        time.time() - start_time)
