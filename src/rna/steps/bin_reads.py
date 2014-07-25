#!/usr/bin/env python
"""
bin_reads.py
Follows Rail-RNA-cointron_fasta / Rail-RNA-align_reads
Precedes Rail-RNA-realign_reads

Associates read sequences that likely overlap introns with genomic partitions.

Input (read from stdin)
----------------------------
Tab-delimited input tuple columns:
(two kinds)

Type 1:
1. Read sequence
2. '\x1c' + FASTA reference name including '>'. The following format is used:
    original RNAME + '+' or '-' indicating which strand is the sense strand
    + '\x1d' + start position of sequence + '\x1d' + comma-separated list of
    subsequence sizes framing introns + '\x1d' + comma-separated list of intron
    sizes) + '\x1d' + 'p' if derived from primary alignment to genome; 's' if
    derived from secondary alignment to genome; 'i' if derived from cointron
    search
3. FASTA sequence + '\x1d' + base 36-encoded integer A such that A & sample
    index != 0 iff sample contains intron combo purportedly overlapped by read
    sequence
4. A random partition in which an exon belongs if a primary alignment; else
   '\x1c'

Type 2:
1. Read sequence
2. QNAME
3. Quality sequence
4. '\x1c'

Type 1 corresponds to a FASTA line to index to which the read sequence is
predicted to align. Type 2 corresponds to a distinct read. Input is partitioned
by field 1.

Hadoop output (written to stdout)
----------------------------
Type 1:
1. Calculated (random) guess as to which partition some exon from read belongs
2. '\x1c'
3. '\x1c' + FASTA reference name including '>'. The following format is used:
    original RNAME + '+' or '-' indicating which strand is the sense strand
    + '\x1d' + start position of sequence + '\x1d' + comma-separated list of
    subsequence sizes framing introns + '\x1d' + comma-separated list of intron
    sizes) + '\x1d' + 'p' if derived from primary alignment to genome; 's' if
    derived from secondary alignment to genome; 'i' if derived from cointron
    search
4. FASTA sequence

Type 2:
1. Calculated (random) guess as to which partition some exon from read belongs
2. Read sequence
3. QNAME
4. Quality sequence
"""

import random
import sys
import time
import os
import site

base_path = os.path.abspath(
                    os.path.dirname(os.path.dirname(os.path.dirname(
                        os.path.realpath(__file__)))
                    )
                )
utils_path = os.path.join(base_path, 'rna', 'utils')
site.addsitedir(utils_path)
site.addsitedir(base_path)

from dooplicity.tools import xstream, dlist

start_time = time.time()
input_line_count, output_line_count = 0, 0

for (seq,), xpartition in xstream(sys.stdin, 1):
    reads = dlist()
    fastas = dlist()
    partitions = set()
    for fasta_or_read in xpartition:
        assert len(fasta_or_read) == 3
        input_line_count += 1
        if fasta_or_read[0][0] == '\x1c':
            # FASTA line
            if fasta_or_read[2] != '\x1c':
                partitions.add(fasta_or_read[-1])
            fastas.append('\t'.join(fasta_or_read[:-1]))
        else:
            # Read line
            reads.append('\t'.join(fasta_or_read[:-1]))
    '''There had better have been at least one soft-clipped alignment of this
    read since all unmapped reads were thrown out in align_reads.'''
    partition_count = len(partitions)
    if partition_count == 0:
        '''This seq never had a soft-clipped alignment; it's a reverse-
        complemented read from cointron_fasta that never had an actual read
        associated with it.'''
        continue
    elif partition_count == 1:
        partition_id = list(partitions)[0]
    else:
        # partition_count > 1: choose one at random
        random.seed(seq)
        partition_id = random.choice(list(partitions))
    for read in reads:
        print partition_id + '\t' + seq + '\t' + read
        output_line_count += 1
    for fasta in fastas:
        print partition_id + '\t\x1c\t' + fasta
        output_line_count += 1

print >>sys.stderr, 'DONE with bin_reads.py; in/out=%d/%d; ' \
        'time=%0.3f s' % (input_line_count, output_line_count,
                                time.time() - start_time)