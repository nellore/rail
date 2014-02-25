"""
Rail-RNA-collapse
Follows Rail-RNA-align
Precedes Rail-RNA-coverage_pre

Reduce step in MapReduce pipelines that improves scalability of coverage
pipeline by collapsing exon_diffs for a sample occurring at the same genomic
position into one exon_diff. (An exon_diff denotes an increment or decrement in
coverage at a genome position.) Before this step, exon_diffs are partitioned
according to genome position but are not sorted. This way, positions in
particularly active partitions, where coverage is high, are (in not-unusual
cases, and at least theoretically) randomly assigned to reducers. This balances
load a bit before Rail-RNA-coverage_pre, where an active partition must be
handled all at once.

Input (read from stdin)
----------------------------
Tab-delimited input tuple columns:
1. Reference name (RNAME in SAM format) + ';' + 
    max(EC start, bin start) (inclusive) on forward strand IFF diff is
    positive and EC end (exclusive) on forward strand IFF diff is negative + 
    ';' + sample label + ('+' or '-' indicating which strand is the sense
    strand if input reads are strand-specific -- that is, --stranded is
    invoked; otherwise, there is no terminal '+' or '-').
2. Bin number
3. A diff -- that is, by how much coverage increases or decreases at the 
    given genomic position: +n or -n for some natural number n.
Input is partitioned by field 1.

Hadoop output (written to stdout)
----------------------------
Tab-delimited output tuple columns (collapsed):
1. Reference name (RNAME in SAM format) + ';' + bin number + ('+' or '-'
    indicating which strand is the sense strand if input reads are
    strand-specific -- that is, --stranded is invoked; otherwise, there is
    no terminal '+' or '-')
2. Sample label
3. max(EC start, bin start) (inclusive) on forward strand IFF diff is
    positive and EC end (exclusive) on forward strand IFF diff is negative.
4. A diff -- that is, by how much coverage increases or decreases at the 
    given genomic position: +n or -n for some natural number n.

IF INPUT COORDINATES ARE (0-) 1-INDEXED, OUTPUT COORDINATES ARE (0-) 1-INDEXED.
"""

import os
import sys
import time
import argparse

parser = argparse.ArgumentParser(description=__doc__, 
            formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument(\
    '--stranded', action='store_const', const=True, default=False,
    help='Assume input reads come from the sense strand; then partitions in '
         'output have terminal + and - indicating sense strand')

args = parser.parse_args()

start_time = time.time()
input_line_count, output_line_count, differential_sum = [0]*3
last_collapse_position = None
last_reverse_strand_string = ''
while True:
    line = sys.stdin.readline()
    if line:
        input_line_count += 1
        tokens = line.rstrip().split('\t')
        assert len(tokens) == 3, 'Bad input line:\n' + line
        (collapse_position, bin_number, differential) = (tokens[0],
            tokens[1], int(tokens[2]))
        (rname, pos, sample_label) = collapse_position.split(';')
        if args.stranded:
            reverse_strand_string = sample_label[-1]
            # Remove terminal +/-
            sample_label = sample_label[:-1]
        pos = int(pos)
    if not line or (last_collapse_position != collapse_position
        and last_collapse_position is not None and differential_sum != 0):
        print 'collapsed\t%s;%s%s\t%s\t%012d\t%d' % \
            (last_rname, last_bin_number, last_reverse_strand_string,
                last_sample_label, last_pos, differential_sum)
        output_line_count += 1
        differential_sum = 0
    if not line: break
    (last_collapse_position, last_rname, last_bin_number, last_sample_label,
        last_pos) = (collapse_position, rname, bin_number, sample_label, pos)
    if args.stranded: last_reverse_strand_string = reverse_strand_string
    differential_sum += differential

print >>sys.stderr, "DONE with collapse.py; in/out = %d/%d; time=%0.3f s" \
    % (input_line_count, output_line_count, start_time - time.time())