#!/usr/bin/env python
"""
Rail-RNA-intron_config
Follows Rail-RNA-intron
Precedes Rail-RNA-intron_fasta

Reduce step in MapReduce pipelines that outputs all possible configurations of
nonoverlapping introns on the same strand that a read of length max_len
(the maximum read length in the data being analyzed) + args.fudge can span.

Input (read from stdin)
----------------------------
Two kinds of input lines, one from Rail-RNA-align (max_len) and the other
from Rail-RNA-intron (intron):

(max_len)
1. Reference name (RNAME in SAM format) + 
    '+' or '-' indicating which strand is the sense strand
2. The character 'a'.
3. A maximum read length output by a Rail-RNA-align mapper.
4. The character '-'.

(intron)
1. Reference name (RNAME in SAM format) + 
    '+' or '-' indicating which strand is the sense strand
2. The character 'i', which will place the row after all 'a's (maximum read
    length rows) in lexicographic sort order.
3. Intron start position (inclusive)
4. Intron end position (exclusive)

Input is partitioned by strand (field 1) and sorted by the remaining fields.

Hadoop output (written to stdout)
----------------------------
Tab-delimited tuple columns:
1. Reference name (RNAME in SAM format) + 
    '+' or '-' indicating which strand is the sense strand
2. Comma-separated list of intron start positions in configuration
3. Comma-separated list of intron end positions in configuration
4. extend_size: by how many bases on either side of an intron the reference
    should extend
"""

import sys
import argparse
import time

# Print file's docstring if -h is invoked
parser = argparse.ArgumentParser(description=__doc__, 
            formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument(\
    '--verbose', action='store_const', const=True, default=False,
    help='Print out extra debugging statements')
parser.add_argument('--fudge', type=int, required=False, default=1, 
    help='A splice junction may be detected at any position along the read '
         'besides directly before or after it; thus, the sequences recorded '
         'in the new index should include max_read_size - 1 bases before and '
         'after an intron in reference space (on the forward strand). These '
         'extensions are extended further by the number of bases --fudge to '
         'accommodate possible small insertions.')

args = parser.parse_args()
start_time = time.time()
input_line_count, output_line_count = 0, 0
last_line, last_line_type, max_read_size= (None,)*3
process_introns = True
while True:
    line = sys.stdin.readline().rstrip()
    if line:
        input_line_count += 1
        if line == last_line:
            if args.verbose:
                print >>sys.stderr, 'Duplicate line encountered; continuing.'
            continue
        tokens = line.rstrip().split('\t')
        assert len(tokens) == 4
        if last_line_type is None and tokens[1] == 'i':
            raise RuntimeError('No max_len records were found.')
        if last_line_type is not None and tokens[1] != last_line_type:
            assert max_read_size is not None
            extend_size = max_read_size - 1 + args.fudge
            extend_size_string = str(extend_size)
            break
        elif tokens[1] == 'a':
            line_type, max_read_size = tokens[1], max(int(tokens[2]),
                                                        max_read_size)
            if args.verbose:
                print >>sys.stderr, 'Read a max_len record from align.'
    else:
        # No intron lines
        process_introns = False
        break
    last_line_type, last_line = line_type, line
if args.verbose:
    print >>sys.stderr, 'Finished processing max_len records from align.'
if process_introns:
    last_line, last_strand = line, tokens[0]
    intron_combos = [[(int(tokens[2]), int(tokens[3]))]]
    output_intron_combos = False
    while True:
        line = sys.stdin.readline().rstrip()
        if line:
            input_line_count += 1
            if line == last_line:
                if args.verbose:
                    print >>sys.stderr, 'Duplicate line encountered; ' \
                        'continuing.'
                continue
            tokens = line.rstrip().split('\t')
            if tokens[1] == 'a':
                # Already retrieved max_len; just continue
                continue
            assert len(tokens) == 4
            strand, pos, end_pos = tokens[0], int(tokens[2]), int(tokens[3])
            if last_strand == strand:
                new_intron_combos, kept_intron_combos = [], []
                for intron_combo in intron_combos:
                    new_intron_combo = []
                    keep_old_combo = False
                    for intron_pos, intron_end_pos in intron_combo:
                        if min(intron_end_pos, end_pos) \
                            - max(intron_pos, pos) <= 0:
                            '''If there's no overlap between the current
                            intron and an intron in an existing
                            configuration, keep that intron when forming
                            a new configuration. Otherwise, toss it.'''
                            new_intron_combo.append(
                                    (intron_pos, intron_end_pos)
                                )
                        else:
                            '''There were overlaps, so the old configuration
                            should be retained as independent.'''
                            keep_old_combo = True
                    if keep_old_combo: kept_intron_combos.append(intron_combo)
                    if len(new_intron_combo):
                        new_intron_combos.append(new_intron_combo)
                '''Output any intron combos from which current intron is
                separated by extend_size.'''
                intron_combos = []
                new_intron_combos.sort()
                for i, intron_combo in enumerate(new_intron_combos):
                    try:
                        if intron_combo == new_intron_combos[i+1]:
                            # Eliminate dupe
                            continue
                    except IndexError: pass
                    if pos - intron_combo[-1][1] >= extend_size:
                        print ('intron\t' + strand + '\t' + 
                                ','.join([str(intron_pos) for 
                                            intron_pos, _ in intron_combo])
                                + '\t' + 
                                ','.join([str(intron_end_pos) for
                                            _, intron_end_pos
                                            in intron_combo])
                                + '\t' + extend_size_string)
                        output_line_count += 1
                    else:
                        # Tack current intron onto combo
                        intron_combos.append(intron_combo + [(pos, end_pos)])
                if not len(intron_combos):
                    # If no intron combos survived, start a new one
                    intron_combos = [[(pos, end_pos)]]
                '''Add kept intron combos---configurations with which the
                current intron overlaps.'''
                intron_combos += kept_intron_combos
            else: output_intron_combos = True
        else: output_intron_combos = True
        if output_intron_combos:
            # Output last intron combos on strand
            intron_combos.sort()
            for i, intron_combo in enumerate(intron_combos):
                try:
                    if intron_combo == intron_combos[i+1]:
                        # Eliminate dupes
                        continue
                except IndexError: pass
                print ('intron\t' + strand + '\t' + 
                        ','.join([str(intron_pos) for 
                                    intron_pos, _ in intron_combo])
                        + '\t' + 
                        ','.join([str(intron_end_pos) for
                                    _, intron_end_pos in intron_combo])
                        + '\t' + extend_size_string)
                output_line_count += 1
        if not line: break
        if output_intron_combos:
            # If there's another strand to handle, set up a new intron_combos
            intron_combos = [[(pos, end_pos)]]
            output_intron_combos = False
        last_line, last_strand = line, strand

print >>sys.stderr, 'DONE with intron_config.py; in/out=%d/%d; ' \
                    'time=%0.3f s' % (input_line_count, output_line_count,
                                        time.time() - start_time)