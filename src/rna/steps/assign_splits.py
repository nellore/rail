"""
Rail-RNA-assign_splits

Follows Rail-RNA-count-inputs
Precedes Rail-RNA-preprocess

Counts up total number of reads across files on the _local_ filesystem and
assigns how they should be divided among workers.
"""

import os
import site
import sys
import argparse

base_path = os.path.abspath(
                    os.path.dirname(os.path.dirname(os.path.dirname(
                        os.path.realpath(__file__)))
                    )
                )
utils_path = os.path.join(base_path, 'rna', 'utils')
site.addsitedir(utils_path)
site.addsitedir(base_path)

from alignment_handlers import running_sum, pairwise
import math
import time

# Print file's docstring if -h is invoked
parser = argparse.ArgumentParser(description=__doc__, 
            formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument(
        '-p', '--num-processes', type=int, required=False, default=1,
        help=('Number of subprocesses that will be opened at once in '
              'future steps')
    )
args = parser.parse_args(sys.argv[1:])

start_time = time.time()
input_line_count, output_line_count = 0, 0
samples = {}
for i, line in enumerate(sys.stdin):
    tokens = line.strip().split('\t')
    token_count = len(tokens)
    if not (token_count >= 6 and tokens[0] == '#!splitload'):
        sys.stdout.write(line)
        output_line_count += 1
        continue
    assert token_count in [6, 8], (
            'Line {} of input has {} fields, but 6 or 8 are expected.'
        ).format(input_line_count+1, token_count)
    if token_count == 6:
        samples[(tokens[3], None)] = (int(tokens[1]), int(tokens[2])) \
                                        + tuple(tokens[3:])
    else:
        # token_count is 8
        samples[(tokens[3], tokens[5])] = (int(tokens[1])*2,
                                            int(tokens[2])*2) \
                                            + tuple(tokens[3:])
input_line_count += 1

critical_sample_values = [
            (critical_value, True) for critical_value in 
            running_sum([sample_data[0] for sample_data in samples.values()])
        ]
total_reads = sum([critical_value[0]
                    for critical_value in critical_sample_values])
reads_per_file = math.ceil(float(total_reads) / args.num_processes)
critical_read_values = [(0, False)]
candidate_critical_value = critical_read_values[-1][0] + reads_per_file
while candidate_critical_value < total_reads:
    critical_read_values.append((candidate_critical_value, False))
    candidate_critical_value = critical_read_values[-1][0] + reads_per_file
critical_values = sorted(
        list(critical_read_values + critical_sample_values),
        key=lambda crit: crit[0]
    )
sample_iter = iter(samples.keys())
try:
    current_sample = sample_iter.next()
except StopIteration:
    pass
else:
    lines_assigned, reads_assigned = [[]], 0
    for critical_pair in pairwise(critical_values):
        if critical_pair[0][-1]:
            current_sample = sample_iter.next()
            sample_index = 0
        if critical_pair[1][0] == critical_pair[0][0]:
            continue
        span = critical_pair[1][0] - critical_pair[0][0]
        if reads_assigned + span > reads_per_file:
            lines_assigned.append([])
            reads_assigned = 0
        lines_assigned[-1].append((current_sample, samples[current_sample][-1],
                                    sample_index, span))
        reads_assigned += span
    for assignment in lines_assigned:
        zipped = zip(*assignment)
        print '\t'.join(['#!splitload'] + [','.join([(element + ';NA'
                    if (type(element) is tuple and element[1] is None)
                    else (';'.join(element) if type(element) is tuple
                                            else str(element)))
                    for element in field]) for field in zipped])
        output_line_count += 1

sys.stdout.flush()
print >>sys.stderr, 'DONE with assign_splits.py; in/out=%d/%d; ' \
        'time=%0.3f s' % (input_line_count, output_line_count,
                            time.time() - start_time)