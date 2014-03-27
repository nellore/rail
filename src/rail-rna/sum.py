"""
Rail-RNA-sum
Used at various points in pipeline

Generic combiner/reducer script for MapReduce pipelines that sums values that
have the same key. A value is an integer/float in the final field, and a key is
the concatenation of fields that precede the final field.

Input (read from stdin)
----------------------------
Tab-delimited input tuple columns
1. Any string
2.     ""
.
.
.
N.     ""
N+1. Integer

Hadoop output (written to stdout)
----------------------------
Tab-delimited output tuple columns:
1. Any string
2.     ""
.
.
.
N.     ""
N+1. Integer

[1, ..., N] now span a unique key, and the value of N+1 is the sum of the
values for all such keys from the input.
"""

import argparse
import sys
import time

start_time = time.time()

# Print file's docstring if -h is invoked
parser = argparse.ArgumentParser(description=__doc__, 
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
parser.add_argument(\
        '--float', action='store_const', const=True, default=False,
        help='Assume final field is a floating-point number rather than '
             'an integer'
    )

args = parser.parse_args()

input_line_count, output_line_count = 0, 0

if args.float:
    last_key, total, write_line = None, 0.0, False
    while True:
        line = sys.stdin.readline()
        if not line:
            if last_key is None:
                # Input is empty
                break
            else:
                # Write final line
                write_line = True
        else:
            input_line_count += 1
            tokens = line.rstrip().split('\t')
            key = tokens[:-1]
            if key != last_key and last_key is not None:
                write_line = True
        if write_line:
            print '\t'.join(last_key + ['%0.11f' % total])
            output_line_count += 1
            total, write_line = 0.0, False
        total += float(tokens[-1])
        last_key = key
else:
    last_key, total, write_line = None, 0, False
    while True:
        line = sys.stdin.readline()
        if not line:
            if last_key is None:
                # Input is empty
                break
            else:
                # Write final line
                write_line = True
        else:
            input_line_count += 1
            tokens = line.rstrip().split('\t')
            key = tokens[:-1]
            if key != last_key and last_key is not None:
                write_line = True
        if write_line:
            print '\t'.join(last_key + ['%d' % total])
            output_line_count += 1
            total, write_line = 0, False
        total += int(tokens[-1])
        last_key = key

print >>sys.stderr, 'DONE with sum.py; in/out=%d/%d; time=%0.3f s' \
                        % (input_line_count, output_line_count, 
                            time.time() - start_time)
