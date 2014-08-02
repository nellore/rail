"""
Rail-RNA-sum

Used at various points in pipeline

Generic combiner/reducer script for MapReduce pipelines that sums values that
have the same key. A value is an integer/float/string from the final K fields,
and a key is the concatenation of fields that precede the final K fields.
Strings are summed by separating them with '\x1d's. An empty string is
represented as '\x1c'.

Input (read from stdin)
----------------------------
Tab-delimited input tuple columns
1. Any string
2.     ""
.
.
.
N.     ""
N+1. Integer/float/string
N+2. Integer/float/string
.
.
.
N+K-1. Integer/float/string

Hadoop output (written to stdout)
----------------------------
Tab-delimited output tuple columns:
1. Any string
2.     ""
.
.
.
N.     ""
N+1. Integer/float/'\x1d'-separated list of strings
N+2. Integer/float/'\x1d'-separated list of strings
.
.
.
N+K-1. Integer/float/comma-separated list of strings

[1, ..., N] now span a unique key, and the value of each of [N+1, ..., N+K-1]
is the "sum" of the values for all such keys from the input.
"""

import argparse
import sys
import time

start_time = time.time()

# Print file's docstring if -h is invoked
parser = argparse.ArgumentParser(description=__doc__, 
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
parser.add_argument('--value-count', type=int, required=False,
        default=1,
        help='Number of value fields; this is K from the docstring')
parser.add_argument(\
        '--type', type=int, required=False, default=1,
        help='If 1, assume final fields are ints. If 2, assume final fields '
             'are floats. Otherwise, assume final fields are strings.'
    )
parser.add_argument(
        '--keep-alive', action='store_const', const=True,
        default=False,
        help='Periodically print Hadoop status messages to stderr to keep ' \
             'job alive'
    )

args = parser.parse_args()

if args.keep_alive:
    import site
    import os

    base_path = os.path.abspath(
                        os.path.dirname(os.path.dirname(os.path.dirname(
                            os.path.realpath(__file__)))
                        )
                    )
    utils_path = os.path.join(base_path, 'rna', 'utils')
    site.addsitedir(utils_path)
    site.addsitedir(base_path)

    from dooplicity.tools import KeepAlive
    keep_alive_thread = KeepAlive(sys.stderr)
    keep_alive_thread.start()

input_line_count, output_line_count = 0, 0

if args.type == 1:
    last_key, totals, write_line = None, [0]*args.value_count, False
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
            key = tokens[:-args.value_count]
            if key != last_key and last_key is not None:
                write_line = True
        if write_line:
            print '\t'.join(last_key + [('%d' % totals[i])
                                            for i in xrange(len(totals))])
            output_line_count += 1
            totals, write_line = [0]*args.value_count, False
        if not line: break
        for i in xrange(1, args.value_count+1):
            totals[-i] += int(tokens[-i])
        last_key = key
elif args.type == 2:
    last_key, totals, write_line = None, [0.0]*args.value_count, False
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
            key = tokens[:-args.value_count]
            if key != last_key and last_key is not None:
                write_line = True
        if write_line:
            print '\t'.join(last_key + [('%0.11f' % totals[i])
                                            for i in xrange(len(totals))])
            output_line_count += 1
            totals, write_line = [0.0]*args.value_count, False
        if not line: break
        for i in xrange(1, args.value_count+1):
            totals[-i] += float(tokens[-i])
        last_key = key
else:
    last_key, totals, write_line \
        = None, [[] for i in xrange(args.value_count)], False
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
            key = tokens[:-args.value_count]
            if key != last_key and last_key is not None:
                write_line = True
        if write_line:
            print '\t'.join(last_key + ['\x1d'.join(totals[i])
                                            if totals[i] != [] else '\x1c'
                                            for i in xrange(len(totals))])
            output_line_count += 1
            totals, write_line = [[] for i in xrange(args.value_count)], False
        if not line: break
        for i in xrange(1, args.value_count+1):
            if tokens[-i] != '\x1c':
                totals[-i].append(tokens[-i])
        last_key = key

print >>sys.stderr, 'DONE with sum.py; in/out=%d/%d; time=%0.3f s' \
                        % (input_line_count, output_line_count, 
                            time.time() - start_time)
