#!/usr/bin/env python
"""
Rail-RNA-junction_collect
Follows Rail-RNA-junction_filter
TERMINUS: no steps follow.

Reduce step in MapReduce pipelines that takes junctions output by
junction_filter and writes/compresses/uploads a single output file containing
junctions and their initially detected coverages by sample
 
Input (read from stdin)
----------------------------
Tab-delimited tuple columns:
1. Reference name (RNAME in SAM format) +
    '+' or '-' indicating which strand is the sense strand
2. Intron start position (inclusive)
3. Intron end position (exclusive)
4. comma-separated list of sample indexes in which junction was found
5. comma-separated list of numbers of reads in which junction was found
    in respective sample specified by field 4

Input is keyed and sorted by the first three fields.

Hadoop output (written to stdout)
----------------------------
None.

Other output (written to directory specified by command-line parameter --out)
----------------------------
File whose name is collected_junctions.tsv.gz by default.
Tab-delimited tuple columns:
1. Reference name (RNAME in SAM format) +
    '+' or '-' indicating which strand is the sense strand
2. Intron start position (inclusive)
3. Intron end position (inclusive)
4. comma-separated list of sample indexes in which junction was found
5. comma-separated list of numbers of reads in which junction was found
    in respective sample specified by field 4
"""
import os
import sys
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

import manifest
from dooplicity.ansibles import Url
from dooplicity.tools import register_cleanup, make_temp_dir, xopen
import filemover
import tempdel

# Print file's docstring if -h is invoked
parser = argparse.ArgumentParser(description=__doc__, 
            formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument(\
    '--out', metavar='URL', type=str, required=False, default=None,
    help='URL to which output should be written. Default is stdout')
parser.add_argument(\
    '--junction-filename', type=str, required=False, 
    default='first_pass_junctions.tsv.gz',
    help='The output filename (excluding path). '
         'Ignored if --out is not specified (that is, if --out is stdout)')
parser.add_argument('--gzip-level', type=int, required=False,
        default=3,
        help=('Level of gzip compression to use for temporary file storing '
              'qnames.'))
parser.add_argument(\
    '--verbose', action='store_const', const=True, default=False,
    help='Print out extra debugging statements')

filemover.add_args(parser)
tempdel.add_args(parser)
args = parser.parse_args()

import time
start_time = time.time()

input_line_count = 0

if args.out is not None:
    '''If --out is a local file, just write directly to that file. Otherwise,
    write to a temporary file that will later be uploaded to the
    destination.'''
    output_url = Url(args.out)
    if output_url.is_local:
        try: os.makedirs(output_url.to_url())
        except: pass
        output_filename = os.path.join(args.out, args.junction_filename)
    else:
        temp_dir_path = make_temp_dir(tempdel.silentexpandvars(args.scratch))
        register_cleanup(tempdel.remove_temporary_directories,
                            [temp_dir_path])
        output_filename = args.junction_filename + '.temp'
        output_filename = os.path.join(temp_dir_path, output_filename)
    with xopen(True, output_filename, 'w', args.gzip_level) as output_stream:
        for line in sys.stdin:
            output_stream.write(line)
            input_line_count += 1
else:
    # Default --out is stdout
    for line in sys.stdin:
        tokens = line.strip().split('\t')
        # Remove leading zeros from ints
        sys.stdout.write('\t'.join([tokens[0], str(int(tokens[1])),
                                    str(int(tokens[2]) - 1), tokens[3],
                                    tokens[4]]))
        input_line_count += 1

if args.out is not None and not output_url.is_local:
    mover = filemover.FileMover(args=args)
    mover.put(output_filename, output_url.plus(args.junction_filename))

print >>sys.stderr, 'DONE with junction_collect.py; in = %d; time=%0.3f s' \
                        % (input_line_count, time.time() - start_time)