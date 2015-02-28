#!/usr/bin/env python
"""
Rail-RNA-coverage_post
Follows Rail-RNA-coverage
TERMINUS: no steps follow.

Reduce step in MapReduce pipelines that takes results from the coverage step as
a single bin and writes out a new manifest file annotated with normalization
factors.
 
Input (read from stdin)
----------------------------
Tab-delimited input tuple columns (just 1 per sample label):
1. Sample label
2. Normalization factor
No binning/sorting before this step.

Hadoop output (written to stdout)
----------------------------
None.

Other output (written to directory specified by command-line parameter --out)
----------------------------
File whose name is normalization_factors.tsv by default.
Tab-delimited tuple columns:
1. The character '-', ensuring there's exactly one partition
2. Sample name
3. Normalization factor
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
from dooplicity.tools import register_cleanup, make_temp_dir
import filemover
import tempdel

# Print file's docstring if -h is invoked
parser = argparse.ArgumentParser(description=__doc__, 
            formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument(\
    '--out', metavar='URL', type=str, required=False, default=None,
    help='URL to which output should be written. Default is stdout')
parser.add_argument(\
    '--normalize-filename', type=str, required=False, 
    default='normalization_factors.tsv',
    help='The output filename (excluding path). '
         'Ignored if --out is not specified (that is, if --out is stdout)')
parser.add_argument(\
    '--verbose', action='store_const', const=True, default=False,
    help='Print out extra debugging statements')

manifest.add_args(parser)
filemover.add_args(parser)
tempdel.add_args(parser)
args = parser.parse_args()

import time
start_time = time.time()

if args.verbose:
    print >>sys.stderr, ' vvv Manifest file vvv'
    with open(args.manifest) as manifest_stream:
        print >>sys.stderr, manifest_stream.read()
    print >>sys.stderr, ' ^^^ Manifest file ^^^ '

'''Get the set of all labels by parsing the manifest file, given on the
filesystem or in the Hadoop file cache.'''
sample_labels = manifest.LabelsAndIndices(args.manifest).label_to_index.keys()
sample_labels.sort()

input_line_count = 0
normalization_factors = {}

for line in sys.stdin:
    tokens = line.rstrip().split('\t')[1:] # Kill partition
    assert len(tokens) == 2, 'Bad input line: %s' % line
    sample_label, normalization_factor = tokens[0], int(tokens[1])
    normalization_factors[sample_label] = normalization_factor
    input_line_count += 1

if args.out is not None:
    '''If --out is a local file, just write directly to that file. Otherwise,
    write to a temporary file that will later be uploaded to the
    destination.'''
    output_url = Url(args.out)
    if output_url.is_local:
        try: os.makedirs(output_url.to_url())
        except: pass
        output_filename = os.path.join(args.out, args.normalize_filename)
    else:
        temp_dir_path = make_temp_dir(args.scratch)
        register_cleanup(tempdel.remove_temporary_directories,
                            [temp_dir_path])
        output_filename = args.normalize_filename + '.temp'
        output_filename = os.path.join(temp_dir_path, output_filename)
    output_stream = open(output_filename, 'w')
else:
    # Default --out is stdout
    output_stream = sys.stdout

# Write normaliation factors where available; otherwise write N/A
for sample_label in sample_labels:
    if sample_label in normalization_factors:
        print >>output_stream, '%s\t%d' % (sample_label,
            normalization_factors[sample_label])
    else:
        print >>output_stream, '%s\tN/A' % sample_label

if args.out is not None:
    output_stream.close()
    if not output_url.is_local:
        mover = filemover.FileMover(args=args)
        mover.put(output_filename, output_url.plus(args.normalize_filename))

print >>sys.stderr, 'DONE with coverage_post.py; in = %d; time=%0.3f s' \
                        % (input_line_count, time.time() - start_time)
