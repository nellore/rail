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
1. Sample name
2. Normalization factor
"""
import os
import sys
import site
import argparse
 
base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
site.addsitedir(os.path.join(base_path, "manifest"))
site.addsitedir(os.path.join(base_path, "util"))

import manifest
import url
import filemover

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

manifest.addArgs(parser)
filemover.addArgs(parser)
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
sample_labels = manifest.labels(args)
sample_labels.sort()

input_line_count = 0
normalization_factors = {}

for line in sys.stdin:
    tokens = line.rstrip().split('\t')
    assert len(tokens) == 2, 'Bad input line: %s' % line
    sample_label, normalization_factor = tokens[0], int(tokens[1])
    normalization_factors[sample_label] = normalization_factor
    input_line_count += 1

if args.out is not None:
    '''If --out is a local file, just write directly to that file. Otherwise,
    write to a temporary file that will later be uploaded to the
    destination.'''
    output_url = url.Url(args.out)
    if output_url.isLocal():
        try: os.makedirs(output_url.toUrl())
        except: pass
        output_filename = os.path.join(args.out, args.normalize_filename)
    else:
        import tempfile
        temp_dir_path = tempfile.mkdtemp()
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
    if not output_url.isLocal():
        mover = filemover.FileMover(args=args)
        mover.put(output_filename, output_url.plus(args.normalize_filename))
        import shutil
        shutil.rmtree(temp_dir_path)

print >>sys.stderr, 'DONE with coverage_post.py; in = %d; time=%0.3f s' \
                        % (input_line_count, time.time() - start_time)
