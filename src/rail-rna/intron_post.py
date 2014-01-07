#!/usr/bin/env python
"""
Rail-RNA-intron_post
Follows Rail-RNA-intron
TERMINUS: no steps follow.

Input (read from stdin)
----------------------------
Tab-delimited input tuple columns are in BED format. (See
https://genome.ucsc.edu/FAQ/FAQformat.html#format1 for a detailed description.)
Input lines mimic the lines of TopHat's junctions.bed; each denotes a different
intron.
1. chrom (chromosome name)
2. chromStart (start position of region)
3. chromEnd (end position of region)
4. name (JUNC[number]; will be renumbered in this script)
5. score (number of reads supporting this intron)
6. strand (forward is sense strand=+; reverse is sense strand=-)
7. thickStart (same as chromStart)
8. thickEnd (same as chromEnd)
9. itemRgb (always 255,0,0; that is, red)
10. blockCount (always 2; each block denotes the maximal left/right overhang of
    any read spanning a junction)
11. blockSizes (comma-separated lengths of the maximal left/right overhangs of
    any read spanning a junction)
12. blockStarts (comma-separated list of the overhang start positions)
----------------------------

Hadoop output (written to stdout)
----------------------------
None.

Other output (written to file specified by command-line parameter --out)
----------------------------
BED file consolidating BED output lines of Rail-RNA-intron.
"""
import os
import sys
import site
import argparse

base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
site.addsitedir(os.path.join(base_path, 'util'))

# Define string version_number
import version
import url
import filemover

# Print file's docstring if -h is invoked
parser = argparse.ArgumentParser(description=__doc__, 
            formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument(\
    '--out', metavar='URL', type=str, required=False, default=None,
    help='URL to which SAM output should be written. Default is stdout')
parser.add_argument(\
    '--bed-filename', type=str, required=False, 
    default='junctions.bed',
    help='Output BED filename. Ignored if --out is not specified (that is, '
         'if --out is stdout)')

filemover.addArgs(parser)
args = parser.parse_args()

import time
start_time = time.time()

output_filename, output_url = None, None

if args.out is not None:
    '''If --out is a local file, just write directly to that file. Otherwise,
    write to a temporary file that will later be uploaded to the
    destination.'''
    output_url = url.Url(args.out)
    if output_url.isLocal():
        try: os.makedirs(output_url.toUrl())
        except: pass
        output_filename = os.path.join(args.out, args.bed_filename)
    else:
        output_filename = args.bed_filename + '.temp'
    output_stream = open(output_filename, 'w')
else:
    # Default --out is stdout
    output_stream = sys.stdout

# Write BED header
print >>output_stream, 'track name=junctions description="Rail-RNA ' \
                       'v%s junctions"' % version.version_number

input_line_count = 0
for line in sys.stdin:
    line = line.rstrip()
    input_line_count += 1
    if len(line) == 0: continue
    tokens = line.split('\t')
    # Renumber junction
    tokens[3] = 'JUNC%08d' % input_line_count
    print >>output_stream, '\t'.join(tokens)

if args.out is not None:
    output_stream.close()
    if not output_url.isLocal():
        mover = filemover.FileMover(args=args)
        mover.put(output_filename, output_url.plus(args.bed_filename))
        os.remove(output_filename)

print >>sys.stderr, 'DONE with intron_post.py; in = %d; time=%0.3f s' \
                        % (input_line_count, time.time() - start_time)