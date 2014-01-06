#!/usr/bin/env python
"""
Rail-RNA-splice_sam
Follows Rail-RNA-align
TERMINUS: no steps follow.

Input (read from stdin)
----------------------------
Tab-delimited input tuple columns are in standard SAM format. (See
http://samtools.sourceforge.net/SAMv1.pdf for detailed descriptions.)
The CIGAR includes only M's and N's, where the N's denote candidate introns.
1. QNAME
2. FLAG
3. RNAME
4. POS
5. MAPQ
6. CIGAR
7. RNEXT
8. PNEXT
9. TLEN
10. SEQ
11. QUAL
----------------------------

Hadoop output (written to stdout)
----------------------------
None.

Other output (written to file specified by command-line parameter --out)
----------------------------
SAM file with valid header consolidating SAM output lines of Rail-RNA-align.
"""
import os
import sys
import site
import argparse

base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
site.addsitedir(os.path.join(base_path, 'util'))
site.addsitedir(os.path.join(base_path, 'fasta'))

import fasta
# Define string version_number
import version
import url
import filemover

# Print file's docstring if -h is invoked
parser = argparse.ArgumentParser(description=__doc__, 
            formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('--refseq', type=str, required=True, 
    help='The fasta sequence of the reference genome. The fasta index of the '
         'reference genome is also required to be built via samtools')
# To be implemented; for now, index is always fasta filename + .fai
parser.add_argument('--faidx', type=str, required=False, 
    help='Fasta index file')
parser.add_argument(\
    '--out', metavar='URL', type=str, required=False, default=None,
    help='URL to which SAM output should be written. Default is stdout')
parser.add_argument(\
    '--sam-filename', type=str, required=False, 
    default='spliced_alignments.sam',
    help='Output SAM filename. Ignored if --out is not specified (that is, '
         'if --out is stdout)')
parser.add_argument(\
    '--verbose', action='store_const', const=True, default=False,
    help='Prints out extra debugging statements')

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
        output_filename = os.path.join(args.out, args.sam_filename)
    else:
        output_filename = args.sam_filename + '.temp'
    output_stream = open(output_filename, 'w')
else:
    # Default --out is stdout
    output_stream = sys.stdout

'''Make RNAME lengths available from reference FASTA so SAM header can be
formed; fasta_object.faidx[RNAME][0] is the length of RNAME.''' 
fasta_object = fasta.fasta(args.refseq)

# Write SAM header
print >>output_stream, '@HD\tVN:1.0\tSO:unsorted'
for rname in fasta_object.faidx:
    print >>output_stream, '@SQ\tSN:%s\tLN:%d' \
                                % (rname, fasta_object.faidx[rname][0])
print >>output_stream, '@PG\tID:Rail-RNA\tVN:%s\tCL:%s %s' \
                                % (version.version_number, sys.executable,
                                    ' '.join(sys.argv))

input_line_count = 0
for line in sys.stdin:
    line = line.rstrip()
    if len(line) == 0: continue
    print >>output_stream, line
    input_line_count += 1

if args.out is not None:
    output_stream.close()
    if not output_url.isLocal():
        mover = filemover.FileMover(args=args)
        mover.put(output_filename, output_url.plus(args.sam_filename))
        os.remove(output_filename)

print >>sys.stderr, 'DONE with splice_sam.py; in = %d; time=%0.3f s' \
                        % (input_line_count, time.time() - start_time)