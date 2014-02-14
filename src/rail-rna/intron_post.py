#!/usr/bin/env python
"""
Rail-RNA-intron_post
Follows Rail-RNA-intron
TERMINUS: no steps follow.

Reduce step in MapReduce pipelines that collects BED output encoding junctions
from Rail-RNA-intron and outputs BED, optionally by chromosome. If not
segregated by chromosome, the single file output is exactly analogous to
TopHat's junctions.bed.

Input (read from stdin)
----------------------------
Tab-delimited input tuple columns are in BED format. (See
https://genome.ucsc.edu/FAQ/FAQformat.html#format1 for a detailed description.)
Input lines mimic the lines of TopHat's junctions.bed; each denotes a different
intron.
Input is partitioned by chrom and sorted in order of ascending
chrom/chromStart/chromEnd.
1. chrom (chromosome name)
2. chromStart (start position of region; 0-BASED)
3. chromEnd (end position of region; 0-BASED)
4. name (includes anchor significance, maximum match rate, 
    and unique displacement count)
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

Other output (written to directory specified by command-line parameter --out)
----------------------------
BED file consolidating BED output lines of Rail-RNA-intron.
If --output-by-chromosome is True, each output file corresponds to a different
chrom. Otherwise, there is only one output file consolidating all BED input.
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
    help='URL to which BED output should be written. Default is stdout. '
         'If left unspecified, --output-by-chromosome is taken to be False.')
parser.add_argument(\
    '--output-by-chromosome', action='store_const', const=True, default=False,
    help='Split BED output files up by chrom (typically chromosome) if '
         'True; otherwise output single BED file consolidating all '
         'spliced alignments')
parser.add_argument(\
    '--bed-basename', type=str, required=False, 
    default='junctions',
    help='The basename (excluding path) of all BED output. Basename is '
         'followed by ".bed" if --by-chromosome=False; otherwise, basename is '
         'followed by ".[chrom].bed" for each [chrom]. Ignored if --out is '
         'not specified (that is, if --out is stdout)')

filemover.addArgs(parser)
args = parser.parse_args()

import time
start_time = time.time()

output_filename, output_stream, output_url, last_chrom = [None]*4
if args.out is not None:
    output_url = url.Url(args.out)
    if output_url.isLocal():
        # Set up destination directory
        try: os.makedirs(output_url.toUrl())
        except: pass
    else:
        # Set up temporary destination
        import tempfile
        temp_dir_path = tempfile.mkdtemp()
input_line_count = 1
move_temporary_file = False
while True:
    line = sys.stdin.readline().rstrip()
    if not line:
        last_output_filename = output_filename
        move_temporary_file = True
    else:
        (chrom, chrom_start, chrom_end, name, score, strand,
            thick_start, thick_end, item_rgb, block_count,
            block_sizes, block_starts) \
            = line.split('\t')
    if move_temporary_file and not output_url.isLocal():
        mover = filemover.FileMover(args=args)
        # Remove .temp in output filename
        mover.put(last_output_filename,
            output_url.plus(last_output_filename[:-5]))
        os.remove(last_output_filename)
        move_temporary_file = False
    if not line: break
    if (chrom != last_chrom and args.output_by_chromosome \
        and args.out is not None) or output_stream is None:
        # Create new output file
        if args.out is not None:
            if output_stream is not None: output_stream.close()
            '''If --out is a local file, just write directly to that file.
            Otherwise, write to a temporary file that will later be uploaded to
            the destination.'''
            last_output_filename = output_filename
            if output_url.isLocal():
                output_filename = os.path.join(args.out,
                    args.bed_basename
                    + (('.' + chrom) if args.output_by_chromosome else '')
                    + '.bed')
            else:
                output_filename = os.path.join(temp_dir_path, args.bed_basename
                    + (('.' + chrom) if args.output_by_chromosome else '')
                    + '.bed.temp')
                # Move last output filename
                move_temporary_file = True
            output_stream = open(output_filename, 'w')
        else:
            # Default --out is stdout
            output_stream = sys.stdout
        # Write BED header
        print >>output_stream, 'track name=%sjunctions ' \
                       'description="Rail-RNA v%s junctions%s"' \
                       % ((chrom + '_') if args.output_by_chromosome else '',
                            version.version_number,
                            (' on chrom ' + chrom)
                            if args.output_by_chromosome else '')
    '''Recall that chrom_start and chrom_end have leading 0's for proper
    sorting; remove them below.'''
    print >>output_stream, ('%s\t'*3 + 'JUNC%08d;%s\t' +
        '%s\t'*7 + '%s') \
        % (chrom, str(int(chrom_start)), str(int(chrom_end)), input_line_count,
            name, score, strand, thick_start, thick_end,
            item_rgb, block_count, block_sizes, block_starts)
    last_chrom = chrom
    input_line_count += 1

if not output_url.isLocal():
    # Clean up
    import shutil
    shutil.rmtree(temp_dir_path)

print >>sys.stderr, 'DONE with intron_post.py; in = %d; time=%0.3f s' \
                        % (input_line_count, time.time() - start_time)