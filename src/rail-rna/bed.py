"""
Rail-RNA-bed
Follows Rail-RNA-bed_pre
TERMINUS: no steps follow.

Consolidates BED output lines from Rail-RNA-bed_pre into TopHat-like
junctions.bed files, one for each sample.

Input (read from stdin)
----------------------------
Tab-delimited input tuple columns are in sample label + 12-column BED format.
(See https://genome.ucsc.edu/FAQ/FAQformat.html#format1 for a detailed
description.)

1. Sample label
2. Number string representing chrom (chromosome name); see BowtieIndexReference
    class in bowtie_index for conversion information
3. chromStart (start position of region; 0-BASED)
4. chromEnd (end position of region; 0-BASED)
5. name (includes maximin anchor size and unique displacement count)
6. score (number of reads supporting this intron)
7. strand (forward is sense strand=+; reverse is sense strand=-)
8. thickStart (same as chromStart)
9. thickEnd (same as chromEnd)
10. itemRgb (always 227,29,118; that is, hot pink)
11. blockCount (always 2; each block denotes the maximal left/right overhang of
    any read spanning a junction)
12. blockSizes (comma-separated lengths of the maximal left/right overhangs of
    any read spanning a junction)
13. blockStarts (comma-separated list of the overhang start positions)

Input is partitioned by column 1 and sorted by columns 2/3/4.

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
for directory_name in ['bowtie', 'util']:
    site.addsitedir(os.path.join(base_path, directory_name))
import bowtie
import bowtie_index

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
    '--bed-basename', type=str, required=False, 
    default='junctions',
    help='The basename (excluding path) of all BED output. Basename is '
         'followed by ".bed"')

bowtie.addArgs(parser)
filemover.addArgs(parser)
args = parser.parse_args()

import time
start_time = time.time()

reference_index = bowtie_index.BowtieIndexReference(args.bowtie_idx)
(output_filename, last_output_filename, output_path, output_stream, output_url,
    last_sample_label) = [None]*6
if args.out is not None:
    output_url = url.Url(args.out)
    if output_url.isLocal():
        # Set up destination directory
        try: os.makedirs(output_url.toUrl())
        except: pass
    else:
        mover = filemover.FileMover(args=args)
        # Set up temporary destination
        import tempfile
        temp_dir_path = tempfile.mkdtemp()
input_line_count = 1
output_line_count = 0
move_temporary_file = False
while True:
    line = sys.stdin.readline().rstrip()
    if not line:
        if output_stream is not None:
            output_stream.close()
        last_output_filename = output_filename
        last_output_path = output_path
        move_temporary_file = True
    else:
        (sample_label, chrom,
            chrom_start, chrom_end, name, score, strand,
            thick_start, thick_end, item_rgb, block_count,
            block_sizes, block_starts) \
            = line.split('\t')
        chrom = reference_index.string_to_rname[chrom]
    if move_temporary_file and last_sample_label is not None \
        and not output_url.isLocal():
        mover.put(last_output_path,
            output_url.plus(last_output_filename))
        os.remove(last_output_path)
        move_temporary_file = False
    if not line: break
    if (sample_label != last_sample_label and args.out is not None) \
        or output_stream is None:
        # Create new output file
        if args.out is not None:
            if output_stream is not None: output_stream.close()
            '''If --out is a local file, just write directly to that file.
            Otherwise, write to a temporary file that will later be uploaded to
            the destination.'''
            last_output_filename = output_filename
            last_output_path = output_path
            output_filename = args.bed_basename  + '.' + sample_label + '.bed'
            if output_url.isLocal():
                output_path = os.path.join(args.out, output_filename)
            else:
                output_path = os.path.join(temp_dir_path, output_filename)
                if last_output_filename is not None:
                    # Move last output filename iff there is one
                    move_temporary_file = True
            output_stream = open(output_path, 'w')
        else:
            # Default --out is stdout
            output_stream = sys.stdout
        # Write BED header
        print >>output_stream, 'track name=%s_junctions ' \
                       'description="Rail-RNA v%s junctions"' \
                       % (sample_label, version.version_number)
        output_line_count += 1
    '''Recall that chrom_start and chrom_end have leading 0's for proper
    sorting; remove them below.'''
    print >>output_stream, ('%s\t'*3 + 'JUNC%08d;%s\t' +
        '%s\t'*7 + '%s') \
        % (chrom, str(int(chrom_start)), str(int(chrom_end)), input_line_count,
            name, score, strand, thick_start, thick_end,
            item_rgb, block_count, block_sizes, block_starts)
    output_line_count += 1
    last_sample_label = sample_label
    input_line_count += 1

if not output_url.isLocal():
    # Clean up
    import shutil
    shutil.rmtree(temp_dir_path)

print >>sys.stderr, 'DONE with bed.py; in/out = %d/%d; time=%0.3f s' \
                        % (input_line_count, output_line_count,
                            time.time() - start_time)