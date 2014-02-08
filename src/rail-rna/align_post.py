#!/usr/bin/env python
"""
Rail-RNA-align_post
Follows Rail-RNA-align
TERMINUS: no steps follow.

Reduce step in MapReduce pipelines that collects spliced-alignment SAM output
of Rail-RNA-align and outputs SAM or BAM, optionally by chromosome. Note that
these alignments correspond to candidate introns, each inferred from exactly
one read; Rail-RNA-intron was not used to make junction calls from clusters of
alignments. This differs from TopHat's accepted_hits.sam, which encodes not
candidate introns, but rather intron calls following an analog to
Rail-RNA-intron.

Input (read from stdin)
----------------------------
Tab-delimited input tuple columns are almost in standard SAM format; the order
of the fields is different. (See http://samtools.sourceforge.net/SAMv1.pdf for
    detailed descriptions.)
Input is partitioned by RNAME and sorted in order of ascending RNAME/POS.
The CIGAR includes only M's and N's, where the N's denote candidate introns.
1. RNAME
2. POS
3. QNAME
4. FLAG
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

Other output (written to directory specified by command-line parameter --out)
----------------------------
BAM or SAM files with valid headers consolidating SAM output lines of
Rail-RNA-align. If --output-by-chromosome is True, each output file corresponds
to a different RNAME. Otherwise, there is only one output file consolidating
all SAM input.
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
    help='URL to which SAM/BAM output should be written. Default is stdout. '
         'If left unspecified, --output-by-chromosome is taken to be False.')
parser.add_argument(\
    '--output-by-chromosome', action='store_const', const=True, default=False,
    help='Split SAM/BAM output files up by RNAME (typically chromosome) if '
         'True; otherwise output single SAM/BAM file consolidating all '
         'spliced alignments')
parser.add_argument(\
    '--bam-basename', type=str, required=False, 
    default='spliced_alignments',
    help='The basename (excluding path) of all SAM/BAM output. Basename is '
         'followed by ".bam" or ".sam" if --by-chromosome=False; otherwise, '
         'basename is followed by ".RNAME.bam" or ".RNAME.sam" for each '
         'RNAME. Ignored if --out is not specified '
         '(that is, if --out is stdout)')
parser.add_argument(\
    '--samtools-exe', metavar='EXE', type=str, required=False,
    default='samtools',
    help='Path to executable for samtools')
parser.add_argument(\
    '--output-sam', action='store_const', const=True, default=False, 
    help='Output SAM files if True; otherwise output BAM files')

filemover.addArgs(parser)
args = parser.parse_args()

if not args.output_sam:
    # Only need subprocess to start samtools if outputting bam
    import subprocess
import time
start_time = time.time()

'''Make RNAME lengths available from reference FASTA so SAM header can be
formed; fasta_object.faidx[RNAME][0] is the length of RNAME.''' 
fasta_object = fasta.fasta(args.refseq)
# Sort RNAMEs in lexicographic order so SAM file is properly sorted
sorted_rnames = sorted(fasta_object.faidx.keys())

output_filename, output_stream, output_url, last_rname = [None]*4
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
input_line_count = 0
move_temporary_file = False # True when temporary file should be uploaded
while True:
    line = sys.stdin.readline()
    if not line:
        last_output_filename = output_filename
        move_temporary_file = True
    else:
        rname, pos, qname, flag, mapq, cigar, rnext, pnext, tlen, seq, qual \
            = line.rstrip().split('\t')
    if move_temporary_file and not output_url.isLocal():
        mover = filemover.FileMover(args=args)
        # Remove .temp in output filename
        mover.put(last_output_filename, 
            output_url.plus(last_output_filename[:-5]))
        os.remove(last_output_filename)
        move_temporary_file = False
    if not line: break
    if (rname != last_rname and args.output_by_chromosome
        and args.out is not None) or output_stream is None:
        # Create new output file
        if args.out is not None:
            if output_stream is not None:
                output_stream.close()
                if not args.output_sam:
                    samtools_return = samtools_process.wait()
                    if samtools_return:
                        raise RuntimeError('samtools returned exitlevel %d' 
                            % samtools_return)
                    subprocess_stdout.close()
            '''If --out is a local file, just write directly to that file.
            Otherwise, write to a temporary file that will later be uploaded to
            the destination.'''
            last_output_filename = output_filename
            if output_url.isLocal():
                output_filename = os.path.join(args.out,
                    args.bam_basename
                    + (('.' + rname) if args.output_by_chromosome else '')
                    + ('.sam' if args.output_sam else '.bam'))
            else:
                output_filename = os.path.join(temp_dir_path, args.bam_basename
                    + (('.' + rname) if args.output_by_chromosome else '')
                    + ('.sam' if args.output_sam else '.bam')
                    + '.temp')
                # Move last output filename
                move_temporary_file = True
            output_stream = open(output_filename, 'w') if args.output_sam \
                else open(output_filename, 'wb')
        else:
            # Default --out is stdout
            output_stream = sys.stdout
            if sys.platform == 'win32' and not args.output_sam:
                # Accommodate Windows newline idiosyncrasy
                import msvcrt
                msvcrt.setmode(sys.stdout.fileno(), os.O_BINARY)
        if not args.output_sam:
            # Start samtools and reroute output_stream to subprocess stdin
            subprocess_stdout = output_stream
            samtools_process = \
                subprocess.Popen([args.samtools_exe, 'view', '-bS', '-'],
                                    stdin=subprocess.PIPE,
                                    stdout=subprocess_stdout
                                )
            output_stream = samtools_process.stdin
        '''Write SAM header; always include all reference sequences to avoid
        confusing users.'''
        print >>output_stream, '@HD\tVN:1.0\tSO:coordinate'
        for header_rname in sorted_rnames:
            print >>output_stream, '@SQ\tSN:%s\tLN:%d' \
                % (header_rname, fasta_object.faidx[header_rname][0])
        print >>output_stream, '@PG\tID:Rail-RNA\tVN:%s\tCL:%s %s' \
                                    % (version.version_number, sys.executable,
                                        ' '.join(sys.argv))
    '''Recall that pos has leading 0's so it is sorted properly; remove them
    below.'''
    print >>output_stream, ('%s\t'*10 + '%s') \
        % (qname, flag, rname, str(int(pos)), mapq, cigar, rnext, pnext, tlen,
            seq, qual)
    last_rname = rname
    input_line_count += 1

if not output_url.isLocal():
    # Clean up
    import shutil
    shutil.rmtree(temp_dir_path)

print >>sys.stderr, 'DONE with align_post.py; in=%d; time=%0.3f s' \
                        % (input_line_count, time.time() - start_time)