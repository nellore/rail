#!/usr/bin/env python
"""
Rail-RNA-bam
Follows Rail-RNA-align / Rail-RNA-realign
TERMINUS: no steps follow.

Reduce step in MapReduce pipelines that collects end-to-end SAM output of 
Rail-RNA-align/spliced alignment SAM output of Rail-RNA-realign and outputs
SAM or BAM for each sample, optionally by chromosome. The output is analogous
to TopHat's accepted_hits.bam and may be used, for example, by Cufflinks.

Input (read from stdin)
----------------------------
Standard SAM except fields are in different order, and the first field 
corresponds to sample label. (Fields are reordered to facilitate partitioning
by sample name/RNAME and sorting by POS.) Each line corresponds to an
end-to-end alignment or an alignment overlapping at least one intron in the 
reference. The order of the fields is as follows.
1. Sample label
2. Number string representing RNAME; see BowtieIndexReference class in
    bowtie_index for conversion information
3. POS
4. QNAME
5. FLAG
6. MAPQ
7. CIGAR
8. RNEXT
9. PNEXT
10. TLEN
11. SEQ
12. QUAL
... + optional fields, including -- for reads overlapping introns --:
XS:A:'+' or '-' depending on which strand is the sense strand

Input is partitioned by sample label (and RNAME if outputing by chromosome)
and sorted in order of ascending RNAME/POS.

Hadoop output (written to stdout)
----------------------------
None.

Other output (written to directory specified by command-line parameter --out)
----------------------------
BAM or SAM files, at least one for each sample, with valid header
consolidating end-to-end SAM output lines of Rail-RNA-align and spliced
alignments of Rail-RNA-realign. If --output-by-chromosome is True, each output
file for a given sample corresponds to a different RNAME. Otherwise, there is
only one output file consolidating all SAM input.
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

import bowtie
import bowtie_index
import manifest
# Define string version_number
import version
import filemover
from dooplicity.ansibles import Url
import subprocess

# Print file's docstring if -h is invoked
parser = argparse.ArgumentParser(description=__doc__, 
            formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument(\
    '--out', metavar='URL', type=str, required=False, default=None,
    help='URL to which SAM/BAM output should be written. Default is stdout. '
         'If left unspecified, --output-by-chromosome is taken to be False.')
parser.add_argument(\
    '--output-by-chromosome', action='store_const', const=True, default=False,
    help='Split SAM/BAM output files up by RNAME (typically chromosome) if '
         'True; otherwise output single SAM/BAM file consolidating all'
         'alignments')
parser.add_argument('--manifest', type=str, required=False,
        default='manifest',
        help='Path to manifest file')
parser.add_argument(\
    '--bam-basename', type=str, required=False, 
    default='alignments',
    help='The basename (excluding path) of all SAM/BAM output. Basename is '
         'followed by ".[sample label].bam" or ".[sample label].sam" if '
         '--by-chromosome=False; otherwise, basename is followed by '
         '".[sample label].RNAME.bam" or ".[sample label].RNAME.sam" for each '
         'RNAME. Ignored if --out is not specified (that is, if --out is '
         'stdout)')
parser.add_argument(\
    '--samtools-exe', metavar='EXE', type=str, required=False,
    default='samtools',
    help='Path to executable for samtools')
parser.add_argument(\
    '--keep-alive', action='store_const', const=True, default=False,
    help='Prints reporter:status:alive messages to stderr to keep EMR '
         'task alive')
parser.add_argument(\
    '--output-sam', action='store_const', const=True, default=False, 
    help='Output SAM files if True; otherwise output BAM files')

filemover.add_args(parser)
bowtie.add_args(parser)
args = parser.parse_args()

if not args.output_sam:
    # Only need subprocess to start samtools if outputting bam
    import subprocess
import time
start_time = time.time()

'''Make RNAME lengths available from reference FASTA so SAM header can be
formed; reference_index.rname_lengths[RNAME] is the length of RNAME.''' 
reference_index = bowtie_index.BowtieIndexReference(args.bowtie_idx)
# Get RNAMEs in order of descending length
sorted_rnames = [reference_index.string_to_rname['%012d' % i]
                    for i in xrange(len(reference_index.string_to_rname) - 1)]
# For mapping sample indices back to original sample labels
manifest_object = manifest.LabelsAndIndices(args.manifest)

keep_alive_interval = 60
keep_alive_next = time.time() + keep_alive_interval
(output_path, output_filename, output_stream, output_url,
    last_rname, last_sample_label) = [None]*6
if args.out is not None:
    output_url = Url(args.out)
    if output_url.is_local:
        # Set up destination directory
        try: os.makedirs(output_url.to_url())
        except: pass
    else:
        mover = filemover.FileMover(args=args)
        # Set up temporary destination
        import tempfile
        temp_dir_path = tempfile.mkdtemp()
        import atexit
        from tempdel import remove_temporary_directories
        atexit.register(remove_temporary_directories,
                            [temp_dir_path])
input_line_count = 0
move_temporary_file = False # True when temporary file should be uploaded
while True:
    line = sys.stdin.readline()
    time_now = time.time()
    if time_now > keep_alive_next:
        print >>sys.stderr, 'reporter:status:alive'
        keep_alive_next += keep_alive_interval
    if not line:
        if output_stream is not None:
            output_stream.close()
            if not args.output_sam:
                samtools_return = samtools_process.wait()
                if samtools_return:
                    raise RuntimeError('samtools returned exitlevel %d' 
                        % samtools_return)
                subprocess_stdout.close()
                # Index bam
                if not output_path.endswith('.unmapped.bam'):
                    subprocess.check_call('samtools index %s'
                                            % output_path, 
                                            shell=True,
                                            bufsize=-1)
        last_output_filename = output_filename
        last_output_path = output_path
        move_temporary_file = True
    else:
        tokens = line.rstrip().split('\t')
        sample_label, rname, pos, qname, flag = tokens[:5]
        sample_label = manifest_object.index_to_label[sample_label]
        rname = reference_index.string_to_rname[rname]
    if move_temporary_file and last_sample_label is not None \
        and not output_url.is_local:
        mover.put(last_output_path, 
            output_url.plus(last_output_filename))
        mover.put(last_output_path + '.bai',
                    output_url.plus(last_output_filename + '.bai'))
        os.remove(last_output_path)
        os.remove(last_output_path + '.bai')
        move_temporary_file = False
    if not line: break
    if ((sample_label != last_sample_label or 
            (rname != last_rname and args.output_by_chromosome)) 
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
                    # Index bam
                    if not output_path.endswith('.unmapped.bam'):
                        subprocess.check_call('samtools index %s'
                                                % output_path, 
                                                shell=True,
                                                bufsize=-1)
            '''If --out is a local file, just write directly to that file.
            Otherwise, write to a temporary file that will later be uploaded to
            the destination.'''
            if rname == '*':
                out_rname = 'unmapped'
            else:
                out_rname = rname
            last_output_filename = output_filename
            last_output_path = output_path
            output_filename = args.bam_basename + '.' + sample_label \
                    + (('.' + out_rname) if args.output_by_chromosome else ''
                        ) + ('.sam' if args.output_sam else '.bam')
            if output_url.is_local:
                output_path = os.path.join(args.out, output_filename)
            else:
                output_path = os.path.join(temp_dir_path, output_filename)
                if last_output_filename is not None:
                    # Move last output filename iff there is one
                    move_temporary_file = True
            output_stream = open(output_path, 'w') if args.output_sam \
                else open(output_path, 'wb')
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
                % (header_rname, reference_index.rname_lengths[header_rname])
        print >>output_stream, '@PG\tID:Rail-RNA\tVN:%s\tCL:%s %s' \
                                    % (version.version_number, sys.executable,
                                        ' '.join(sys.argv))
    '''Recall that pos has leading 0's so it is sorted properly; remove them
    below.'''
    print >>output_stream, ((('%s\t'*4) % (qname[:254], flag, rname,
                                            str(int(pos))))
                                            + '\t'.join(tokens[5:]))
    last_rname = rname
    last_sample_label = sample_label
    input_line_count += 1

if not output_url.is_local:
    # Clean up
    import shutil
    shutil.rmtree(temp_dir_path)

print >>sys.stderr, 'DONE with bam.py; in=%d; time=%0.3f s' \
                        % (input_line_count, time.time() - start_time)
