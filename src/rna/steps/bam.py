#!/usr/bin/env python
"""
Rail-RNA-bam
Follows Rail-RNA-align / Rail-RNA-realign
Precedes Rail-RNA-collect_read_stats

Reduce step in MapReduce pipelines that collects end-to-end SAM output of 
Rail-RNA-align/spliced alignment SAM output of other steps and outputs
SAM or BAM for each sample, optionally by chromosome. The output is analogous
to TopHat's accepted_hits.bam and may be used, for example, by Cufflinks.

Input (read from stdin)
----------------------------
Standard SAM except fields are in different order, and the first field 
corresponds to sample label. (Fields are reordered to facilitate partitioning
by sample name/RNAME and sorting by POS.) Each line corresponds to an
end-to-end alignment or an alignment overlapping at least one intron in the 
reference. The order of the fields is as follows.
1. Sample index if outputting BAMs by sample OR sample-rname index if
    outputting BAMs by chr
2. (Number string representing RNAME; see BowtieIndexReference class in
    bowtie_index for conversion information) OR '0' if outputting BAMs by chr
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
Tab-delimited output columns (collect):
1. '-' to enforce that all records are placed in the same partition
2. Sample index
3. RNAME index
4. Number of primary alignments overlapping contig with RNAME
5. Number of unique alignments overlapping contig with RNAME

If RNAME corresponds to unmapped reads, unique total = 0 and primary total =
number of unmapped reads.

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
from dooplicity.tools import register_cleanup, make_temp_dir, xstream
from alignment_handlers import SampleAndRnameIndexes
import subprocess
import tempdel

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
parser.add_argument(\
    '--suppress-bam', action='store_const', const=True, default=False, 
    help='Do not write any alignments; overrides all other output parameters')

filemover.add_args(parser)
bowtie.add_args(parser)
tempdel.add_args(parser)
from alignment_handlers import add_args as alignment_handlers_add_args
alignment_handlers_add_args(parser)
args = parser.parse_args()

# Start keep_alive thread immediately
if args.keep_alive:
    from dooplicity.tools import KeepAlive
    keep_alive_thread = KeepAlive(sys.stderr)
    keep_alive_thread.start()

'''Make RNAME lengths available from reference FASTA so SAM header can be
formed; reference_index.rname_lengths[RNAME] is the length of RNAME.''' 
reference_index = bowtie_index.BowtieIndexReference(args.bowtie_idx)
# For mapping sample indices back to original sample labels
manifest_object = manifest.LabelsAndIndices(args.manifest)
# To convert sample-rname index to sample index-rname index tuple
sample_and_rname_indexes = SampleAndRnameIndexes(
                                                    manifest_object,
                                                    args.output_by_chromosome
                                                )

import time
start_time = time.time()
from alignment_handlers import AlignmentPrinter
alignment_printer = AlignmentPrinter(manifest_object, reference_index,
                                        tie_margin=args.tie_margin)
input_line_count = 0
if args.suppress_bam:
    # Just grab stats
    if args.output_by_chromosome:
        for (index, _), xpartition in xstream(sys.stdin, 2):
            sample_index, rname_index = (
                    sample_and_rname_indexes.sample_and_rname_indexes(index)
                )
            unique_count, total_count = 0, 0
            for record in xpartition:
                if not (int(record[2]) & 256):
                    total_count += 1
                    try:
                        # seq is at position 8
                        if alignment_printer.unique(record, seq_index=8):
                            unique_count += 1
                    except IndexError:
                        # Unmapped read; it's unique
                        unique_count += 1
                input_line_count += 1
            # Only primary alignments (flag & 256 != 1)
            print 'counts\t-\t%s\t%s\t%d\t%d' % (sample_index, rname_index,
                                                 total_count, unique_count)
    else:
        for (sample_index, rname_index), xpartition in xstream(sys.stdin, 2):
            unique_count, total_count = 0, 0
            for record in xpartition:
                if not (int(record[2]) & 256):
                    total_count += 1
                    try:
                        if alignment_printer.unique(record, seq_index=8):
                            unique_count += 1
                    except IndexError:
                        # Unmapped read; it's unique
                        unique_count += 1
            # Only primary alignments (flag & 256 != 1)
            print 'counts\t-\t%s\t%s\t%d\t%d' % (sample_index, rname_index,
                                                 total_count, unique_count)
else:
    # Grab stats _and_ output SAM/BAMs
    if not args.output_sam:
        # Only need subprocess to start samtools if outputting bam
        import subprocess

    # Get RNAMEs in order of descending length
    sorted_rnames = [reference_index.string_to_rname['%012d' % i]
                        for i in xrange(
                                    len(reference_index.string_to_rname) - 1
                                )]
    (output_path, output_filename, output_stream, output_url,
        last_rname, last_sample_label) = [None]*6
    total_count, unique_count = 0, 0
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
            temp_dir_path = make_temp_dir(args.scratch)
            register_cleanup(tempdel.remove_temporary_directories,
                                [temp_dir_path])
    move_temporary_file = False # True when temporary file should be uploaded
    while True:
        line = sys.stdin.readline()
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
                        subprocess.check_call([args.samtools_exe, 'index',
                                                output_path],
                                                bufsize=-1)
            last_output_filename = output_filename
            last_output_path = output_path
            move_temporary_file = True
        else:
            tokens = line.rstrip().split('\t')
            sample_index, rname_index, pos, qname, flag = tokens[:5]
            if args.output_by_chromosome:
                (sample_index, rname_index) \
                    = sample_and_rname_indexes.sample_and_rname_indexes(
                            sample_index
                        )
            sample_label = manifest_object.index_to_label[sample_index]
            rname = reference_index.string_to_rname[rname_index]
        if move_temporary_file and last_sample_label is not None \
            and not output_url.is_local:
            mover.put(last_output_path, 
                output_url.plus(last_output_filename))
            os.remove(last_output_path)
            if not last_output_path.endswith('.unmapped.bam'):
                mover.put(
                    ''.join([last_output_path, '.bai']),
                    output_url.plus(''.join([last_output_filename, '.bai']))
                )
                os.remove(''.join([last_output_path, '.bai']))
            move_temporary_file = False
        try:
            if (sample_label != last_sample_label or rname != last_rname
                or not line):
                print 'counts\t-\t%s\t%s\t%d\t%d' % (last_sample_index,
                                                     last_rname_index,
                                                     total_count, unique_count)
                total_count, unique_count = 0, 0
        except NameError:
            # First record
            pass
        if not line:
            break
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
                            subprocess.check_call([args.samtools_exe, 'index',
                                                    output_path],
                                                    bufsize=-1)
                '''If --out is a local file, just write directly to that file.
                Otherwise, write to a temporary file that will later be
                uploaded to the destination.'''
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
            '''Write SAM header; always include all reference sequences to
            avoid confusing users.'''
            print >>output_stream, '@HD\tVN:1.0\tSO:coordinate'
            for header_rname in sorted_rnames:
                print >>output_stream, '@SQ\tSN:%s\tLN:%d' \
                    % (header_rname,
                        reference_index.rname_lengths[header_rname])
            print >>output_stream, '@PG\tID:Rail-RNA\tVN:%s\tCL:%s %s' % (
                                    version.version_number,
                                    sys.executable,
                                    ' '.join(sys.argv)
                                )
        '''Recall that pos has leading 0's so it is sorted properly; remove
        them below.'''
        sam_line_to_print = [qname[:254], flag, rname,
                                str(int(pos))] + [('ZS:i:' + token[5:]
                                                    if token[:5] == 'XS:i:'
                                                    else token)
                                                    for token in tokens[5:]]
        try:
            print >>output_stream, '\t'.join(sam_line_to_print)
        except IOError:
            raise IOError('Error writing line "%s".' % sam_line_to_print)
        if not (int(flag) & 256):
            total_count += 1
            try:
                if alignment_printer.unique(sam_line_to_print):
                    unique_count += 1
            except IndexError:
                # Unmapped read; it's unique
                unique_count += 1
        last_rname = rname
        last_sample_label = sample_label
        last_rname_index = rname_index
        last_sample_index = sample_index
        input_line_count += 1

print >>sys.stderr, 'DONE with bam.py; in=%d; time=%0.3f s' % (
                                input_line_count, time.time() - start_time
                            )