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
end-to-end alignment or an alignment overlapping at least one junction in the 
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
... + optional fields, including -- for reads overlapping junctions --:
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

if '--test' in sys.argv:
    print("No unit tests")
    #unittest.main(argv=[sys.argv[0]])
    sys.exit(0)

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
from dooplicity.counters import Counter
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
reference_index = bowtie_index.BowtieIndexReference(
                            os.path.expandvars(args.bowtie_idx)
                        )
# For mapping sample indices back to original sample labels
manifest_object = manifest.LabelsAndIndices(os.path.expandvars(args.manifest))
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
counter = Counter('bam')
register_cleanup(counter.flush)

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
            counter.add('print_stats')
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
            counter.add('print_stats')
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
    total_count, unique_count = 0, 0
    if args.out is not None:
        output_url = Url(args.out)
        if output_url.is_local:
            # Set up destination directory
            try: os.makedirs(output_url.to_url())
            except: pass
            output_dir = args.out
        else:
            mover = filemover.FileMover(args=args)
            # Set up temporary destination
            import tempfile
            temp_dir_path = make_temp_dir(
                                tempdel.silentexpandvars(args.scratch)
                            )
            register_cleanup(tempdel.remove_temporary_directories,
                                [temp_dir_path])
            output_dir = temp_dir_path

    from contextlib import contextmanager
    @contextmanager
    def stream_and_upload(rnames, filename=None, mover=None, output_url=None,
                            sam=False, samtools_exe=None):
        """ Yields output stream to write to and uploads as necessary

            sorted_rnames: list of rnames in order of descending length
            filename: full path to file to write or None if writing to stdout
            mover: FileMover object or None if no moving should be performed
            output_url: url to which to write or None if no moving should be
                performed
            sam: True iff sam should be output
            samtools_exe: path to samtools exe if necessary; else None

            Yield value: stream
        """
        '''Write SAM header; always include all reference sequences to
        avoid confusing users.'''
        counter.add('stream_and_upload')
        header = '\n'.join(
                ['@HD\tVN:1.0\tSO:coordinate']
                + [('@SQ\tSN:%s\tLN:%d' % (
                        header_rname,
                        reference_index.rname_lengths[header_rname]
                    )) for header_rname in sorted_rnames]
                + ['@PG\tID:Rail-RNA\tVN:%s\tCL:%s %s' % (
                                version.version_number,
                                sys.executable,
                                ' '.join(sys.argv)
                            )]
            )
        if filename is None:
            try:
                print header
                yield sys.stdout
            finally:
                pass
        elif sam:
            try:
                output_stream = open(filename, 'w')
                print >>output_stream, header
                yield output_stream
            finally:
                output_stream.close()
                if not output_url.is_local:
                    mover.put(filename, 
                              output_url.plus(os.path.basename(filename)))
                    os.remove(filename)
        else:
            try:
                samtools_stream = open(filename, 'wb')
                samtools_process = subprocess.Popen(
                            [samtools_exe, 'view', '-bS', '-'],
                            stdin=subprocess.PIPE,
                            stdout=samtools_stream
                        )
                output_stream = samtools_process.stdin
                print >>output_stream, header
                yield output_stream
            finally:
                output_stream.close()
                samtools_return = samtools_process.wait()
                if samtools_return:
                    raise RuntimeError(
                            'samtools returned exit code %d' 
                                % samtools_return
                        )
                samtools_stream.close()
                # Index bam
                unmapped = filename.endswith('.unmapped.bam')
                if not unmapped:
                    subprocess.check_call([samtools_exe, 'index',
                                            filename],
                                            bufsize=-1)
                if not output_url.is_local:
                    mover.put(filename, 
                              output_url.plus(os.path.basename(filename)))
                    if not unmapped:
                        bai = filename + '.bai'
                        mover.put(
                                filename, 
                                output_url.plus(os.path.basename(bai))
                            )
                        os.remove(bai)
                    os.remove(filename)

    counter.flush()
    if args.output_by_chromosome:
        for (index, _), xpartition in xstream(sys.stdin, 2):
            sample_index, rname_index = (
                    sample_and_rname_indexes.sample_and_rname_indexes(index)
                )
            sample_label = manifest_object.index_to_label[sample_index]
            rname = reference_index.string_to_rname[rname_index]
            unique_count, total_count = 0, 0
            with stream_and_upload(
                        sorted_rnames,
                        filename=(
                            os.path.join(
                                    output_dir,
                                    args.bam_basename + '.' + sample_label
                                        + ('.unmapped' if rname == '*'
                                             else ('.' + rname))
                                        + ('.sam' if args.output_sam
                                             else '.bam')
                                ) if args.out is not None else None
                            ),
                        mover=(mover if (args.out is not None
                                         and not output_url.is_local)
                               else None),
                        output_url=(None if args.out is None else output_url),
                        sam=args.output_sam,
                        samtools_exe=args.samtools_exe
                    ) as output_stream:
                for record in xpartition:
                    sam_line_to_print = [record[1][:254], record[2], rname,
                                         str(int(record[0]))] + [
                                            ('ZS:i:' + token[5:]
                                                if token[:5] == 'XS:i:'
                                                else token)
                                                for token in record[3:]]
                    try:
                        print >>output_stream, '\t'.join(sam_line_to_print)
                        counter.add('sam_line')
                    except IOError:
                        raise IOError(
                                'Error writing line "%s".' % sam_line_to_print
                            )
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
            sample_label = manifest_object.index_to_label[sample_index]
            rname = reference_index.string_to_rname[rname_index]
            unique_count, total_count = 0, 0
            with stream_and_upload(
                        sorted_rnames,
                        filename=(
                            os.path.join(
                                    output_dir,
                                    args.bam_basename + '.' + sample_label
                                        + ('.sam' if args.output_sam
                                             else '.bam')
                                ) if args.out is not None else None
                            ),
                        mover=(mover if (args.out is not None
                                         and not output_url.is_local)
                               else None),
                        output_url=(None if args.out is None else output_url),
                        sam=args.output_sam,
                        samtools_exe=args.samtools_exe
                    ) as output_stream:
                for record in xpartition:
                    sam_line_to_print = [record[1][:254], record[2], rname,
                                         str(int(record[0]))] + [
                                            ('ZS:i:' + token[5:]
                                                if token[:5] == 'XS:i:'
                                                else token)
                                                for token in record[3:]]
                    try:
                        print >>output_stream, '\t'.join(sam_line_to_print)
                        counter.add('sam_line')
                    except IOError:
                        raise IOError(
                                'Error writing line "%s".' % sam_line_to_print
                            )
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

print >>sys.stderr, 'DONE with bam.py; in=%d; time=%0.3f s' % (
                                input_line_count, time.time() - start_time
                            )
