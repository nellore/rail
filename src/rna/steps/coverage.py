#!/usr/bin/env python
"""
Rail-RNA-coverage
Follows Rail-RNA-coverage_pre
Precedes Rail-RNA-coverage_post

Reduce step in MapReduce pipelines that outputs normalization factors for
sample coverages. The normalization factor is computed from the histogram of 
base coverage (horizontal axis: number of exonic chunks covering given base;
    vertical axis: number of bases covered) as the (k*100)-th coverage
percentile, where k is input by the user via the command-line parameter
--percentile. bigwig files encoding coverage per sample are also written to a
specified destination, local or remote. Rail-RNA-coverage_post merely collects
the normalization factors and writes them to a file.

Input (read from stdin)
----------------------------
Tab-delimited input tuple columns:
1. Sample label
2. Number string representing reference name (RNAME in SAM format; see 
    BowtieIndexReference class in bowtie_index for conversion information)
3. Position
4. Coverage (that is, the number of called ECs in the sample
    overlapping the position)
Input is partitioned first by sample label, then sorted by fields 2-3.

Hadoop output (written to stdout)
----------------------------
Tab-delimited output tuple columns (only 1 per sample):
1. The character '-', ensuring there's exactly one partition
2. Original sample label
3. Normalization factor

Other output (written to directory specified by command-line parameter --out)
----------------------------
One bigWig file per sample encoding coverage of genome by exonic alignments
"""
import os
import sys
import site
import argparse
import subprocess
import threading

base_path = os.path.abspath(
                    os.path.dirname(os.path.dirname(os.path.dirname(
                        os.path.realpath(__file__)))
                    )
                )
utils_path = os.path.join(base_path, 'rna', 'utils')
site.addsitedir(utils_path)
site.addsitedir(base_path)

import manifest
import bowtie
import bowtie_index
import filemover
import itertools
from collections import defaultdict
from dooplicity.tools import xstream, register_cleanup, make_temp_dir
from dooplicity.ansibles import Url
import tempdel

# Print file's docstring if -h is invoked
parser = argparse.ArgumentParser(description=__doc__, 
            formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument(\
    '--percentile', metavar='FRACTION', type=float, required=False,
    default=0.75,
    help='For a given sample, the per-position percentile to extract as the '
         'normalization factor')
parser.add_argument(\
    '--out', metavar='URL', type=str, required=False, default='.',
    help='URL to which bigwig coverage output should be written. '
         'DEFAULT IS CURRENT WORKING DIRECTORY, NOT STDOUT')
parser.add_argument('--manifest', type=str, required=False,
        default='manifest',
        help='Path to manifest file')
parser.add_argument(\
    '--bigwig-exe', type=str, required=False, default='bedGraphToBigWig',
    help='Location of the Kent Tools bedGraphToBigWig executable')
parser.add_argument('--bigwig-basename', type=str, required=False, default='',
    help='The basename (excluding path) of all bigwig output. Basename is'
         'followed by ".[sample label].bw"; if basename is an empty string, '
         'a sample\'s bigwig filename is simply [sample label].bw')
parser.add_argument(\
    '--keep-alive', action='store_const', const=True, default=False,
    help='Prints reporter:status:alive messages to stderr to keep EMR '
         'task alive')
parser.add_argument(\
    '--verbose', action='store_const', const=True, default=False,
    help='Print out extra debugging statements')

filemover.add_args(parser)
bowtie.add_args(parser)
tempdel.add_args(parser)
args = parser.parse_args()

# Start keep_alive thread immediately
if args.keep_alive:
    from dooplicity.tools import KeepAlive
    keep_alive_thread = KeepAlive(sys.stderr)
    keep_alive_thread.start()

if args.keep_alive:
    class BedTobigwigThread(threading.Thread):
        """ Wrapper class for bedtobigwig that permits polling for completion.
        """
        def __init__(self, command_list):
            super(BedTobigwigThread, self).__init__()
            self.command_list = command_list
            self.bedtobigwig_process = None
        def run(self):
            self.bedtobigwig_process = subprocess.Popen(self.command_list,
                                            stdout=sys.stderr,
                                            stderr=sys.stderr).wait()

def percentile(histogram, percentile=0.75):
    """ Given histogram, computes desired percentile.

        histogram: a dictionary whose keys are integers
            and whose values are frequencies.
        percentile: a value k on [0, 100] specifying that the (k*100)-th
            percentile should be returned

        Return value: Integer key closest to desired percentile.
    """
    covered = 0
    normalization = sum(histogram.values())
    for key, frequency in sorted(histogram.items(), reverse=True):
        covered += frequency
        assert covered <= normalization
        if covered > ((1.0 - percentile) * normalization):
            return key
    raise RuntimeError('Percentile computation should have terminated '
                       'mid-loop.')

import time
start_time = time.time()

temp_dir_path = make_temp_dir(args.scratch)
# Clean up after script
register_cleanup(tempdel.remove_temporary_directories, [temp_dir_path])
bed_filename = os.path.join(temp_dir_path, 'temp.bed')
if args.verbose:
    print >>sys.stderr, 'Writing to temporary bed %s .' % bed_filename
output_filename, output_url = None, None

'''Make RNAME lengths available from reference FASTA so SAM header can be
formed; reference_index.rname_lengths[RNAME] is the length of RNAME.''' 
reference_index = bowtie_index.BowtieIndexReference(args.bowtie_idx)
# For mapping sample indices back to original sample labels
manifest_object = manifest.LabelsAndIndices(args.manifest)
# Create file with chromosome sizes for bedTobigwig
sizes_filename = os.path.join(temp_dir_path, 'chrom.sizes')
if args.verbose:
    print >>sys.stderr, 'Sizes file: %s .' % sizes_filename
with open(sizes_filename, 'w') as sizes_stream:
    for rname in reference_index.rname_lengths:
        print >>sizes_stream, '%s %d' % (rname, 
            reference_index.rname_lengths[rname])

input_line_count, output_line_count = 0, 0
output_url = Url(args.out)
if output_url.is_local:
    # Set up destination directory
    try: os.makedirs(output_url.to_url())
    except: pass
mover = filemover.FileMover(args=args)
for (sample_label,), xpartition in xstream(sys.stdin, 1):
    try:
        sample_label = manifest_object.index_to_label[sample_label]
    except KeyError:
        raise RuntimeError('Sample label index "%s" was not recorded.'
                                % sample_label)
    '''Dictionary for which each key is a coverage (i.e., number of ECs
    covering a given base). Its corresponding value is the number of bases with
    that coverage.'''
    coverage_histogram = defaultdict(int)
    with open(bed_filename, 'w') as bed_stream:
        print >>bed_stream, ('track type=bedGraph name="%s" '
         'description="base coverage by reads" visibility=full '
         'color=227,29,118 altColor=0,179,220 priority=400') % sample_label
        for rname, coverages in itertools.groupby(xpartition, 
                                                    key=lambda val: val[0]):
            try:
                rname = reference_index.string_to_rname[rname]
            except KeyError:
                raise RuntimeError(
                        'RNAME number string "%s" not in Bowtie index.' 
                        % rname
                    )
            last_pos, last_coverage = 0, 0
            for _, pos, coverage in coverages:
                # BED is zero-indexed, while input is 1-indexed
                pos, coverage = int(pos) - 1, int(coverage)
                input_line_count += 1
                print >>bed_stream, '%s\t%d\t%d\t%d' % (rname,
                    last_pos, pos, last_coverage)
                if last_coverage != 0:
                    # Only care about nonzero-coverage regions
                    coverage_histogram[last_coverage] += pos - last_pos
                last_pos, last_coverage = pos, coverage
            if last_pos != reference_index.rname_lengths[rname]:
                # Print coverage up to end of strand
                print >>bed_stream, '%s\t%d\t%d\t%d' % (rname,
                    last_pos, reference_index.rname_lengths[rname], coverage)
    # Output normalization factor
    print '-\t%s\t%d' % (sample_label, percentile(coverage_histogram,
                                                    args.percentile))
    output_line_count += 1
    # Write bigwig
    assert os.path.exists(sizes_filename)
    bigwig_filename = ((args.bigwig_basename + '.') 
        if args.bigwig_basename != '' else '') + sample_label + '.bw'
    if output_url.is_local:
        # Write directly to local destination
        bigwig_file_path = os.path.join(args.out, bigwig_filename)
    else:
        # Write to temporary directory, and later upload to URL
        bigwig_file_path = os.path.join(temp_dir_path, bigwig_filename)
    bigwig_command = [args.bigwig_exe, bed_filename, sizes_filename,
                        bigwig_file_path]
    if args.verbose:
        print >>sys.stderr, 'Writing bigwig with command %s .' \
            % ' '.join(bigwig_command)
    bedtobigwig_process = subprocess.Popen(
                                bigwig_command,
                                stderr=sys.stderr,
                                stdout=sys.stderr,
                                bufsize=-1
                            )
    bedtobigwig_process.wait()
    if bedtobigwig_process.returncode:
        raise RuntimeError('bedgraphtobigwig process failed w/ '
                           'exitlevel %d.'
                            % bedtobigwig_process.returncode)
    if args.verbose:
        print >>sys.stderr, ('bedTobigwig command '
                             + ' '.join([args.bigwig_exe, bed_filename,
                                         sizes_filename, bigwig_file_path])
                             + ' succeeded.')
    if not output_url.is_local:
        # bigwig must be uploaded to URL and deleted
        mover.put(bigwig_file_path, output_url.plus(bigwig_filename))
        os.remove(bigwig_file_path)

print >>sys.stderr, 'DONE with coverage.py; in/out=%d/%d; time=%0.3f s' \
                        % (input_line_count, output_line_count,
                            time.time() - start_time)
