"""
Rail-RNA-coverage
Follows Rail-RNA-coverage_pre
Precedes Rail-RNA-coverage_post

Reduce step in MapReduce pipelines that outputs normalization factors for
sample coverages. The normalization factor is computed from the histogram of 
base coverage (horizontal axis: number of exonic chunks covering given base;
    vertical axis: number of bases covered) as the (k*100)-th coverage
percentile, where k is input by the user via the command-line parameter
--percentile. bigBed files encoding coverage per sample are also written to a
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
One bigBed file per sample encoding coverage of genome by exonic alignments
"""
import os
import sys
import site
import argparse
import subprocess

base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for directory_name in ['util', 'fasta', 'bowtie', 'manifest']:
    site.addsitedir(os.path.join(base_path, directory_name))

import manifest
import bowtie
import bowtie_index
import url
import path
import filemover
import fasta

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
    help='URL to which bigBed coverage output should be written. '
         'DEFAULT IS CURRENT WORKING DIRECTORY, NOT STDOUT')
parser.add_argument('--manifest', type=str, required=False,
        default='manifest',
        help='Path to manifest file')
parser.add_argument(\
    '--bigbed-exe', type=str, required=False, default='bedToBigBed',
    help='Location of the Kent Tools bedToBigBed executable')
parser.add_argument('--bigbed-basename', type=str, required=False, default='',
    help='The basename (excluding path) of all bigBed output. Basename is'
         'followed by ".[sample label].bb"; if basename is an empty string, '
         'a sample\'s bigBed filename is simply [sample label].bb')
parser.add_argument(\
    '--verbose', action='store_const', const=True, default=False,
    help='Print out extra debugging statements')

filemover.addArgs(parser)
bowtie.addArgs(parser)
args = parser.parse_args()

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

# For storing BED files before conversion to bigBed
import tempfile
temp_dir_path = tempfile.mkdtemp()
bed_filename = os.path.join(temp_dir_path, 'temp.bed')
output_filename, output_url = None, None

'''Make RNAME lengths available from reference FASTA so SAM header can be
formed; reference_index.rname_lengths[RNAME] is the length of RNAME.''' 
reference_index = bowtie_index.BowtieIndexReference(args.bowtie_idx)
# For mapping sample indices back to original sample labels
manifest_object = manifest.LabelsAndIndices(args.manifest)
# Create file with chromosome sizes for bedToBigBed
sizes_filename = os.path.join(temp_dir_path, 'chrom.sizes')
with open(sizes_filename, 'w') as sizes_stream:
    for rname in reference_index.rname_lengths:
        print >>sizes_stream, '%s %d' % (rname, 
            reference_index.rname_lengths[rname])

input_line_count, output_line_count = 0, 0
last_sample_label, last_rname, last_pos, last_coverage = [None]*4
'''Dictionary for which each key is a coverage (i.e., number of ECs covering
a given base). Its corresponding value is the number of bases with that
coverage.'''
coverage_histogram = {}
output_url = url.Url(args.out)
if output_url.isLocal():
    # Set up destination directory
    try: os.makedirs(output_url.toUrl())
    except: pass
bed_stream = open(bed_filename, 'w')
mover = filemover.FileMover(args=args)
while True:
    line = sys.stdin.readline()
    if line:
        input_line_count += 1
        tokens = line.rstrip().split()
        assert len(tokens) == 4, 'Bad input line:\n' + line
        sample_label, rname, pos, coverage = (tokens[0], tokens[1],
                                                int(tokens[2]), int(tokens[3]))
        assert rname in reference_index.string_to_rname, \
            'RNAME number string "%s" not in Bowtie index.' % rname
        rname = reference_index.string_to_rname[rname]
        sample_label = manifest_object.index_to_label[sample_label]
    if (not line or sample_label != last_sample_label) \
        and last_sample_label is not None:
        # All of a sample's coverage entries have been read
        if last_coverage != 0 \
            and last_pos < reference_index.rname_lengths[last_rname]:
            # Output final coverage entry for sample
            coverage_start_pos, coverage_end_pos = (last_pos - 1,
                reference_index.rname_lengths[last_rname])
            print >>bed_stream, '%s\t%d\t%d\t%d' % (last_rname,
                coverage_start_pos, coverage_end_pos, last_coverage)
            coverage_histogram[last_coverage] = \
                coverage_histogram.get(last_coverage, 0) + coverage_end_pos \
                - coverage_start_pos
        print '-\t%s\t%d' % (last_sample_label, percentile(coverage_histogram,
                                                              args.percentile))
        output_line_count += 1
        coverage_histogram = {}
        bed_stream.close()
        # Write bigBed
        assert os.path.exists(sizes_filename)
        assert path.is_exe(args.bigbed_exe)
        bigbed_filename = ((args.bigbed_basename + '.') 
            if args.bigbed_basename != '' else '') + last_sample_label + '.bb'
        if output_url.isLocal():
            # Write directly to local destination
            bigbed_file_path = os.path.join(args.out, bigbed_filename)
        else:
            # Write to temporary directory, and later upload to URL
            bigbed_file_path = os.path.join(temp_dir_path, bigbed_filename)
        bedtobigbed_command = ' '.join([args.bigbed_exe, bed_filename,
            sizes_filename, bigbed_file_path])
        bedtobigbed_process = subprocess.Popen(bedtobigbed_command, shell=True,
            bufsize=-1, stdout=(sys.stderr if args.verbose else os.devnull))
        bigbed_return = bedtobigbed_process.wait()
        if bigbed_return:
            raise RuntimeError('bedToBigBed command ' + bedtobigbed_command
                + (' returned with exitlevel %d' % bigbed_return))
        if args.verbose:
            print >>sys.stderr, ('bedToBigBed command ' + bedtobigbed_command
                + ' succeeded.' )
        if not output_url.isLocal():
            # bigBed must be uploaded to URL and deleted
            mover.put(bigbed_file_path, output_url.plus(bigbed_filename))
            os.remove(bigbed_file_path)
        bed_stream = open(bed_filename, 'w')
    elif last_sample_label is not None and sample_label == last_sample_label:
        if last_rname == rname:
            if last_coverage != 0:
                '''Add to histogram only if coverage > 0 to minimize
                dictionary size.'''
                coverage_histogram[last_coverage] \
                    = coverage_histogram.get(last_coverage, 0) \
                    + (pos - last_pos)
            if coverage != last_coverage and last_coverage != 0:
                print >>bed_stream, '%s\t%d\t%d\t%d' % (last_rname,
                    last_pos - 1, pos - 1, last_coverage)
            elif last_coverage != 0:
                assert coverage == last_coverage
                # So next output interval extends back to previous pos
                pos = last_pos
        elif last_coverage != 0 \
            and last_pos < reference_index.rname_lengths[last_rname]:
            # Output final coverage entry for RNAME
            coverage_start_pos, coverage_end_pos = (last_pos - 1,
                reference_index.rname_lengths[last_rname])
            print >>bed_stream, '%s\t%d\t%d\t%d' % (last_rname,
                coverage_start_pos, coverage_end_pos, last_coverage)
            coverage_histogram[last_coverage] = \
                coverage_histogram.get(last_coverage, 0) + coverage_end_pos \
                - coverage_start_pos
    if not line: break
    last_sample_label, last_rname, last_pos, last_coverage = (sample_label,
        rname, pos, coverage)

bed_stream.close()

if not output_url.isLocal():
    # Clean up
    import shutil
    shutil.rmtree(temp_dir_path)

print >>sys.stderr, 'DONE with coverage.py; in=%d; time=%0.3f s' \
                        % (input_line_count, time.time() - start_time)
