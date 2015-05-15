"""
collect_read_stats.py
Follows Rail-RNA-bam.py
TERMINUS: no steps follow

Reduce step in MapReduce pipelines that collects numbers of primary alignments
and unique alignments overlapping contigs.

Input (read from stdin)
----------------------------
Tab-delimited columns:
1. '-' to enforce that all records are placed in the same partition
2. Sample index
3. RNAME index
4. Number of primary alignments overlapping contig with RNAME
5. Number of unique alignments overlapping contig with RNAME

Input is partitioned by the first field and sorted by the second one.

Hadoop output (written to stdout)
----------------------------
None.

Other output (written to directory specified by command-line parameter --out)
----------------------------
Matrix stored in TSV whose (i, j)th element is the comma-delimited pair
(X, Y), where X is the number of primary alignments to contig j and Y
is the number of unique alignments to contig j. X = number of unmapped reads
and Y = 0 if the RNAME is "*". When contig j is "mapped totals", X = total
number of mapped reads, and Y = total number of uniquely mapped reads.
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

from dooplicity.ansibles import Url
from dooplicity.tools import xstream, register_cleanup, make_temp_dir, xopen
import bowtie
import bowtie_index
import manifest
import filemover
import tempdel
from collections import defaultdict

# Print file's docstring if -h is invoked
parser = argparse.ArgumentParser(description=__doc__, 
            formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument(
    '--out', metavar='URL', type=str, required=False, default=None,
    help='URL to which BED output should be written. Default is current '
         'working directory.')
parser.add_argument('--manifest', type=str, required=False,
        default='manifest',
        help='Path to manifest file')
parser.add_argument(
    '--tsv-basename', type=str, required=False, 
    default='',
    help='The basename (excluding path) of TSV output. Basename is '
         'followed by ".counts.tsv.gz"')
parser.add_argument('--gzip-level', type=int, required=False,
        default=3,
        help='Level of gzip compression to use for temporary files')
parser.add_argument(
    '--keep-alive', action='store_const', const=True, default=False,
    help='Prints reporter:status:alive messages to stderr to keep EMR '
         'task alive')

bowtie.add_args(parser)
filemover.add_args(parser)
tempdel.add_args(parser)
args = parser.parse_args()

# Start keep_alive thread immediately
if args.keep_alive:
    from dooplicity.tools import KeepAlive
    keep_alive_thread = KeepAlive(sys.stderr)
    keep_alive_thread.start()

import time
start_time = time.time()

reference_index = bowtie_index.BowtieIndexReference(args.bowtie_idx)
# For mapping sample indices back to original sample labels
manifest_object = manifest.LabelsAndIndices(args.manifest)
output_url = Url(args.out) if args.out is not None \
    else Url(os.getcwd())
input_line_count = 0
if output_url.is_local:
    # Set up destination directory
    try: os.makedirs(output_url.to_url())
    except: pass
else:
    mover = filemover.FileMover(args=args)
    # Set up temporary destination
    import tempfile
    temp_dir_path = make_temp_dir(args.scratch)
    register_cleanup(tempdel.remove_temporary_directories, [temp_dir_path])

output_filename = ((args.tsv_basename + '.' 
                          if args.tsv_basename != '' else '')
                          + type_string + '.tsv.gz')
if output_url.is_local:
    output_path = os.path.join(args.out, output_filename)
else:
    output_path = os.path.join(temp_dir_path, output_filename)

input_line_count = 0
# Get RNAMEs in order of descending length
sorted_rnames = [reference_index.string_to_rname['%012d' % i]
                    for i in xrange(
                                len(reference_index.string_to_rname) - 1
                            )]
with xopen(True, output_path, 'w', args.gzip_level) as output_stream:
    print >>output_stream, '\t'.join([''] + sorted_rnames + 'mapped totals')
    for (_, sample_index), xpartition in xstream(sys.stdin, 2):
        sample_label = manifest_object.index_to_label[sample_index]
        total_counts, unique_counts = defaultdict(int), defaultdict(int)
        for rname_index, total_count, unique_count in xpartition:
            rname = reference_index.string_to_rname[rname_index]
            total_counts[rname] = total_count
            unique_counts[rname] = unique_count
        total_mapped_reads = sum(total_counts.values()) - total_counts['*']
        total_unique_alignments = (
                    sum(unique_counts.values()) - unique_counts['*']
                )
        print >>output_stream, '\t'.join(
                [sample_label] + ['%d,%d' % (total_counts[rname],
                                                unique_counts[rname])]
                    + ['%d,%d' % (total_mapped_reads, total_unique_alignments)]
            )

if not output_url.is_local:
    mover.put(output_path, output_url.plus(output_filename))
    os.remove(output_path)

print >>sys.stderr, 'DONE with collect_read_stats.py; in=%d; time=%0.3f s' % (
        input_line_count, time.time() - start_time
    )