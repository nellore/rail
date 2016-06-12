#!/usr/bin/env python
"""
Rail-RNA-tsv
Follows Rail-RNA-bed_pre / Rail-RNA-coverage
TERMINUS: no steps follow.

Reducer for MapReduce pipelines that writes cross-sample tables storing
coverages of junctions, insertions, and deletions across samples.

Input (read from stdin)
----------------------------
Tab-delimited output tuple columns (collect)
1. '0' if insertion, '1' if deletion, '2' if junction line, '3'
    if normalization factor
2. Number string representing RNAME (+ '+ or -' if junction; same as field 6)
    or sample index if field 1 is '3'
3. Start position (Last base before insertion, first base of deletion,
                    or first base of intron)
    or '\x1c' if field 1 is '3'
4. End position (Last base before insertion, last base of deletion (exclusive),
                    or last base of intron (INCLUSIVE))
    or '\x1c' if field 1 is '3'
5. '+' or '-' indicating which strand is the sense strand for junctions,
   inserted sequence for insertions, or deleted sequence for deletions
   or '\x1c' if field 1 is '3'
6. Coverage of feature for sample with index N
    or normalization factor
...
N + 6. Coverage of feature in sample with index N

Input is partitioned by field 1 and sorted by fields 2-5.

Hadoop output (written to stdout)
----------------------------
None.

Other output (written to directory specified by command-line parameter --out)
----------------------------
1) Three coverage matrices, each stored in a different TSV file: one is for
junctions, one is for insertions, and one is for deletions. A given element
(i, j) of a given matrix specifies the number of reads in which the feature
(junction, insertion, or deletion) i was found in sample j.

A feature i is specified as
RNAME + ';' + strand for junctions or inserted sequence for insertions or
deleted sequence for deletions + ';' + start position (last base before
insertion, first base of deletion, or first base of intron) + ';' +
end position (last base before insertion, last base of deletion (exclusive), or
last base of intron (exclusive))

2) Normalization factors for sample read coverage distributions
Tab-delimited tuple columns:
1. Sample name
2. Normalization factor for coverage vector counting all primary alignments
3. Normalization factor for coverage vector counting only unique alignments
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

from dooplicity.ansibles import Url
from dooplicity.tools import xstream, register_cleanup, make_temp_dir, xopen
import bowtie
import bowtie_index
import manifest
import filemover
import tempdel

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
         'followed by ".[junctions/insertions/deletions].tsv.gz"')
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

reference_index = bowtie_index.BowtieIndexReference(
                        os.path.expandvars(args.bowtie_idx)
                    )
# For mapping sample indices back to original sample labels
manifest_object = manifest.LabelsAndIndices(
                        os.path.expandvars(args.manifest)
                    )
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
    temp_dir_path = make_temp_dir(tempdel.silentexpandvars(args.scratch))
    register_cleanup(tempdel.remove_temporary_directories, [temp_dir_path])

input_line_count = 0
for (line_type,), xpartition in xstream(sys.stdin, 1):
    type_string = ('insertions' if line_type == '0' else
                    ('deletions' if line_type == '1' else
                      ('junctions' if line_type == '2' else
                        'normalization')))
    output_filename = ((args.tsv_basename + '.' 
                          if args.tsv_basename != '' else '')
                          + type_string + '.tsv.gz')
    if output_url.is_local:
        output_path = os.path.join(args.out, output_filename)
    else:
        output_path = os.path.join(temp_dir_path, output_filename)
    with xopen(True, output_path, 'w', args.gzip_level) as output_stream:
        if line_type != '3':
            '''Print all labels in the order in which they appear in the
            manifest file.'''
            sample_count = len(manifest_object.index_to_label)
            for i in xrange(sample_count):
                output_stream.write(
                        '\t' + manifest_object.index_to_label[str(i)]
                    )
            output_stream.write('\n')
            for coverage_line in xpartition:
                input_line_count += 1
                (rname, pos, end_pos, strand_or_seq) = coverage_line[:4]
                '''Handle missing zeros at end of line here; in previous step,
                the total number of samples was unknown, so this was not
                done.'''
                print >>output_stream, '\t'.join(
                    (';'.join(
                        [reference_index.string_to_rname[rname],
                            strand_or_seq, str(int(pos)), str(int(end_pos))]
                    ),) + coverage_line[4:] + ('0',)*(-len(coverage_line) + 4
                                                        + sample_count)
                )
        else:
            for (sample_index, _, _, _, factor, unique_factor) in xpartition:
                input_line_count += 1
                print >>output_stream, '\t'.join([
                        manifest_object.index_to_label[sample_index],
                        factor, unique_factor
                    ])

    if not output_url.is_local:
        mover.put(output_path, output_url.plus(output_filename))
        os.remove(output_path)

print >>sys.stderr, 'DONE with tsv.py; in=%d; time=%0.3f s' \
                        % (input_line_count, time.time() - start_time)