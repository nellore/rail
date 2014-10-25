"""
Rail-RNA-bed
Follows Rail-RNA-bed_pre
TERMINUS: no steps follow.

Reduce step in MapReduce pipelines that turns output of Rail-RNA-bed_pre
into TopHat-like indel/intron BED output.

Input (read from stdin)
----------------------------
Tab-delimited input tuple columns:
1. 'I', 'D', or 'N' for insertion, deletion, or intron line
2. Sample label
3. RNAME
4. Start position (Last base before insertion, first base of deletion,
                    or first base of intron)
5. End position (Last base before insertion, last base of deletion (exclusive),
                    or last base of intron (exclusive))
6. '+' or '-' indicating which strand is the sense strand for introns,
   inserted sequence for insertions, or deleted sequence for deletions
----Next fields are for introns only; they are '\x1c' for indels----
7. MAX number of nucleotides between 5' end of intron and 5' end of read from
    which it was inferred, ASSUMING THE SENSE STRAND IS THE FORWARD STRAND.
    That is, if the sense strand is the reverse strand, this is the distance
    between the 3' end of the read and the 3' end of the intron.
8. MAX number of nucleotides between 3' end of intron and 3' end of read from
    which it was inferred, ASSUMING THE SENSE STRAND IS THE FORWARD STRAND.
--------------------------------------------------------------------
9. Number of instances of intron, insertion, or deletion in sample; this is
    always +1 before bed_pre combiner/reducer

Input is partitioned by fields 1-2 and sorted by fields 3-5

Hadoop output (written to stdout)
----------------------------
None.

Other output (written to directory specified by command-line parameter --out)
----------------------------
TopHat-like junctions.bed, insertions.bed, and deletions.bed.
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
from dooplicity.tools import xstream
import bowtie
import bowtie_index
import manifest
import atexit
# Define string version_number
from version import version_number
import filemover

# Print file's docstring if -h is invoked
parser = argparse.ArgumentParser(description=__doc__, 
            formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument(\
    '--out', metavar='URL', type=str, required=False, default=None,
    help='URL to which BED output should be written. Default is current '
         'working directory.')
parser.add_argument('--manifest', type=str, required=False,
        default='manifest',
        help='Path to manifest file')
parser.add_argument(\
    '--bed-basename', type=str, required=False, 
    default='',
    help='The basename (excluding path) of all BED output. Basename is '
         'followed by ".[junctions/insertions/deletions].[sample_label].bed"')

bowtie.add_args(parser)
filemover.add_args(parser)
args = parser.parse_args()

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
    temp_dir_path = tempfile.mkdtemp()
    from tempdel import remove_temporary_directories
    atexit.register(remove_temporary_directories, [temp_dir_path])
for (line_type, sample_label), xpartition in xstream(sys.stdin, 2):
    assert line_type in 'NID'
    sample_label = manifest_object.index_to_label[sample_label]
    type_string = ('insertions' if line_type == 'I' else
                    ('deletions' if line_type == 'D' else 'junctions'))
    output_filename = ((args.bed_basename + '.' 
                          if args.bed_basename != '' else '')
                          + type_string + '.' + sample_label + '.bed')
    if output_url.is_local:
        output_path = os.path.join(args.out, output_filename)
    else:
        output_path = os.path.join(temp_dir_path, output_filename)
    with open(output_path, 'w') as output_stream:
        print >>output_stream, 'track name="%s_%s" description="' \
                                   'Rail-RNA v%s %s for sample %s"' \
                                                      % (sample_label,
                                                         type_string,
                                                         version_number,
                                                         type_string,
                                                         sample_label)
        if line_type == 'N':
            for i, (rname, pos, end_pos, reverse_strand_string,
                    max_left_overhang, max_right_overhang, 
                    maximin_overhang, coverage) \
                in enumerate(xpartition):
                pos, end_pos = int(pos) - 1, int(end_pos) - 1
                max_left_overhang, max_right_overhang \
                    = int(max_left_overhang), int(max_right_overhang)
                start_position = pos - max_left_overhang
                end_position = end_pos + max_right_overhang
                print >>output_stream, '%s\t%d\t%d\tJUNC%08d' \
                                       ';maximin_overhang=%s\t' \
                                       '%s\t%s\t%d\t%d\t227,29,118\t2\t' \
                                       '%d,%d\t0,%d' % (
                                            reference_index.string_to_rname[
                                                    rname
                                                ],
                                            start_position,
                                            end_position, i+1,
                                            maximin_overhang,
                                            coverage,
                                            reverse_strand_string,
                                            start_position, end_position,
                                            max_left_overhang,
                                            max_right_overhang,
                                            max_left_overhang + end_pos - pos
                                        )
            input_line_count += i
        else:
            for i, (rname, pos, end_pos, seq, _, _, _, coverage) \
                in enumerate(xpartition):
                pos, end_pos = int(pos) - 1, int(end_pos) - 1
                print >>output_stream, '%s\t%d\t%d\t%s\t%s' \
                                        % (reference_index.string_to_rname[
                                                rname
                                            ], pos, end_pos, seq, coverage)
            input_line_count += i
    if not output_url.is_local:
        mover.put(output_path, output_url.plus(output_filename))
        os.remove(output_path)

print >>sys.stderr, 'DONE with bed.py; in=%d; time=%0.3f s' \
                        % (input_line_count, time.time() - start_time)
