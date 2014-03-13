#!/usr/bin/env python
"""
Rail-RNA-intron_post
Follows Rail-RNA-intron
Precedes Rail-RNA-realign

Reduce step in MapReduce pipelines that builds a new Bowtie index from only
the introns called in Rail-RNA-intron. Reference names are designed for
simple reconstruction of intron positions from alignments.

Input (read from stdin)
----------------------------
Two kinds of input lines, one from Rail-RNA-align (max_len) and the other
from Rail-RNA-intron (intron):

(max_len)
1. The character 'a'.
2. A maximum read length output by a Rail-RNA-align mapper.
3. The character '0'.
4. The character '0'.

(intron)
1. The character 'i', which will place the row after all 'a's (maximum read
    length rows) in lexicographic sort order.
2. Reference name (RNAME in SAM format) + 
    '+' or '-' indicating which strand is the sense strand
3. Intron start position (inclusive)
4. Intron end position (exclusive)

Input is sorted by all four columns.

Hadoop output (written to stdout)
----------------------------
None.

Other output (written to directory specified by command-line parameter --out)
----------------------------
Bowtie index files for realignment only to regions framing introns of kept
unmapped reads from Rail-RNA-align.

A given reference name in the index is in the following format:    
    original RNAME + '+' or '-' indicating which strand is the sense
    strand + ';' + start position of sequence + ';' + comma-separated
    list of subsequence sizes framing introns + ';' + comma-separated
    list of intron sizes.
"""
import os
import sys
import site
import subprocess
import argparse
import atexit

base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for directory_name in ['bowtie', 'util']:
    site.addsitedir(os.path.join(base_path, directory_name))

import bowtie
import bowtie_index
import url
import filemover

# Print file's docstring if -h is invoked
parser = argparse.ArgumentParser(description=__doc__, 
            formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument(\
    '--out', metavar='URL', type=str, required=False,
    default='.',
    help='Bowtie index files are written to this URL. DEFAULT IS CURRENT '
         ' WORKING DIRECTORY.')
parser.add_argument('--fudge', type=int, required=False, default=1, 
    help='A splice junction may be detected at any position along the read '
         'besides directly before or after it; thus, the sequences recorded '
         'in the new index should include max_read_size - 1 bases before and '
         'after an intron in reference space (on the forward strand). These '
         'extensions are extended further by the number of bases --fudge to '
         'accommodate possible small insertions.')

filemover.addArgs(parser)
bowtie.addArgs(parser)
args = parser.parse_args()

def handle_temporary_directory(temp_dir_path):
    """ Deletes temporary directory.

        temp_dir_path: path of temporary directory for storing intermediate
            alignments; archived if archive is not None.

        No return value.
    """
    # Kill temporary directory
    import shutil
    shutil.rmtree(temp_dir_path)

import time
start_time = time.time()
_input_line_count = 0

output_filename, output_stream, output_url = [None]*3
reference_index = bowtie_index.BowtieIndexReference(args.bowtie_idx)
output_url = url.Url(args.out)
# Set up temporary destination
import tempfile
temp_dir_path = tempfile.mkdtemp()
# For deleting temporary directory, even on unexpected exit
atexit.register(handle_temporary_directory, temp_dir_path)
if output_url.isLocal():
    # Set up final destination
    try: os.makedirs(output_url.toUrl())
    except: pass
else:
    # Set up temporary destination
    try: os.makedirs(os.path.join(temp_dir_path, 'index'))
    except: pass
if output_url.isLocal():
    # Write directly to local destination
    index_basename = os.path.join(output_url.toUrl(), 'intron')
else:
    # Write to temporary directory, and later upload to URL
    index_basename = os.path.join(temp_dir_path, 'index/intron')
last_line, last_line_type, max_read_size, last_rname = [None]*4
write_sequence = False
fasta_file = os.path.join(output_url.toUrl(), 'temp.fa')
with open(fasta_file, 'w') as fasta_stream:
    while True:
        line = sys.stdin.readline().rstrip()
        if line:
            _input_line_count += 1
            if line == last_line: continue
            tokens = line.rstrip().split('\t')
            token_count = len(tokens)
            assert token_count == 4
            if last_line_type is not None and tokens[0] != last_line_type:
                extend_size = last_max_read_size - 1 + args.fudge
            elif tokens[0] == 'a': 
                line_type, max_read_size = tokens[0], int(tokens[1])
            if tokens[0] == 'i':
                rname, pos, end_pos = tokens[1:]
                # Grab reverse_strand_string attached to rname
                reverse_strand_string = rname[-1]
                rname = rname[:-1]
                pos, end_pos = int(pos), int(end_pos)
                if last_rname is None or not (
                        last_rname == rname 
                        and last_reverse_strand_string == reverse_strand_string
                        and pos - last_end_pos < extend_size
                    ):
                    if last_rname is not None: write_sequence = True
                    intron_combos = set([frozenset([(pos, end_pos)])])
                else:
                    new_intron_combos = set()
                    for intron_combo in intron_combos:
                        new_intron_combo = []
                        add_old_combo = False
                        for intron_pos, intron_end_pos in intron_combo:
                            if min(intron_end_pos, end_pos) \
                                - max(intron_pos, pos) <= 0:
                                '''If there's no overlap between the current
                                intron and an intron in an existing
                                configuration, keep that intron when forming
                                a new configuration. Otherwise, toss it.'''
                                new_intron_combo.append(
                                        (intron_pos, intron_end_pos)
                                    )
                            else: 
                                '''If there's any overlap between an existing
                                intron configuration and the current intron,
                                that configuration should be preserved as
                                independent.'''
                                add_old_combo = True
                        if add_old_combo:
                            new_intron_combos.add(intron_combo)
                        new_intron_combos.add(
                                frozenset(new_intron_combo + [(pos, end_pos)])
                            )
                    intron_combos = new_intron_combos
        else: write_sequence = True
        if write_sequence:
            reference_length = reference_index.rname_lengths[last_rname]
            '''Some intron combos may still contain introns separated by a
            distance > extend_size. Correct this here.''' 
            new_intron_combos = set()
            for intron_combo in last_intron_combos:
                intron_combo_list = sorted(list(intron_combo))
                start_index = 0
                for i in xrange(1, len(intron_combo_list)):
                    if intron_combo_list[i][0] - intron_combo_list[i-1][1] \
                        >= extend_size:
                        # Form new intron combo
                        new_intron_combos.add(
                                frozenset(intron_combo_list[start_index:i])
                            )
                        start_index = i
                # Add last intron combo
                new_intron_combos.add(
                        frozenset(intron_combo_list[start_index:])
                    )
            for intron_combo in new_intron_combos:
                subseqs = []
                intron_combo_list = sorted(list(intron_combo))
                left_start = max(intron_combo_list[0][0] - extend_size, 1)
                # Add sequence before first intron
                subseqs.append(
                        reference_index.get_stretch(last_rname, left_start - 1, 
                            intron_combo_list[0][0] - left_start)
                    )
                # Add sequences between introns
                for i in xrange(1, len(intron_combo_list)):
                    subseqs.append(
                            reference_index.get_stretch(last_rname, 
                                intron_combo_list[i-1][1] - 1,
                                intron_combo_list[i][0]
                                - intron_combo_list[i-1][1]
                            )
                        )
                # Add final sequence
                subseqs.append(
                        reference_index.get_stretch(last_rname,
                            intron_combo_list[-1][1] - 1,
                            min(extend_size, reference_length - 
                                                intron_combo_list[-1][1] + 1))
                    )
                '''A given reference name in the index will be in the following
                format:
                original RNAME + '+' or '-' indicating which strand is the
                sense strand + ';' + start position of sequence + ';' +
                comma-separated list of subsequence sizes framing introns + ';'
                + comma-separated list of intron sizes.'''
                print >>fasta_stream, ('>' + last_rname 
                    + last_reverse_strand_string + ';' + str(left_start) + ';'
                    + ','.join([str(len(subseq)) for subseq in subseqs]) + ';'
                    + ','.join([str(intron_end_pos - intron_pos)
                                    for intron_pos, intron_end_pos
                                    in intron_combo_list])
                )
                full_seq = ''.join(subseqs)
                fasta_stream.write(
                        '\n'.join([full_seq[i:i+80] for i 
                                    in xrange(0, len(full_seq), 80)])
                    )
                fasta_stream.write('\n')
            write_sequence = False
        if not line: break
        if tokens[0] == 'a':
            last_line_type, last_max_read_size = line_type, max_read_size
        else:
            # Line corresponds to an intron
            last_rname, last_reverse_strand_string, last_pos, last_end_pos \
                = rname, reverse_strand_string, pos, end_pos
            last_intron_combos = intron_combos
            last_line = line
# Build index
bowtie_build_process = subprocess.Popen(
                                [args.bowtie_build_exe,
                                    fasta_file,
                                    index_basename],
                                stderr=subprocess.STDOUT,
                                stdout=subprocess.PIPE
                            )
bowtie_build_process.wait()
sys.stderr.write(bowtie_build_process.stdout.read())
sys.stderr.flush()
if bowtie_build_process.returncode:
    RuntimeError('Bowtie index construction failed.')
if not output_url.isLocal():
    # Upload index files and clean up
    mover = filemover.FileMover(args=args)
    for index_file in os.listdir(os.path.join(temp_dir_path, 'index')):
        index_file_path = os.path.join(temp_dir_path, 'index', index_file)
        mover.put(index_file_path, output_url.plus(index_file))
        os.remove(index_file_path)

print >>sys.stderr, 'DONE with intron_post.py; in = %d; time=%0.3f s' \
                        % (_input_line_count, time.time() - start_time)