#!/usr/bin/env python
"""
Rail-RNA-intron_index
Follows Rail-RNA-intron_fasta
Precedes Rail-RNA-realign

Reduce step in MapReduce pipelines that builds a new Bowtie index from the
FASTA output by RailRNA-intron_fasta.

Input (read from stdin)
----------------------------
Tab-delimited tuple columns:
1. '-' to enforce that all lines end up in the same partition
2. FASTA reference name including '>'. The following format is used:
    original RNAME + '+' or '-' indicating which strand is the sense strand
    + ';' + start position of sequence + ';' + comma-separated list of
    subsequence sizes framing introns + ';' + comma-separated list of intron
    sizes
3. Sequence

Input is partitioned by the first column and needs no sort.

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
import tarfile
import atexit
import threading

base_path = os.path.abspath(
                    os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
                )
utils_path = os.path.join(base_path, 'rna/utils')
dp_path = os.path.join(base_path, 'dooplicity')
site.addsitedir(utils_path)
site.addsitedir(dp_path)

import bowtie
import dooplicity as dp
import filemover

# Print file's docstring if -h is invoked
parser = argparse.ArgumentParser(description=__doc__, 
            formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument(\
    '--out', metavar='URL', type=str, required=False,
    default='None',
    help='Bowtie index files are written to this URL. DEFAULT IS CURRENT '
         'WORKING DIRECTORY.')
parser.add_argument(\
    '--basename', type=str, required=False,
    default='intron',
    help='Basename for index to be written')
parser.add_argument(\
    '--keep-alive', action='store_const', const=True, default=False,
    help='Prints reporter:status:alive messages to stderr to keep EMR '
         'task alive')

filemover.add_args(parser)
bowtie.add_args(parser)
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

output_filename, output_stream, output_url = [None]*3
output_url = dp.Url(args.out) if args.out is not None \
    else dp.Url(os.getcwd())
# Set up temporary destination
import tempfile
temp_dir_path = tempfile.mkdtemp()
# For deleting temporary directory, even on unexpected exit
atexit.register(handle_temporary_directory, temp_dir_path)
if output_url.is_local():
    # Set up final destination
    try: os.makedirs(output_url.to_url())
    except: pass
else:
    # Set up temporary destination
    try: os.makedirs(os.path.join(temp_dir_path, 'index'))
    except: pass
if output_url.is_local():
    # Write directly to local destination
    index_basename = os.path.join(output_url.to_url(), args.basename)
else:
    # Write to temporary directory, and later upload to URL
    index_basename = os.path.join(temp_dir_path, 'index/' + args.basename)
fasta_file = os.path.join(temp_dir_path, 'temp.fa')
print >>sys.stderr, 'Opened %s for writing....' % fasta_file
with open(fasta_file, 'w') as fasta_stream:
    for input_line_count, line in enumerate(sys.stdin):
        if args.keep_alive: print >>sys.stderr, 'reporter:status:alive'
        tokens = line.rstrip().split('\t')
        assert len(tokens) == 3
        rname, seq = tokens[1:]
        '''A given reference name in the index will be in the following
        format:
        original RNAME + '+' or '-' indicating which strand is the
        sense strand + ';' + start position of sequence + ';' +
        comma-separated list of subsequence sizes framing introns + ';'
        + comma-separated list of intron sizes.'''
        print >>fasta_stream, rname
        fasta_stream.write(
                '\n'.join([seq[i:i+80] for i 
                            in xrange(0, len(seq), 80)])
            )
        fasta_stream.write('\n')

# Build index
print >>sys.stderr, 'Running bowtie-build....'

if args.keep_alive:
    class BowtieBuildThread(threading.Thread):
        """ Wrapper class for bowtie-build that permits polling for completion.
        """
        def __init__(self, command_list):
            super(BowtieBuildThread, self).__init__()
            self.command_list = command_list
            self.bowtie_build_process = None
        def run(self):
            self.bowtie_build_process = subprocess.Popen(self.command_list,
                                            stdout=sys.stderr).wait()
    bowtie_build_thread = BowtieBuildThread([args.bowtie_build_exe,
                                                fasta_file,
                                                index_basename])
    bowtie_build_thread.start()
    while bowtie_build_thread.is_alive():
        print >>sys.stderr, 'reporter:status:alive'
        sys.stderr.flush()
        time.sleep(5)
    if bowtie_build_thread.bowtie_build_process:
        raise RuntimeError('Bowtie index construction failed w/ exitlevel %d.'
                                % bowtie_build_thread.bowtie_build_process)
else:
    bowtie_build_process = subprocess.Popen(
                                [args.bowtie_build_exe,
                                    fasta_file,
                                    index_basename],
                                stderr=sys.stderr,
                                stdout=sys.stderr
                            )
    bowtie_build_process.wait()
    if bowtie_build_process.returncode:
        raise RuntimeError('Bowtie index construction failed w/ exitlevel %d.'
                                % bowtie_build_process.returncode)

if not output_url.is_local():
    # Compress index files
    intron_index_filename = args.basename + '.tar.gz'
    intron_index_path = os.path.join(temp_dir_path, intron_index_filename)
    index_path = os.path.join(temp_dir_path, 'index')
    tar = tarfile.TarFile.gzopen(intron_index_path, mode='w', compresslevel=3)
    for index_file in os.listdir(index_path):
        tar.add(os.path.join(index_path, index_file), arcname=index_file)
    tar.close()
    # Upload compressed index
    mover = filemover.FileMover(args=args)
    mover.put(intron_index_path, output_url.plus(intron_index_filename))

print >>sys.stderr, 'DONE with intron_index.py; in=%d; time=%0.3f s' \
                        % (input_line_count, time.time() - start_time)
