#!/usr/bin/env python
"""
Rail-RNA-cojunction_enum

Follows Rail-RNA-junction_index
Precedes Rail-RNA-cojunction_fasta

Alignment script for MapReduce pipelines that wraps Bowtie 2. Finds junctions
that cooccur on reads by local alignments to transcriptome elements.

Input (read from stdin)
----------------------------
Single input tuple column:
1. SEQ or its reversed complement -- must be unique
    (but not necessarily in alphabetical order)

Hadoop output (written to stdout)
----------------------------
Tab-delimited tuple columns:
1. Reference name (RNAME in SAM format) + 
    '+' or '-' indicating which strand is the sense strand
2. Comma-separated list of intron start positions in configuration
3. Comma-separated list of intron end positions in configuration
4. left_extend_size: by how many bases on the left side of an intron the
    reference should extend
5. right_extend_size: by how many bases on the right side of an intron the
    reference should extend
6. Read sequence or reversed complement, whatever's first in alphabetical order

ALL OUTPUT COORDINATES ARE 1-INDEXED.
"""
import sys
import os
import site
import subprocess

base_path = os.path.abspath(
                    os.path.dirname(os.path.dirname(os.path.dirname(
                        os.path.realpath(__file__)))
                    )
                )
utils_path = os.path.join(base_path, 'rna', 'utils')
site.addsitedir(utils_path)
site.addsitedir(base_path)

import bowtie
from dooplicity.tools import xstream, register_cleanup, xopen, \
    make_temp_dir
from dooplicity.counters import Counter
from dooplicity.ansibles import Url
import tempdel
import filemover

# Initialize global variable for tracking number of input lines
_input_line_count = 0
counter = Counter('cojunction_enum')
register_cleanup(counter.flush)

def go(input_stream=sys.stdin, output_stream=sys.stdout, bowtie2_exe='bowtie2',
    bowtie2_index_base='genome', bowtie2_args='', verbose=False,
    report_multiplier=1.2, stranded=False, fudge=5, max_refs=300, score_min=60,
    gzip_level=3, mover=filemover.FileMover(), intermediate_dir='.',
    scratch=None):
    """ Runs Rail-RNA-cojunction_enum 

        Alignment script for MapReduce pipelines that wraps Bowtie 2. Finds
        junctions that cooccur on reads by local alignments to transcriptome
        elements from Bowtie 2.

        Input (read from stdin)
        ----------------------------
        Tab-delimited output tuple columns (readletize)
        1. SEQ or its reversed complement, whichever is first in alphabetical
            order
        2. Comma-separated list of sample labels if field 1 is the read
            sequence; '\x1c' if empty
        3. Comma-separated list of sample labels if field 1 is the reversed
            complement of the read sequence; '\x1c' if empty

        Hadoop output (written to stdout)
        ----------------------------
        Tab-delimited tuple columns:
        1. Reference name (RNAME in SAM format) + 
            '+' or '-' indicating which strand is the sense strand
        2. Comma-separated list of intron start positions in configuration
        3. Comma-separated list of intron end positions in configuration
        4. left_extend_size: by how many bases on the left side of an intron
            the reference should extend
        5. right_extend_size: by how many bases on the right side of an intron
            the reference should extend
        6. Read sequence

        input_stream: where to find input reads.
        output_stream: where to emit exonic chunks and junctions.
        bowtie2_exe: filename of Bowtie 2 executable; include path if not in
            $PATH.
        bowtie2_index_base: the basename of the Bowtie index files associated
            with the reference.
        bowtie2_args: string containing precisely extra command-line arguments
            to pass to Bowtie 2, e.g., "--tryhard --best"; or None.
        verbose: True iff more informative messages should be written to
            stderr.
        report_multiplier: if verbose is True, the line number of an alignment
            written to stderr increases exponentially with base
            report_multiplier.
        stranded: True iff input reads are strand-specific; this affects
            whether an output partition has a terminal '+' or '-' indicating
            the sense strand. Further, if stranded is True, an alignment is
            returned only if its strand agrees with the junction's strand.
        fudge: by how many bases to extend left and right extend sizes
                to accommodate potential indels
        max_refs: hard limit on number of reference seqs to enumerate per
            read per strand
        score_min: Bowtie2 CONSTANT minimum alignment score
        gzip_level: compression level to use for temporary files
        mover: FileMover object, for use in case Bowtie2 idx needs to be
            pulled from S3
        intermediate_dir: where intermediates are stored; for temporarily
            storing transcript index if it needs to be pulled from S3
        scratch: scratch directory for storing temporary files or None if 
            securely created temporary directory

        No return value.
    """
    bowtie2_index_base_url = Url(bowtie2_index_base)
    if bowtie2_index_base_url.is_s3:
        index_basename = os.path.basename(bowtie2_index_base)
        index_directory = os.path.join(intermediate_dir, 'transcript_index')
        if not os.path.exists(os.path.join(index_directory, '_STARTED')):
            # Download index
            counter.add('index_download')
            with open(os.path.join(index_directory, '_STARTED'), 'w') \
                as started_stream:
                print >>started_stream, 'STARTED'
            for extension in ['.1.bt2', '.2.bt2', '.3.bt2', '.4.bt2', 
                                '.rev.1.bt2', '.rev.2.bt2']:
                mover.get(bowtie2_index_base_url + extension, index_directory)
            with open(os.path.join(index_directory, '_SUCCESS'), 'w') \
                as success_stream:
                print >>success_stream, 'SUCCESS'
        while not os.path.exists(os.path.join(index_directory, '_SUCCESS')):
            time.sleep(0.5)
        bowtie2_index_base = os.path.join(index_directory, index_basename)  
    global _input_line_count
    temp_dir_path = make_temp_dir(scratch)
    register_cleanup(tempdel.remove_temporary_directories, [temp_dir_path])
    reads_file = os.path.join(temp_dir_path, 'reads.temp.gz')
    with xopen(True, reads_file, 'w', gzip_level) as reads_stream:
        for _input_line_count, line in enumerate(input_stream):
            seq = line.strip()
            counter.add('reads_to_temp')
            print >>reads_stream, '\t'.join([seq, seq, 'I'*len(seq)])
    gzip_command = 'gzip -cd {}'.format(reads_file)
    bowtie_command = ' '.join([bowtie2_exe,
        bowtie2_args if bowtie2_args is not None else '',
        ' --local -t --no-hd --mm -x', bowtie2_index_base, '--12 -',
        '--score-min L,%d,0' % score_min,
        '-D 24 -R 3 -N 1 -L 20 -i L,4,0'])
    delegate_command = ''.join(
            [sys.executable, ' ', os.path.realpath(__file__)[:-3],
                ('_delegate.py --report-multiplier %08f --fudge %d '
                 '--max-refs %d %s %s') % (report_multiplier, fudge, max_refs,
                                            '--stranded' if stranded else '',
                                            '--verbose' if verbose else '')]
        )
    full_command = ' | '.join([gzip_command,
                                bowtie_command, delegate_command])
    print >>sys.stderr, 'Starting Bowtie2 with command: ' + full_command
    bowtie_process = subprocess.Popen(' '.join(
                ['set -exo pipefail;', full_command]
            ), bufsize=-1, stdout=sys.stdout, stderr=sys.stderr,
        shell=True, executable='/bin/bash')
    _input_line_count += 1
    return_code = bowtie_process.wait()
    counter.add('bowtie2_subprocess_done')
    if return_code:
        raise RuntimeError('Error occurred while reading Bowtie 2 output; '
                           'exitlevel was %d.' % return_code)

if __name__ == '__main__':
    import argparse
    # Print file's docstring if -h is invoked
    parser = argparse.ArgumentParser(description=__doc__, 
                formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--report_multiplier', type=float, required=False,
        default=1.2,
        help='When --verbose is also invoked, the only lines of lengthy '
             'intermediate output written to stderr have line number that '
             'increases exponentially with this base')
    parser.add_argument('--verbose', action='store_const', const=True,
        default=False,
        help='Print out extra debugging statements')
    parser.add_argument('--test', action='store_const', const=True,
        default=False,
        help='Run unit tests; DOES NOT NEED INPUT FROM STDIN')
    parser.add_argument('--keep-alive', action='store_const', const=True,
        default=False,
        help='Periodically print Hadoop status messages to stderr to keep ' \
             'job alive')
    parser.add_argument('--fudge', type=int, required=False,
        default=5,
        help='Permits a sum of exonic bases for a junction combo to be '
             'within the specified number of bases of a read sequence\'s '
             'size; this allows for indels with respect to the reference')
    parser.add_argument(
        '--stranded', action='store_const', const=True, default=False,
        help='Assume input reads come from the sense strand; then partitions '
             'in output have terminal + and - indicating sense strand')
    parser.add_argument('--score-min', type=int, required=False,
        default=48,
        help='Bowtie2 minimum CONSTANT score to use')
    parser.add_argument('--max-refs', type=int, required=False,
        default=300,
        help='Hard limit on the number of reference sequences to emit '
             'per read per strand. Prioritizes reference sequences that '
             'overlap the fewest junctions')
    parser.add_argument('--gzip-level', type=int, required=False,
        default=3,
        help='Gzip compression level to use for temporary Bowtie input file')
    parser.add_argument('--intermediate-dir', type=str, required=False,
        default='./',
        help='Where to put transcript index if it needs to be downloaded')

    # Add command-line arguments for dependencies
    bowtie.add_args(parser)
    filemover.add_args(parser)
    tempdel.add_args(parser)

    # Collect Bowtie arguments, supplied in command line after the -- token
    argv = sys.argv
    bowtie2_args = ''
    in_args = False
    for i, argument in enumerate(sys.argv[1:]):
        if in_args:
            bowtie2_args += argument + ' '
        if argument == '--':
            argv = sys.argv[:i + 1]
            in_args = True

    '''Now collect other arguments. While the variable args declared below is
    global, properties of args are also arguments of the go() function so
    different command-line arguments can be passed to it for unit tests.'''
    args = parser.parse_args(argv[1:])
    mover = filemover.FileMover(args=args)

    # Start keep_alive thread immediately
    if args.keep_alive:
        from dooplicity.tools import KeepAlive
        keep_alive_thread = KeepAlive(sys.stderr)
        keep_alive_thread.start()
    
if __name__ == '__main__' and not args.test:
    import time
    start_time = time.time()
    go(bowtie2_exe=os.path.expandvars(args.bowtie2_exe),
        bowtie2_index_base=os.path.expandvars(args.bowtie2_idx),
        bowtie2_args=bowtie2_args, 
        verbose=args.verbose,
        report_multiplier=args.report_multiplier,
        stranded=args.stranded,
        fudge=args.fudge,
        max_refs=args.max_refs,
        score_min=args.score_min,
        mover=mover,
        intermediate_dir=args.intermediate_dir,
        scratch=tempdel.silentexpandvars(args.scratch))
    print >>sys.stderr, ('DONE with cojunction_enum.py; in=%d; '
                            'time=%0.3f s') % (_input_line_count,
                                                time.time() - start_time)
elif __name__ == '__main__':
    # Test units
    del sys.argv[1:] # Don't choke on extra command-line parameters
    import unittest
    # Add unit tests here
    unittest.main()
