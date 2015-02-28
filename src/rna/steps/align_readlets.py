#!/usr/bin/env python
"""
Rail-RNA-align_readlets

Follows Rail-RNA-align_reads
Precedes Rail-RNA-intron_search

Alignment script for MapReduce pipelines that wraps Bowtie. Aligns input
readlet sequences and writes a single output line per readlet belonging to
a distinct read sequence. THIS CODE ASSUMES THAT BOWTIE RUNS ON JUST ONE
THREAD; it exploits that Bowtie returns reads in the order in which they were
sent, which is guaranteed only when on a single thread.

Input (read from stdin)
----------------------------
Tab-delimited input tuple columns:
1. Readlet sequence or its reversed complement, whichever is first in
    alphabetical order
2. '\x1d'-separated list of [read sequence ID + ('-' if readlet sequence is
    reverse-complemented; else '+') + '\x1e' + displacement of readlet's 5' end
    from read's 5' end + '\x1e' + displacement of readlet's 3' end from read's
    3' end (+, for EXACTLY one readlet of a read sequence, '\x1e' +
    read sequence + '\x1e' + (an '\x1f'-separated list A of unique sample
    labels with read sequences that match the original read sequence) + '\x1e'
    + (an '\x1f'-separated list  of unique sample labels B with read sequences
    that match the reversed complement of the original read sequence))] + 
    '\x1e' + (an '\x1f'-separated list of the number of instances of the read
    sequence for each respective sample in list A) + '\x1e' + (an
    '\x1f'-separated list of the number of instances of the read
    sequence's reversed complement for each respective sample in list B)].
    Here, a read sequence ID takes the form X:Y, where X is the
    "mapred_task_partition" environment variable -- a unique index for a task
    within a job -- and Y is the index of the read sequence relative to the
    beginning of the input stream.

Input is partitioned by field 1, the readlet sequence or its reversed
complement.

Hadoop output (written to stdout)
----------------------------
Tab-delimited output tuple columns, where each line corresponds to a readlet
from a distinct read rather than a unique readlet sequence:
1. Read sequence ID
2. Displacement of readlet's 5' end from read's 5' end + '\x1e' +
    displacement of readlet's 3' end from read's 3' end (+, for EXACTLY one
    readlet of a read sequence, '\x1e' + read sequence + '\x1e' + number of
    instances of read sequence + '\x1e' + number of instances of read
    sequence's reversed complement + '\x1e' (+, for EXACTLY one readlet of a
    read sequence, '\x1e' + read sequence + '\x1e' + (an '\x1f'-separated list
    A of unique sample labels with read sequences that match the original read
    sequence) + '\x1e' + (an '\x1f'-separated list  of unique sample labels B
    with read sequences that match the reversed complement of the original read
    sequence))] + '\x1e' + (an '\x1f'-separated list of the number of instances
    of the read sequence for each respective sample in list A) + '\x1e' + (an
    '\x1f'-separated list of the number of instances of the read sequence's
    reversed complement for each respective sample in list B)
3. '\x1f'-separated list of alignment RNAMEs or '\x1c' if no alignments found
4. '\x1f'-separated list of alignment FLAGs or '\x1c' if no alignments found
5. '\x1f'-separated list of alignment POSes or '\x1c' if no alignments found

ALL OUTPUT COORDINATES ARE 1-INDEXED.
"""
import sys
import os
import site
import subprocess
import tempfile

base_path = os.path.abspath(
                    os.path.dirname(os.path.dirname(os.path.dirname(
                        os.path.realpath(__file__)))
                    )
                )
utils_path = os.path.join(base_path, 'rna', 'utils')
site.addsitedir(utils_path)
site.addsitedir(base_path)

import bowtie
from dooplicity.tools import xstream, register_cleanup, xopen, make_temp_dir
import tempdel

# Initialize global variable for tracking number of input lines
_input_line_count = 0

def go(input_stream=sys.stdin, output_stream=sys.stdout, bowtie_exe='bowtie',
    bowtie_index_base='genome', bowtie_args='', gzip_level=3, verbose=False,
    report_multiplier=1.2, scratch=None):
    """ Runs Rail-RNA-align_readlets.

        Aligns input readlet sequences and writes a single output line per
        readlet belonging to a distinct read sequence.

        Input (read from stdin)
        ----------------------------
        Tab-delimited input tuple columns:
        1. Readlet sequence or its reversed complement, whichever is first in
            alphabetical order
        2. '\x1d'-separated list of [read sequence ID + ('-' if readlet
            sequence is reverse-complemented; else '+') + '\x1e' + displacement
            of readlet's 5' end from read's 5' end + '\x1e' + displacement of
            readlet's 3' end from read's 3' end (+, for EXACTLY one readlet of
            a read sequence, '\x1e' + read sequence + '\x1e' +
            (an '\x1f'-separated list A of unique sample labels with read
            sequences that match the original read sequence) + '\x1e' +
            (an '\x1f'-separated list  of unique sample labels B with read
            sequences that match the reversed complement of the original read
            sequence)) + '\x1e' + (an '\x1f'-separated list of the number of
            instances of the read sequence for each respective sample in list
            A) + '\x1e' + (an '\x1f'-separated list of the number of instances
            of the read sequence's reversed complement for each respective
            sample in list B)]. Here, a read sequence ID takes the form X:Y,
            where X is the "mapred_task_partition" environment variable -- a
            unique index for a task within a job -- and Y is the index of the
            read sequence relative to the beginning of the input stream.

        Input is partitioned by field 1, the readlet sequence or its reversed
        complement.

        Hadoop output (written to stdout)
        ----------------------------
        Tab-delimited output tuple columns, where each line corresponds to a
        readlet from a distinct read rather than a unique readlet sequence:
        1. Read sequence ID
        2. Displacement of readlet's 5' end from read's 5' end + '\x1e' +
            displacement of readlet's 3' end from read's 3' end (+, for EXACTLY
            one readlet of a read sequence, '\x1e' + read sequence + '\x1e' +
            number of instances of read sequence + '\x1e' + number of instances
            of read sequence's reversed complement + '\x1e' (+, for EXACTLY one
            readlet of a read sequence, '\x1e' + read sequence + '\x1e' +
            (an '\x1f'-separated list A of unique sample labels with read
            sequences that match the original read sequence) + '\x1e' +
            (an '\x1f'-separated list  of unique sample labels B with read
            sequences that match the reversed complement of the original read
            sequence))] + '\x1e' + (an '\x1f'-separated list of the number of
            instances of the read sequence for each respective
            sample in list A) + '\x1e' + (an '\x1f'-separated list of the
            number of instances of the read sequence's reversed complement for
            each respective sample in list B)
        3. '\x1f'-separated list of alignment RNAMEs or '\x1c' if no alignments
        4. '\x1f'-separated list of alignment FLAGs or '\x1c' if no alignments
        5. '\x1f-separated list of alignment POSes or '\x1c' if no alignments

        ALL OUTPUT COORDINATES ARE 1-INDEXED.

        input_stream: where to find input reads.
        output_stream: where to emit exonic chunks and introns.
        bowtie_exe: filename of Bowtie executable; include path if not in
            $PATH.
        bowtie_index_base: the basename of the Bowtie index files associated
            with the reference.
        bowtie_args: string containing precisely extra command-line arguments
            to pass to first-pass Bowtie, e.g., "--tryhard --best"; or None.
        gzip_level: level of gzip compression to use for qname file
        verbose: True iff more informative messages should be written to
            stderr.
        report_multiplier: if verbose is True, the line number of an alignment
            written to stderr increases exponentially with base
            report_multiplier.
        scratch: scratch directory for storing temporary files or None if 
            securely created temporary directory

        No return value.
    """
    global _input_line_count
    # For storing long qnames
    temp_dir = make_temp_dir(scratch)
    register_cleanup(tempdel.remove_temporary_directories, [temp_dir])
    qnames_file = os.path.join(temp_dir, 'qnames.temp.gz')
    readlet_file = os.path.join(temp_dir, 'readlets.temp.gz')
    with xopen(True, qnames_file, 'w', gzip_level) as qname_stream:
        with xopen(True, readlet_file, 'w', gzip_level) as readlet_stream:
            for (seq_count, ((seq,), xpartition)) \
                in enumerate(xstream(input_stream, 1)):
                print >>readlet_stream, \
                    '\t'.join([str(seq_count), seq, 'I'*len(seq)])
                qname_stream.write(next(iter(xpartition))[0])
                for (qname,) in xpartition:
                    _input_line_count += 1
                    qname_stream.write('\x1d' + qname)
                qname_stream.write('\n')
    input_command = 'gzip -cd %s' % readlet_file
    bowtie_command = ' '.join([bowtie_exe, bowtie_args,
        '-S -t --sam-nohead --mm', bowtie_index_base, '--12 -'])
    delegate_command = ''.join(
                [sys.executable, ' ', os.path.realpath(__file__)[:-3],
                    '_delegate.py --report-multiplier %08f --qnames-file %s %s'
                        % (report_multiplier, qnames_file,
                            '--verbose' if verbose else '')]
            )
    full_command = ' | '.join([input_command, 
                                bowtie_command, delegate_command])
    print >>sys.stderr, 'Starting Bowtie with command: ' + full_command
    bowtie_process = subprocess.Popen(' '.join(
                    ['set -exo pipefail;', full_command]
                ),
            bufsize=-1, stdout=sys.stdout, stderr=sys.stderr, shell=True,
            executable='/bin/bash')
    return_code = bowtie_process.wait()
    if return_code:
        raise RuntimeError('Error occurred while reading Bowtie output; '
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
    parser.add_argument('--gzip-level', type=int, required=False,
        default=3,
        help=('Level of gzip compression to use for temporary file storing '
              'qnames.'))
    parser.add_argument('--keep-alive', action='store_const', const=True,
        default=False,
        help='Periodically print Hadoop status messages to stderr to keep ' \
             'job alive')

    # Add command-line arguments for dependencies
    bowtie.add_args(parser)
    tempdel.add_args(parser)

    # Collect Bowtie arguments, supplied in command line after the -- token
    argv = sys.argv
    bowtie_args = ''
    in_args = False
    for i, argument in enumerate(sys.argv[1:]):
        if in_args:
            bowtie_args += argument + ' '
        if argument == '--':
            argv = sys.argv[:i + 1]
            in_args = True

    '''Now collect other arguments. While the variable args declared below is
    global, properties of args are also arguments of the go() function so
    different command-line arguments can be passed to it for unit tests.'''
    args = parser.parse_args(argv[1:])

    # Start keep_alive thread immediately
    if args.keep_alive:
        from dooplicity.tools import KeepAlive
        keep_alive_thread = KeepAlive(sys.stderr)
        keep_alive_thread.start()

if __name__ == '__main__':
    import time
    start_time = time.time()
    go(bowtie_exe=args.bowtie_exe,
        bowtie_index_base=args.bowtie_idx,
        bowtie_args=bowtie_args,
        gzip_level=args.gzip_level,
        verbose=args.verbose,
        report_multiplier=args.report_multiplier,
        scratch=args.scratch)
    print >>sys.stderr, 'DONE with align_readlets.py; in=%d; ' \
        'time=%0.3f s' % (_input_line_count, time.time() - start_time)
