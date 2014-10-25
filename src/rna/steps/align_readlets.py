#!/usr/bin/env python
"""
Rail-RNA-align_readlets
Follows Rail-RNA-readletize after an optional Rail-RNA-sum combine/reduce step
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
    read sequence + '\x1e' + number of instances of read sequence + '\x1e' +
    number of instances of read sequence's reversed complement + '\x1e' + (an
    '\x1f'-separated set of unique sample labels with read sequences that match
    the original read sequence) + '\x1e' + (an '\x1f'-separated set of unique
    sample labels with read sequences that match the reversed complement of the
    original read sequence))]. Here, a read sequence ID takes the form X:Y,
    where X is the "mapred_task_partition" environment variable -- a unique
    index for a task within a job -- and Y is the index of the read sequence
    relative to the beginning of the input stream.

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
    sequence's reversed complement + '\x1e' + (an '\x1f'-separated set of
    unique sample labels with read sequences that match the original read
    sequence) + '\x1e' + (an '\x1f'-separated set of unique sample labels with
    read sequences that match the reversed complement of the original read
    sequence))
3. '\x1f'-separated list of alignment RNAMEs or '\x1c' if no alignments found
4. '\x1f'-separated list of alignment FLAGs or '\x1c' if no alignments found
5. '\x1f'-separated list of alignment POSes or '\x1c' if no alignments found

ALL OUTPUT COORDINATES ARE 1-INDEXED.
"""
import sys
import os
import site
import threading
import subprocess
import tempfile
import atexit

base_path = os.path.abspath(
                    os.path.dirname(os.path.dirname(os.path.dirname(
                        os.path.realpath(__file__)))
                    )
                )
utils_path = os.path.join(base_path, 'rna', 'utils')
site.addsitedir(utils_path)
site.addsitedir(base_path)

import bowtie
from dooplicity.tools import xstream

# Initialize global variables for tracking number of input/output lines
_input_line_count = 0
_output_line_count = 0

def handle_temporary_directory(temp_dir_path):
    """ Deletes temporary directory.

        temp_dir_paths: path to temporary directory to delete

        No return value.
    """
    import shutil
    shutil.rmtree(temp_dir_path)

import re
_expression = re.compile('\x1d')

def itersplit(a_string):
    """ Generator for splitting a string between '\x1d's.

        Based on http://stackoverflow.com/questions/3862010/
        is-there-a-generator-version-of-string-split-in-python

        a_string: string

        Yield value: next part of string between \x1ds.
    """
    pos = 0
    while True:
        next = _expression.search(a_string, pos)
        if not next:
            yield a_string[pos:]
            break
        yield a_string[pos:next.start()]
        pos = next.end()

class BowtieOutputThread(threading.Thread):
    """ Processes Bowtie alignments, emitting tuples for exons and introns. """
    
    def __init__(self, input_stream, qname_stream, return_set,
        output_stream=sys.stdout, verbose=False, report_multiplier=1.2):
        """ Constructor for BowtieOutputThread.

            input_stream: where to retrieve Bowtie's SAM output, typically a
                Bowtie process's stdout.
            output_stream: where to emit exon and intron tuples; typically,
                this is sys.stdout.
            qname_stream: where to find long names containing read information
                associated with readlets
            return_set: 0 is added to this set if a thread finishes
                successfully
            verbose: True if alignments should occasionally be written 
                to stderr.
            bin_size: genome is partitioned in units of bin_size for later load
                balancing.
            report_multiplier: if verbose is True, the line number of an
                alignment written to stderr increases exponentially with base
                report_multiplier.
        """
        super(BowtieOutputThread, self).__init__()
        self.daemon = True
        self.input_stream = input_stream
        self.output_stream = output_stream
        self.verbose = verbose
        self.report_multiplier = report_multiplier
        self.qname_stream = qname_stream
        self.return_set = return_set

    def run(self):
        """ Prints exons for end-to-end alignments.

            Overrides default method containing thread activity.

            No return value.
        """
        global _output_line_count
        next_report_line, i = 0, 0
        for (qname,), xpartition in xstream(self.input_stream, 1):
            '''While labeled multireadlet, this list may end up simply a
            unireadlet.'''
            multireadlet = []
            for tokens in xpartition:
                (flag, rname, pos, mapq, cigar,
                    rnext, pnext, tlen, seq, qual) = tokens[:10]
                flag = int(flag)
                multireadlet.append((rname, flag, pos))
                if self.verbose and next_report_line == i:
                    print >>sys.stderr, \
                        'SAM output record %d: rdname="%s", flag=%d' % (i,
                                                                        qname,
                                                                        flag)
                    next_report_line = int((next_report_line + 1)
                                            * self.report_multiplier + 1) - 1
                i += 1
            '''If the next qname doesn't match the last qname or there are no
            more lines, all of a multireadlet's alignments have been
            collected.'''
            reads = itersplit(self.qname_stream.readline().rstrip())
            if not flag & 4:
                '''Last readlet has at least one alignment; print all
                alignments for each read from which readlet sequence is
                derived.'''
                rnames, flags, poses = zip(*multireadlet)
                reverse_flags = [a_flag ^ 16 for a_flag in flags]
                flags = '\x1f'.join([str(a_flag) for a_flag in flags])
                reverse_flags = '\x1f'.join(
                                        [str(a_flag) for a_flag
                                            in reverse_flags]
                                    )
                rnames = '\x1f'.join(rnames)
                poses = '\x1f'.join(poses)
                for read in reads:
                    read_id, _, read_rest = read.partition('\x1e')
                    if read_id[-1] == '-':
                        current_flags = reverse_flags
                    else:
                        current_flags = flags
                    print >>self.output_stream, '%s\t%s\t%s\t%s\t%s' % \
                        (read_id[:-1], read_rest, rnames,
                            current_flags, poses)
                    _output_line_count += 1
            else:
                '''Readlet had no reported alignments; print ONLY when readlet
                contains general info about read.'''
                for read in reads:
                    read_id, _, read_rest = read.partition('\x1e')
                    if len(read_rest.split('\x1e')) > 2:
                        print >>self.output_stream, \
                            '%s\t%s\t\x1c\t\x1c\t\x1c' % (read_id[:-1],
                                                            read_rest)
                    _output_line_count += 1
        self.output_stream.flush()
        self.return_set.add(0)

def go(input_stream=sys.stdin, output_stream=sys.stdout, bowtie_exe='bowtie',
    bowtie_index_base='genome', bowtie_args='', verbose=False,
    report_multiplier=1.2):
    """ Runs Rail-RNA-align_readlets.

        Aligns input readlet sequences and writes a single output line per
        readlet belonging to a distinct read sequence.

        Input (read from stdin)
        ----------------------------
        Tab-delimited input tuple columns:
        1. Readlet sequence or its reversed complement, whichever is first in
            alphabetical order
        2. '\x1d'-separated list of [read sequence ID + 
            ('-' if readlet sequence is reverse-complemented; else '+')
            + '\x1e' + displacement of readlet's 5' end from read's
            5' end + '\x1e' + displacement of readlet's 3' end from read's 3'
            end (+, for EXACTLY one readlet of a read sequence, '\x1e' + 
            read sequence + '\x1e' + number of instances of read sequence +
            '\x1e' + number of instances of read sequence's reversed complement
            + '\x1e' + (an '\x1f'-separated set of unique sample labels with
            read sequences that match the original read sequence) + '\x1e' +
            (an '\x1f'-separated set of unique sample labels with read
            sequences that match the reversed complement of the original read
            sequence)]. Here, a read sequence ID takes the form X:Y, where X is
            the "mapred_task_partition" environment variable -- a unique index
            for a task within a job -- and Y is the index of the read sequence
            relative to the beginning of the input stream.

        Input is partitioned by field 1, the readlet sequence or its reversed
        complement.

        Hadoop output (written to stdout)
        ----------------------------
        Tab-delimited output tuple columns, where each line corresponds to a
        readlet from a distinct read rather than a unique readlet sequence:
        1. Read sequence ID
        2. Displacement of readlet's 5' end from read's 5' end + ';' +
            displacement of readlet's 3' end from read's 3' end (+, for EXACTLY
            one readlet of a read sequence, ';' + number of readlets in read
            sequence + ';' + read sequence + ';' + number of instances of read
            sequence + ';' + number of instances of read sequence's reversed
            complement).
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
        verbose: True iff more informative messages should be written to
            stderr.
        report_multiplier: if verbose is True, the line number of an alignment
            written to stderr increases exponentially with base
            report_multiplier.    

        No return value.
    """
    global _input_line_count
    # For storing long qnames
    temp_dir = tempfile.mkdtemp()
    atexit.register(handle_temporary_directory, temp_dir)
    qname_file = os.path.join(temp_dir, 'qnames.temp')
    readlet_file = os.path.join(temp_dir, 'readlets.temp')
    output_file = os.path.join(temp_dir, 'out.sam')
    with open(qname_file, 'w') as qname_stream:
        with open(readlet_file, 'w') as readlet_stream:
            for (seq_count, ((seq,), xpartition)) \
                in enumerate(xstream(input_stream, 1)):
                print >>readlet_stream, \
                    '\t'.join([str(seq_count), seq, 'I'*len(seq)])
                qname_stream.write(next(iter(xpartition))[0])
                for (qname,) in xpartition:
                    _input_line_count += 1
                    qname_stream.write('\x1d' + qname)
                qname_stream.write('\n')
                qname_stream.flush()
    bowtie_command = ' '.join([bowtie_exe,
        bowtie_args,
        '-S -t --sam-nohead --mm', bowtie_index_base, '--12', readlet_file,
        output_file])
    print >>sys.stderr, 'Starting Bowtie with command: ' + bowtie_command
    bowtie_process = subprocess.Popen(bowtie_command, bufsize=-1, shell=True,
        stdout=subprocess.PIPE, stderr=sys.stderr)
    bowtie_process.wait()
    if os.path.exists(output_file):
        return_set = set()
        with open(qname_file) as qname_stream:
            output_thread = BowtieOutputThread(
                    open(output_file),
                    qname_stream,
                    return_set,
                    verbose=verbose, 
                    output_stream=output_stream,
                    report_multiplier=report_multiplier
                )
            output_thread.start()
            # Join thread to pause execution in main thread
            if verbose: print >>sys.stderr, 'Joining thread...'
            output_thread.join()
        if not return_set:
            raise RuntimeError('Error occurred in BowtieOutputThread.')

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

    # Add command-line arguments for dependencies
    bowtie.add_args(parser)

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

if __name__ == '__main__' and not args.test:
    import time
    start_time = time.time()
    go(bowtie_exe=args.bowtie_exe,
        bowtie_index_base=args.bowtie_idx,
        bowtie_args=bowtie_args, 
        verbose=args.verbose,
        report_multiplier=args.report_multiplier)
    print >>sys.stderr, 'DONE with align_readlets.py; in/out=%d/%d; ' \
        'time=%0.3f s' % (_input_line_count, _output_line_count,
                            time.time() - start_time)
elif __name__ == '__main__':
    # Test units
    del sys.argv[1:] # Don't choke on extra command-line parameters
    import unittest
    # Add unit tests here
    unittest.main()
