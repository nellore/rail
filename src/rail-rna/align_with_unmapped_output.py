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
(align)
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

(unmapped) -- unmapped readlets saved for later alignment
Tab-delimited output tuple columns, where each line 
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

ALL OUTPUT COORDINATES ARE 1-INDEXED.
"""
import sys
import os
import site
import threading
import subprocess
import string
import tempfile
import atexit

base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
site.addsitedir(os.path.join(base_path, 'bowtie'))

import bowtie

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

class BowtieOutputThread(threading.Thread):
    """ Processes Bowtie alignments, emitting tuples for exons and introns. """
    
    def __init__(self, input_stream, seq_and_qname_stream, 
        output_unmapped_readlets=False, readlet_size=25,
        output_stream=sys.stdout, verbose=False, report_multiplier=1.2):
        """ Constructor for BowtieOutputThread.

            input_stream: where to retrieve Bowtie's SAM output, typically a
                Bowtie process's stdout.
            output_stream: where to emit exon and intron tuples; typically,
                this is sys.stdout.
            seq_and_qname_stream: where to find seqs and long names containing
                read information associated with readlets
            output_unmapped_readlets: output unmapped reads for later alignment
                to intron reference.
            readlet_size: size of unmapped readlets to retain.
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
        self.output_unmapped_readlets = output_unmapped_readlets
        self.readlet_size = readlet_size
        self.verbose = verbose
        self.report_multiplier = report_multiplier
        self.seq_and_qname_stream = seq_and_qname_stream

    def run(self):
        """ Prints exons for end-to-end alignments.

            Overrides default method containing thread activity.

            No return value.
        """
        global _output_line_count
        next_report_line = 0
        import random
        '''Next readlet must be known to tell if a readlet mapped to multiple
        locations, so always work with previous read.'''
        while True:
            line = self.input_stream.readline()
            if not line:
                return # Bowtie output nothing
            # Skip header line
            if line[0] == '@': continue
            last_tokens = line.rstrip().split('\t')
            (last_qname, last_flag, last_rname, last_pos, last_mapq,
                last_cigar, last_rnext, last_pnext,
                last_tlen, last_seq, last_qual) = last_tokens[:11]
            last_flag = int(last_flag)
            break
        # While labeled multireadlet, this list may end up simply a unireadlet
        multireadlet = []
        while True:
            line = self.input_stream.readline()
            if line:
                tokens = line.rstrip().split('\t')
                (qname, flag, rname, pos, mapq, cigar, rnext,
                    pnext, tlen, seq, qual) = tokens[:11]
                flag = int(flag)
            if self.verbose and next_report_line == i:
                print >>sys.stderr, \
                    'SAM output record %d: rdname="%s", flag=%d' % (i,
                                                                    last_qname,
                                                                    last_flag)
                next_report_line = int((next_report_line + 1)
                    * self.report_multiplier + 1) - 1
            multireadlet.append((last_rname, last_flag, last_pos))
            if not line or qname != last_qname:
                '''If the next qname doesn't match the last qname or there are
                no more lines, all of a multireadlet's alignments have been
                collected.'''
                reads_line = self.seq_and_qname_stream.readline()
                reads = reads_line.rstrip().split('\t')[1].split('\x1d')
                if not (last_flag & 4):
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
                        print >>self.output_stream, ('align\t%s\t%s\t%s\t%s'
                            + '\t%s') % (read_id[:-1], read_rest, rnames,
                                        current_flags, poses)
                        _output_line_count += 1
                elif self.output_unmapped_readlets:
                    '''Readlet had no reported alignments; print ONLY when 
                    readlet contains general info about read and unmapped
                    readlets aren't being output.'''
                    for read in reads:
                        read_id, _, read_rest = read.partition('\x1e')
                        if len(read_rest.split('\x1e')) > 2:
                            print >>self.output_stream, \
                                'align\t%s\t%s\t\x1c\t\x1c\t\x1c' \
                                % (read_id[:-1], read_rest)
                        _output_line_count += 1
                    if len(seq) == self.readlet_size:
                        sys.stdout.write('unmapped\t' + reads_line)
                multireadlet = []
            if not line: break
            last_tokens = tokens
            (last_qname, last_flag, last_rname, last_pos) = (qname, flag,
                rname, pos)

def go(input_stream=sys.stdin, output_stream=sys.stdout, bowtie_exe='bowtie',
    bowtie_index_base='genome', bowtie_args=None,
    output_unmapped_readlets=False, readlet_size=25, verbose=False,
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
        (align)
        Tab-delimited output tuple columns, where each line corresponds to a
        readlet from a distinct read rather than a unique readlet sequence:
        1. Read sequence ID
        2. Displacement of readlet's 5' end from read's 5' end + '\x1e' +
            displacement of readlet's 3' end from read's 3' end (+, for EXACTLY
            one readlet of a read sequence, '\x1e' + read sequence + '\x1e'
            + number of instances of read sequence + '\x1e' + number of
            instances of read sequence's reversed complement + '\x1e'
            + (an '\x1f'-separated set of unique sample labels with read
            sequences that match the original read sequence) + '\x1e'
            + (an '\x1f'-separated set of unique sample labels with
            read sequences that match the reversed complement of the original
            read sequence))
        3. '\x1f'-separated list of alignment RNAMEs or '\x1c' if no alignments
            found
        4. '\x1f'-separated list of alignment FLAGs or '\x1c' if no alignments
            found
        5. '\x1f'-separated list of alignment POSes or '\x1c' if no alignments
            found

        (unmapped) -- unmapped readlets saved for later alignment
        Tab-delimited output tuple columns, where each line 
        1. Readlet sequence or its reversed complement, whichever is first in
            alphabetical order
        2. '\x1d'-separated list of [read sequence ID + ('-' if readlet
            sequence is reverse-complemented; else '+') + '\x1e' + displacement
            of readlet's 5' end from read's 5' end + '\x1e' + displacement of
            readlet's 3' end from read's 3' end (+, for EXACTLY one readlet of
            a read sequence, '\x1e' + read sequence + '\x1e' + number of
            instances of read sequence + '\x1e' + number of instances of read
            sequence's reversed complement + '\x1e' + (an '\x1f'-separated set 
            of unique sample labels with read sequences that match the original
            read sequence) + '\x1e' + (an '\x1f'-separated set of unique
            sample labels with read sequences that match the reversed
            complement of the original read sequence))]. Here, a read sequence
            ID takes the form X:Y, where X is the "mapred_task_partition"
            environment variable -- a unique index for a task within a job --
            and Y is the index of the read sequence relative to the beginning
            of the input stream.

        ALL OUTPUT COORDINATES ARE 1-INDEXED.

        input_stream: where to find input reads.
        output_stream: where to emit exonic chunks and introns.
        bowtie_exe: filename of Bowtie executable; include path if not in
            $PATH.
        bowtie_index_base: the basename of the Bowtie index files associated
            with the reference.
        bowtie_args: string containing precisely extra command-line arguments
            to pass to first-pass Bowtie, e.g., "--tryhard --best"; or None.
        output_unmapped_readlets: output unmapped reads for later alignment to
            intron reference.
        readlet_size: size of unmapped readlets to retain.
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
    seq_and_qname_file = os.path.join(temp_dir, 'seqs_and_qnames.temp')
    readlet_file = os.path.join(temp_dir, 'readlets.temp')
    with open(seq_and_qname_file, 'w') as seq_and_qname_stream:
        with open(readlet_file, 'w') as readlet_stream:
            for _input_line_count, line in enumerate(input_stream):
                tokens = line.rstrip().split('\t')
                assert len(tokens) == 2
                seq, qname = tokens
                print >>readlet_stream, \
                    '\t'.join([str(_input_line_count), seq, 'I'*len(seq)])
                print >>seq_and_qname_stream, '%s\t%s' % (seq, qname)
    bowtie_process, bowtie_command, threads = bowtie.proc(
            bowtieExe=bowtie_exe, bowtieIdx=bowtie_index_base,
            readFn=readlet_file, bowtieArgs=bowtie_args, sam=True,
            stdoutPipe=True, stdinPipe=False
        )
    with open(seq_and_qname_file) as seq_and_qname_stream:
        output_thread = BowtieOutputThread(
                bowtie_process.stdout,
                seq_and_qname_stream,
                readlet_size=readlet_size,
                output_unmapped_readlets=output_unmapped_readlets,
                verbose=verbose, 
                output_stream=output_stream,
                report_multiplier=report_multiplier
            )
        threads.append(output_thread)
        output_thread.start()
        # Join threads to pause execution in main thread
        for thread in threads:
            if verbose: print >>sys.stderr, 'Joining thread...'
            thread.join()
    output_stream.flush()

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
    parser.add_argument('--readlet-size', type=int, required=False,
        default=25,
        help='If --output-unmapped-readlets is True, unmapped readlets are '
             'filtered for this readlet size; for suppressing capping '
             'readlets')

    # Add command-line arguments for dependencies
    bowtie.addArgs(parser)

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

if __name__ == '__main__' and not args.test:
    import time
    start_time = time.time()
    go(bowtie_exe=args.bowtie_exe,
        bowtie_index_base=args.bowtie_idx,
        bowtie_args=bowtie_args,
        readlet_size=args.readlet_size,
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