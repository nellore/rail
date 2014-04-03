#!/usr/bin/env python
"""
Rail-RNA-align_readlets
Follows Rail-RNA-readletize after an optional Rail-RNA-sum combine/reduce step
Precedes Rail-RNA-intron_search

Alignment script for MapReduce pipelines that wraps Bowtie. Aligns input
readlet sequences and writes a single output line per readlet belonging to
a distinct read sequence.

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
import string
import tempfile
import atexit

base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for directory_name in ['bowtie', 'util']:
    site.addsitedir(os.path.join(base_path, directory_name))

'''For creating a persistent dictionary with an SQLite3 backend; see
https://pypi.python.org/pypi/sqlitedict'''
import sqlitedict
import bowtie

# Initialize global variables for tracking number of input/output lines
_input_line_count = 0
_output_line_count = 0

def handle_temporary_directory(temp_dir_path):
    """ Deletes temporary directory.

        temp_dir_path: path of temporary directory for storing intermediate
            alignments; archived if archive is not None.

        No return value.
    """
    import shutil
    shutil.rmtree(temp_dir_path)

class BowtieOutputThread(threading.Thread):
    """ Processes Bowtie alignments, emitting tuples for exons and introns. """
    
    def __init__(self, input_stream, qname_to_reads, output_stream=sys.stdout,
        verbose=False, report_multiplier=1.2):
        """ Constructor for BowtieOutputThread.

            input_stream: where to retrieve Bowtie's SAM output, typically a
                Bowtie process's stdout.
            output_stream: where to emit exon and intron tuples; typically,
                this is sys.stdout.
            qname_to_reads: maps Bowtie's output qname to seqs associated
                with aligned readlet; that is, the second column in the input
                to this file
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
        self.qname_to_reads = qname_to_reads

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
            last_qname = self.qname_to_reads[last_qname]
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
                qname = self.qname_to_reads[qname]
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
                    reads = last_qname.split('\x1d')
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
                    '''Readlet had no reported alignments; print ONLY when 
                    readlet contains general info about read.'''
                    reads = last_qname.split('\x1d')
                    for read in reads:
                        read_id, _, read_rest = read.partition('\x1e')
                        if len(read_rest.split('\x1e')) > 2:
                            print >>self.output_stream, \
                                '%s\t%s\t\x1c\t\x1c\t\x1c' % (read_id[:-1],
                                                                read_rest)
                        _output_line_count += 1
                multireadlet = []
            if not line: break
            last_tokens = tokens
            (last_qname, last_flag, last_rname, last_pos) = (qname, flag,
                rname, pos)

def go(input_stream=sys.stdin, output_stream=sys.stdout, bowtie_exe='bowtie',
    bowtie_index_base='genome', bowtie_args=None, verbose=False,
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
    # For creating persistent dictionary
    db_directory = tempfile.mkdtemp()
    atexit.register(handle_temporary_directory, db_directory)
    db_filename = os.path.join(db_directory, 'readinfo.db')
    qname_to_reads = sqlitedict.SqliteDict(db_filename)
    bowtie_process, bowtie_command, threads = bowtie.proc(
            bowtieExe=bowtie_exe, bowtieIdx=bowtie_index_base,
            readFn=None, bowtieArgs=bowtie_args, sam=True,
            stdoutPipe=True, stdinPipe=True
        )
    output_thread = BowtieOutputThread(
            bowtie_process.stdout,
            qname_to_reads,
            verbose=verbose, 
            output_stream=output_stream,
            report_multiplier=report_multiplier
        )
    threads.append(output_thread)
    output_thread.start()
    line = None
    try:
        for _input_line_count, line in enumerate(input_stream):
            tokens = line.rstrip().split('\t')
            assert len(tokens) == 2
            seq, qname = tokens
            qual_seq = 'I'*len(seq)
            index = str(_input_line_count)
            qname_to_reads[index] = qname
            print >>bowtie_process.stdin, '\t'.join([index, seq, qual_seq])
        bowtie_process.stdin.close()
    except Exception:
        print >>sys.stderr, 'Error. Stats: input line count=%d, ' \
            'output line count=%d, line=%s' % (_input_line_count + 1, 
                                                _output_line_count,
                                                line if line is not None else
                                                'line not defined')
        raise
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