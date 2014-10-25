#!/usr/bin/env python
"""
Rail-RNA-realign_reads

Follows Rail-RNA-cointron_search
Precedes Rail-RNA-compare_alignments

Realignment script for MapReduce pipelines that wraps Bowtie2. Creates Bowtie2
index including only sequences framing introns to align only those reads for
which Bowtie2 did not report alignments in Rail-RNA-align. Reference names in
the index encode intron sizes and locations in the (presmably) exonic
sequences it records. Infers exonic chunks and introns from alignments.

THIS CODE ASSUMES THAT BOWTIE RUNS ON JUST ONE THREAD; it exploits that Bowtie
returns reads in the order in which they were sent, which is guaranteed only
when on a single thread.

Input (read from stdin)
----------------------------
Tab-delimited input tuple columns:
(two kinds)

Type 1:
1. Read sequence
2. '\x1c' + FASTA reference name including '>'. The following format is used:
    original RNAME + '+' or '-' indicating which strand is the sense strand
    + '\x1d' + start position of sequence + '\x1d' + comma-separated list of
    subsequence sizes framing introns + '\x1d' + comma-separated list of intron
    sizes) + '\x1d' + 'p' if derived from primary alignment to genome; 's' if
    derived from secondary alignment to genome; 'i' if derived from cointron
    search
3. FASTA sequence

Type 2:
1. Read sequence
2. QNAME
3. Quality sequence

Type 1 corresponds to a FASTA line to index to which the read sequence is
predicted to align. Type 2 corresponds to a distinct read. Input is partitioned
by field 1.

Hadoop output (written to stdout)
----------------------------
Tab-delimited output tuple columns:
Standard 11+ -column SAM output.

ALL OUTPUT COORDINATES ARE 1-INDEXED.
"""
import sys
import os
import site
import threading
import tempfile
import atexit
import subprocess
import time
import string

base_path = os.path.abspath(
                    os.path.dirname(os.path.dirname(os.path.dirname(
                        os.path.realpath(__file__)))
                    )
                )
utils_path = os.path.join(base_path, 'rna', 'utils')
site.addsitedir(utils_path)
site.addsitedir(base_path)

from dooplicity.tools import xstream
import bowtie
import argparse

# Initialize global variables for tracking number of input/output lines
_input_line_count = 0
_output_line_count = 0

def input_files_from_input_stream(input_stream,
                                    temp_dir_path=tempfile.mkdtemp(),
                                    verbose=False):
    """ Creates FASTA reference to index, rname file, and file with reads.

        Each line of the read file is in the following format:

        read number + TAB + SEQ + TAB + QUAL

        Each line of the rname file is in the following format:

        (number of qnames associated with rnames) + '\x1e'
        + ('\x1e'-separated list of valid rnames) .

        The RNAME file contains those RNAMEs to which a given sequence is
        "allowed" to align. When Bowtie2 is run in -a mode and only these
        alignments are permitted, this ensures the reproducibility of
        results regardless of partitioning (and hence # cores).

        input_stream: where to find Hadoop input
        temp_dir_path: where to store files
        verbose: output extra debugging messages

        Return value: tuple (path to FASTA reference file, path to read file)
    """
    global _input_line_count
    prefasta_filename = os.path.join(temp_dir_path, 'temp.prefa')
    deduped_fasta_filename = os.path.join(temp_dir_path, 'temp.deduped.prefa')
    final_fasta_filename = os.path.join(temp_dir_path, 'temp.fa')
    reads_filename = os.path.join(temp_dir_path, 'reads.temp')
    rname_filename = os.path.join(temp_dir_path, 'rnames.temp')
    if verbose:
        print >>sys.stderr, 'Writing prefasta and input reads...'
    with open(prefasta_filename, 'w') as fasta_stream:
        with open(reads_filename, 'w') as read_stream:
            with open(rname_filename, 'w') as rname_stream:
                for (read_seq,), xpartition in xstream(input_stream, 1):
                    rnames = []
                    fasta_lines = []
                    read_count = 0
                    for value in xpartition:
                        _input_line_count += 1
                        if value[0][0] == '\x1c':
                            # Print FASTA line
                            fasta_lines.append('\t'.join([value[0][1:-2],
                                                             value[1]]))
                            rnames.append(value[0][2:-2])
                        else:
                            # Add to temporary seq stream
                            print >>read_stream, '\t'.join([value[0], read_seq,
                                                                value[1]])
                            read_count += 1
                    if read_count:
                        # Print FASTA line iff there's a read to align it to
                        rname_stream.write(str(read_count) + '\x1e')
                        print >>rname_stream, '\x1e'.join(rnames)
                        if fasta_lines:
                            # Print FASTA line iff it exists
                            print >>fasta_stream, '\n'.join(fasta_lines)

    if verbose:
        print >>sys.stderr, 'Done! Sorting and deduplicating prefasta...'
    # Sort prefasta and eliminate duplicate lines
    dedup_process_return = subprocess.call(
            r'''sort %s | uniq > %s'''
            % (prefasta_filename, deduped_fasta_filename), shell=True
        )
    if dedup_process_return != 0:
        raise RuntimeError('Problem encountered deduplicating FASTA reference')
    if verbose:
        print >>sys.stderr, 'Done! Writing final FASTA.'
    with open(final_fasta_filename, 'w') as output_stream:
        with open(deduped_fasta_filename) as fasta_stream:
            for line in fasta_stream:
                rname, seq = line.strip().split('\t')
                print >>output_stream, rname
                output_stream.write(
                    '\n'.join([seq[i:i+80] for i 
                                in xrange(0, len(seq), 80)])
                )
                output_stream.write('\n')
    return final_fasta_filename, reads_filename, rname_filename

def create_index_from_reference_fasta(bowtie2_build_exe, fasta_file,
        index_basename):
    """ Creates Bowtie2 index from reference fasta.

        bowtie2_build_exe: Path to Bowtie2
        fasta_file: Path to reference FASTA to index
        index_basename: Path to index basename to be created

        Return value: return value of bowtie-build process
    """
    bowtie_build_process = subprocess.Popen(
                                [args.bowtie2_build_exe,
                                    fasta_file,
                                    index_basename],
                                stderr=sys.stderr,
                                stdout=sys.stderr
                            )
    bowtie_build_process.wait()
    return bowtie_build_process.returncode

class BowtieOutputThread(threading.Thread):
    """ Processes Bowtie alignments, emitting tuples for exons and introns. """
    
    def __init__(self, input_stream, rname_stream, return_set,
        output_stream=sys.stdout, verbose=False, report_multiplier=1.2):
        """ Constructor for BowtieOutputThread.

            input_stream: where to retrieve Bowtie's SAM output, typically a
                Bowtie process's stdout.
            rname_stream: where to retrieve RNAMES, each of which is in
                the form (number of qnames associated with rnames) + '\x1e'
                + ('\x1e'-separated list of valid rnames)
            return_set: 0 is added to this set if the thread finishes
                successfully
            output_stream: where to emit exon and intron tuples; typically,
                this is sys.stdout.
            verbose: True if alignments should occasionally be written 
                to stderr.
            report_multiplier: if verbose is True, the line number of an
                alignment written to stderr increases exponentially with base
                report_multiplier.
        """
        super(BowtieOutputThread, self).__init__()
        self.daemon = True
        self.input_stream = input_stream
        self.rname_stream = rname_stream
        self.return_set = return_set
        self.output_stream = output_stream
        self.verbose = verbose
        self.report_multiplier = report_multiplier

    def run(self):
        """ Prints raw SAM output.

            Overrides default method containing thread activity.

            No return value.
        """
        global _output_line_count
        reversed_complement_translation_table \
            = string.maketrans('ATCG', 'TAGC')
        next_report_line = 0
        i = 0
        qname_count = 0
        tokens = self.rname_stream.readline().strip().split('\x1e')
        qname_total, rnames = int(tokens[0]), set(tokens[1:])
        done = False
        for (qname,), xpartition in xstream(self.input_stream, 1):
            assert not done
            printed = False
            for rest_of_line in xpartition:
                i += 1
                flag = int(rest_of_line[0])
                if self.verbose and next_report_line == i:
                    print >>sys.stderr, \
                        'SAM output record %d: rdname="%s", flag=%d' \
                        % (i, qname, flag)
                    next_report_line = max(int(next_report_line
                        * self.report_multiplier), next_report_line + 1)
                rname = rest_of_line[1]
                if rname in rnames:
                    print >>self.output_stream, \
                        '\t'.join((qname,) + rest_of_line)
                    printed = True
                    _output_line_count += 1
            if flag & 4 or not printed:
                # This is an unmapped read
                if flag & 16:
                    seq_to_write = rest_of_line[8][::-1].translate(
                                    reversed_complement_translation_table
                                )
                    qual_to_write = rest_of_line[9][::-1]
                else:
                    seq_to_write = rest_of_line[8]
                    qual_to_write = rest_of_line[9]
                # Write only essentials; handle "formal" writing in next step
                print >>self.output_stream, ('%s\t4\t\x1c\t\x1c\t\x1c\t\x1c'
                                             '\t\x1c\t\x1c\t\x1c\t%s\t%s') % (
                                                    qname,
                                                    seq_to_write,
                                                    qual_to_write
                                                )
            qname_count += 1
            if qname_count == qname_total:
                tokens = self.rname_stream.readline().strip().split('\x1e')
                try:
                    qname_total, rnames = int(tokens[0]), set(tokens[1:])
                    qname_count = 0
                except ValueError:
                    done = True
        self.output_stream.flush()
        self.return_set.add(0)

def handle_temporary_directory(archive, temp_dir_path):
    """ Archives or deletes temporary directory.

        archive: directory name, including path, to which temp_dir_path should
            be renamed when script is complete; or None if temp_dir_path
            should be trashed.
        temp_dir_path: path of temporary directory for storing intermediate
            alignments; archived if archive is not None.

        No return value.
    """
    if archive is not None:
        # Rename temporary directory for archiving
        os.rename(temp_dir_path, archive)
    else:
        # Kill temporary directory
        import shutil
        shutil.rmtree(temp_dir_path)

def go(input_stream=sys.stdin, output_stream=sys.stdout, bowtie2_exe='bowtie2',
    bowtie2_build_exe='bowtie2-build', bowtie2_args=None,
    temp_dir_path=tempfile.mkdtemp(), verbose=False, report_multiplier=1.2,
    replicable=False, count_multiplier=6):
    """ Runs Rail-RNA-realign.

        Uses Bowtie index including only sequences framing introns to align
        only those reads for which Bowtie did not report at least one alignment
        in Rail-RNA-align. Referencenames in this index encode intron sizes and
        locations in the (presmably) exonic sequences it records. Infers exonic
        chunks and introns from alignments.

        THIS CODE ASSUMES THAT BOWTIE RUNS ON JUST ONE THREAD; it exploits that
        Bowtie returns reads in the order in which they were sent, which is
        guaranteed only when on a single thread.

        Input (read from stdin)
        ----------------------------
        Tab-delimited input tuple columns:
        (two kinds)

        Type 1:
        1. Read sequence
        2. '\x1c' + FASTA reference name including '>'. The following format is
            used: original RNAME + '+' or '-' indicating which strand is the
            sense strand + '\x1d' + start position of sequence + '\x1d'
            + comma-separated list of subsequence sizes framing introns
            + '\x1d' + comma-separated list of intron sizes) + '\x1d'
            + 'p' if derived from primary alignment to genome; 's' if derived
            from secondary alignment to genome; 'i' if derived from cointron
            search
        3. FASTA sequence

        Type 2:
        1. Read sequence
        2. QNAME
        3. Quality sequence

        Type 1 corresponds to a FASTA line to index to which the read sequence
        is predicted to align. Type 2 corresponds to a distinct read. Input is
        partitioned by field 1.

        Hadoop output (written to stdout)
        ----------------------------

        Tab-delimited output tuple columns:
        Standard 11+ -column SAM output.

        ALL OUTPUT COORDINATES ARE 1-INDEXED.

        input_stream: where to find input reads.
        output_stream: where to emit exonic chunks and introns.
        bowtie2_exe: filename of Bowtie executable; include path if not in
            $PATH.
        bowtie2_build_exe: path to bowtie2-build executable
        bowtie2_args: string containing precisely extra command-line arguments
            to pass to Bowtie2
        temp_dir_path: path of temporary directory for storing intermediate
            alignments
        verbose: True iff more informative messages should be written to
            stderr.
        report_multiplier: if verbose is True, the line number of an alignment,
            read, or first readlet of a read written to stderr increases
            exponentially with base report_multiplier.
        replicable: use bowtie2's -a to ensure results are reproducible.
        count_multiplier: the bowtie2 -k parameter used is
            alignment_count_to_report * count_multiplier, where
            alignment_count_to_report is the user-specified bowtie2 -k arg.
            Ignored if replicable = True.

        No return value.
    """
    start_time = time.time()
    fasta_file, reads_file, rnames_file = input_files_from_input_stream(
                                                input_stream,
                                                verbose=verbose,
                                                temp_dir_path=temp_dir_path
                                            )
    bowtie2_index_base = os.path.join(temp_dir_path, 'tempidx')
    bowtie_build_return_code = create_index_from_reference_fasta(
                                    bowtie2_build_exe,
                                    fasta_file,
                                    bowtie2_index_base)
    if bowtie_build_return_code == 0:
        output_file = os.path.join(temp_dir_path, 'out.sam')
        alignment_count_to_report, _, _ \
            = bowtie.parsed_bowtie_args(bowtie2_args)
        bowtie_command = ' ' .join([bowtie2_exe,
            bowtie2_args if bowtie2_args is not None else '',
            '{0} --local -t --no-hd --mm -x'.format(
            ('-a' if replicable else 
            ('-k {0}'.format(alignment_count_to_report * count_multiplier)))),
            bowtie2_index_base, '--12', reads_file, '-S', output_file])
        print >>sys.stderr, 'Starting Bowtie2 with command: ' + bowtie_command
        bowtie_process = subprocess.Popen(bowtie_command, bufsize=-1,
            stdout=subprocess.PIPE, stderr=sys.stderr, shell=True)
        bowtie_process.wait()
        if os.path.exists(output_file):
            return_set = set()
            output_thread = BowtieOutputThread(
                                open(output_file),
                                open(rnames_file),
                                return_set=return_set,
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
    elif bowtie_build_return_code == 1:
        print >>sys.stderr, ('Bowtie build failed, but probably because '
                             'FASTA file was empty. Continuing...')
    else:
        raise RuntimeError('Bowtie build process failed with exitlevel %d.'
                            % bowtie_build_return_code)
    print >> sys.stderr, 'DONE with realign_reads.py; in/out=%d/%d; ' \
        'time=%0.3f s' % (_input_line_count, _output_line_count,
                                time.time() - start_time)

if __name__ == '__main__':
    # Print file's docstring if -h is invoked
    parser = argparse.ArgumentParser(description=__doc__, 
                formatter_class=argparse.RawDescriptionHelpFormatter)
    # Add command-line arguments
    parser.add_argument('--report_multiplier', type=float, required=False,
        default=1.2,
        help='When --verbose is also invoked, the only lines of lengthy '
             'intermediate output written to stderr have line number that '
             'increases exponentially with this base')
    parser.add_argument(\
        '--keep-alive', action='store_const', const=True, default=False,
        help='Prints reporter:status:alive messages to stderr to keep EMR '
             'task alive')
    parser.add_argument('--verbose', action='store_const', const=True,
        default=False,
        help='Print out extra debugging statements')
    parser.add_argument('--replicable', action='store_const',
        const=True,
        default=False, 
        help='Ensures that results are completely reproducible across '
             'different cluster configurations.')
    parser.add_argument('--count-multiplier', type=int, required=False,
        default=4,
        help='User-specified bowtie2 -k parameter is multiplied by this '
             'value when enumerating alignments')
    parser.add_argument(\
        '--stranded', action='store_const', const=True, default=False,
        help='Assume input reads come from the sense strand; then partitions '
             'in output have terminal + and - indicating sense strand')
    parser.add_argument('--test', action='store_const', const=True,
        default=False,
        help='Run unit tests; DOES NOT NEED INPUT FROM STDIN, AND DOES NOT '
             'WRITE EXONS AND INTRONS TO STDOUT')
    parser.add_argument('--archive', metavar='PATH', type=str, 
        default=None,
        help='Save output and Bowtie command to a subdirectory (named using ' 
             'this process\'s PID) of PATH')

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
    temp_dir_path = tempfile.mkdtemp()
    archive = os.path.join(args.archive,
        str(os.getpid())) if args.archive is not None else None
    # Handle temporary directory if CTRL+C'd
    atexit.register(handle_temporary_directory, archive, temp_dir_path)
    if args.verbose:
        print >>sys.stderr, 'Creating temporary directory %s' \
            % temp_dir_path
    go(bowtie2_exe=args.bowtie2_exe,
        bowtie2_build_exe=args.bowtie2_build_exe,
        bowtie2_args=bowtie_args,
        temp_dir_path=temp_dir_path,
        verbose=args.verbose, 
        report_multiplier=args.report_multiplier,
        replicable=args.replicable,
        count_multiplier=args.count_multiplier)
elif __name__ == '__main__':
    # Test units
    del sys.argv[1:] # Don't choke on extra command-line parameters
    import unittest
    unittest.main()