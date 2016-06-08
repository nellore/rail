#!/usr/bin/env python
"""
Rail-RNA-realign_reads

Follows Rail-RNA-cojunction_fasta, Rail-RNA-align_reads
Precedes Rail-RNA-compare_alignments

Realignment script for MapReduce pipelines that wraps Bowtie2. Creates Bowtie2
indexes including only sequences framing introns to align only those reads for
which Bowtie2 did not report alignments in Rail-RNA-align. Each group of
input reads (specified by first field below) is associated with a distinct
set of transcript fragments to which they are aligned. Reference names in
the index encode intron sizes and locations in the (presmably) exonic
sequences it records. Exonic chunks and junctions are inferred from alignments
in the next step.

Input (read from stdin)
----------------------------
Tab-delimited input tuple columns:
(two kinds)

Type 1:
1. Transcriptome Bowtie 2 index group number
2. Read sequence
3. '0' + FASTA reference name including '>'. The following format is used:
    original RNAME + '+' or '-' indicating which strand is the sense strand
    + '\x1d' + start position of sequence + '\x1d' + comma-separated list of
    subsequence sizes framing introns + '\x1d' + comma-separated list of intron
    sizes) + '\x1d' + 'p' if derived from primary alignment to genome; 's' if
    derived from secondary alignment to genome; 'i' if derived from cojunction
    search
4. FASTA sequence

Type 2:
1. Transcriptome Bowtie 2 index group number
2. Read sequence
3. 2 if SEQ is reverse-complemented, else 1
4. QNAME
5. QUAL

Type 1 corresponds to a FASTA line to index to which the read sequence is
predicted to align. Type 2 corresponds to a distinct read. Input is partitioned
by field 1 and sorted by field 2.

Hadoop output (written to stdout)
----------------------------
Tab-delimited output tuple columns:
Standard 11+ -column SAM output.

ALL OUTPUT COORDINATES ARE 1-INDEXED.
"""
import sys
import os
import site
import tempfile
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

from dooplicity.tools import xstream, register_cleanup, xopen, \
    make_temp_dir
import bowtie
import argparse
import tempdel
import itertools

# Initialize global variable for tracking number of input lines
_input_line_count = 0

_reversed_complement_translation_table = string.maketrans('ATCG', 'TAGC')

def input_files_from_input_stream(input_stream,
                                    output_stream,
                                    temp_dir_path=None,
                                    verbose=False,
                                    gzip_level=3):
    """ Generates FASTA reference to index and file with reads.

        Each line of the read file is in the following format:

        read number <TAB> SEQ <TAB> QUAL

        input_stream: where to find Hadoop input
        output_stream: where to write unmapped reads
        temp_dir_path: where to store files
        verbose: output extra debugging messages
        gzip_level: gzip compression level (0-9)

        Yield value: tuple (path to FASTA reference file, path to read file)
    """
    global _input_line_count
    if temp_dir_path is None: temp_dir_path = tempfile.mkdtemp()
    prefasta_filename = os.path.join(temp_dir_path, 'temp.prefa')
    deduped_fasta_filename = os.path.join(temp_dir_path, 'temp.deduped.prefa')
    final_fasta_filename = os.path.join(temp_dir_path, 'temp.fa')
    reads_filename = os.path.join(temp_dir_path, 'reads.temp.gz')
    for (counter, ((index_group,), xpartition)) in enumerate(
                                                    xstream(input_stream, 1)
                                                ):
        if verbose:
            print >>sys.stderr, (
                        'Group %d: Writing prefasta and input reads...'
                        % counter
                    )
        with open(prefasta_filename, 'w') as fasta_stream:
            with xopen(True, reads_filename, 'w') as read_stream:
                for read_seq, values in itertools.groupby(xpartition, 
                                                    key=lambda val: val[0]):
                    fasta_printed = False
                    for value in values:
                        _input_line_count += 1
                        if value[1][0] == '0':
                            # Print FASTA line
                            print >>fasta_stream, '\t'.join([value[1][1:-2],
                                                                value[2]])
                            fasta_printed = True
                        elif fasta_printed:
                            '''Add to temporary seq stream only if an
                            associated FASTA line was found.'''
                            if value[1] == '1':
                                print >>read_stream, '\t'.join([value[2],
                                                                    read_seq,
                                                                    value[3]])
                            else:
                                print >>read_stream, '\t'.join([
                                            value[2],
                                            read_seq[::-1].translate(
                                        _reversed_complement_translation_table
                                    ),
                                            value[3][::-1]])
                        else:
                            # Print unmapped read
                            if value[1] == '1':
                                seq_to_write = read_seq
                                qual_to_write = value[3]
                            else:
                                seq_to_write = read_seq[::-1].translate(
                                        _reversed_complement_translation_table
                                    )
                                qual_to_write = value[3][::-1]
                            '''Write only essentials; handle "formal" writing
                            in next step.'''
                            output_stream.write(
                                        '%s\t4\t\x1c\t\x1c\t\x1c\t\x1c'
                                        '\t\x1c\t\x1c\t\x1c\t%s\t%s\n' % (
                                                                value[2],
                                                                seq_to_write,
                                                                qual_to_write
                                                            )
                                    )
        if verbose:
            print >>sys.stderr, (
                        'Group %d: Done! Sorting and deduplicating prefasta...'
                        % counter
                    )
        # Sort prefasta and eliminate duplicate lines
        dedup_process_return = subprocess.call(
                r'''sort %s | uniq >%s'''
                % (prefasta_filename, deduped_fasta_filename), shell=True,
                executable='/bin/bash'
            )
        if dedup_process_return != 0:
            raise RuntimeError(
                    'Problem encountered deduplicating FASTA reference'
                )
        if verbose:
            print >>sys.stderr, (
                    'Group %d Done! Writing final FASTA.' % counter
                )
        with open(final_fasta_filename, 'w') as final_fasta_stream:
            with open(deduped_fasta_filename) as fasta_stream:
                for line in fasta_stream:
                    rname, seq = line.strip().split('\t')
                    print >>final_fasta_stream, rname
                    final_fasta_stream.write(
                        '\n'.join([seq[i:i+80] for i 
                                    in xrange(0, len(seq), 80)])
                    )
                    final_fasta_stream.write('\n')
        os.remove(deduped_fasta_filename)
        os.remove(prefasta_filename)
        yield final_fasta_filename, reads_filename

def create_index_from_reference_fasta(bowtie2_build_exe, fasta_file,
        index_basename):
    """ Creates Bowtie2 index from reference fasta.

        bowtie2_build_exe: Path to Bowtie2
        fasta_file: Path to reference FASTA to index
        index_basename: Path to index basename to be created

        Return value: return value of bowtie-build process
    """
    with open(os.devnull) as null_stream:
        bowtie_build_process = subprocess.Popen(
                                    [args.bowtie2_build_exe,
                                        fasta_file,
                                        index_basename],
                                    stderr=null_stream,
                                    stdout=null_stream
                                )
    bowtie_build_process.wait()
    return bowtie_build_process.returncode

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
    temp_dir_path=None, verbose=False, report_multiplier=1.2, gzip_level=3,
    count_multiplier=4, tie_margin=0):
    """ Runs Rail-RNA-realign.

        Realignment script for MapReduce pipelines that wraps Bowtie2. Creates
        Bowtie2 indexes including only sequences framing introns to align only
        those reads for which Bowtie2 did not report alignments in
        Rail-RNA-align. Each group of input reads (specified by first field
        below) is associated with a distinct set of transcript fragments to
        which they are aligned. Reference names in the index encode intron
        sizes and locations in the (presmably) exonic sequences it records.
        Exonic chunks and junctions are inferred from alignments in the next 
        step.

        Input (read from stdin)
        ----------------------------
        Tab-delimited input tuple columns:
        (two kinds)

        Type 1:
        1. Transcriptome Bowtie 2 index group number
        2. Read sequence
        3. '\x1c' + FASTA reference name including '>'. The following format is
            used: original RNAME + '+' or '-' indicating which strand is the
            sense strand + '\x1d' + start position of sequence + '\x1d'
            + comma-separated list of subsequence sizes framing junctions
            + '\x1d' + comma-separated list of intron sizes) + '\x1d'
            + 'p' if derived from primary alignment to genome; 's' if derived
            from secondary alignment to genome; 'i' if derived from cojunction
            search
        4. FASTA sequence

        Type 2:
        1. Transcriptome Bowtie 2 index group number
        2. Read sequence
        3. 1 if SEQ is reverse-complemented, else 0
        4. QNAME
        5. QUAL

        Type 1 corresponds to a FASTA line to index to which the read sequence
        is predicted to align. Type 2 corresponds to a distinct read. Input is
        partitioned by field 1 and sorted by field 2.

        Hadoop output (written to stdout)
        ----------------------------

        Tab-delimited output tuple columns:
        Standard 11+ -column SAM output.

        ALL OUTPUT COORDINATES ARE 1-INDEXED.

        input_stream: where to find input reads.
        output_stream: where to emit exonic chunks and junctions.
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
        gzip_level: level of gzip compression to use for some temporary files
        count_multiplier: the bowtie2 -k parameter used is
            alignment_count_to_report * count_multiplier, where
            alignment_count_to_report is the user-specified bowtie2 -k arg
        tie_margin: allowed score difference per 100 bases among ties in 
             max alignment score.

        No return value.
    """
    start_time = time.time()
    if temp_dir_path is None: temp_dir_path = tempfile.mkdtemp()
    bowtie2_index_base = os.path.join(temp_dir_path, 'tempidx')
    alignment_count_to_report, _, _ \
            = bowtie.parsed_bowtie_args(bowtie2_args)
    reads_filename = os.path.join(temp_dir_path, 'reads.temp')
    input_command = 'gzip -cd %s' % reads_filename
    bowtie_command = ' ' .join([bowtie2_exe,
        bowtie2_args if bowtie2_args is not None else '',
        '{0} --local -t --no-hd --mm -x'.format(
                '-k {0}'.format(alignment_count_to_report * count_multiplier)
            ),
        bowtie2_index_base, '--12 -'])
    delegate_command = ''.join(
            [sys.executable, ' ', os.path.realpath(__file__)[:-3],
                ('_delegate.py --report-multiplier %08f '
                 '--alignment-count-to-report %d '
                 '--tie-margin %d %s')
                    % (report_multiplier, alignment_count_to_report,
                        tie_margin, '--verbose' if verbose else '')]
        )
    # Use grep to kill empty lines terminating python script
    full_command = ' | '.join([input_command, 
                                bowtie_command, delegate_command])
    print >>sys.stderr, 'Bowtie2 command to execute: ' + full_command
    for fasta_file, reads_file in input_files_from_input_stream(
                                                input_stream,
                                                output_stream,
                                                verbose=verbose,
                                                temp_dir_path=temp_dir_path,
                                                gzip_level=gzip_level
                                            ):
        bowtie_build_return_code = create_index_from_reference_fasta(
                                        bowtie2_build_exe,
                                        fasta_file,
                                        bowtie2_index_base
                                    )
        if bowtie_build_return_code == 0:
            try:
                os.remove(fasta_file)
            except OSError:
                pass
            bowtie_process = subprocess.Popen(' '.join(
                        ['set -exo pipefail;', full_command]
                    ), bufsize=-1,
                stdout=sys.stdout, stderr=sys.stderr, shell=True,
                executable='/bin/bash')
            return_code = bowtie_process.wait()
            if return_code:
                raise RuntimeError(
                            'Error occurred while reading Bowtie 2 output; '
                            'exitlevel was %d.' % return_code
                        )
        elif bowtie_build_return_code == 1:
            print >>sys.stderr, ('Bowtie build failed, but probably because '
                                 'FASTA file was empty. Continuing...')
        else:
            raise RuntimeError('Bowtie build process failed with exitlevel %d.'
                                % bowtie_build_return_code)

    print >>sys.stderr, 'DONE with realign_reads.py; in=%d; ' \
        'time=%0.3f s' % (_input_line_count, time.time() - start_time)

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
    parser.add_argument(
        '--keep-alive', action='store_const', const=True, default=False,
        help='Prints reporter:status:alive messages to stderr to keep EMR '
             'task alive')
    parser.add_argument('--verbose', action='store_const', const=True,
        default=False,
        help='Print out extra debugging statements')
    parser.add_argument('--count-multiplier', type=int, required=False,
        default=4,
        help='User-specified bowtie2 -k parameter is multiplied by this '
             'value when enumerating alignments')
    parser.add_argument(
        '--stranded', action='store_const', const=True, default=False,
        help='Assume input reads come from the sense strand; then partitions '
             'in output have terminal + and - indicating sense strand')
    parser.add_argument('--test', action='store_const', const=True,
        default=False,
        help='Run unit tests; DOES NOT NEED INPUT FROM STDIN, AND DOES NOT '
             'WRITE EXONS AND JUNCTIONS TO STDOUT')
    parser.add_argument('--archive', metavar='PATH', type=str, 
        default=None,
        help='Save output and Bowtie command to a subdirectory (named using ' 
             'this process\'s PID) of PATH')
    parser.add_argument('--gzip-level', type=int, required=False,
        default=3,
        help='Level of gzip compression to use, if applicable')

    # Add command-line arguments for dependencies
    bowtie.add_args(parser)
    tempdel.add_args(parser)
    from alignment_handlers import add_args as alignment_handlers_add_args
    alignment_handlers_add_args(parser)

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
    temp_dir_path = make_temp_dir(tempdel.silentexpandvars(args.scratch))
    archive = os.path.join(args.archive,
        str(os.getpid())) if args.archive is not None else None
    # Handle temporary directory if CTRL+C'd
    register_cleanup(handle_temporary_directory, archive, temp_dir_path)
    if args.verbose:
        print >>sys.stderr, 'Creating temporary directory %s' \
            % temp_dir_path
    go(bowtie2_exe=os.path.expandvars(args.bowtie2_exe),
        bowtie2_build_exe=os.path.expandvars(args.bowtie2_build_exe),
        bowtie2_args=bowtie_args,
        temp_dir_path=temp_dir_path,
        verbose=args.verbose, 
        report_multiplier=args.report_multiplier,
        gzip_level=args.gzip_level,
        count_multiplier=args.count_multiplier,
        tie_margin=args.tie_margin)
elif __name__ == '__main__':
    # Test units
    del sys.argv[1:] # Don't choke on extra command-line parameters
    import unittest
    unittest.main()