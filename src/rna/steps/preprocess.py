#!/usr/bin/env python
"""
Rail-RNA-preprocess
START: follows no step
Precedes Rail-RNA-align_reads

Map step in MapReduce pipelines converts paired- and single-end RNA-seq
input data in FASTA or FASTQ format to single-line format parseable in
future steps.

Input (read from stdin)
----------------------------
Myrna-style manifest file, each of whose lines is in one of the following
formats:

(for single-end reads)
<URL>(tab)<Optional MD5>(tab)<Sample label>

(for paired-end reads)
<URL 1>(tab)<Optional MD5 1>(tab)<URL 2>(tab)<Optional MD5 2>(tab)
<Sample label>

Hadoop output (written to stdout)
----------------------------
None.

Other output (written to directory specified by command-line parameter --push)
____________________________
Files containing input data in one of the following formats:
Format 1 (single-end, 3-column):
  1. Name
  2. Nucleotide sequence
  3. Quality sequence
Format 2 (paired, 6-column):
  1. Name for mate 1
  2. Nucleotide sequence for mate 1
  3. Quality sequence for mate 1
  4. Name for mate 2
  5. Nucleotide sequence for mate 2
  6. Quality sequence for mate 2

Quality sequences are strings of Is for FASTA input.

Files are gzipped by default.
"""

import os
import sys
import atexit
import tempfile
import site

base_path = os.path.abspath(
                    os.path.dirname(os.path.dirname(os.path.dirname(
                        os.path.realpath(__file__)))
                    )
                )
utils_path = os.path.join(base_path, 'rna', 'utils')
site.addsitedir(utils_path)
site.addsitedir(base_path)

from dooplicity.ansibles import Url
from dooplicity.tools import xopen
import filemover
from tempdel import remove_temporary_directories

_input_line_count, _output_line_count = 0, 0

def qname_from_read(original_qname, seq, sample_label):
    """ Returns QNAME including sample label and ID formed from hash.

        New QNAME takes the form:
            <original_qname> + '\x1d' + 
            <short hash of original_qname + seq + sample label> + '\x1d' +
            <sample_label>

        The purpose of the ID is to distinguish reads with the same name in
        the original FASTA/FASTQ. This should never happen in well-formed
        inputs, but can happen in, for example, Flux-generated data.

        original_qname: original name of read
        seq: read sequence (+ qual sequence; doesn't matter)
        sample_label: sample label

        Return value: new QNAME string
    """
    return ''.join([original_qname,
                    '\x1d',
                    hex(hash(original_qname + seq + sample_label))[-12:],
                    '\x1d', sample_label])

def go(nucleotides_per_input=8000000, gzip_output=True, gzip_level=3,
        push='.', mover=filemover.FileMover()):
    """ Runs Rail-RNA-preprocess

        Input (read from stdin)
        ----------------------------
        Myrna-style manifest file, each of whose lines is in one of the
        following formats:

        (for single-end reads)
        <URL>(tab)<Optional MD5>(tab)<Sample label>

        (for paired-end reads)
        <URL 1>(tab)<Optional MD5 1>(tab)<URL 2>(tab)<Optional MD5 2>(tab)
        <Sample label>

        Hadoop output (written to stdout)
        ----------------------------
        None.

        Other output (written to directory specified by command-line parameter
            --push)
        ____________________________
        Files containing input data in one of the following formats:
        Format 1 (single-end, 3-column):
          1. Name
          2. Nucleotide sequence
          3. Quality sequence
        Format 2 (paired, 6-column):
          1. Name for mate 1
          2. Nucleotide sequence for mate 1
          3. Quality sequence for mate 1
          4. Name for mate 2
          5. Nucleotide sequence for mate 2
          6. Quality sequence for mate 2

        Quality sequences are strings of Is for FASTA input.

        nucleotides_per_input: maximum number of nucleotides to put in a given
            input file
        gzip_output: True iff preprocessed input should be gzipped
        gzip_level: level of gzip compression to use
        push: where to send output
        verbose: True iff extra debugging statements should be printed to
            stderr

        No return value
    """
    global _input_line_count, _output_line_count
    temp_dir = tempfile.mkdtemp()
    print >>sys.stderr, 'Created local destination directory "%s".' % temp_dir
    atexit.register(remove_temporary_directories, [temp_dir])
    input_line_count, output_line_count = 0, 0
    token_lengths = set([3, 5])
    push_url = Url(push)
    if push_url.is_local:
        destination = push
    elif push_url.is_s3:
        destination = temp_dir
    else:
        raise RuntimeError('Push destination must be on S3 or local.')
    fastq_cues = set(['@'])
    fasta_cues = set(['>', ';'])
    for k, line in enumerate(sys.stdin):
        if line[0] == '#' or not line.strip():
            # Allow comments and blank lines
            continue
        _input_line_count += 1
        # Kill offset from start of manifest file
        tokens = line.strip().split('\t')[1:]
        token_count = len(tokens)
        assert token_count in token_lengths
        if token_count == 3:
            # Single-end reads
            source_urls = [Url(tokens[0])]
        else:
            # token_count == 5
            source_urls = [Url(tokens[0]), Url(tokens[2])]
        sample_label = tokens[-1]
        downloaded = set()
        sources = []
        for source_url in source_urls:
            if not source_url.is_local:
                # Download
                print >>sys.stderr, 'Retrieving URL "%s"...' \
                    % source_url.to_url()
                mover.get(source_url, temp_dir)
                downloaded = list(
                        set(os.listdir(temp_dir)).difference(downloaded)
                    )
                sources.append(os.path.join(temp_dir, list(downloaded)[0]))
            else:
                sources.append(source_url.to_url())
        '''Use os.devnull so single- and paired-end data can be handled in one
        loop.'''
        if len(sources) == 1:
            sources.append(os.devnull)
        with xopen(None, sources[0]) as source_stream_1, \
            xopen(None, sources[1]) as source_stream_2:
            source_streams = [source_stream_1, source_stream_2]
            file_number = 0
            break_outer_loop = False
            while True:
                '''Name files using Hadoop task environment property
                mapred.task.partition.'''
                if gzip_output:
                    output_file = os.path.join(
                                destination, 
                                '.'.join([
                                    os.environ['mapred_task_partition'],
                                    str(k), str(file_number), 'gz'
                                ])
                            )
                    open_args = [output_file, 'w', gzip_level]
                else:
                    output_file = os.path.join(
                                destination, 
                                '.'.join([
                                    os.environ['mapred_task_partition'],
                                    str(k), str(file_number)
                                ])
                            )
                    open_args = [output_file, 'w']
                try:
                    os.makedirs(os.path.dirname(output_file))
                except OSError:
                    pass
                '''Use xopen to handle compressed streams and normal streams
                generally.'''
                with xopen(gzip_output, *open_args) as output_stream:
                    line_numbers = [0, 0]
                    read_next_line = True
                    nucs_read = 0
                    while True:
                        if read_next_line:
                            # Read next line only if FASTA mode didn't already
                            lines = []
                            for source_stream in source_streams:
                                lines.append(source_stream.readline())
                        read_next_line = True
                        if not lines[0]:
                            break_outer_loop = True
                            break
                        line_numbers = [i + 1 for i in line_numbers]
                        lines = [line.strip() for line in lines]
                        if lines[0][0] in fastq_cues:
                            seqs = [source_stream.readline().strip()
                                        for source_stream in source_streams]
                            line_numbers = [i + 1 for i in line_numbers]
                            plus_lines = [source_stream.readline().strip()
                                            for source_stream
                                            in source_streams]
                            line_numbers = [i + 1 for i in line_numbers]
                            assert plus_lines[0][0] == '+', (
                                    'Malformed read "%s" at line %d of '
                                    'file "%s".'
                                ) % (lines[0], line_numbers[0], sources[0])
                            if plus_lines[1]:
                                assert plus_lines[1][0] == '+', (
                                        'Malformed read "%s" at line %d of '
                                        'file "%s".'
                                    ) % (lines[1], line_numbers[1], sources[1])
                            try:
                                # Kill spaces in name
                                original_qnames = [line[1:].replace(' ', '_')
                                                    for line in lines]
                            except IndexError:
                                print >>sys.stderr, 'Error finding QNAME at ' \
                                    ('line %d of either %s or %s' % (
                                            sources[0],
                                            sources[1]
                                        ))
                                raise
                            quals = [source_stream.readline().strip()
                                        for source_stream in source_streams]
                            line_numbers = [i + 1 for i in line_numbers]
                            for i in xrange(2):
                                assert len(seqs[i]) == len(quals[i]), (
                                        'Length of read sequence does not '
                                        'match length of quality string at '
                                        'line %d of file "%s".'
                                    ) % (line_numbers[i], sources[i])
                        elif lines[0][0] in fasta_cues:
                            original_qnames, seqs, quals, lines = []*4
                            for i, source_stream in enumerate(source_streams):
                                next_line = source_stream.readline()
                                if not next_line: break
                                line_numbers[i] += 1
                                while next_line[0] == ';':
                                    # Skip comment lines
                                    next_line = source_stream.readline()
                                    line_numbers[i] += 1
                                try:
                                    # Kill spaces in name
                                    original_qnames.append(
                                            line[1:].replace(' ', '_')
                                        )
                                except IndexError:
                                    raise RuntimeError(
                                            ('No QNAME for read '
                                             'above line %d in file "%s".') % (
                                                            line_numbers[i],
                                                            sources[i]
                                                        ) 
                                        )
                                assert next_line not in fasta_cues, (
                                        'No read sequence for read named "%s" '
                                        'above line %d in file "%s".'
                                    ) % (
                                        original_qnames[i],
                                        line_numbers[i],
                                        sources[i]
                                    )
                                read_lines = []
                                while next_line[0] not in fasta_cues:
                                    read_lines.append(next_line.strip())
                                    next_line = source_stream.readline()
                                    line_numbers[i] += 1
                                seqs.append(''.join(read_lines))
                                quals.append('I'*len(read_lines))
                                lines.append(next_line)
                                read_next_line = False
                        if len(original_qnames) == 2 and original_qnames[1]:
                            # Paired-end write
                            if original_qnames[0] == original_qnames[1]:
                                # Add paired-end identifiers
                                original_qnames[0] += '/1'
                                original_qnames[1] += '/2'
                            assert seqs[1]
                            assert quals[1]
                            print >>output_stream, '\t'.join(
                                        [qname_from_read(
                                                original_qnames[0],
                                                seqs[0] + quals[0], 
                                                sample_label
                                            ),
                                        seqs[0],
                                        quals[0],
                                        qname_from_read(
                                                original_qnames[1],
                                                seqs[1] + quals[1], 
                                                sample_label
                                            ),
                                        seqs[1],
                                        quals[1]]
                                    )
                        else:
                            # Single-end write
                            print >>output_stream, '\t'.join(
                                        [qname_from_read(
                                                original_qnames[0],
                                                seqs[0] + quals[0], 
                                                sample_label
                                            ),
                                        seqs[0],
                                        quals[0]]
                                    )
                        _output_line_count += 1
                        for seq in seqs:
                            nucs_read += len(seq)
                        if nucs_read > nucleotides_per_input:
                            break
                if push_url.is_s3:
                    print >>sys.stderr, 'Pushing "%s" to "%s" ...' % (
                                                            output_file,
                                                            push_url.to_url()
                                                        )
                    mover.put(output_file, push_url.plus(os.path.basename(
                                                                output_file
                                                            )))
                    try:
                        os.remove(output_file)
                    except OSError:
                        pass
                if break_outer_loop: break
                file_number += 1
        # Clear temporary directory
        for input_file in os.listdir(temp_dir):
            try:
                os.remove(os.path.join(temp_dir, input_file))
            except OSError:
                pass

if __name__ == '__main__':
    import argparse
    # Print file's docstring if -h is invoked
    parser = argparse.ArgumentParser(description=__doc__, 
                formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--nucs-per-file', type=int, required=False,
        default=8000000,
        help='Write to next file if more than this many nucleotides are ' \
             'found to have been written to given output file')
    parser.add_argument('--gzip-output', action='store_const', const=True,
        default=False,
        help='Compress output files with gzip')
    parser.add_argument('--gzip-level', type=int, required=False,
        default=3,
        help='Level of gzip compression to use, if applicable')
    parser.add_argument('--push', type=str, required=False,
        default='.',
        help='Directory in which to push output files')
    filemover.add_args(parser)
    args = parser.parse_args()

if __name__ == '__main__':
    import time
    start_time = time.time()
    mover = filemover.FileMover(args=args)
    go(nucleotides_per_input=args.nucs_per_file,
        gzip_output=args.gzip_output,
        gzip_level=args.gzip_level,
        push=args.push,
        mover=mover)
    print >>sys.stderr, 'DONE with preprocess.py; in/out=%d/%d; ' \
        'time=%0.3f s' % (_input_line_count, _output_line_count,
                            time.time() - start_time)