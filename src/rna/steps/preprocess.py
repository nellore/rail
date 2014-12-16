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
Tab-separated fields:
---If URL is local:
1. #!splitload
2. \x1d-separated list of 0-based indexes of reads at which to start
    each new file
3. \x1d-separated list of numbers of reads to include in gzipped files
4. \x1d-separated list of manifest lines whose tabs are replaced by \x1es

---Otherwise:
manifest line

A manifest line has the following format

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
  1. Nucleotide sequence or its reversed complement, whichever is first in 
    alphabetical order
  2. 1 if sequence was reverse-complemented else 0
  3. Name
  4. Quality sequence or its reverse, whichever corresponds to field 1

Format 2 (paired, 2 lines, 3 columns each)
(so this is the same as single-end)
  1. Nucleotide sequence for mate 1 or its reversed complement, whichever is
    first in alphabetical order
  2. 1 if sequence was reverse-complemented else 0
  3. Name for mate 1
  4. Quality sequence for mate 1 or its reverse, whichever corresponds to
    field 1
    
    (new line)

  1. Nucleotide sequence for mate 2 or its reversed complement, whichever is
    first in alphabetical order
  2. 1 if sequence was reverse complemented else 0
  3. Name for mate 2
  4. Quality sequence for mate 2 or its reverse, whichever corresponds to
    field 1

Quality sequences are strings of Is for FASTA input.

Files are gzipped by default.
"""

import os
import sys
import atexit
import tempfile
import site
import string

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

_reversed_complement_translation_table = string.maketrans('ATCG', 'TAGC')

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
        to_stdout=False, push='.', mover=filemover.FileMover()):
    """ Runs Rail-RNA-preprocess

        Input (read from stdin)
        ----------------------------
        Tab-separated fields:
        ---If URL is local:
        1. #!splitload
        2. \x1d-separated list of 0-based indexes of reads at which to start
            each new file
        3. \x1d-separated list of numbers of reads to include in gzipped files
        4. \x1d-separated list of manifest lines whose tabs are replaced by
            \x1es

        ---Otherwise:
        manifest line

        A manifest line has the following format

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
          1. Nucleotide sequence or its reversed complement, whichever is first
            in alphabetical order
          2. 1 if sequence was reverse-complemented else 0
          3. Name
          4. Quality sequence or its reverse, whichever corresponds to field 1

        Format 2 (paired, 2 lines, 3 columns each)
        (so this is the same as single-end)
          1. Nucleotide sequence for mate 1 or its reversed complement,
            whichever is first in alphabetical order
          2. 1 if sequence was reverse-complemented else 0
          3. Name for mate 1
          4. Quality sequence for mate 1 or its reverse, whichever corresponds
            to field 1
            
            (new line)

          1. Nucleotide sequence for mate 2 or its reversed complement,
            whichever is first in alphabetical order
          2. 1 if sequence was reverse complemented else 0
          3. Name for mate 2
          4. Quality sequence for mate 2 or its reverse, whichever corresponds
            to field 1

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
    if not to_stdout:
        push_url = Url(push)
        if push_url.is_local:
            destination = push
        elif push_url.is_s3 or push_url.is_hdfs:
            destination = temp_dir
        else:
            raise RuntimeError('Push destination must be '
                               'on S3, HDFS, or local.')
    fastq_cues = set(['@'])
    fasta_cues = set(['>', ';'])
    source_dict = {}
    for line in sys.stdin:
        if not line.strip() or line[0] == '#': continue
        _input_line_count += 1
        # Kill offset from start of manifest file
        tokens = line.strip().split('\t')[1:]
        token_count = len(tokens)
        if token_count == 4 and tokens[0] == '#!splitload':
            '''Line specifies precisely how records from files should be
            placed.'''
            assert not to_stdout, ('Split manifest line inconsistent with '
                                   'writing to stdout.')
            indexes = tokens[1].split('\x1d')
            read_counts = tokens[2].split('\x1d')
            manifest_lines = [token.split('\x1e')
                                for token in tokens[3].split('\x1d')]
            assert len(indexes) == len(read_counts) == len(manifest_lines)
            for i, manifest_line in enumerate(manifest_lines):
                manifest_line_field_count = len(manifest_line)
                if manifest_line_field_count == 3:
                    source_dict[(Url(manifest_line[0]),)] = \
                        (manifest_line[-1],
                            int(indexes[i]), int(read_counts[i]))
                else:
                    assert manifest_line_field_count == 5
                    source_dict[(Url(manifest_line[0]),
                                 Url(manifest_line[2]))] = (
                                                        manifest_line[-1],
                                                        int(indexes[i]),
                                                        int(read_counts[i])
                                                    )
        elif token_count == 3:
            # Single-end reads
            source_dict[(Url(tokens[0]),)] = (tokens[-1],)
        elif token_count == 5:
            # token_count == 5
            source_dict[(Url(tokens[0]), Url(tokens[2]))] = (tokens[-1],)
        else:
            # Not a valid line, but continue for robustness
            continue
    file_number = 0
    for source_urls in source_dict:
        sample_label = source_dict[source_urls][0]
        downloaded = set()
        sources = []
        records_printed = 0
        if len(source_dict[source_urls]) == 3:
            skip_count = source_dict[source_urls][1]
            if len(source_urls) == 2:
                records_to_consume = source_dict[source_urls][2]
                if skip_count % 2:
                    skip_count -= 1
                    records_to_consume += 1
                if records_to_consume % 2:
                    records_to_consume -= 1
            else:
                records_to_consume = source_dict[source_urls][2]
        else:
            skip_count = 0
            records_to_consume = None # Consume all records
        assert (records_to_consume >= 0 or records_to_consume is None), (
                'Negative value %d of records to consume encountered.'
            ) % records_to_consume
        if records_to_consume == 0: continue
        skipped = False
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
            break_outer_loop = False
            while True:
                if not to_stdout:
                    '''Name files using Hadoop task environment property
                    mapred.task.partition.'''
                    if gzip_output:
                        try:
                            output_file = os.path.join(
                                    destination, 
                                    '.'.join([
                                        os.environ['mapred_task_partition'],
                                        str(file_number), 'gz'
                                    ])
                                )
                        except KeyError:
                            '''Hadoop 2.x: mapreduce.task.partition; see 
                            http://hadoop.apache.org/docs/r2.0.3-alpha/
                            hadoop-project-dist/hadoop-common/
                            DeprecatedProperties.html.'''
                            output_file = os.path.join(
                                    destination, 
                                    '.'.join([
                                        os.environ['mapreduce_task_partition'],
                                        str(file_number), 'gz'
                                    ])
                                )
                        open_args = [output_file, 'a', gzip_level]
                    else:
                        try:
                            output_file = os.path.join(
                                    destination, 
                                    '.'.join([
                                        os.environ['mapred_task_partition'],
                                        str(file_number)
                                    ])
                                )
                        except KeyError:
                            output_file = os.path.join(
                                    destination, 
                                    '.'.join([
                                        os.environ['mapreduce_task_partition'],
                                        str(k), str(file_number)
                                    ])
                                )
                        open_args = [output_file, 'a']
                    try:
                        os.makedirs(os.path.dirname(output_file))
                    except OSError:
                        pass
                else:
                    open_args = []
                '''Use xopen to handle compressed streams and normal streams
                generally.'''
                with xopen(gzip_output if not to_stdout else '-', *open_args) \
                    as output_stream:
                    perform_push = False
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
                            if records_to_consume and not skipped:
                                '''Skip lines as necessary; for paired-end
                                reads skip the largest even number of records 
                                less than records_to_consume.'''
                                if len(source_urls) == 1:
                                    # single-end
                                    line_skip_count = max(
                                            skip_count * 4 - 1, 0
                                        )
                                else:
                                    # paired-end
                                    line_skip_count = max(
                                            ((skip_count / 2) * 4 - 1), 0
                                        )
                                    for _ in xrange(line_skip_count):
                                        next(source_stream_2)
                                for _ in xrange(line_skip_count):
                                    next(source_stream_1)
                                if skip_count:
                                    lines = []
                                    for source_stream in source_streams:
                                        lines.append(source_stream.readline())
                                    if not lines[0]:
                                        break_outer_loop = True
                                        break
                                    lines = [line.strip() for line in lines]
                                skipped = True
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
                            if records_to_consume and not skipped:
                                '''Skip lines as necessary; for paired-end
                                reads skip the largest even number of records 
                                less than records_to_consume.'''
                                if len(source_urls) == 1:
                                    # single-end
                                    line_skip_count = max(
                                            skip_count * 2 - 1, 0
                                        )
                                else:
                                    # paired-end
                                    line_skip_count = max(
                                            ((skip_count / 2) * 2 - 1), 0
                                        )
                                    for _ in xrange(line_skip_count):
                                        next(source_stream_2)
                                for _ in xrange(line_skip_count):
                                    next(source_stream_1)
                                if skip_count:
                                    lines = []
                                    for source_stream in source_streams:
                                        lines.append(source_stream.readline())
                                    if not lines[0]:
                                        break_outer_loop = True
                                        break
                                    lines = [line.strip() for line in lines]
                                skipped = True
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
                            seqs = [seq.upper() for seq in seqs]
                            reversed_complement_seqs = [
                                    seqs[0][::-1].translate(
                                        _reversed_complement_translation_table
                                    ),
                                    seqs[1][::-1].translate(
                                        _reversed_complement_translation_table
                                    )
                                ]
                            if seqs[0] < reversed_complement_seqs[0]:
                                left_seq = seqs[0]
                                left_qual = quals[0]
                                left_reversed = '0'
                            else:
                                left_seq = reversed_complement_seqs[0]
                                left_qual = quals[0][::-1]
                                left_reversed = '1'
                            if seqs[1] < reversed_complement_seqs[1]:
                                right_seq = seqs[1]
                                right_qual = quals[1]
                                right_reversed = '0'
                            else:
                                right_seq = reversed_complement_seqs[1]
                                right_qual = quals[1][::-1]
                                right_reversed = '1'
                            print >>output_stream, '\t'.join(
                                        [
                                            left_seq,
                                            left_reversed,
                                            qname_from_read(
                                                    original_qnames[0],
                                                    seqs[0] + quals[0], 
                                                    sample_label
                                                ),
                                            '\n'.join([left_qual, right_seq]),
                                            right_reversed,
                                            qname_from_read(
                                                    original_qnames[1],
                                                    seqs[1] + quals[1], 
                                                    sample_label
                                                ),
                                            right_qual
                                        ]
                                    )
                            records_printed += 2
                        else:
                            seqs[0] = seqs[0].upper()
                            reversed_complement_seqs = [
                                    seqs[0][::-1].translate(
                                        _reversed_complement_translation_table
                                    )
                                ]
                            # Single-end write
                            if seqs[0] < reversed_complement_seqs[0]:
                                seq = seqs[0]
                                qual = quals[0]
                                is_reversed = '0'
                            else:
                                seq = reversed_complement_seqs[0]
                                qual = quals[0][::-1]
                                is_reversed = '1'
                            print >>output_stream, '\t'.join(
                                        [
                                            seq,
                                            is_reversed,
                                            qname_from_read(
                                                original_qnames[0],
                                                seqs[0] + quals[0], 
                                                sample_label
                                            ),
                                            qual
                                        ]
                                    )
                            records_printed += 1
                        _output_line_count += 1
                        for seq in seqs:
                            nucs_read += len(seq)
                        if records_printed == records_to_consume:
                            break_outer_loop = True
                            perform_push = True
                            break
                        if not to_stdout and not records_to_consume and \
                            nucs_read > nucleotides_per_input:
                            file_number += 1
                            break
                if verbose:
                    print >>sys.stderr, (
                            'Exited with statement; line numbers are %s' 
                            % line_numbers
                        )
                if (not to_stdout) and (push_url.is_s3 or push_url.is_hdfs) \
                    and ((not records_to_consume) or
                         (records_to_consume and perform_push)):
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
            if verbose:
                print >>sys.stderr, 'Exiting source streams...'
        if verbose:
            print >>sys.stderr, 'Exited source streams.'
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
    parser.add_argument('--stdout', action='store_const', const=True,
        default=False,
        help='Write to stdout rather than to gzip files; useful in hadoop ' \
             'modes when using LZO output formats')
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
        to_stdout=args.stdout,
        push=args.push,
        mover=mover)
    print >>sys.stderr, 'DONE with preprocess.py; in/out=%d/%d; ' \
        'time=%0.3f s' % (_input_line_count, _output_line_count,
                            time.time() - start_time)