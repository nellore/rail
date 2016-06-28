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
from dooplicity.tools import xopen, register_cleanup, make_temp_dir
from dooplicity.counters import Counter
import filemover
import tempdel
import subprocess
from guess import phred_converter
from encode import encode, encode_sequence

_reversed_complement_translation_table = string.maketrans('ATCG', 'TAGC')

_input_line_count, _output_line_count = 0, 0
counter = Counter('preprocessing')
register_cleanup(counter.flush)

# Maximum length of read to ignore if dealing with barcodes
_max_stubby_read_length = 10

'''Set Bowtie 2's default quality-based mismatch penalties; see
http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#bowtie2-options-mp
for more information.'''
_MX, _MN = 6, 2
'''This table depends on _MX, _MN above; it maps mismatch penalties to quality
scores, each an "exemplar" from a different bin.'''
_mismatch_penalties_to_quality_scores = string.maketrans('23456', '#05Hh')

def qname_from_read(qname, seq, sample_label, mate=None):
    """ Returns QNAME including sample label and ID formed from hash.

        New QNAME takes the form:
            <qname> + '\x1d'
            + <short hash of original_qname + seq + sample label> + 
            colon character + (base 36-encoded mate sequence if present else
            blank)
            '\x1d'
            + <sample_label>

        The purpose of the ID is to distinguish reads with the same name in
        the original FASTA/FASTQ. This should never happen in well-formed
        inputs, but can happen in, for example, Flux-generated data.

        qname: original name of read OR a shortened name
        seq: read sequence (+ qual sequence; doesn't matter)
        short_name: short name to ultimately give read if qname shortening
            is turned on
        sample_label: sample label
        mate: mate sequence if paired-end read; else None

        Return value: new QNAME string
    """
    return ''.join([qname, '\x1d',
                    hex(hash(qname + seq + sample_label))[-12:],
                    ':', encode_sequence(mate) if mate is not None else '',
                    '\x1d', sample_label])

def max_min_read_lengths_from_fastq_stream(fastq_stream):
    """ Grabs maximum and minimum read lengths from a FASTQ stream.

        fastq_stream: input fastq stream

        Return value: tuple (max read length, min read length)
    """
    max_len, min_len = 0, 1000000000
    for i, line in enumerate(fastq_stream):
        if (i - 1) % 4 == 0:
            max_len = max(max_len, len(line.strip()))
            min_len = min(min_len, len(line.strip()))
    return max_len, min_len

def go(nucleotides_per_input=8000000, gzip_output=True, gzip_level=3,
        to_stdout=False, push='.', mover=filemover.FileMover(),
        verbose=False, scratch=None, bin_qualities=True, short_qnames=False,
        skip_bad_records=False, workspace_dir=None,
        fastq_dump_exe='fastq-dump', ignore_missing_sra_samples=False):
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
        scratch: scratch directory for storing temporary files or None if 
            securely created temporary directory
        bin_qualities: True iff quality string should be binned according to
            rules in _mismatch_penalties_to_quality_scores
            and round_quality_string() defined in go()
        short_qnames: True iff original qname should be killed and a new qname
            should be written in a short base64-encoded format
        skip_bad_records: True iff bad records should be skipped; otherwise,
            raises exception if bad record is encountered
        workspace_dir: where to use fastq-dump -- needed for working with
            dbGaP data. None if temporary dir should be used.
        fastq_dump_exe: path to fastq-dump executable
        ignore_missing_sra_samples: does not return error if fastq-dump doesn't
            find a sample

        No return value
    """
    if bin_qualities:
        import math
        def round_quality_string(qual):
            """ Bins phred+33 quality string to improve compression.

                Uses 5-bin scheme that does not affect Bowtie 2 alignments

                qual: quality string

                Return value: "binned" quality string.
            """
            return ''.join(
                [str(int(
                    _MN + math.floor((_MX - _MN) * min(
                                                    ord(qual_char) - 33.0, 40.0
                                                ) / 40.0)
                        )) for qual_char in qual]).translate(
                                _mismatch_penalties_to_quality_scores
                            )
    else:
        def round_quality_string(qual):
            """ Leaves quality string unbinned and untouched.

                qual: quality string

                Return value: qual
            """
            return qual
    global _input_line_count, _output_line_count
    skip_stubs = False
    temp_dir = make_temp_dir(scratch)
    print >>sys.stderr, 'Created local destination directory "%s".' % temp_dir
    register_cleanup(tempdel.remove_temporary_directories, [temp_dir])
    input_line_count, output_line_count = 0, 0
    if not to_stdout:
        push_url = Url(push)
        if push_url.is_local:
            destination = push
        elif push_url.is_s3 or push_url.is_hdfs or push_url.is_nfs:
            destination = temp_dir
        else:
            raise RuntimeError('Push destination must be '
                               'on S3, HDFS, NFS, or local.')
    fastq_cues = set(['@'])
    fasta_cues = set(['>', ';'])
    source_dict = {}
    onward = False
    for line in sys.stdin:
        _input_line_count += 1
        if not line.strip(): continue
        # Kill offset from start of manifest file
        try:
            tokens = line.strip().split('\t')[1:]
            if tokens[0][0] == '#' and tokens[0] != '#!splitload':
                # Comment line
                continue
        except IndexError:
            # Be robust to bad lines
            continue
        token_count = len(tokens)
        qual_getter = None
        if tokens[0] == '#!splitload':
            '''Line specifies precisely how records from files should be
            placed.'''
            counter.add('splitload_lines')
            assert not to_stdout, ('Split manifest line inconsistent with '
                                   'writing to stdout.')
            qual_getter = phred_converter(phred_format=tokens[-1])
            indexes = tokens[1].split('\x1d')
            read_counts = tokens[2].split('\x1d')
            manifest_lines = [token.split('\x1e')
                                for token in tokens[3].split('\x1d')]
            assert len(indexes) == len(read_counts) == len(manifest_lines)
            for i, manifest_line in enumerate(manifest_lines):
                manifest_line_field_count = len(manifest_line)
                if manifest_line_field_count == 3:
                    source_dict[(Url(manifest_line[0]),)] = (
                            manifest_line[-1],
                            int(indexes[i]),
                            int(read_counts[i])
                        )
                else:
                    assert manifest_line_field_count == 5
                    source_dict[(Url(manifest_line[0]),
                                 Url(manifest_line[2]))] = (
                                                        manifest_line[-1],
                                                        int(indexes[i]),
                                                        int(read_counts[i])
                                                    )
        elif token_count == 3:
            # SRA or single-end reads
            counter.add('single_end_lines')
            source_dict[(Url(tokens[0]),)] = (tokens[-1],)
        elif token_count == 5:
            # Paired-end reads
            counter.add('paired_end_lines')
            source_dict[(Url(tokens[0]), Url(tokens[2]))] = (tokens[-1],)
        else:
            # Not a valid line, but continue for robustness
            continue
    file_number = 0
    counter.flush()
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
                # Index reads according to order in input to shorten read names
                read_index = skip_count / 2 # Index reads in pairs
            else:
                records_to_consume = source_dict[source_urls][2]
                read_index = skip_count
        else:
            skip_count = 0
            records_to_consume = None # Consume all records
            read_index = 0
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
                if source_url.is_dbgap:
                    download_dir = workspace_dir
                elif source_url.is_sra:
                    download_dir = temp_dir
                if source_url.is_sra:
                    sra_accession = source_url.to_url()
                    fastq_dump_command = (
                            'set -exo pipefail; cd {download_dir}; '
                            '{fastq_dump_exe} -I -X 10000 --split-files '
                            '{sra_accession}'
                        ).format(download_dir=download_dir,
                                    fastq_dump_exe=fastq_dump_exe,
                                    sra_accession=sra_accession)
                    try:
                        subprocess.check_call(
                            fastq_dump_command, shell=True, 
                            executable='/bin/bash',
                            stdout=sys.stderr
                        )
                    except subprocess.CalledProcessError as e:
                        if e.returncode == 3 and ignore_missing_sra_samples:
                            onward = True
                            break
                        else:
                            raise RuntimeError(
                                ('Error "%s" encountered executing '
                                 'command "%s".') % (e.output,
                                                        fastq_dump_command))
                    import glob
                    sra_fastq_files = sorted(
                                        glob.glob(os.path.join(download_dir,
                                            '%s[_.]*' % sra_accession))
                                        ) # ensure 1 before 2 if paired-end
                    # Schedule for deletion
                    def silent_remove(filename):
                        try:
                            os.remove(filename)
                        except OSError as e:
                            pass
                    for sra_fastq_file in sra_fastq_files:
                        register_cleanup(silent_remove, sra_fastq_file)
                    sra_file_count = len(sra_fastq_files)
                    check_for_paired = False
                    if sra_file_count == 1:
                        sra_paired_end = False
                        print >>sys.stderr, 'Detected single-end SRA sample.'
                    elif sra_file_count in [2, 3]:
                        print >>sys.stderr, ('2 or 3 FASTQ files detected. '
                                             'Checking for barcodes...')
                        check_for_paired = True
                    else:
                        raise RuntimeError(
                                ('Unexpected number of files "%d" output '
                                 'by fastq-dump command "%s".')
                                    % (sra_file_count, fastq_dump_command)
                            )
                    if check_for_paired:
                        # Get max/min read lengths from FASTQ
                        with open(
                                    sra_fastq_files[sra_file_count - 2]
                                ) as fastq_stream:
                            max_len, min_len = (
                                    max_min_read_lengths_from_fastq_stream(
                                        fastq_stream
                                    )
                                )
                            print >>sys.stderr, (
                                    'Max/min read length found in candidate '
                                    'barcode FASTQ was {}/{}.'
                                ).format(max_len, min_len)
                            if max_len <= _max_stubby_read_length:
                                print >>sys.stderr, (
                                        'Assumed barcode FASTQ.'
                                    )
                                skip_stubs = True
                                if sra_file_count == 2:
                                    sra_paired_end = False
                                else:
                                    sra_paired_end = True
                            else:
                                if sra_file_count == 2:
                                    sra_paired_end = True
                                else:
                                    raise RuntimeError(
                                        '3 FASTQs detected, but one of them '
                                        'was not recognized as containing '
                                        'barcodes.'
                                    )
                    # Guess quality from first 10k lines
                    with xopen(None, sra_fastq_files[0]) as source_stream:
                        qual_getter = phred_converter(
                                            fastq_stream=source_stream
                                        )
                    for sra_fastq_file in sra_fastq_files:
                        os.remove(sra_fastq_file)
                    sources.append(os.devnull)
                    fastq_dump_command = (
                            'set -exo pipefail; cd {download_dir}; '
                            '{fastq_dump_exe} --split-spot -I --stdout '
                            '{sra_accession}'
                        ).format(download_dir=download_dir,
                                    fastq_dump_exe=fastq_dump_exe,
                                    sra_accession=sra_accession)
                    if skip_stubs:
                        fastq_dump_command += (
                                ' | awk \'BEGIN {{OFS = "\\n"}} '
                                '{{header = $0; '
                                'getline seq; getline qheader; getline qseq; '
                                'if (length(seq) > {min_len}) {{print header, '
                                'seq, qheader, qseq}}}}\''
                            ).format(min_len=_max_stubby_read_length)
                    print >>sys.stderr, fastq_dump_command
                    sra_process = subprocess.Popen(fastq_dump_command,
                                                    shell=True,
                                                    executable='/bin/bash',
                                                    stdout=subprocess.PIPE,
                                                    bufsize=-1)
                else:
                    mover.get(source_url, temp_dir)
                    downloaded = list(
                            set(os.listdir(temp_dir)).difference(downloaded)
                        )
                    sources.append(os.path.join(temp_dir, list(downloaded)[0]))
            else:
                sources.append(source_url.to_url())
        if onward: continue
        '''Use os.devnull so single- and paired-end data can be handled in one
        loop.'''
        if len(sources) == 1:
            sources.append(os.devnull)
        if qual_getter is None:
            # Figure out Phred format
            with xopen(None, sources[0]) as source_stream:
                qual_getter = phred_converter(fastq_stream=source_stream)
        with xopen(None, sources[0]) as source_stream_1, xopen(
                None, sources[1]
            ) as source_stream_2:
            source_streams = [source_stream_1, source_stream_2]
            reorganize = all([source == os.devnull for source in sources])
            if reorganize:
                # SRA data is live
                if sra_paired_end:
                    source_streams = [sra_process.stdout, sra_process.stdout]
                else:
                    source_streams = [sra_process.stdout, open(os.devnull)]
            break_outer_loop = False
            while True:
                if not to_stdout:
                    '''Name files using Hadoop task environment property
                    mapred.task.partition, and tack ".gz" onto output filename
                    and gzip compression level onto arguments passed to opener
                    as appropriate.'''
                    if gzip_output:
                        filename_suffix = '.gz'
                        open_args_suffix = [gzip_level]
                    else:
                        filename_suffix = ''
                        open_args_suffix = []
                    try:
                        output_file = os.path.join(
                                destination, 
                                '.'.join([
                                    os.environ['mapred_task_partition'],
                                    str(file_number)
                                ]) + filename_suffix
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
                                    str(file_number)
                                ]) + filename_suffix
                            )
                    open_args = [output_file, 'a'] + open_args_suffix
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
                    pairs_read = 0
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
                        bad_record_skip = False
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
                            quals = [source_stream.readline().strip()
                                        for source_stream in source_streams]
                            if reorganize and sra_paired_end:
                                # Fix order!
                                lines, seqs, plus_lines, quals = (
                                        [lines[0], plus_lines[0]],
                                        [lines[1], plus_lines[1]],
                                        [seqs[0], quals[0]],
                                        [seqs[1], quals[1]]
                                    )
                            try:
                                assert plus_lines[0][0] == '+', (
                                        'Malformed read "%s" at line %d of '
                                        'file "%s".'
                                    ) % (lines[0], line_numbers[0], sources[0])
                                if plus_lines[1]:
                                    assert plus_lines[1][0] == '+', (
                                            'Malformed read "%s" at line %d '
                                            'of file "%s".'
                                        ) % (
                                        lines[1], line_numbers[1], sources[1]
                                    )
                                try:
                                    # Kill spaces in name
                                    original_qnames = \
                                        [line[1:].replace(' ', '_')
                                            for line in lines]
                                except IndexError:
                                    raise RuntimeError(
                                            'Read name not found at ' 
                                            'line %d of %s.' % (
                                                        line_numbers[0],
                                                        sources
                                                    )
                                        )
                            except (AssertionError,
                                    IndexError, RuntimeError) as e:
                                if skip_bad_records:
                                    print >>sys.stderr, ('Error "%s" '
                                            'encountered; skipping bad record.'
                                        ) % e.message
                                    for source_stream in source_streams:
                                        source_stream.readline()
                                    line_numbers = [
                                            i + 1 for i in line_numbers
                                        ]
                                    bad_record_skip = True
                                else:
                                    raise
                            else:
                                try:
                                    quals = [
                                            qual_getter(qual) for qual in quals
                                        ]
                                except Exception as e:
                                    if skip_bad_records:
                                        print >>sys.stderr, (
                                                'Error "%s" encountered '
                                                'trying to convert quality '
                                                'string to Sanger format; '
                                                'skipping bad record.'
                                            ) % e.message
                                        bad_record_skip = True
                                    else:
                                        raise
                                line_numbers = [i + 1 for i in line_numbers]
                                try: 
                                    for i in xrange(2):
                                        assert len(seqs[i]) == len(quals[i]), (
                                            'Length of read sequence does not '
                                            'match length of quality string '
                                            'at line %d of file "%s".'
                                        ) % (line_numbers[i], sources[i])
                                except (AssertionError, IndexError) as e:
                                    if skip_bad_records:
                                        print >>sys.stderr, (
                                                'Error "%s" encountered; '
                                                'skipping bad record.'
                                            ) % e.message
                                        bad_record_skip = True
                                    else:
                                        raise
                        elif lines[0][0] in fasta_cues:
                            seqs = [[], []]
                            next_lines = []
                            for p, source_stream in enumerate(source_streams):
                                while True:
                                    next_line \
                                        = source_stream.readline().strip()
                                    try:
                                        if next_line[0] in fasta_cues:
                                            break
                                        else:
                                            try:
                                                seqs[p].append(next_line)
                                            except IndexError:
                                                raise
                                    except IndexError:
                                        break
                                next_lines.append(next_line)
                            seqs = [''.join(seq) for seq in seqs]
                            line_numbers = [i + 1 for i in line_numbers]
                            try:
                                try:
                                    # Kill spaces in name
                                    original_qnames = \
                                        [line[1:].replace(' ', '_')
                                            for line in lines]
                                except IndexError:
                                    raise RuntimeError(
                                            'Error finding QNAME at ' 
                                            'line %d of either %s or %s' % (
                                                        sources[0],
                                                        sources[1]
                                                    )
                                        )
                            except (AssertionError,
                                    IndexError, RuntimeError) as e:
                                if skip_bad_records:
                                    print >>sys.stderr, ('Error "%s" '
                                            'encountered; skipping bad record.'
                                        ) % e.message
                                    for source_stream in source_streams:
                                        source_stream.readline()
                                    line_numbers = [
                                            i + 1 for i in line_numbers
                                        ]
                                    bad_record_skip = True
                                else:
                                    raise
                            else:
                                try:
                                    quals = [
                                        'h'*len(seq) for seq in seqs
                                        ]
                                except Exception as e:
                                    if skip_bad_records:
                                        print >>sys.stderr, (
                                                'Error "%s" encountered '
                                                'trying to convert quality '
                                                'string to Sanger format; '
                                                'skipping bad record.'
                                            ) % e.message
                                        bad_record_skip = True
                                    else:
                                        raise
                                line_numbers = [i + 1 for i in line_numbers]
                            lines = next_lines
                            read_next_line = False
                        if bad_record_skip:
                            seqs = []
                            # Fake record-printing to get to records_to_consume
                            if source_streams[-1].name == os.devnull:
                                records_printed += 1
                            else:
                                records_printed += 2
                        elif len(original_qnames) == 2 and original_qnames[1]:
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
                            if short_qnames:
                                left_qname_to_write = encode(read_index) + '/1'
                                right_qname_to_write = encode(
                                                            read_index
                                                        ) + '/2'
                            else:
                                left_qname_to_write = original_qnames[0]
                                right_qname_to_write = original_qnames[1]
                            print >>output_stream, '\t'.join(
                                        [
                                            left_seq,
                                            left_reversed,
                                            qname_from_read(
                                                    left_qname_to_write,
                                                    seqs[0] + quals[0], 
                                                    sample_label,
                                                    mate=seqs[1]
                                                ),
                                            '\n'.join([
                                                round_quality_string(
                                                    left_qual
                                                ), right_seq
                                            ]),
                                            right_reversed,
                                            qname_from_read(
                                                    right_qname_to_write,
                                                    seqs[1] + quals[1], 
                                                    sample_label,
                                                    mate=seqs[0]
                                                ),
                                            round_quality_string(right_qual)
                                        ]
                                    )
                            records_printed += 2
                            _output_line_count += 1
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
                            if short_qnames:
                                qname_to_write = encode(read_index)
                            else:
                                qname_to_write = original_qnames[0]
                            print >>output_stream, '\t'.join(
                                        [
                                            seq,
                                            is_reversed,
                                            qname_from_read(
                                                qname_to_write,
                                                seqs[0] + quals[0], 
                                                sample_label
                                            ),
                                            round_quality_string(qual)
                                        ]
                                    )
                            records_printed += 1
                            _output_line_count += 1
                        read_index += 1
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
                if (not to_stdout) and (push_url.is_nfs or
                    push_url.is_s3 or push_url.is_hdfs) \
                    and ((not records_to_consume) or
                         (records_to_consume and perform_push)):
                    print >>sys.stderr, 'Pushing "%s" to "%s" ...' % (
                                                            output_file,
                                                            push_url.to_url()
                                                        )
                    print >>sys.stderr, 'reporter:status:alive'
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
        if 'sra_process' in locals():
            sra_process.stdout.close()
            sra_return_code = sra_process.wait()
            if sra_return_code > 0:
                raise RuntimeError(('fastq-dump terminated with exit '
                                    'code %d. Command run was "%s".')
                                        % (sra_return_code,
                                            fastq_dump_command))
            del sra_process

if __name__ == '__main__':
    import unittest

    class TestEmpty(unittest.TestCase):

        def test_empty(self):
            pass

    if '--test' in sys.argv:
        unittest.main(argv=[sys.argv[0]])
        sys.exit(0)

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
    parser.add_argument('--bin-qualities', action='store_const', const=True,
        default=False,
        help='Rounds base quality score, placing it in one of five bins')
    parser.add_argument('--keep-alive', action='store_const', const=True,
        default=False,
        help='Periodically print Hadoop status messages to stderr to keep ' \
             'job alive')
    parser.add_argument('--shorten-read-names', action='store_const',
        const=True, default=False,
        help='Gives each read (pair) a unique short read name to save space')
    parser.add_argument('--skip-bad-records', action='store_const',
        const=True, default=False,
        help=('Skips bad records rather than raising exception; here, '
              '"bad" could mean that the quality sequence length does '
              'not match the read length or the record is not in the '
              'proper format'))
    parser.add_argument(
            '--workspace-dir', required=False, 
            default=None,
            help='Path to SRA Tools workspace directory; needed if '
                 'downloading dbGaP data'
        )
    parser.add_argument(
            '--fastq-dump-exe', required=False, 
            default='fastq-dump',
            help='Path to fastq-dump executable'
        )
    parser.add_argument('--ignore-missing-sra-samples', action='store_const',
        const=True, default=False,
        help='Does not raise exception if fastq-dump doesn\'t find an SRA '
             'sample; instead, sample is skipped')
    parser.add_argument('--verbose', action='store_const', const=True,
        default=False,
        help='Print out extra debugging statements')
    filemover.add_args(parser)
    tempdel.add_args(parser)
    args = parser.parse_args()
    # Start keep_alive thread immediately
    if args.keep_alive:
        from dooplicity.tools import KeepAlive
        keep_alive_thread = KeepAlive(sys.stderr)
        keep_alive_thread.start()
    import time
    start_time = time.time()
    mover = filemover.FileMover(args=args)
    go(nucleotides_per_input=args.nucs_per_file,
        gzip_output=args.gzip_output,
        gzip_level=args.gzip_level,
        to_stdout=args.stdout,
        push=args.push,
        bin_qualities=args.bin_qualities,
        short_qnames=args.shorten_read_names,
        skip_bad_records=args.skip_bad_records,
        scratch=tempdel.silentexpandvars(args.scratch),
        verbose=args.verbose,
        mover=mover,
        workspace_dir=args.workspace_dir,
        fastq_dump_exe=args.fastq_dump_exe,
        ignore_missing_sra_samples=args.ignore_missing_sra_samples)
    print >>sys.stderr, 'DONE with preprocess.py; in/out=%d/%d; ' \
        'time=%0.3f s' % (_input_line_count, _output_line_count,
                            time.time() - start_time)
