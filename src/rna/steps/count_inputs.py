"""
Rail-RNA-count_inputs

START: follows no step
Precedes Rail-RNA-assign-splits

Runs "wc -cl" on each line of input (which is a manifest file) and prepends
the result to the line if URLs are local.

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
Tab-separated output fields:
---If URL is local:
1. #!splitload
2. number of read(s) (pairs) in sample; number of pairs if paired-end and
    number of reads if single-end
3. number of uncompressed bytes in left (or right) reads file
4 ... end. same as manifest line

---Otherwise:
same as manifest line
"""

import os
import sys
import site
import time

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
import subprocess
import argparse
# Print file's docstring if -h is invoked
parser = argparse.ArgumentParser(description=__doc__, 
            formatter_class=argparse.RawDescriptionHelpFormatter)

args = parser.parse_args(sys.argv[1:])

start_time = time.time()
# For telling apart FASTQs and FASTAs
fastq_cues = set(['@'])
fasta_cues = set(['>', ';'])

input_line_count, output_line_count = 0, 0

for input_line_count, line in enumerate(sys.stdin):
    # Kill offset from start of manifest file
    tokens = line.strip().split('\t')[1:]
    try:
        stripped = tokens[0].strip()
        if stripped[0] == '#' or not line.strip():
            continue
    except IndexError:
        continue
    token_count = len(tokens)
    assert token_count in [3, 5], (
            'Line {} of input has {} fields, but 3 or 5 are expected.'
        ).format(input_line_count + 1, token_count)
    file_to_count = tokens[0]
    if (not ((token_count == 3 and Url(tokens[0]).is_local) or
        (token_count == 5 and Url(tokens[0]).is_local
            and Url(tokens[2]).is_local))):
        sys.stdout.write(line)
        output_line_count += 1
        continue
    with open(file_to_count, 'rb') as binary_input_stream:
        if binary_input_stream.read(2) == '\x1f\x8b':
            # Magic number of gzip'd file found
            command_to_run = 'wc -lc $(gzip -cd {})'.format(file_to_count)
        else:
            command_to_run = 'wc -lc {}'.format(file_to_count)
    with xopen(None, file_to_count) as input_stream:
        first_char = input_stream.readline()[0]
        if first_char in fastq_cues:
            # 4 lines per record
            line_divider = 4
        elif first_char in fastq_cues:
            line_divider = 2
        else:
            raise RuntimeError(
                    'File "{}" is neither a FASTA nor a FASTQ file.'.format(
                            file_to_count
                        )
                )
    try:
        lines_and_bytes = subprocess.check_output(
                                        command_to_run,
                                        shell=True,
                                        executable='/bin/bash',
                                        bufsize=-1
                                    ).split()[:-1]
    except subprocess.CalledProcessError:
        from traceback import format_exc
        print >>sys.stderr, \
                'Error\n\n{}\nencountered counting lines with {}.'.format(
                    format_exc(), command_to_run
                )
    lines_and_bytes[0] = str(int(lines_and_bytes[0]) / line_divider)
    sys.stdout.write(
        '\t'.join(
            ['#!splitload'] + lines_and_bytes + [line.partition('\t')[2]]
        )
    )
    output_line_count += 1

sys.stdout.flush()
print >>sys.stderr, 'DONE with count_inputs.py; in/out=%d/%d; ' \
        'time=%0.3f s' % (input_line_count + 1, output_line_count,
                            time.time() - start_time)