#!/usr/bin/env python
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
3 ... end. same as manifest line
4. Phred format (Sanger or Phred64)

---Otherwise:
same as manifest line
"""

import os
import sys
import site
import time

if '--test' in sys.argv:
    print("No unit tests")
    #unittest.main(argv=[sys.argv[0]])
    sys.exit(0)

base_path = os.path.abspath(
                    os.path.dirname(os.path.dirname(os.path.dirname(
                        os.path.realpath(__file__)))
                    )
                )
utils_path = os.path.join(base_path, 'rna', 'utils')
site.addsitedir(utils_path)
site.addsitedir(base_path)

from dooplicity.ansibles import Url
from dooplicity.tools import xopen, register_cleanup
from dooplicity.counters import Counter
import argparse
from guess import inferred_phred_format
# Print file's docstring if -h is invoked
parser = argparse.ArgumentParser(description=__doc__, 
            formatter_class=argparse.RawDescriptionHelpFormatter)

args = parser.parse_args(sys.argv[1:])

start_time = time.time()
# For telling apart FASTQs and FASTAs
fastq_cues = set(['@'])
fasta_cues = set(['>', ';'])

input_line_count, output_line_count = 0, 0
counter = Counter('count_inputs')
register_cleanup(counter.flush)
line_divider = 4

for input_line_count, line in enumerate(sys.stdin):
    # Kill offset from start of manifest file
    tokens = line.strip().split('\t')[1:]
    try:
        stripped = tokens[0].strip()
        if stripped[0] == '#' or not line.strip():
            counter.add('comment_lines')
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
            and Url(tokens[2]).is_local)) or 
            (Url(tokens[0]).is_local and 
                (tokens[0].endswith('.tar.gz')
                    or tokens[0].endswith('.tar.bz2')
                    or tokens[0].endswith('.tar')))
            ):
        sys.stdout.write(line)
        output_line_count += 1
        continue
    with xopen(None, file_to_count) as input_stream:
        first_char = input_stream.readline()[0]
        if first_char in fasta_cues:
            sys.stdout.write(line)
            counter.add('fasta_line')
            output_line_count += 1
        elif first_char in fastq_cues:
            counter.add('fastq_line')
        else:
            raise RuntimeError(
                    'File "{}" is neither a FASTA nor a FASTQ file.'.format(
                            file_to_count
                        )
                )
    with xopen(None, file_to_count) as input_stream:
        phred_format, line_count = inferred_phred_format(input_stream)
        counter.add('inferred_' + phred_format)
    lines_and_bytes = str((int(line_count) + 1) / line_divider)
    print '\t'.join(
            ['#!splitload', lines_and_bytes, line.partition('\t')[2].strip(),
                phred_format]
        )
    output_line_count += 1
    counter.flush()

sys.stdout.flush()
print >>sys.stderr, 'DONE with count_inputs.py; in/out=%d/%d; ' \
        'time=%0.3f s' % (input_line_count + 1, output_line_count,
                            time.time() - start_time)
