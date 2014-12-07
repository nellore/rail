"""
realign_reads_delegate.py 

Output of Bowtie2 from realign_reads.py is streamed to this script to obtain
final output. See realign_reads.py for output format information.
"""

import sys
import os
import site

base_path = os.path.abspath(
                    os.path.dirname(os.path.dirname(os.path.dirname(
                        os.path.realpath(__file__)))
                    )
                )
utils_path = os.path.join(base_path, 'rna', 'utils')
site.addsitedir(utils_path)
site.addsitedir(base_path)

from dooplicity.tools import xstream, xopen

import string
_reversed_complement_translation_table = string.maketrans('ATCG', 'TAGC')

def go(rname_stream, output_stream=sys.stdout, input_stream=sys.stdin,
        verbose=False, report_multiplier=1.2):
    """ Processes Bowtie 2 alignments, emitting filtered SAM output.

        Only those alignments whose RNAMEs appear in rname_stream are output
        to ensure replicability of SAM output independently of which contigs
        appear in the Bowtie 2 index.

        rname_stream: where to retrieve RNAMES, each of which is in the form
            (number of qnames associated with rnames) + '\x1e' +
            ('\x1e'-separated list of valid rnames)
        output_stream: where to emit exon and intron tuples; typically, this is
            sys.stdout.
        verbose: True if alignments should occasionally be written to stderr.
        report_multiplier: if verbose is True, the line number of an
            alignment written to stderr increases exponentially with base
            report_multiplier.
    """
    output_line_count, next_report_line, i, qname_count = 0, 0, 0, 0
    tokens = rname_stream.readline().strip().split('\x1e')
    qname_total, rnames = int(tokens[0]), set(tokens[1:])
    done = False
    for (qname,), xpartition in xstream(input_stream, 1):
        assert not done
        printed = False
        for rest_of_line in xpartition:
            i += 1
            flag = int(rest_of_line[0])
            if verbose and next_report_line == i:
                print >>sys.stderr, \
                    'SAM output record %d: rdname="%s", flag=%d' \
                    % (i, qname, flag)
                next_report_line = max(int(next_report_line
                    * report_multiplier), next_report_line + 1)
            rname = rest_of_line[1]
            if rname in rnames:
                print >>output_stream, \
                    '\t'.join((qname,) + rest_of_line)
                printed = True
                output_line_count += 1
        if flag & 4 or not printed:
            # This is an unmapped read
            if flag & 16:
                seq_to_write = rest_of_line[8][::-1].translate(
                                _reversed_complement_translation_table
                            )
                qual_to_write = rest_of_line[9][::-1]
            else:
                seq_to_write = rest_of_line[8]
                qual_to_write = rest_of_line[9]
            # Write only essentials; handle "formal" writing in next step
            print >>output_stream, ('%s\t4\t\x1c\t\x1c\t\x1c\t\x1c'
                                         '\t\x1c\t\x1c\t\x1c\t%s\t%s') % (
                                                qname,
                                                seq_to_write,
                                                qual_to_write
                                            )
        qname_count += 1
        if qname_count == qname_total:
            tokens = rname_stream.readline().strip().split('\x1e')
            try:
                qname_total, rnames = int(tokens[0]), set(tokens[1:])
                qname_count = 0
            except ValueError:
                done = True
    output_stream.flush()
    print >>sys.stderr, ('realign_reads_delegate.py reports %d output lines.'
                            % output_line_count)

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__, 
                formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--verbose', action='store_const', const=True,
        default=False,
        help='Print out extra debugging statements')
    parser.add_argument('--rnames-file', type=str, required=True,
        help=('Where extended RNAMEs storing those RNAMEs to which a given '
              'read sequence is "allowed" to align are stored'))
    parser.add_argument('--report-multiplier', type=float, required=False,
        default=1.2,
        help='When --verbose is also invoked, the only lines of lengthy '
             'intermediate output written to stderr have line number that '
             'increases exponentially with this base')
    args = parser.parse_args()

    with xopen(None, args.rnames_file) as rname_stream:
        go(rname_stream, verbose=args.verbose,
            report_multiplier=args.report_multiplier)