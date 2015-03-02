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

def go(output_stream=sys.stdout, input_stream=sys.stdin,
        verbose=False, report_multiplier=1.2):
    """ Processes Bowtie 2 alignments, emitting filtered SAM output.

        output_stream: where to emit exon and intron tuples; typically, this is
            sys.stdout.
        input_stream: where to find input to process
        verbose: True if alignments should occasionally be written to stderr.
        report_multiplier: if verbose is True, the line number of an
            alignment written to stderr increases exponentially with base
            report_multiplier.
    """
    output_line_count, next_report_line, i = 0, 0, 0
    for (qname,), xpartition in xstream(input_stream, 1):
        for rest_of_line in xpartition:
            i += 1
            flag = int(rest_of_line[0])
            if verbose and next_report_line == i:
                print >>sys.stderr, \
                    'SAM output record %d: rdname="%s", flag=%d' \
                    % (i, qname, flag)
                next_report_line = max(int(next_report_line
                    * report_multiplier), next_report_line + 1)
            print >>output_stream, \
                '\t'.join((qname,) + rest_of_line)
            output_line_count += 1
        if flag & 4:
            '''Unmapped read; write only essentials, and handle "formal"
            writing in next step.'''
            print >>output_stream, ('%s\t4\t\x1c\t\x1c\t\x1c\t\x1c'
                                         '\t\x1c\t\x1c\t\x1c\t%s\t%s') % (
                                                qname,
                                                rest_of_line[8],
                                                rest_of_line[9]
                                            )
            output_line_count += 1
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
    parser.add_argument('--report-multiplier', type=float, required=False,
        default=1.2,
        help='When --verbose is also invoked, the only lines of lengthy '
             'intermediate output written to stderr have line number that '
             'increases exponentially with this base')
    args = parser.parse_args()

    go(verbose=args.verbose,
        report_multiplier=args.report_multiplier)