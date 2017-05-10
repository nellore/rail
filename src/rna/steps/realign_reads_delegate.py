"""
realign_reads_delegate.py 

Output of Bowtie2 from realign_reads.py is streamed to this script to obtain
final output. See realign_reads.py for output format information.
"""

import sys
import os
import site

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

from dooplicity.tools import xstream, xopen, register_cleanup
from dooplicity.counters import Counter

counter = Counter('realign_reads_delegate')
register_cleanup(counter.flush)

import string
_reversed_complement_translation_table = string.maketrans('ATCG', 'TAGC')

def go(output_stream=sys.stdout, input_stream=sys.stdin,
        verbose=False, report_multiplier=1.2,
        alignment_count_to_report=1, tie_margin=0):
    """ Processes Bowtie 2 alignments, emitting filtered SAM output.

        Only max(# tied alignments, alignment_count_to_report) alignments
        are printed. This way, the compare_alignments step always has enough
        information to fill the XS field.

        output_stream: where to emit exon and junction tuples; typically, this
            is sys.stdout.
        input_stream: where to find input to process
        verbose: True if alignments should occasionally be written to stderr.
        report_multiplier: if verbose is True, the line number of an
            alignment written to stderr increases exponentially with base
            report_multiplier.
        alignment_count_to_report: argument of Bowtie 2's -k field
        tie_margin: allowed score difference per 100 bases among ties in 
             max alignment score.
    """
    output_line_count, next_report_line = 0, 0
    threshold_alignment_count = max(2, alignment_count_to_report)
    for (qname,), xpartition in xstream(input_stream, 1):
        counter.add('partitions')
        max_score, alignments_output, current_tie_margin = None, 0, None
        for rest_of_line in xpartition:
            counter.add('inputs')
            # Note Bowtie 2 outputs alignments in order of descending score
            try:
                score = int([field[5:] for field in rest_of_line
                                if field[:5] == 'AS:i:'][0])
            except IndexError:
                # Unmapped read; flag should be 4. Print only essentials.
                counter.add('unaligned_records')
                assert int(rest_of_line[0]) == 4
                print >>output_stream, ('%s\t4\t\x1c\t\x1c\t\x1c\t\x1c'
                                         '\t\x1c\t\x1c\t\x1c\t%s\t%s') % (
                                                qname,
                                                rest_of_line[8],
                                                rest_of_line[9]
                                            )
                output_line_count += 1
            else:
                if current_tie_margin is None:
                    current_tie_margin = round(
                            tie_margin * float(len(rest_of_line[8])) / 100
                        )
                if score + current_tie_margin >= max_score:
                    max_score = max(max_score, score)
                elif alignments_output >= threshold_alignment_count:
                    break
                counter.add('aligned_records')
                print >>output_stream, '\t'.join((qname,) + rest_of_line)
                alignments_output += 1
                output_line_count += 1
                if verbose and next_report_line == output_line_count:
                    print >>sys.stderr, \
                        'SAM output record %d: rdname="%s", flag=%d' \
                        % (output_line_count, qname, int(rest_of_line[0]))
                    next_report_line = max(int(next_report_line
                        * report_multiplier), next_report_line + 1)
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
    parser.add_argument('--alignment-count-to-report', type=int,
        required=False, default=1,
        help='Argument of Bowtie 2\'s -k parameter')
    from alignment_handlers import add_args as alignment_handlers_add_args
    alignment_handlers_add_args(parser)
    args = parser.parse_args()

    go(verbose=args.verbose,
        report_multiplier=args.report_multiplier,
        alignment_count_to_report=args.alignment_count_to_report,
        tie_margin=args.tie_margin)
