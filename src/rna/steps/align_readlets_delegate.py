"""
align_readlets_delegate.py 

Output of Bowtie from align_readlets.py is streamed to this script to obtain
final output. See align_readlets.py for output format information.
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

def go(qname_stream, output_stream=sys.stdout, input_stream=sys.stdin,
        verbose=False, report_multiplier=1.2):
    """ Emits readlet alignments.

        qname_stream contains long QNAMEs in the order in which readlets passed
        to Bowtie appeared. These names would have been truncated. Each QNAME
        takes the form

            '\x1d'-separated list of [read sequence ID + ('-' if readlet
            sequence is reverse-complemented; else '+') + '\x1e' + displacement
            of readlet's 5' end from read's 5' end + '\x1e' + displacement of
            readlet's 3' end from read's 3' end (+, for EXACTLY one readlet of
            a read sequence, '\x1e' + read sequence + '\x1e' +
            (an '\x1f'-separated list A of unique sample labels with read
            sequences that match the original read sequence) + '\x1e' +
            (an '\x1f'-separated list  of unique sample labels B with read
            sequences that match the reversed complement of the original read
            sequence)) + '\x1e' + (an '\x1f'-separated list of the number of
            instances of the read sequence for each respective sample in list
            A) + '\x1e' + (an '\x1f'-separated list of the number of instances
            of the read sequence's reversed complement for each respective
            sample in list B)]

        A line is written per readlet per associated read sequence. So if a
        given readlet can be found on 3 reads, 3 lines are written, each
        containing the readlet's alignments.

        qname_stream: where to retrieve extended qnames
        input_stream: where to retrieve Bowtie output
        output_stream: where to emit exon and junction tuples; typically, this
            is sys.stdout.
        verbose: True if alignments should occasionally be written to stderr.
        report_multiplier: if verbose is True, the line number of an
            alignment written to stderr increases exponentially with base
            report_multiplier.
    """
    output_line_count, next_report_line, i = 0, 0, 0
    for (qname,), xpartition in xstream(input_stream, 1):
        '''While labeled multireadlet, this list may end up simply a
        unireadlet.'''
        multireadlet = []
        flag = 4
        for tokens in xpartition:
            (flag, rname, pos, mapq, cigar,
                rnext, pnext, tlen, seq, qual) = tokens[:10]
            flag = int(flag)
            multireadlet.append((rname, flag, pos))
            if verbose and next_report_line == i:
                print >>sys.stderr, \
                    'SAM output record %d: rdname="%s", flag=%d' % (i,
                                                                    qname,
                                                                    flag)
                next_report_line = int((next_report_line + 1)
                                        * report_multiplier + 1) - 1
            i += 1
        '''If the next qname doesn't match the last qname or there are no
        more lines, all of a multireadlet's alignments have been
        collected.'''
        if not flag & 4:
            '''Last readlet has at least one alignment; print all
            alignments for each read from which readlet sequence is
            derived.'''
            rnames, flags, poses = zip(*multireadlet)
            reverse_flags = [a_flag ^ 16 for a_flag in flags]
            flags = '\x1f'.join([str(a_flag) for a_flag in flags])
            reverse_flags = '\x1f'.join(
                                    [str(a_flag) for a_flag
                                        in reverse_flags]
                                )
            rnames = '\x1f'.join(rnames)
            poses = '\x1f'.join(poses)
            read = qname_stream.readline().strip()
            while read != '+':
                read_id, _, read_rest = read.partition('\x1e')
                if read_id[-1] == '-':
                    current_flags = reverse_flags
                else:
                    current_flags = flags
                print >>output_stream, '%s\t%s\t%s\t%s\t%s' % \
                    (read_id[:-1], read_rest, rnames,
                        current_flags, poses)
                output_line_count += 1
                read = qname_stream.readline().strip()
        else:
            '''Readlet had no reported alignments; print ONLY when readlet
            contains general info about read.'''
            read = qname_stream.readline().strip()
            while read != '+':
                read_id, _, read_rest = read.partition('\x1e')
                if len(read_rest.split('\x1e')) > 2:
                    print >>output_stream, \
                        '%s\t%s\t\x1c\t\x1c\t\x1c' % (read_id[:-1],
                                                        read_rest)
                output_line_count += 1
                read = qname_stream.readline().strip()
    output_stream.flush()
    print >>sys.stderr, ('align_readlets_delegate.py reports %d output lines.'
                            % output_line_count)

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__, 
                formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--verbose', action='store_const', const=True,
        default=False,
        help='Print out extra debugging statements')
    parser.add_argument('--qnames-file', type=str, required=True,
        help=('Where to find extended QNAMEs storing read sequence IDs to '
              'which readlets belong and other pertinent information'))
    parser.add_argument('--report-multiplier', type=float, required=False,
        default=1.2,
        help='When --verbose is also invoked, the only lines of lengthy '
             'intermediate output written to stderr have line number that '
             'increases exponentially with this base')
    args = parser.parse_args()

    with xopen(None, args.qnames_file) as qname_stream:
        go(qname_stream, verbose=args.verbose,
            report_multiplier=args.report_multiplier)