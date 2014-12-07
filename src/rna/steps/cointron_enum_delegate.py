"""
cointron_enum_delegate.py 

Output of Bowtie 2 from cointron_enum.py is streamed to this script to obtain
final output. See cointron_enum.py for output format information.
"""

import sys
import os
import site
from collections import defaultdict

base_path = os.path.abspath(
                    os.path.dirname(os.path.dirname(os.path.dirname(
                        os.path.realpath(__file__)))
                    )
                )
utils_path = os.path.join(base_path, 'rna', 'utils')
site.addsitedir(utils_path)
site.addsitedir(base_path)

from dooplicity.tools import xstream
from alignment_handlers import multiread_with_introns, indels_introns_and_exons

import string
_reversed_complement_translation_table = string.maketrans('ATCG', 'TAGC')

def maximal_cliques(cointrons):
    """ Finds maximal cliques of graph of intron combinations.

        Consider an undirected graph where each node corresponds to
        cointron, as specified by an element of the list "cointrons". Place
        an edge between two nodes for which no two _distinct_ introns overlap
        where one intron is associated with one node, and the other intron is
        associated with the other node. Now enumerate maximal cliques. This
        gives all possible valid clusters of _consistent_ introns.

        This code is adapted from NetworkX's find_cliques(), which requires
        inclusion of the following copyright notice.

        --------
        Copyright (C) 2004-2012, NetworkX Developers
        Aric Hagberg <hagberg@lanl.gov>
        Dan Schult <dschult@colgate.edu>
        Pieter Swart <swart@lanl.gov>
        All rights reserved.

        Redistribution and use in source and binary forms, with or without
        modification, are permitted provided that the following conditions are
        met:

          * Redistributions of source code must retain the above copyright
            notice, this list of conditions and the following disclaimer.

          * Redistributions in binary form must reproduce the above
            copyright notice, this list of conditions and the following
            disclaimer in the documentation and/or other materials provided
            with the distribution.

          * Neither the name of the NetworkX Developers nor the names of its
            contributors may be used to endorse or promote products derived
            from this software without specific prior written permission.


        THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
        "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
        LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
        A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
        OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
        SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
        LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
        DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
        THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
        (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
        OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
        --------

        cointrons: list of introns from the same strand
        
        Yield value: A maximal clique -- a list of introns.
    """
    # Cache nbrs and find first pivot (highest degree)
    maxconn=-1
    nnbrs={}
    pivotnbrs=set() # handle empty graph
    for n,nbrs in [(intron, [an_intron for an_intron in cointrons
                                if not ((min(an_intron[1], intron[1])
                                           - max(an_intron[0], intron[0]))
                                            >= 0)])
                    for intron in cointrons]:
        nbrs=set(nbrs)
        nbrs.discard(n)
        conn = len(nbrs)
        if conn > maxconn:
            nnbrs[n] = pivotnbrs = nbrs
            maxconn = conn
        else:
            nnbrs[n] = nbrs
    # Initial setup
    cand=set(nnbrs)
    smallcand = set(cand - pivotnbrs)
    done=set()
    stack=[]
    clique_so_far=[]
    # Start main loop
    while smallcand or stack:
        try:
            # Any nodes left to check?
            n=smallcand.pop()
        except KeyError:
            # back out clique_so_far
            cand,done,smallcand = stack.pop()
            clique_so_far.pop()
            continue
        # Add next node to clique
        clique_so_far.append(n)
        cand.remove(n)
        done.add(n)
        nn=nnbrs[n]
        new_cand = cand & nn
        new_done = done & nn
        # check if we have more to search
        if not new_cand:
            if not new_done:
                # Found a clique!
                yield clique_so_far[:]
            clique_so_far.pop()
            continue
        # Shortcut--only one node left!
        if not new_done and len(new_cand)==1:
            yield clique_so_far + list(new_cand)
            clique_so_far.pop()
            continue
        # find pivot node (max connected in cand)
        # look in done nodes first
        numb_cand=len(new_cand)
        maxconndone=-1
        for n in new_done:
            cn = new_cand & nnbrs[n]
            conn=len(cn)
            if conn > maxconndone:
                pivotdonenbrs=cn
                maxconndone=conn
                if maxconndone==numb_cand:
                    break
        # Shortcut--this part of tree already searched
        if maxconndone == numb_cand:
            clique_so_far.pop()
            continue
        # still finding pivot node
        # look in cand nodes second
        maxconn=-1
        for n in new_cand:
            cn = new_cand & nnbrs[n]
            conn=len(cn)
            if conn > maxconn:
                pivotnbrs=cn
                maxconn=conn
                if maxconn == numb_cand-1:
                    break
        # pivot node is max connected in cand from done or cand
        if maxconndone > maxconn:
            pivotnbrs = pivotdonenbrs
        # save search status for later backout
        stack.append( (cand, done, smallcand) )
        cand=new_cand
        done=new_done
        smallcand = cand - pivotnbrs

def separated_introns(introns, separation):
    """ Splits introns up if successive introns are separated by > separation

        Two introns overlapped by readlet alignments may be on the same 
        strand but much further apart than a read can overlap. This function
        splits up alignments if successive introns are separated by more
        than separation. IT ALSO ENUMERATES ALL INTRON COMBINATIONS WITH AND
        WITHOUT 'EDGE' INTRONS.

        alignment_collections: list of intron tuples (pos, end pos)
        separation: number of bases at or above which introns should be
            separated

        Return value: list of lists of introns.
    """
    if not introns: return []
    introns.sort()
    prereturn = [[introns[0]]]
    for i in xrange(1, len(introns)):
        if introns[i][0] - prereturn[-1][-1][1] >= separation:
            prereturn.append([introns[i]])
        else:
            prereturn[-1].append(introns[i])
    to_return = set()
    for intron_combo in prereturn:
        to_return.add(tuple(intron_combo[:]))
        to_return.add(tuple(intron_combo[1:-1]))
        to_return.add(tuple(intron_combo[1:]))
        to_return.add(tuple(intron_combo[:-1]))
    return [list(intron_combo) for intron_combo in to_return if intron_combo]

def go(input_stream=sys.stdin, output_stream=sys.stdout, fudge=5,
        stranded=False, verbose=False, report_multiplier=1.2):
    """ Emits intron combinations associated with reads.

        Soft-clipped Bowtie 2 alignments of read sequences to the transcript
        fragment index are used infer which cointrons could possibly be
        overlapped by reads. Then maximal cliques of the graph described in
        the maximal_cliques() function are enumerated to obtain which
        intron combinations could possibly be overlapped by reads.

        input_stream: where to retrieve Bowtie 2 output
        output_stream: where to emit exon and intron tuples; typically, this is
            sys.stdout.
        verbose: True if alignments should occasionally be written to stderr.
        stranded: True iff input reads are strand-specific; this affects
            whether an output partition has a terminal '+' or '-' indicating
            the sense strand. Further, if stranded is True, an alignment is
            returned only if its strand agrees with the intron's strand.
        fudge: by how many bases to extend left and right extend sizes
            to accommodate potential indels
        report_multiplier: if verbose is True, the line number of an
            alignment written to stderr increases exponentially with base
            report_multiplier.
    """
    output_line_count, next_report_line, i = 0, 0, 0
    for (qname,), xpartition in xstream(input_stream, 1):
        '''While labeled multireadlet, this list may end up simply a
        unireadlet.'''
        multiread = []
        for tokens in xpartition:
            flag = int(tokens[0])
            if verbose and next_report_line == i:
                print >>sys.stderr, \
                    'SAM output record %d: rdname="%s", flag=%d' % (i,
                                                                    qname,
                                                                    flag)
                next_report_line = int((next_report_line + 1)
                                        * report_multiplier + 1) - 1
            i += 1
            multiread.append((qname,) + tokens)
        if flag & 4: continue
        corrected_multiread = multiread_with_introns(multiread,
                                                        stranded)
        all_introns = {}
        for alignment in multiread_with_introns(multiread, stranded):
            cigar = alignment[5]
            md = [field for field in alignment
                    if field[:5] == 'MD:Z:'][0][5:]
            pos = int(alignment[3])
            seq = alignment[9]
            reversed_complement_seq = seq[::-1].translate(
                    _reversed_complement_translation_table
                )
            if seq < reversed_complement_seq:
                seq_to_print = seq
            else:
                seq_to_print = reversed_complement_seq
            seq_size = len(seq)
            rname = alignment[2]
            sense = [field for field in alignment
                        if field[:5] == 'XS:A:'][0][5:]
            if (rname, sense) not in all_introns:
                all_introns[(rname, sense)] = defaultdict(list)
            _, _, introns, _ = indels_introns_and_exons(
                                                cigar, md, pos, seq
                                            )
            for intron in introns:
                if (intron[0], intron[1]) \
                    not in all_introns[(rname, sense)]:
                    all_introns[(rname, sense)][(intron[0], intron[1])] \
                        = [intron[2], intron[3]]
                else:
                    all_introns[(rname, sense)][
                            (intron[0], intron[1])
                        ][0] = max(all_introns[(rname, sense)][
                                (intron[0], intron[1])
                            ][0], intron[2])
                    all_introns[(rname, sense)][
                            (intron[0], intron[1])
                        ][1] = max(all_introns[(rname, sense)][
                                (intron[0], intron[1])
                            ][1], intron[3])
        for rname, sense in all_introns:
            to_write = set()
            # Grab maximal cliques
            for clique in \
                maximal_cliques(all_introns[(rname, sense)].keys()):
                for cointrons in separated_introns(
                            clique,
                            separation=(seq_size + fudge)
                        ):
                    cointrons.sort()
                    left_extend_size = all_introns[(rname, sense)][
                                            (cointrons[0][0],
                                                cointrons[0][1])
                                        ][0]
                    right_extend_size = all_introns[(rname, sense)][
                                             (cointrons[-1][0],
                                                cointrons[-1][1])
                                         ][1]
                    to_write.add(('{rname}{sense}\t{start}'
                           '\t{other_starts}'
                           '\t{ends}\t{left_size}'
                           '\t{right_size}\t{seq}').format(
                                rname=rname,
                                sense=sense,
                                start=cointrons[0][0],
                                other_starts=(
                                        ','.join(
                                        [str(intron[0]) for intron
                                            in cointrons[1:]]
                                    ) if len(cointrons) > 1 else '\x1c'
                                ),
                                ends=','.join(
                                        [str(intron[1])
                                            for intron in cointrons]
                                    ),
                                left_size=(left_extend_size
                                            + fudge),
                                right_size=(right_extend_size
                                            + fudge),
                                seq=seq_to_print
                           ))
            for line_to_write in to_write:
                print line_to_write
                output_line_count += 1
    output_stream.flush()
    print >>sys.stderr, ('cointron_enum_delegate.py reports %d output lines.'
                            % output_line_count)

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__, 
                formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--verbose', action='store_const', const=True,
        default=False,
        help='Print out extra debugging statements')
    parser.add_argument('--fudge', type=int, required=False,
        default=5,
        help='Permits a sum of exonic bases for an intron combo to be within '
             'the specified number of bases of a read sequence\'s size; '
             'this allows for indels with respect to the reference')
    parser.add_argument(
        '--stranded', action='store_const', const=True, default=False,
        help='Assume input reads come from the sense strand; then partitions '
             'in output have terminal + and - indicating sense strand')
    parser.add_argument('--report-multiplier', type=float, required=False,
        default=1.2,
        help='When --verbose is also invoked, the only lines of lengthy '
             'intermediate output written to stderr have line number that '
             'increases exponentially with this base')
    args = parser.parse_args()

    go(stranded=args.stranded, fudge=args.fudge,
        verbose=args.verbose, report_multiplier=args.report_multiplier)