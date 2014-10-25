#!/usr/bin/env python
"""
Rail-RNA-cointron_enum

Follows Rail-RNA-intron_index
Precedes Rail-RNA-cointron_fasta

Alignment script for MapReduce pipelines that wraps Bowtie. Finds introns
that cooccur on reads by local alignments to transcriptome elements from
Bowtie 2.

Input (read from stdin)
----------------------------
Single input tuple column:
1. SEQ or its reversed complement, whichever is first in alphabetical order
    --- should be unique for efficient processing, but this isn't required

Hadoop output (written to stdout)
----------------------------
Tab-delimited tuple columns:
1. Reference name (RNAME in SAM format) + 
    '+' or '-' indicating which strand is the sense strand
2. First intron start position in configuration
3. Rest of intron start positions in configuration or '\x1c' if there are none
4. Comma-separated list of intron end positions in configuration
5. left_extend_size: by how many bases on the left side of an intron the
    reference should extend
6. right_extend_size: by how many bases on the right side of an intron the
    reference should extend
7. Read sequence

ALL OUTPUT COORDINATES ARE 1-INDEXED.
"""
import sys
import os
import site
import threading
import subprocess
import tempfile
import atexit
from collections import defaultdict

base_path = os.path.abspath(
                    os.path.dirname(os.path.dirname(os.path.dirname(
                        os.path.realpath(__file__)))
                    )
                )
utils_path = os.path.join(base_path, 'rna', 'utils')
site.addsitedir(utils_path)
site.addsitedir(base_path)

import bowtie
from dooplicity.tools import xstream
from alignment_handlers import multiread_with_introns, indels_introns_and_exons

# Initialize global variables for tracking number of input/output lines
_input_line_count = 0
_output_line_count = 0

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

def handle_temporary_directory(temp_dir_path):
    """ Deletes temporary directory.

        temp_dir_paths: path to temporary directory to delete

        No return value.
    """
    import shutil
    shutil.rmtree(temp_dir_path)

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

class BowtieOutputThread(threading.Thread):
    """ Processes Bowtie alignments, emitting tuples for exons and introns. """
    
    def __init__(self, input_stream, return_set, output_stream=sys.stdout,
        verbose=False, report_multiplier=1.2, stranded=False, fudge=5):
        """ Constructor for BowtieOutputThread.

            input_stream: where to retrieve Bowtie's SAM output, typically a
                Bowtie process's stdout.
            output_stream: where to emit exon and intron tuples; typically,
                this is sys.stdout.
            return_set: 0 is added to this set if a thread finishes
                successfully
            verbose: True if alignments should occasionally be written 
                to stderr.
            bin_size: genome is partitioned in units of bin_size for later load
                balancing.
            report_multiplier: if verbose is True, the line number of an
                alignment written to stderr increases exponentially with base
                report_multiplier.
            fudge: by how many bases to extend left and right extend sizes
                to accommodate potential indels
        """
        super(BowtieOutputThread, self).__init__()
        self.daemon = True
        self.input_stream = input_stream
        self.output_stream = output_stream
        self.verbose = verbose
        self.report_multiplier = report_multiplier
        self.return_set = return_set
        self.stranded = stranded
        self.fudge = fudge

    def run(self):
        """ Prints exons for end-to-end alignments.

            Overrides default method containing thread activity.

            No return value.
        """
        global _output_line_count
        next_report_line, i = 0, 0
        for (qname,), xpartition in xstream(self.input_stream, 1):
            '''While labeled multireadlet, this list may end up simply a
            unireadlet.'''
            multiread = []
            for tokens in xpartition:
                flag = int(tokens[0])
                if self.verbose and next_report_line == i:
                    print >>sys.stderr, \
                        'SAM output record %d: rdname="%s", flag=%d' % (i,
                                                                        qname,
                                                                        flag)
                    next_report_line = int((next_report_line + 1)
                                            * self.report_multiplier + 1) - 1
                i += 1
                multiread.append((qname,) + tokens)
            if flag & 4: continue
            corrected_multiread = multiread_with_introns(multiread,
                                                            self.stranded)
            all_introns = {}
            for alignment in multiread_with_introns(multiread, self.stranded):
                cigar = alignment[5]
                md = [field for field in alignment
                        if field[:5] == 'MD:Z:'][0][5:]
                pos = int(alignment[3])
                seq = alignment[9]
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
                                separation=(seq_size + self.fudge)
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
                                                + self.fudge),
                                    right_size=(right_extend_size
                                                + self.fudge),
                                    seq=seq
                               ))
                for line_to_write in to_write:
                    print line_to_write
                    _output_line_count += 1
        self.output_stream.flush()
        self.return_set.add(0)

def go(input_stream=sys.stdin, output_stream=sys.stdout, bowtie2_exe='bowtie2',
    bowtie2_index_base='genome', bowtie2_args='', verbose=False,
    report_multiplier=1.2, stranded=False, fudge=5, score_min=60):
    """ Runs Rail-RNA-cointron_enum 

        Alignment script for MapReduce pipelines that wraps Bowtie. Finds
        introns that cooccur on reads by local alignments to transcriptome
        elements from Bowtie 2.

        Input (read from stdin)
        ----------------------------
        Tab-delimited output tuple columns (readletize)
        1. SEQ or its reversed complement, whichever is first in alphabetical
            order
        2. Comma-separated list of sample labels if field 1 is the read
            sequence; '\x1c' if empty
        3. Comma-separated list of sample labels if field 1 is the reversed
            complement of the read sequence; '\x1c' if empty

        Hadoop output (written to stdout)
        ----------------------------
        Tab-delimited tuple columns:
        1. Reference name (RNAME in SAM format) + 
            '+' or '-' indicating which strand is the sense strand
        2. First intron start position in configuration
        3. Rest of intron start positions in configuration or '\x1c' if there
            are none
        4. Comma-separated list of intron end positions in configuration
        5. left_extend_size: by how many bases on the left side of an intron
            the reference should extend
        6. right_extend_size: by how many bases on the right side of an intron
            the reference should extend
        7. Read sequence

        input_stream: where to find input reads.
        output_stream: where to emit exonic chunks and introns.
        bowtie2_exe: filename of Bowtie 2 executable; include path if not in
            $PATH.
        bowtie2_index_base: the basename of the Bowtie index files associated
            with the reference.
        bowtie2_args: string containing precisely extra command-line arguments
            to pass to Bowtie 2, e.g., "--tryhard --best"; or None.
        verbose: True iff more informative messages should be written to
            stderr.
        report_multiplier: if verbose is True, the line number of an alignment
            written to stderr increases exponentially with base
            report_multiplier.
        stranded: True iff input reads are strand-specific; this affects
            whether an output partition has a terminal '+' or '-' indicating
            the sense strand. Further, if stranded is True, an alignment is
            returned only if its strand agrees with the intron's strand.
        fudge: by how many bases to extend left and right extend sizes
                to accommodate potential indels
        score_min: Bowtie2 CONSTANT minimum alignment score

        No return value.
    """
    global _input_line_count
    # For storing long qnames
    temp_dir = tempfile.mkdtemp()
    atexit.register(handle_temporary_directory, temp_dir)
    reads_file = os.path.join(temp_dir, 'reads.temp')
    output_file = os.path.join(temp_dir, 'out.sam')
    with open(reads_file, 'w') as read_stream:
        for _input_line_count, line in enumerate(input_stream):
            seq = line.strip()
            print >>read_stream, \
                '\t'.join([str(_input_line_count), seq, 'I'*len(seq)])
    bowtie_command = ' '.join([bowtie2_exe,
        bowtie2_args if bowtie2_args is not None else '',
        ' --local -t --no-hd --mm -x', bowtie2_index_base, '--12',
        reads_file, '-S', output_file, '--score-min L,%d,0' % score_min, 
        '-D 24 -R 3 -N 1 -L 20 -i L,4,0'])
    print >>sys.stderr, 'Starting Bowtie2 with command: ' + bowtie_command
    # Because of problems with buffering, write output to file
    bowtie_process = subprocess.Popen(bowtie_command, bufsize=-1,
        stdout=subprocess.PIPE, stderr=sys.stderr, shell=True)
    bowtie_process.wait()
    if os.path.exists(output_file):
        return_set = set()
        output_thread = BowtieOutputThread(
                open(output_file),
                return_set,
                verbose=verbose, 
                output_stream=output_stream,
                report_multiplier=report_multiplier,
                stranded=stranded,
                fudge=fudge
            )
        output_thread.start()
        # Join thread to pause execution in main thread
        if verbose: print >>sys.stderr, 'Joining thread...'
        output_thread.join()
        if not return_set:
            raise RuntimeError('Error occurred in BowtieOutputThread.')

if __name__ == '__main__':
    import argparse
    # Print file's docstring if -h is invoked
    parser = argparse.ArgumentParser(description=__doc__, 
                formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--report_multiplier', type=float, required=False,
        default=1.2,
        help='When --verbose is also invoked, the only lines of lengthy '
             'intermediate output written to stderr have line number that '
             'increases exponentially with this base')
    parser.add_argument('--verbose', action='store_const', const=True,
        default=False,
        help='Print out extra debugging statements')
    parser.add_argument('--test', action='store_const', const=True,
        default=False,
        help='Run unit tests; DOES NOT NEED INPUT FROM STDIN')
    parser.add_argument('--keep-alive', action='store_const', const=True,
        default=False,
        help='Periodically print Hadoop status messages to stderr to keep ' \
             'job alive')
    parser.add_argument('--fudge', type=int, required=False,
        default=5,
        help='Permits a sum of exonic bases for an intron combo to be within '
             'the specified number of bases of a read sequence\'s size; '
             'this allows for indels with respect to the reference')
    parser.add_argument(
        '--stranded', action='store_const', const=True, default=False,
        help='Assume input reads come from the sense strand; then partitions '
             'in output have terminal + and - indicating sense strand')
    parser.add_argument('--score-min', type=int, required=False,
        default=48,
        help='Bowtie2 minimum CONSTANT score to use')

    # Add command-line arguments for dependencies
    bowtie.add_args(parser)

    # Collect Bowtie arguments, supplied in command line after the -- token
    argv = sys.argv
    bowtie2_args = ''
    in_args = False
    for i, argument in enumerate(sys.argv[1:]):
        if in_args:
            bowtie2_args += argument + ' '
        if argument == '--':
            argv = sys.argv[:i + 1]
            in_args = True

    '''Now collect other arguments. While the variable args declared below is
    global, properties of args are also arguments of the go() function so
    different command-line arguments can be passed to it for unit tests.'''
    args = parser.parse_args(argv[1:])

    # Start keep_alive thread immediately
    if args.keep_alive:
        from dooplicity.tools import KeepAlive
        keep_alive_thread = KeepAlive(sys.stderr)
        keep_alive_thread.start()
    
if __name__ == '__main__' and not args.test:
    import time
    start_time = time.time()
    go(bowtie2_exe=args.bowtie2_exe,
        bowtie2_index_base=args.bowtie2_idx,
        bowtie2_args=bowtie2_args, 
        verbose=args.verbose,
        report_multiplier=args.report_multiplier,
        stranded=args.stranded,
        fudge=args.fudge,
        score_min=args.score_min)
    print >>sys.stderr, 'DONE with cointron_enum.py; in/out=%d/%d; ' \
        'time=%0.3f s' % (_input_line_count, _output_line_count,
                            time.time() - start_time)
elif __name__ == '__main__':
    # Test units
    del sys.argv[1:] # Don't choke on extra command-line parameters
    import unittest
    # Add unit tests here
    unittest.main()