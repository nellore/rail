#!/usr/bin/env python
"""
Rail-RNA-intron

Follows Rail-RNA-align
Precedes Rail-RNA-intron_post

Reduce step in MapReduce pipelines that infers intron start/end positions from 
spliced alignments output by Rail-RNA-align. Each worker operates on a set of
genome partitions, clustering alignments that overlap and associating clusters
with donor/acceptor motifs.

If the input reads are NOT stranded, then all relevant splice junction
(donor/acceptor) motifs must be permitted, including, e.g. both GT-AG and
CT-AC. If the input reads are stranded, motifs that are inconsistent with the
strand of the flanking aligned sequences can be ignored. The --stranded option
handles this.

Input (read from stdin)
----------------------------
Tab-delimited columns:
1. Reference name (RNAME in SAM format) + ';' + partition number +  
    ('+' or '-' indicating which strand is the sense strand if input reads are
        strand-specific -- that is, --stranded in Rail-RNA-align was invoked;
        otherwise, there is no terminal '+' or '-')
2. Candidate intron start (inclusive) on forward strand (1-indexed)
3. Candidate intron end (exclusive) on forward strand (1-indexed)

Input is partitioned by column 1.

Hadoop output (written to stdout)
----------------------------
Tab-delimited columns:
1. The character '-', enforcing a single partition.
2. The character 'i', which will place the row after 'a' (maximum read length
    rows) in lexicographic sort order.
3. Reference name (RNAME in SAM format) +
    '+' or '-' indicating which strand is the sense strand
4. Intron start position (inclusive)
5. Intron end position (exclusive)

OUTPUT COORDINATES ARE 1-INDEXED.
"""
import os
import sys
# Regular expressions are used to identify splice-site motifs
import re
import itertools
from collections import defaultdict
import numpy as np
import site
import random
import time

base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for directory_name in ['interval', 'bowtie']:
    site.addsitedir(os.path.join(base_path, directory_name))

import partition
import bowtie
import bowtie_index

'''Initialize lists of donor/acceptor motifs in order of descending priority;
the unstranded motif list is used if --stranded is False; otherwise, 
_forward_strand_motifs are used if the sense strand is the forward strand and
_reverse_strand_motifs are used if the sense strand is the reverse strand.
Each tuple in a given sublist denotes (donor motif, acceptor motif, True iff
    sense strand is reverse strand). The tuples in a sublist have the same
priority (rank).'''
_unstranded_motifs = [[('GT', 'AG', False), ('CT', 'AC', True)],
                        [('GC', 'AG', False), ('CT', 'GC', True)],
                        [('AT', 'AC', False), ('GT', 'AT', True)]]
_forward_strand_motifs = [('GT', 'AG', False), ('GC', 'AG', False),
                            ('AT', 'AC', False)]
_reverse_strand_motifs = [('CT', 'AC', True), ('CT', 'GC', True),
                            ('GT', 'AT', True)]

def intron_clusters_in_partition(candidate_introns, partition_start, 
    partition_end, cluster_radius=5, verbose=False):
    """ Clusters candidate introns from a strand in a genomic partition.

        Each candidate intron is specified by its start position start_position
        and its end position end_position (exclusive) and is associated with a
        count read_count of the number of reads supporting the candidate.
        The algorithm iterates through a dynamic lineup of candidates. The
        lineup is initially the list of candidates in order of descending
        read_count. A cluster is formed by associating a given candidate C
        under examination with candidates {C_i} of the same size that lie
        within min(intron_size, cluster_radius) bases of C. The {C_i} are then
        removed from the lineup, and the next candidate is examined.
        
        candidate_introns: a dictionary. Each key is a tuple
            (start_position, end_position) and its corresponding value is
            supporting_read_count, the number of reads (across samples)
            supporting the candidate intron.
        partition_start: start position (inclusive) of partition.
        partition_end: end position (exclusive) of partition.
        cluster_radius: distance from a candidate intron under examination
            within which to search for candidate introns in the same cluster.
            See above for a detailed description.
        verbose: True iff counts of possible splice junctions and clusters 
            should be written to stderr.

        Return value: a list of lists, each of which corresponds to a cluster
            of candidate introns. Each item in a cluster is a tuple
            (start_position, end_position, supporting_read_count) corresponding
            to an intron call agreed on by a set of reads.
    """
    '''Construct list of candidate introns sorted in order of descending
    read_count.'''
    candidate_intron_list = [(supporting_read_count, 
                                end_position - start_position, end_position)
                                for ((start_position, end_position),
                                supporting_read_count)
                                in candidate_introns.items()]
    candidate_intron_list.sort()
    # Construct set of candidate introns for fast searching
    candidate_intron_set = set()
    for _, intron_size, end_position in candidate_intron_list:
        candidate_intron_set.add((intron_size, end_position))
    candidate_intron_count = len(candidate_intron_set)
    # Initialize set for storing candidate introns already clustered
    clustered_introns = set()
    intron_clusters = []
    total_cluster_count = 0
    for _, intron_size, end_position in candidate_intron_list:
        if (intron_size, end_position) not in clustered_introns:
            current_cluster_radius = min(intron_size, cluster_radius)
            total_cluster_count += 1
            clustered_introns.add((intron_size, end_position))
            intron_cluster = [(end_position - intron_size, end_position)]
            for an_end_position in xrange(end_position 
                                            - current_cluster_radius,
                                            end_position
                                            + current_cluster_radius + 1):
                if (intron_size, an_end_position) in candidate_intron_set \
                    and (intron_size,
                            an_end_position) not in clustered_introns:
                    '''If a nearby candidate hasn't been clustered, absorb it
                    into the current cluster.'''
                    intron_cluster.append((an_end_position - intron_size,
                                                an_end_position))
                    clustered_introns.add((intron_size, an_end_position))
            intron_clusters.append(intron_cluster)
            assert len(clustered_introns) <= candidate_intron_count
            if len(clustered_introns) == candidate_intron_count:
                if verbose:
                    print >> sys.stderr, '%d possible splice junction(s) ' \
                        'clustered down to %d.' % (candidate_intron_count, 
                            total_cluster_count)
                return [[an_intron + (candidate_introns[an_intron],) 
                            for an_intron in an_intron_cluster] 
                            for an_intron_cluster in intron_clusters]
    raise RuntimeError('For loop should never be completed.')

def ranked_splice_sites_from_cluster(reference_index, intron_cluster,
    rname, motifs, motif_radius=1, verbose=False):
    """ Ranks possible splice sites for a cluster using donor/acceptor motifs.

        Consider the following cluster of three candidate introns for which the
        sense strand is the forward strand, used here as an illustrative
        example.

                     5'                                                 3'
        I1              ============================================
        I2                ============================================
        I3                  =========================================
        Reference       GCGGGTAATAG................................AGAATA

        Call the length of the full reference strand reference_length. Consider 
        the minimum and maximum start positions (inclusive) S_min and S_max of 
        the candidate introns. Above, I1 sets S_min and I3 sets S_max. Consider 
        also the minimum and maximum end positions (exclusive) E_min and E_max
        of the candidate introns. Above, I1 sets E_min and I2 sets E_max. The
        algo scans for splice-site motifs GT..AG, GC..AG, and AT..AC on the
        corresponding intervals
        [S_min, min(S_max + 2, reference_length))
        ..
        [min(E_min - 2, 0), E_max) . (These splice-site motifs assume the
        sense strand is the forward strand. If instead the sense strand is the
        reverse strand, the motifs are CT..AC, CT..GC, and GT..AT.) Splice
        sites are then ranked according to motif:
                        sense strand = +        sense strand = -
                1.          GT..AG                  CT..AC
                2.          GC..AG                  CT..GC
                3.          AT..AC                  GT..AT
        (References: 
            -Turunen JJ, Niemela EH, Verma B, Frilander MJ. The significant
            other: splicing by the minor spliceosome. Wiley Interdiscip Rev
            RNA. 2013 Jan-Feb;4(1):61-76. doi: 10.1002/wrna.1141.

            -http://onlinelibrary.wiley.com/doi/10.1002/wrna.1141/full)
        The only motifs considered appear at the ends of an intron of exactly
        the (unique) length of the introns in a given cluster. If more than one
        possible motif (pair) of the same rank is found, the tie is broken as
        follows. The means and standard deviations of the start and end
        positions of the candidate introns are computed. Then the z-score of
        the start motif position is added to the z-score of the end motif
        position, and the motif pairs are ranked in order of ascending z-score.
        If no motifs are found, the return value of this function is an empty
        list. So above, ranked first would be a GT..AG, which is spanned by I3,
        and second would be a GC..AG, spanned by
        [I1's start position, I3's end position].

        ALL INPUT COORDINATES ARE ASSUMED TO BE 1-INDEXED.

        reference_index: object of class bowtie_index.BowtieIndexReference
            that permits access to reference; used to find splice-site motifs.
        intron_cluster: a list of lists, each of which corresponds to a cluster
            of candidate introns. Each item in a cluster is a tuple
            (start_position, end_position, supporting_read_count) corresponding
            to an intron call supported by several reads.
        rname: SAM-format RNAME indicating the chromosome.
        motifs: List of lists of motif tuples (donor motif, acceptor motif, 
            reverse strand boolean) in order of descending priority. The
            reverse strand boolean is False iff the forward strand is the 
            sense strand for the motif. Each sublist contains motif tuples
            whose priority (rank) is the same.
            Example tuple: ('GT', 'AG', False).
        motif_radius: distance (in bp) from each of the start and end positions
            of a cluster within which to search for motifs; that is, the
            intervals [min(start_positions) - motif_radius,
                max(start_positions) + motif_radius + 2) and
            [min(end_positions - 2 - motif_radius),
                max(end_positions) + motif_radius) are scanned for motifs,
            where start_positions (end_positions) aggregates the start (end)
            positions of the candidate introns in a cluster.
        verbose: If True, writes to stderr when no splice sites are identified
            for a cluster.

        Return value: List of tuples representing possible final introns in
            order of descending rank. Each tuple is of the form
            (start position, end position, summed z-score,
                left motif, right motif, reverse strand boolean). Example:
            [(142451, 143128, 2.8324, 'GT', 'AG', False),
                (142449, 143128, 3.1124, 'CT', 'AC', True)].
    """
    start_positions, end_positions, weights \
        = zip(*intron_cluster)
    denominator = sum(weights)
    weights = [float(count) / denominator for count in weights]
    reference_length = reference_index.rname_lengths[rname]
    min_start_position = max(min(start_positions) - motif_radius, 1)
    max_start_position = min(max(start_positions) + motif_radius + 2,
                                reference_length)
    min_end_position = max(min(end_positions) - 2 - motif_radius, 1)
    max_end_position = min(max(end_positions) + motif_radius, reference_length)
    left_sequence = reference_index.get_stretch(rname, min_start_position - 1,
        max_start_position - min_start_position)
    right_sequence = reference_index.get_stretch(rname, min_end_position - 1,
        max_end_position - min_end_position)
    assert max_end_position >= min_end_position and \
        max_start_position >= min_start_position
    # For computing z-scores
    mean_start_position = np.average(start_positions, weights=weights)
    mean_end_position = np.average(end_positions, weights=weights)
    # Maxes below avoid future ZeroDivisionError exceptions
    stdev_start_position = max(np.sqrt(
                                np.average(
                                    (np.array(start_positions)
                                        -mean_start_position)**2, 
                                    weights=weights)
                                ), 1e-6)
    stdev_end_position = max(np.sqrt(
                                np.average(
                                    (np.array(end_positions)
                                        -mean_end_position)**2, 
                                    weights=weights)
                                ), 1e-6)
    # Initialize list for storing ranked intron start/end positions
    ranked_introns = []
    for motif_priority_class in motifs:
        positions_and_z_scores = []
        for motif in motif_priority_class:
            '''Use regex lookahead to identify possibly overlapping motifs.
            Each *_offset record offset from beginning of left_sequence or
            right_sequence.'''
            left_motif_offsets = [a_match.start() for a_match in 
                                    re.finditer(r'(?=(%s))' % motif[0],
                                                    left_sequence)]
            right_motif_offsets = [a_match.start() for a_match in 
                                    re.finditer(r'(?=(%s))' % motif[1],
                                                    right_sequence)]
            '''Find all possible combinations of left and right offsets for a
            given motif (pair).'''
            motif_pairs = \
                itertools.product(*[left_motif_offsets, right_motif_offsets])
            for left_motif_offset, right_motif_offset in motif_pairs:
                if left_motif_offset != right_motif_offset: 
                    '''Eliminate motifs that do not yield introns of the 
                    proper length.'''
                    continue
                left_motif_start_position = min_start_position \
                    + left_motif_offset
                right_motif_end_position = min_end_position \
                    + right_motif_offset + 2
                if right_motif_end_position \
                    < left_motif_start_position + 4:
                    # Filter out overlapping donor/acceptor motifs
                    continue
                z_score_sum = abs(left_motif_start_position 
                    - mean_start_position) / float(stdev_start_position) \
                    + abs(right_motif_end_position - mean_end_position) / \
                    float(stdev_end_position)
                positions_and_z_scores.append((left_motif_start_position,
                                                right_motif_end_position,
                                                z_score_sum) + motif)
        # Sort all matches of same priority (rank) by z-score sum
        positions_and_z_scores.sort(
                key=lambda positions_and_z_score: positions_and_z_score[2]
            )
        ranked_introns += positions_and_z_scores
    if len(ranked_introns) == 0 and verbose:
        print >>sys.stderr, \
            'Warning: No splice site found for cluster with ' \
            '%d candidate introns' % len(intron_cluster)
    return ranked_introns

def go(bowtie_index_base="genome", input_stream=sys.stdin,
        output_stream=sys.stdout, bin_size=10000, cluster_radius=5,
        stranded=False, intron_partition_overlap=20, motif_radius=1, 
        min_anchor_significance=9, seed=0, verbose=False):
    """ Runs Rail-RNA-intron.

        Input lines are binned, so they are examined two at a time. When the
        first column of the input --- the partition --- differs between the
        two lines, all the candidate introns starting in a given partition
        have been collected. The algo then clusters candidate introns within
        a partition (see intron_clusters_in_partition()), ranks each cluster's
        possible splice sites by searching for nearby donor/acceptor motifs
        (see ranked_splice_sites_from_cluster()), and writes output in any
        or all of three formats.

        If the input reads are NOT stranded, then all relevant splice junction
        (donor/acceptor) motifs must be permitted, including, e.g. both GT-AG
        and CT-AC. If the input reads are stranded, motifs that are
        inconsistent with the strand of the flanking aligned sequences can be
        ignored. The --stranded option handles this.

        Input (read from stdin)
        ----------------------------
        Tab-delimited columns:
        1. Reference name (RNAME in SAM format) + ';' + partition number +  
            ('+' or '-' indicating which strand is the sense strand if input
            reads are strand-specific -- that is, --stranded in Rail-RNA-align
            was invoked; otherwise, there is no terminal '+' or '-')
        2. Candidate intron start (inclusive) on forward strand (1-indexed)
        3. Candidate intron end (exclusive) on forward strand (1-indexed)

        Input is partitioned by column 1.

        Hadoop output (written to stdout)
        ----------------------------
        Tab-delimited columns:
        1. The character '-', enforcing a single partition.
        2. The character 'i', which will place the row after 'a' (maximum read length
            rows) in lexicographic sort order.
        3. Reference name (RNAME in SAM format) +
            '+' or '-' indicating which strand is the sense strand
        4. Intron start position (inclusive)
        5. Intron end position (exclusive)

        OUTPUT COORDINATES ARE 1-INDEXED.

        bowtie_index_base: the basename of the Bowtie index files associated
            with the reference.
        input_stream: where to find input.
        output_stream: where to emit splice sites.
        bin_size: genome is partitioned in units of bin_size, and candidate
            introns from input are in first-column partitions numbered
            accordingly.
        cluster_radius: distance from a candidate intron under examination
            within which to search for candidate introns in the same cluster.
            See intron_clusters_in_partition()'s docstring for more
            information.
        stranded: True iff input reads are strand-specific, and lines read
            from input_stream have first column with terminal '+' or '-'
            indicating which strand is the sense strand.
        intron_partition_overlap: number of bases to subtract from
            reference start position of candidate intron when determining
            genome partition it is in.
        motif_radius: distance (in bp) from each of the start and end
            positions of a cluster within which to search for motifs.
        verbose: True iff informative messages should be written to stderr.

        No return value.
    """
    start_time = time.time()
    reference_index = bowtie_index.BowtieIndexReference(bowtie_index_base)
    input_line_count, output_line_count = 0, 0
    handle_partition = False
    last_partition_id = None
    candidate_introns = {}
    # Make results reproducible
    random.seed(seed)
    while True:
        line = input_stream.readline()
        if line:
            input_line_count += 1
            tokens = line.rstrip().split('\t')
            assert len(tokens) == 3
            partition_id, pos, end_pos = tokens
            pos, end_pos = int(pos), int(end_pos)
            assert end_pos > pos, line
            if stranded:
                '''If input reads are strand-specific, the partition_id has a
                terminal '+' or '-' indicating which strand is the sense
                strand. Make it the reverse_strand'''
                rname, partition_start, partition_end = \
                    partition.parse(partition_id[:-1], bin_size)
                reverse_strand = (True if partition_id[-1] == '-' else False)
            else:
                rname, partition_start, partition_end = \
                    partition.parse(partition_id, bin_size)
            assert pos >= partition_start - intron_partition_overlap and \
                pos < partition_end + intron_partition_overlap, \
                'Intron start %d is not in partition [%d, %d)' \
                ', partition id=%s' \
                % (pos, partition_start, partition_end, partition_id)
            if last_partition_id is not None and \
                last_partition_id != partition_id:
                handle_partition = True
        else:
            # If there's no next line, handle the final partition
            handle_partition = True
        if handle_partition and len(candidate_introns):
            if verbose:
                print >> sys.stderr, 'For partition %s:[%d, %d)' \
                    % (last_partition_id, last_partition_start,
                        last_partition_end)
            intron_clusters = intron_clusters_in_partition(candidate_introns,
                last_partition_start, last_partition_end, 
                cluster_radius=cluster_radius,
                verbose=verbose)
            cluster_splice_sites = {}
            if stranded:
                # The sense strand is known, so narrow motif set used
                if last_reverse_strand:
                    for intron_cluster in intron_clusters:
                        ranked_splice_sites = ranked_splice_sites_from_cluster(
                                reference_index, intron_cluster, last_rname,
                                _reverse_strand_motifs,
                                motif_radius=motif_radius, verbose=verbose
                            )
                        if len(ranked_splice_sites) != 0:
                            # Break any ties in top-ranked intron at random
                            cluster_splice_site = \
                                random.sample([splice_site for splice_site in
                                                    ranked_splice_sites if
                                                    ranked_splice_sites[0][2:]
                                                    == splice_site[2:]], 1)[0]
                            intron_strand = ('-' if cluster_splice_site[-1]
                                                else '+')
                            print >>output_stream, \
                                'intron\t-\ti\t%s%s\t%012d\t%012d' \
                                % (last_rname, intron_strand, 
                                    cluster_splice_site[0],
                                    cluster_splice_site[1])
                            output_line_count += 1
                else:
                    for intron_cluster in intron_clusters:
                        ranked_splice_sites = ranked_splice_sites_from_cluster(
                                reference_index, intron_cluster, last_rname,
                                _forward_strand_motifs,
                                motif_radius=motif_radius, verbose=verbose
                            )
                        if len(ranked_splice_sites) != 0:
                            # Break any ties in top-ranked intron at random
                            cluster_splice_site = \
                                random.sample([splice_site for splice_site in
                                                    ranked_splice_sites if
                                                    ranked_splice_sites[0][2:]
                                                    == splice_site[2:]], 1)[0]
                            intron_strand = ('-' if cluster_splice_site[-1]
                                                else '+')
                            print >>output_stream, \
                                'intron\t-\ti\t%s%s\t%012d\t%012d' \
                                % (last_rname, intron_strand, 
                                    cluster_splice_site[0],
                                    cluster_splice_site[1])
                            output_line_count += 1
            else:
                # The sense strand is unknown, so use a general motif set
                for intron_cluster in intron_clusters:
                    ranked_splice_sites = ranked_splice_sites_from_cluster(
                            reference_index, intron_cluster, last_rname,
                            _unstranded_motifs,
                            motif_radius=motif_radius, verbose=verbose
                        )
                    if len(ranked_splice_sites) != 0:
                        # Break any ties in top-ranked intron at random
                        cluster_splice_site = \
                            random.sample([splice_site for splice_site in
                                                ranked_splice_sites if
                                                ranked_splice_sites[0][2:]
                                                 == splice_site[2:]], 1)[0]
                        intron_strand = ('-' if cluster_splice_site[-1] 
                                            else '+')
                        print >>output_stream, \
                                    'intron\t-\ti\t%s%s\t%012d\t%012d' \
                                    % (last_rname, intron_strand, 
                                        cluster_splice_site[0],
                                        cluster_splice_site[1])
                        output_line_count += 1
            candidate_introns = {}
            handle_partition = False
        if line:
            candidate_introns[(pos, end_pos)] = \
                candidate_introns.get((pos, end_pos), 0) + 1
            (last_partition_id, last_partition_start, last_partition_end, 
                last_rname) \
            = (partition_id, partition_start, partition_end, rname)
            if stranded:
                last_reverse_strand = reverse_strand
        else: break

    print >> sys.stderr, 'DONE with intron.py; in/out=%d/%d; ' \
        'time=%0.3f s' % (input_line_count, output_line_count,
                                time.time() - start_time)

if __name__ == '__main__':
    import argparse
    # Print file's docstring if -h is invoked
    parser = argparse.ArgumentParser(description=__doc__, 
                formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--test', action='store_const', const=True,
        default=False,
        help='Run unit tests; DOES NOT NEED INPUT FROM STDIN, AND DOES NOT '
             'WRITE INTRONS TO STDOUT')
    parser.add_argument('--verbose', action='store_const', const=True,
        default=False,
        help='Print out extra debugging statements')
    parser.add_argument(\
        '--cluster-radius', type=int, required=False, default=20,
        help='The maximum radius of a cluster of candidate introns for which '
             'splice sites are called')
    parser.add_argument(\
        '--motif-radius', type=int, required=False, default=4,
        help='Distance (in bp) from each of the start and end positions '
             'of a cluster within which to search for motifs')
    parser.add_argument('--intron-partition-overlap', type=int, required=False,
        default=20, 
        help='Amount by which partitions overlap their left and right '
             'neighbors')
    parser.add_argument(\
        '--stranded', action='store_const', const=True, default=False,
        help='Assume input reads come from the sense strand')

    # Add command-line arguments for dependencies
    partition.addArgs(parser)
    bowtie.addArgs(parser)

    args = parser.parse_args(sys.argv[1:])

if __name__ == '__main__' and not args.test:
    go(bowtie_index_base=args.bowtie_idx,
        bin_size=args.partition_length,
        cluster_radius=args.cluster_radius,
        intron_partition_overlap=args.intron_partition_overlap,
        motif_radius=args.motif_radius,
        verbose=args.verbose)
elif __name__ == '__main__':
    # Test units
    del sys.argv[1:] # Don't choke on extra command-line parameters
    import unittest
    import shutil
    import tempfile
    import subprocess

    class TestIntronClustersInPartition(unittest.TestCase):
        """ Tests intron_clusters_in_partition(). """

        def test_that_overlapping_introns_of_equal_size_cluster(self):
            """ Fails if all introns are not placed in the same cluster
            """
            candidate_introns = {
                    (5, 25)  : 1,
                    (4, 24)  : 1,
                    (10, 30) : 1
                }
            clusters = intron_clusters_in_partition(candidate_introns,
                1, 40, cluster_radius=10, verbose=False)
            self.assertEquals(len(clusters), 1)
            self.assertEquals(
                    sorted(clusters[0]), 
                    [(4, 24),
                     (5, 25),
                     (10, 30)]
                )

        def test_that_overlapping_introns_of_different_size_separate(self):
            """ Fails if introns do not separate properly into two clusters.
            """
            candidate_introns = {
                    (5, 25)  : 1,
                    (4, 24)  : 1,
                    (10, 40) : 1
                }
            '''The third intron above has a different size than the first two,
            so it should be placed in its own cluster.'''
            clusters = sorted(intron_clusters_in_partition(candidate_introns,
                1, 40, cluster_radius=10, verbose=False), key=len)
            self.assertEquals(len(clusters), 2)
            self.assertEquals(
                    sorted(clusters[0]),
                    [(10, 40)]
                )
            self.assertEquals(
                    sorted(clusters[1]),
                    [(4, 24), 
                     (5, 25)]
                )

    class TestRankedSpliceSitesFromCluster(unittest.TestCase):
        """ Tests ranked_splice_sites_from_cluster(). """

        def setUp(self):
            """ Creates temporary directory and Bowtie index. """
            reference_seq = 'ATGGCATACGATACGTCAGACCATGCAggACctTTacCTACATACTG' \
                            'GTTACATAGTACATCTAGGCATACTACGTgaCATACGgaCTACGTAG' \
                            'GTCCAGATTACGATACAAaTACGAAcTCccATAGCAaCATaCTAGac' \
                            'CAttAaaGACTAGACTAACAGACAaAACTAGCATacGATCATGACaA' \
                            'ACGAGATCCATATAtTTAGCAaGACTAaACGATACGATACAGTACaA' \
                            'ATACAGaaATCAGaGCAGAAaATACAGATCAaAGCTAGCAaAAtAtA'
            self.temp_dir_path = tempfile.mkdtemp()
            fasta_file = os.path.join(self.temp_dir_path, 'test.fa')
            self.bowtie_build_base = os.path.join(self.temp_dir_path, 'test')
            fasta_stream = open(fasta_file, 'w')
            print >>fasta_stream, '>chr1'
            print >>fasta_stream, \
                '\n'.join([reference_seq[i:i+50] for i in range(0, 251, 50)])
            fasta_stream.close()
            bowtie_build_process = subprocess.call(
                    [args.bowtie_build_exe,
                    fasta_file,
                    self.bowtie_build_base],
                    stdout=open(os.devnull, 'w'),
                    stderr=subprocess.STDOUT
                )
            self.reference_index = bowtie_index.BowtieIndexReference(
                                    self.bowtie_build_base
                                   )

        def test_that_second_line_of_reference_is_called_as_an_intron(self):
            """ Fails if reference_seq's second line isn't called as an intron.
            """
            cluster = [(49, 96),
                       (48, 95),
                       (47, 94)]
            splice_sites = ranked_splice_sites_from_cluster(
                    self.reference_index, cluster,
                    'chr1', _unstranded_motifs, motif_radius=1, verbose=False
                )
            self.assertEquals((48, 95), splice_sites[0][:2])
            self.assertEquals(('GT', 'AG'), splice_sites[0][3:5])

        def tearDown(self):
            # Kill temporary directory
            shutil.rmtree(self.temp_dir_path)

    unittest.main()