#!/usr/bin/env python
"""
Rail-RNA-intron

Follows Rail-RNA-align
Precedes Rail-RNA-intron_post

Reduce step in MapReduce pipelines that infers intron/exon boundaries from 
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
2. Sample label
3. Candidate intron start (inclusive) on forward strand (1-indexed)
4. Candidate intron end (exclusive) on forward strand (1-indexed)
5. Number of nucleotides between 5' end of candidate intron and 5' end of read
from which it was inferred, ASSUMING THE SENSE STRAND IS THE FORWARD STRAND.
That is, if the sense strand is the reverse strand, this is the distance
between the 5' end of the reverse-complemented read and the 5' end of the
reverse-complemented intron.
6. Number of nucleotides between 3' end of candidate intron and 3' end of read 
from which it was inferred, ASSUMING THE SENSE STRAND IS THE FORWARD STRAND.
7. Number of nucleotides spanned by EC on the left (that is, towards the 5'
end of the read) of the intron, ASSUMING THE SENSE STRAND IS THE FORWARD
STRAND.
8. Number of nucleotides spanned by EC on the right (that is, towards the 3'
end of the read) of the intron, ASSUMING THE SENSE STRAND IS THE FORWARD
STRAND.
9. Match rate: the number of bases that match the reference per base of the
aligned readlets comprising the read from which the intron was inferred.
10. '-' if reversed complements of readlets from which intron was inferred
aligned to forward strand; else '+'

Hadoop output (written to stdout)
----------------------------
Format 1 (span): tab-delimited columns, one line per read spanning a splice
site:
1. SAM-format RNAME + ';' + ('+' or '-' indicating which strand is the sense
    strand)
2. Intron start (inclusive) on forward strand (1-indexed)
3. Intron end (exclusive) on forward strand (1-index)
4. Left motif (donor motif if sense strand is '+'; acceptor motif otherwise)
5. Right motif (acceptor motif if sense strand is '+'; donor motif otherwise)
6. Sample label from which read was derived

Format 2 (site): tab-delimited columns, one line per splice site spanned by
the reads in a sample
1. SAM-format RNAME + ';' + ('+' or '-' indicating which strand is the sense
    strand)
2. Intron start (inclusive) on forward strand (1-indexed)
3. Intron end (exclusive) on forward strand (1-index)
4. Left motif (donor motif if sense strand is '+'; acceptor motif otherwise)
5. Right motif (acceptor motif if sense strand is '+'; donor motif otherwise)
6. Sample label with nonzero number of reads spanning splice site
7. Number of reads in sample spanning splice site

Format 3 (junction):
12-column (3 required fields + 9 optional fields) BED output mimicking TopHat's
junctions.bed, except anchor significance, maximum match rate, and unique
displacement count are also included in the name field. From the TopHat manual
http://tophat.cbcb.umd.edu/manual.shtml: "Each junction consists of two
connected BED blocks, where each block is as long as the maximal overhang of
any read spanning the junction. The score is the number
of alignments spanning the junction." Anchor significance of a junction is
NOT defined as in the MapSplice paper
http://nar.oxfordjournals.org/content/38/18/e178.long. Consider the set {R_i}
of reads that span a junction J. For a given read R_i, consider the two anchors
A_i and B_i on either side of J, and let L_(A_i) and L_(B_i) be their lengths
(in bp). The anchor significance of J is max_i min(L(A_i), L(B_i)).) Maximum
match rate is the highest match rate from among the reads spanning a junction,
where match rate is defined in the description of field 9 from the input above.
Unique displacement count is the number of unique displacements of J from the
5' end of the original read in the set {R_i}. The 5' end of the read is
determined from the strand to which contituent readlets align: if the readlets'
reversed complements align to the forward strand, the 5' end of the original
read sequence is the right end; if the unaltered readlet sequences align to
the forward strand, the 5' end of the original read sequence is the left end.

OUTPUT BED COORDINATES ARE 0-INDEXED; HADOOP OUTPUT COORDINATES ARE 1-INDEXED.
"""
import os
import sys
# Regular expressions are used to identify splice-site motifs
import re
import itertools
from collections import defaultdict
import numpy as np
import site

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
        removed from the lineup, and the next candidate is examined. The
        algorithm also filters out every cluster whose leftmost intron is not
        wholly within [partition_start, partition_end).
        
        candidate_introns: a dictionary. Each key is a tuple
            (start_position, end_position) and its corresponding value is a
            list (of length read_count), each of whose items is a tuple
            (sample_label, five_prime_displacement, three_prime_displacement,
                left_EC_size, right_EC_size, match_rate, reverse_strand)
            associated with a read supporting the candidate. Here,
            five_prime_displacement is the displacement of the 5' end of the
            candidate from the 5' end of the read, while
            three_prime_displacement is the displacement of the 3' end of the
            candidate from the 3' end of the read. left_EC_size (right_EC_size)
            is the number of nucleotides spanned by the EC on the left (right)
            of the intron. Here, "left" ("right") means "towards the 5' (3')
            end of the read." match_rate is the number of bases that match the
            reference per base of the aligned readlets comprising the read from
            which the intron was inferred. reverse_strand is is True iff the
            readlets' reversed complements aligned to the forward strand.
        partition_start: start position (inclusive) of partition.
        partition_end: end position (exclusive) of partition.
        cluster_radius: distance from a candidate intron under examination
            within which to search for candidate introns in the same cluster.
            See above for a detailed description.
        verbose: True iff counts of possible splice junctions, clusters, and
            filtered clusters should be written to stderr.

        Return value: a list of lists, each of which corresponds to a cluster
        of candidate introns. Each item in a cluster is a tuple
        (start_position, end_position, sample_label, five_prime_displacement,
            three_prime_displacement, left_EC_size, right_EC_size,
            match_rate, reverse_strand),
        which corresponds to a read supporting a candidate intron spanning
        [start_position, end_position) in the cluster.
    """
    '''Construct list of candidate introns sorted in order of descending
    read_count.'''
    candidate_intron_list = [(len(sample_labels), 
                                end_position - start_position, end_position)
                                for (start_position, end_position),
                                sample_labels in candidate_introns.items()]
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
    filtered_cluster_count = 0
    for _, intron_size, end_position in candidate_intron_list:
        if (intron_size, end_position) not in clustered_introns:
            current_cluster_radius = min(intron_size, cluster_radius)
            total_cluster_count += 1
            clustered_introns.add((intron_size, end_position))
            intron_cluster = [(intron_size, end_position)]
            for an_end_position in xrange(end_position 
                                            - current_cluster_radius,
                                            end_position
                                            + current_cluster_radius + 1):
                if (intron_size, an_end_position) in candidate_intron_set \
                    and (intron_size,
                            an_end_position) not in clustered_introns:
                    '''If a nearby candidate hasn't been clustered, absorb it
                    into the current cluster.'''
                    intron_cluster.append((intron_size, an_end_position))
                    clustered_introns.add((intron_size, an_end_position))
            if partition_start <= min([an_end_position 
                for _, an_end_position in intron_cluster]) \
                <= partition_end:
                '''Add a cluster iff its leftmost element lies in
                [start_position, end_position)'''
                intron_clusters.append(intron_cluster)
            else:
                filtered_cluster_count += 1
            assert len(clustered_introns) <= candidate_intron_count
            if len(clustered_introns) == candidate_intron_count:
                if verbose:
                    print >> sys.stderr, '%d possible splice junction(s) ' \
                        'clustered down to %d; then %d cluster(s) ' \
                        'filtered out.' % (candidate_intron_count, 
                            total_cluster_count, filtered_cluster_count)
                # Collect reads supporting cluster
                reads_for_intron_clusters = []
                for intron_cluster in intron_clusters:
                    reads_for_intron_cluster = []
                    for an_intron_size, an_end_position in intron_cluster:
                        a_start_position = an_end_position - an_intron_size
                        reads_for_intron_cluster.append((a_start_position, 
                            an_end_position) \
                            + candidate_introns[(a_start_position,
                                                    an_end_position)])
                    reads_for_intron_clusters.append(reads_for_intron_cluster)
                return reads_for_intron_clusters
    raise RuntimeError('For loop should never be completed.')

def ranked_splice_sites_from_cluster(reference_index, intron_cluster,
    rname, motifs, motif_radius=3, verbose=False):
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
        If more than one possible motif (pair) of the same rank is found, the
        tie is broken as follows. The means and standard deviations of the
        start and end positions of the candidate introns are computed. Then the
        z-score of the start motif position is added to the z-score of the end
        end motif position, and the motif pairs are ranked in order of
        ascending z-score. If no motifs are found, the return value of this
        function is an empty list. So above, ranked first would be a GT..AG,
        which is spanned by I3, and second would be a GC..AG, spanned by
        [I1's start position, I3's end position].

        ALL INPUT COORDINATES ARE ASSUMED TO BE 1-INDEXED.

        reference_index: object of class bowtie_index.BowtieIndexReference
                that permits access to reference; used to find splice-site
                motifs.
        intron_cluster: a list of lists, each of which corresponds to a cluster
            of candidate introns. Each item in a cluster is a tuple
            (start_position, end_position, sample_label, 
                five_prime_displacement, three_prime_displacement,
                left_EC_size, right_EC_size, match_rate, reverse_strand), which
            corresponds to a read supporting a candidate intron spanning
            [start_position, end_position) in the cluster.
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
    start_positions, end_positions, _, _, _, _, _, _, _ = zip(*intron_cluster)
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
    mean_start_position = np.mean(start_positions)
    mean_end_position = np.mean(end_positions)
    # Maxes below avoid future ZeroDivisionError exceptions
    stdev_start_position = max(np.std(start_positions), 1e-6)
    stdev_end_position = max(np.std(end_positions), 1e-6)
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
                left_motif_start_position = min_start_position \
                    + left_motif_offset
                right_motif_end_position = min_end_position \
                    + right_motif_offset + 2
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
        ranked_introns += [positions_and_z_score 
                            for positions_and_z_score in 
                            positions_and_z_scores]
    if len(ranked_introns) == 0 and verbose:
        print >>sys.stderr, \
            'Warning: No splice site found for cluster with ' \
            '%d candidate introns' % len(intron_cluster)
    return ranked_introns

def go(bowtie_index_base="genome", input_stream=sys.stdin,
        output_stream=sys.stdout, bin_size=10000, cluster_radius=5,
        per_site=False, per_span=True, output_bed=True, stranded=False,
        intron_partition_overlap=20, verbose=False):
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

        Input (read from input_stream)---
        Tab-delimited columns:
        1. Reference name (RNAME in SAM format) + ';' + partition number +  
            ('+' or '-' indicating which strand is the sense strand if input
                reads are strand-specific -- that is, --stranded in
                Rail-RNA-align was invoked; otherwise, there is no terminal
                '+' or '-')
        2. Sample label
        3. Candidate intron start (inclusive) on forward strand (1-indexed)
        4. Candidate intron end (exclusive) on forward strand (1-indexed)
        5. Number of nucleotides between 5' end of candidate intron and 5' end
        of read from which it was inferred, ASSUMING THE SENSE STRAND IS THE
        FORWARD STRAND. That is, if the sense strand is the reverse strand,
        this is the distance between the 5' end of the reverse-complemented
        read and the 5' end of the reverse-complemented intron.
        6. Number of nucleotides between 3' end of candidate intron and 3' end
        of read from which it was inferred, ASSUMING THE SENSE STRAND IS THE
        FORWARD STRAND.
        7. Number of nucleotides spanned by EC on the left (that is, towards
        the 5' end of the read) of the intron, ASSUMING THE SENSE STRAND IS THE
        FORWARD STRAND.
        8. Number of nucleotides spanned by EC on the right (that is, towards
        the 3' end of the read) of the intron, ASSUMING THE SENSE STRAND IS THE
        FORWARD STRAND.
        9. Match rate: the number of bases that match the reference per base of
        the aligned readlets comprising the read from which the intron was
        inferred.
        10. '-' if reversed complements of readlets from which intron was
        inferred aligned to forward strand; else '+'

        Hadoop output (written to output_stream)---
        Format 1 (span): tab-delimited columns, one line per read spanning a
        splice site:
        1. SAM-format RNAME + ';' + ('+' or '-' indicating which strand is the
            sense strand)
        2. Intron start (inclusive) on forward strand (1-indexed)
        3. Intron end (exclusive) on forward strand (1-index)
        4. Left motif (donor motif if sense strand is '+'; acceptor motif 
            otherwise)
        5. Right motif (acceptor motif if sense strand is '+'; donor motif
            otherwise)
        6. Sample label from which read was derived

        Format 2 (site): tab-delimited columns, one line per splice site
        spanned by the reads in a sample
        1. SAM-format RNAME + ';' + ('+' or '-' indicating which strand is the
            sense strand)
        2. Intron start (inclusive) on forward strand (1-indexed)
        3. Intron end (exclusive) on forward strand (1-index)
        4. Left motif (donor motif if sense strand is '+'; acceptor motif
            otherwise)
        5. Right motif (acceptor motif if sense strand is '+'; donor motif
            otherwise)
        6. Sample label with nonzero number of reads spanning splice site
        7. Number of reads in sample spanning splice site

        Format 3 (junction):
        12-column (3 required fields + 9 optional fields) BED output mimicking
        TopHat's junctions.bed, except anchor significance, maximum match rate,
        and unique displacement count are also included in the name field. From
        the TopHat manual http://tophat.cbcb.umd.edu/manual.shtml: "Each
        junction consists of two connected BED blocks, where each block is as
        long as the maximal overhang of any read spanning the junction. The
        score is the number of alignments spanning the junction." Anchor
        significance of a junction is NOT defined as in the MapSplice paper
        http://nar.oxfordjournals.org/content/38/18/e178.long. Consider the set
        {R_i} of reads that span a junction J. For a given read R_i, consider
        the two anchors A_i and B_i on either side of J, and let L_(A_i) and
        L_(B_i) be their lengths (in bp). The anchor significance of J is
        max_i min(L(A_i), L(B_i)).) Maximum match rate is the highest match
        rate from among the reads spanning a junction, where match rate is
        defined in the description of field 9 from the input above. Unique
        displacement count is the number of unique displacements of J from the
        5' end of the original read in the set {R_i}. The 5' end of the read is
        determined from the strand to which contituent readlets align: if the
        readlets' reversed complements align to the forward strand, the 5' end
        of the original read sequence is the right end; if the unaltered
        readlet sequences align to the forward strand, the 5' end of the
        original read sequence is the left end.

        OUTPUT BED COORDINATES ARE 0-INDEXED; HADOOP OUTPUT COORDINATES ARE
        1-INDEXED.

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
        per_site: False iff Format 2 (site) output lines from above should be
            suppressed.
        per_span: False iff Format 1 (span) output lines from above should be
            suppressed.
        output_bed: False iff Format 3 (junction) output lines from above
            should be suppressed.
        stranded: True iff input reads are strand-specific, and lines read
            from input_stream have first column with terminal '+' or '-'
            indicating which strand is the sense strand.
        intron_partition_overlap: number of bases to subtract from
            reference start position of candidate intron when determining
            genome partition it is in.
        verbose: True iff informative messages should be written to stderr.

        No return value.
"""
    reference_index = bowtie_index.BowtieIndexReference(bowtie_index_base)
    input_line_count = 0
    junction_number = 0
    handle_partition = False
    last_partition_id = None
    candidate_introns = {}
    while True:
        line = input_stream.readline()
        if line:
            input_line_count += 1
            tokens = line.rstrip().split('\t')
            assert len(tokens) == 10
            (partition_id, sample_label, pos, end_pos, five_prime_displacement,
                three_prime_displacement, left_EC_size, right_EC_size, 
                match_rate, reverse_strand) = tokens
            (pos, end_pos, five_prime_displacement,
                three_prime_displacement, left_EC_size, right_EC_size,
                match_rate, reverse_strand) \
            = (int(pos), int(end_pos), int(five_prime_displacement),
                int(three_prime_displacement), int(left_EC_size), 
                int(right_EC_size), float(match_rate), True if 
                reverse_strand == '-' else False)
            assert end_pos > pos, line
            if stranded:
                '''If input reads are strand-specific, the partition_id has a
                terminal '+' or '-' indicating which strand is the sense
                strand. Remove it.'''
                rname, partition_start, partition_end = \
                    partition.parse(partition_id[:-1], bin_size)
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
        if handle_partition:
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
                    for i, intron_cluster in enumerate(intron_clusters):
                        ranked_splice_sites = ranked_splice_sites_from_cluster(
                                reference_index, intron_cluster, last_rname,
                                _reverse_strand_motifs, verbose=verbose
                            )
                        if len(ranked_splice_sites) != 0:
                            # Pick top-ranked intron
                            cluster_splice_sites[i] = ranked_splice_sites[0]
                else:
                    for i, intron_cluster in enumerate(intron_clusters):
                        ranked_splice_sites = ranked_splice_sites_from_cluster(
                                reference_index, intron_cluster, last_rname,
                                _forward_strand_motifs, verbose=verbose
                            )
                        if len(ranked_splice_sites) != 0:
                            # Pick top-ranked intron
                            cluster_splice_sites[i] = ranked_splice_sites[0]
            else:
                # The sense strand is unknown, so use a general motif set
                for i, intron_cluster in enumerate(intron_clusters):
                    ranked_splice_sites = ranked_splice_sites_from_cluster(
                            reference_index, intron_cluster, last_rname,
                            _unstranded_motifs, verbose=verbose
                        )
                    if len(ranked_splice_sites) != 0:
                        # Pick top-ranked intron
                        cluster_splice_sites[i] = ranked_splice_sites[0]
            if per_span:
                for i, (start_position, end_position, z_score_sum, left_motif,
                    right_motif, motif_reverse_strand) \
                    in cluster_splice_sites.items():
                    motif_reverse_strand_string = '-' if motif_reverse_strand \
                        else '+'
                    for _, _, sample_label, _, _, _, _, _, _ \
                        in intron_clusters[i]:
                        print >>output_stream, 'span\t%s\t%012d\t%012d\t' \
                            '%s\t%s\t%s' \
                            % (last_rname + ';' + motif_reverse_strand_string,
                                start_position, end_position, 
                                left_motif, right_motif, sample_label)
            if per_site:
                for i, (start_position, end_position, z_score_sum, left_motif,
                    right_motif, motif_reverse_strand) \
                    in cluster_splice_sites.items():
                    motif_reverse_strand_string = '-' if motif_reverse_strand \
                        else '+'
                    sample_label_counts = defaultdict(int)
                    for _, _, sample_label, _, _, _, _, _, _\
                        in intron_clusters[i]:
                        sample_label_counts[sample_label] += 1
                    for sample_label in sample_label_counts:
                        print >>output_stream, \
                            'site\t%s\t%012d\t%012d\t%s\t%s\t%s\t%s' \
                            % (last_rname + ';' + motif_reverse_strand_string,
                                start_position, end_position, left_motif,
                                right_motif, sample_label, 
                                sample_label_counts[sample_label])
            if output_bed:
                '''The output bed mimics TopHat's junctions.bed, except
                anchor significance, maximum match rate, and unique
                displacement count are in the name field. Anchor significance
                of a junction is NOT defined as in the MapSplice paper
                http://nar.oxfordjournals.org/content/38/18/e178.long.
                Consider the set {R_i} of reads spanning a junction J.
                For a given read R_i, consider the two anchors A_i and B_i on
                either side of J, and let L_(A_i) and L_(B_i) be their lengths
                (in bp). The anchor significance of J is
                max_i min(L(A_i), L(B_i)).) Maximum match rate is the highest
                match rate from among the reads spanning a junction. Unique
                displacement count is the number of unique displacements of J
                from the 5' end of the original read in the set {R_i}. The 5'
                end of the read is determined from the strand to which
                contituent readlets align: if the readlets' reversed
                complements align to the forward strand, the 5' end of the
                original read sequence is the right end; if the unaltered
                readlet sequences align to the forward strand, the 5' end of
                the original read sequence is the left end. See TopHat
                documentation for more information on junctions.bed.'''
                for i, (start_position, end_position, z_score_sum, left_motif,
                    right_motif, motif_reverse_strand) \
                    in cluster_splice_sites.items():
                    motif_reverse_strand_string = '-' if motif_reverse_strand \
                        else '+'
                    '''Identify longest read overhangs on either side of splice
                    junction.'''
                    left_overhang, right_overhang = 0, 0
                    read_anchor_significance = None
                    # For counting number of unique displacements
                    displacement_set = set()
                    maximum_match_rate = 0
                    for (candidate_start_position, candidate_end_position,
                            _, candidate_five_prime_displacement,
                            candidate_three_prime_displacement,
                            candidate_left_EC_size,
                            candidate_right_EC_size, candidate_match_rate,
                            candidate_reverse_strand) in intron_clusters[i]:
                        if candidate_reverse_strand:
                            displacement_set.add(
                                    candidate_three_prime_displacement
                                )
                        else:
                            displacement_set.add(
                                    candidate_five_prime_displacement
                                )
                        if candidate_match_rate > maximum_match_rate:
                            maximum_match_rate = candidate_match_rate
                        read_anchor_significance = max(min(
                                    candidate_left_EC_size,
                                    candidate_right_EC_size
                                ), read_anchor_significance)
                        left_overhang = max(start_position
                                        - candidate_start_position
                                        + candidate_five_prime_displacement, 
                                        left_overhang)
                        right_overhang = max(candidate_end_position
                            - end_position
                            + candidate_three_prime_displacement,
                            right_overhang)
                    junction_number += 1
                    '''Print line of bed file; where a first column 'junction'
                    is inserted so the line can be distinguished from other
                    types of output lines. Subtract 1 from every position to
                    accommodate how bed is 0-based.'''
                    left_pos = start_position - left_overhang - 1
                    right_pos = end_position + right_overhang - 1
                    print >>output_stream, \
                        'junction\t%s\t%012d\t%012d\t' \
                        'read_anchor_significance=%d;' \
                        'maximum_match_rate=%.12f;' \
                        'unique_displacement_count=%d\t%d\t%s\t%d\t' \
                        '%d\t255,0,0\t2\t%d,%d\t0,%d' \
                        % (last_rname, left_pos, right_pos,
                            read_anchor_significance,
                            maximum_match_rate,
                            len(displacement_set),
                            len(intron_clusters[i]),
                            motif_reverse_strand_string, left_pos, right_pos,
                            left_overhang, right_overhang, 
                            end_position - left_pos - 1)
            candidate_introns = {}
            handle_partition = False
        if line:
            candidate_introns[(pos, end_pos)] = (sample_label,
                five_prime_displacement, three_prime_displacement, 
                left_EC_size, right_EC_size, match_rate, reverse_strand)
            (last_partition_id, last_partition_start, last_partition_end, 
                last_rname, last_reverse_strand) \
            = (partition_id, partition_start, partition_end, rname, 
                reverse_strand)
        else: break

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
    parser.add_argument('--intron-partition-overlap', type=int, required=False,
        default=20, 
        help='Amount by which partitions overlap their left and right '
             'neighbors')
    parser.add_argument(\
        '--per-site', action='store_const', const=True, default=False,
        help='Output one record for every splice site, giving information '
             'about the number of times a read from each label spans the site')
    parser.add_argument(\
        '--per-span', action='store_const', const=True, default=False,
        help='Output one record for every instance where a read '
             'spans a splice site')
    parser.add_argument(\
        '--output-bed', action='store_const', const=True, default=False,
        help='Output BED lines denoting splice junctions analogous to '
             'TopHat\'s junctions.bed')
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
        per_site=args.per_site,
        per_span=args.per_span,
        output_bed=args.output_bed,
        intron_partition_overlap=args.intron_partition_overlap,
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
                    (5, 25)  : ('sample', 25, 75, 25, 75, .9, True),
                    (4, 24)  : ('sample', 25, 75, 25, 75, .95, True),
                    (10, 30) : ('sample2', 23, 77, 23, 77, .92, False)
                }
            clusters = intron_clusters_in_partition(candidate_introns,
                1, 40, cluster_radius=10, verbose=False)
            self.assertEquals(len(clusters), 1)
            self.assertEquals(
                    sorted(clusters[0]), 
                    [(4, 24, 'sample', 25, 75, 25, 75, .95, True),
                     (5, 25, 'sample', 25, 75, 25, 75, .9, True),
                     (10, 30, 'sample2', 23, 77, 23, 77, .92, False)]
                )

        def test_that_overlapping_introns_of_different_size_separate(self):
            """ Fails if introns do not separate properly into two clusters.
            """
            candidate_introns = {
                    (5, 25)  : ('sample', 25, 75, 25, 75, .9, True),
                    (4, 24)  : ('sample', 25, 75, 25, 75, .95, True),
                    (10, 40) : ('sample2', 23, 77, 23, 77, .92, False)
                }
            '''The third intron above has a different size than the first two,
            so it should be placed in its own cluster.'''
            clusters = sorted(intron_clusters_in_partition(candidate_introns,
                1, 40, cluster_radius=10, verbose=False), key=len)
            self.assertEquals(len(clusters), 2)
            self.assertEquals(
                    sorted(clusters[0]),
                    [(10, 40, 'sample2', 23, 77, 23, 77, .92, False)]
                )
            self.assertEquals(
                    sorted(clusters[1]),
                    [(4, 24, 'sample', 25, 75, 25, 75, .95, True), 
                     (5, 25, 'sample', 25, 75, 25, 75, .9, True)]
                )

        def test_that_cluster_is_filtered(self):
            """ Fails if cluster is not filtered out. 

                A cluster is filtered out if its leftmost member does not lie
                entirely in the partition. It would hopefully be caught
                by the next partition given the partition overlap.
            """
            candidate_introns = {
                    (5, 25)  : ('sample', 25, 75, 25, 75, .9, True),
                    (4, 24)  : ('sample', 25, 75, 25, 75, .95, True),
                    (10, 40) : ('sample2', 23, 77, 23, 77, .92, False)
                }
            '''Note that the partition end is set to 30 below, so the last
            candidate intron above should be filtered out.'''
            clusters = intron_clusters_in_partition(candidate_introns,
                1, 30, cluster_radius=10, verbose=False)
            self.assertEquals(len(clusters), 1)
            self.assertEquals(
                    sorted(clusters[0]),
                    [(4, 24, 'sample', 25, 75, 25, 75, .95, True), 
                     (5, 25, 'sample', 25, 75, 25, 75, .9, True)]
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
            cluster = [(45, 97, 'sample', 25, 75, 25, 75, .99, True),
                       (47, 99, 'sample', 25, 75, 25, 75, .99, True),
                       (44, 96, 'sample', 25, 75, 25, 75, .87, True)]
            splice_sites = ranked_splice_sites_from_cluster(
                    self.reference_index, cluster,
                    'chr1', _unstranded_motifs, motif_radius=3, verbose=False
                )
            self.assertEquals((48, 95), splice_sites[0][:2])
            self.assertEquals(('GT', 'AG'), splice_sites[0][3:5])

        def tearDown(self):
            # Kill temporary directory
            shutil.rmtree(self.temp_dir_path)

    unittest.main()