#!/usr/bin/env python
"""
alignment_handlers.py
Part of Rail-RNA

Includes various functions for handling alignments output by Bowtie2 -- most
importantly:
-a function that outputs indels, introns, and exons from a genome
position, CIGAR string, and MD string (indels_introns_and_exons)
-a function that inserts introns in a CIGAR string (multread_with_introns).
"""

import re
import random
import sys
import bisect
import partition

def running_sum(iterable):
    """ Generates a running sum of the numbers in an iterable

        iterable: some iterable with numbers

        Yield value: next value in running sum
    """
    total = 0
    for number in iterable:
        total += number
        yield total

def multiread_with_introns(multiread, stranded=False):
    """ Modifies read alignments to fix CIGARs/positions/primary alignments.

        An alignment takes the form of a line of SAM.
        Introns are encoded in a multiread's RNAME as follows:

        original RNAME + '+' or '-' indicating which strand is the sense
        strand + '\x1d' + start position of sequence + '\x1d' + comma-separated
        list of subsequence sizes framing introns + '\x1d' + comma-separated
        list of intron sizes.

        If there are no introns present, the multiread's RNAME takes the
        following form:
        original RNAME + '\x1d' + start position of sequence + '\x1d\x1d'

        More than one alignment output by Bowtie may correspond to the same
        alignment in reference space because at least two reference names in
        the intron Bowtie index contained overlapping bases. In this case, the
        alignments are collapsed into a single alignment. If an alignment is
        found to overlap introns, the XS:A:('+' or '-') field is appended to
        indicate which strand is the sense strand. An NH:i:(integer) field is
        also added to each alignment to indicate the number of alignments.
        These extra fields are used by Cufflinks.

        multiread: a list of lists, each of whose elements are the tokens
            from a line of SAM representing an alignment.
        stranded: if input reads are stranded, an alignment is returned only
            if the strand of any introns agrees with the strand of the 
            alignment.

        Return value: alignments modified according to the rules given above;
            a list of tuples, each of whose elements are the tokens from a line
            of SAM representing an alignment. NOTE: ALL ALIGNMENTS ARE
            MARKED AS SECONDARY, and multiread_to_report must decide which
            is primary!
    """
    if not multiread:
        return []
    seq = multiread[0][9]
    qual = multiread[0][10]
    new_multiread = []
    for alignment in multiread:
        qname = alignment[0]
        tokens = alignment[2].split('\x1d')
        offset = int(alignment[3]) - 1
        cigar = re.split(r'([MINDS])', alignment[5])[:-1]
        flag = int(alignment[1])
        if not tokens[-1] or len(tokens) == 1:
            # No introns can be found
            try:
                pos = offset + int(tokens[1])
            except IndexError:
                # Ordinary alignment without augmented RNAME
                pos = offset + 1
            rname = tokens[0]
            '''Use second field in each element of new_multiread to store which
            items should be tested to find whether two alignments are
            "identical".'''
            new_multiread.append(
                    ([qname, str(flag | 256), rname,
                        str(pos), alignment[4],
                        alignment[5]] + list(alignment[6:]),
                    (qname, (flag & 16 != 0),
                        rname,
                        pos, alignment[5],
                        [field for field in alignment
                         if field[:5] == 'MD:Z:'][0][5:])
                        )
                )
            continue
        reverse_strand_string = tokens[0][-1]
        assert reverse_strand_string in '+-'
        reverse_strand = (True if reverse_strand_string == '-' else False)
        if stranded and (flag & 16 != 0) == reverse_strand:
            # Strand of alignment doesn't agree with strand of intron
            continue
        rname = tokens[0][:-1]
        exon_sizes = map(int, tokens[2].split(','))
        intron_sizes = map(int, tokens[3].split(','))
        for i, exon_sum in enumerate(running_sum(exon_sizes)):
            if exon_sum > offset: break
        # Compute start position of alignment
        pos = offset + sum(intron_sizes[:i]) + int(tokens[1])
        # Adjust exon/intron lists so they start where alignment starts
        exon_sizes = exon_sizes[i:]
        exon_sizes[0] = exon_sum - offset
        intron_sizes = intron_sizes[i:]
        new_cigar = []
        for i in xrange(0, len(cigar), 2):
            char_type = cigar[i+1]
            base_count = int(cigar[i])
            if char_type in 'MD':
                for j, exon_sum in enumerate(running_sum(exon_sizes)):
                    if exon_sum >= base_count: break
                for k in xrange(j):
                    new_cigar.extend(
                            [(str(exon_sizes[k]) + char_type)
                                if exon_sizes[k] != 0 else '',
                             str(intron_sizes[k]), 'N']
                        )
                last_size = base_count - (exon_sum - exon_sizes[j])
                new_cigar.extend([str(last_size), char_type])
                new_size = exon_sum - base_count
                exon_sizes = exon_sizes[j:]
                exon_sizes[0] = new_size
                intron_sizes = intron_sizes[j:]
            elif char_type in 'IS':
                new_cigar.extend([cigar[i], char_type])
            else:
                raise RuntimeError('Bowtie2 CIGAR chars are expected to be '
                                   'in set (DIMS).')
        new_cigar = ''.join(new_cigar)
        if 'N' not in new_cigar:
            '''Alignment to transcriptome was purely exonic; this case should
            be ignored.'''
            continue
        new_multiread.append(
                    ([alignment[0], str(flag | 256),
                        rname, str(pos), alignment[4], new_cigar]
                     + list(alignment[6:])
                     + ['XS:A:' + reverse_strand_string],
                     (alignment[0],
                        (flag & 16 != 0),
                        rname,
                        pos, new_cigar,
                        [field for field in alignment[::-1]
                         if field[:5] == 'MD:Z:'][0][5:]))
                )
    if not new_multiread:
        return []
    new_multiread.sort(key=lambda alignment: alignment[1])
    # Eliminate duplicate alignments and set primary alignment
    multiread_to_return = [new_multiread[0][0]]
    for i in xrange(1, len(new_multiread)):
        if new_multiread[i][1] == new_multiread[i-1][1]:
            continue
        multiread_to_return.append(new_multiread[i][0])
    if len(multiread_to_return) == 1:
        '''Only one alignment; remove XS:i field if it's present.'''
        for i in xrange(10, len(multiread_to_return[0])):
            if multiread_to_return[0][i][:5] == 'XS:i:':
                multiread_to_return[0].remove(multiread_to_return[0][i])
                break
        return [tuple(alignment) for alignment in multiread_to_return]
    # Correct XS:i fields
    alignment_scores = [[int(token[5:]) for token in alignment[::-1]
                            if token[:5] == 'AS:i:'][0]
                            for alignment in multiread_to_return]
    sorted_alignment_scores = sorted(alignment_scores, reverse=True)
    XS_field = 'XS:i:%d' % sorted_alignment_scores[1]
    for i in xrange(len(multiread_to_return)):
        for j in xrange(10, len(multiread_to_return[i])):
            if multiread_to_return[i][j][:5] == 'XS:i:':
                multiread_to_return[i][j] = XS_field
    return [tuple(alignment) for alignment in multiread_to_return]

def multiread_to_report(multiread, alignment_count_to_report=1, seed=0,
    non_deterministic=False, weights=[], tie_margin=6):
    """ Parses Bowtie2 args and returns which alignments should be reported.

        Primary alignment is decided here. flag & 256 should == 1 for
        every incoming alignment.

        multiread: list of tuples, each of whose elements are the tokens from
            a line of SAM representing an alignment
        alignment_count_to_report: number of alignments to report or -1 if
            all alignments are to be reported
        seed: Bowtie2's --seed parameter
        non_deterministic: Bowtie2's --non-deterministic parameter
        weights: integer weights to use when breaking ties in alignments
        tie_margin: allowed score difference per 100 bases among ties in
            max score. For example, 150 and 144 are tied alignment scores
            for a 100-bp read when tie_margin is 6.

        Return value: either:
            1) 2-tuple whose second element is a list of "tied" alignments
                and whose first element is a list of "resolved" alignments
                that are ready to be output as secondary SAM lines; tied
                alignments will be resolved by determining primary in a
                subsequent step. No alignment is primary yet.
            2) Tuple whose sole element is a list of resolved alignments,
                where one alignment is a primary
    """
    if not multiread:
        return ([],)
    # Set random seed
    if non_deterministic:
        # Use system time to set seed, as in Bowtie2
        random.seed()
    else:
        '''Seed is sum of qname, seq, qual, and seed, approximating
        Bowtie2's behavior.'''
        random.seed(multiread[0][0] + multiread[0][9] + multiread[0][10]
                     + str(seed))
    if weights:
        # Choose primary alignment using weights
        assert len(multiread) == len(weights)
        total = sum(weights)
        weight_bounds = [0.] + list(running_sum([float(weight) / total
                                            for weight in weights]))
        primary_index = bisect.bisect_left(weight_bounds, random.random()) - 1
        prereturn_multiread = (
                [(multiread[primary_index][0],
                    str(int(multiread[primary_index][1]) & ~256))
                    + multiread[primary_index][2:]] +
                sorted([(alignment[0], str(int(alignment[1]) | 256))
                    + alignment[2:] for i, alignment in enumerate(multiread)
                    if i != primary_index],
                    key=lambda alignment: [int(token[5:]) for token
                                            in alignment[::-1]
                                            if token[:5] == 'AS:i:'][0],
                    reverse=True)
            ,)
    else:
        seq_size = len(multiread[0][9])
        tie_margin = round(tie_margin * float(seq_size) / 100)
        # Determine number of "ties" including margin
        alignment_count = len(multiread)
        alignments_and_scores = [(alignment,
                                    [int(token[5:]) for token
                                     in alignment[::-1]
                                     if token[:5] == 'AS:i:'][0])
                                     for alignment in multiread]
        alignments_and_scores.sort(key=lambda alignment_and_score:
                                        alignment_and_score[1],
                                        reverse=True)
        min_score = alignments_and_scores[0][1] - tie_margin
        for tie_count in xrange(alignment_count):
            if alignments_and_scores[tie_count][1] < min_score: break
        if alignments_and_scores[-1][1] >= min_score: tie_count += 1
        if tie_count == 1:
            # Primary alignment is resolved
            prereturn_multiread = (
                [(alignments_and_scores[0][0][0],
                    str(int(alignments_and_scores[0][0][1]) & ~256))
                    + alignments_and_scores[0][0][2:]] +
                [(alignment[0], str(int(alignment[1]) | 256))
                    + alignment[2:] for alignment, _
                    in alignments_and_scores[1:]]
            ,)
        else:
            # Primary alignment must be resolved in a later step
            if tie_count >= alignment_count_to_report:
                '''Nothing to report yet; all alignments to be reported are
                among ties.'''
                return ([], [(alignment[0], str(int(alignment[1]) | 256))
                                + alignment[2:] + ('NH:i:%d'
                                    % alignment_count_to_report,)
                                for alignment, _
                                in alignments_and_scores[:tie_count]])
            prereturn_multiread = (
                [(alignment[0], str(int(alignment[1]) | 256))
                            + alignment[2:] for alignment, _
                            in alignments_and_scores[tie_count:]],
                [(alignment[0], str(int(alignment[1]) | 256))
                            + alignment[2:] for alignment, _
                            in alignments_and_scores[:tie_count]]
            )
    prereturn_multiread_count = len(prereturn_multiread)
    if alignment_count_to_report != -1:
        if prereturn_multiread_count == 1:
            '''All alignments have been resolved; retain as many as there are
            to be reported.'''
            reports_to_retain_count = alignment_count_to_report
        else:
            assert prereturn_multiread_count == 2
            reports_to_retain_count = alignment_count_to_report - len(
                    prereturn_multiread[1]
                )
        assert reports_to_retain_count > 0
        scores = [[int(token[5:]) for token in alignment[::-1]
                   if token[:5] == 'AS:i:'][0] 
                   for alignment in prereturn_multiread[0]]
        try:
            min_permitted_score = scores[reports_to_retain_count - 1]
        except IndexError:
            # k value larger than number of alignments; retain all
            reports_to_return = prereturn_multiread[0]
        else:
            left_count = scores[:reports_to_retain_count - 1].count(
                                min_permitted_score
                            )
            min_permitted_count = scores.count(min_permitted_score)
            if min_permitted_score == scores[0] \
                and not (int(prereturn_multiread[0][0][1]) & 256):
                '''Primary alignment has been decided, but its score is the
                same as the min score. Exclude it from random sample of 
                min scores.'''
                left_count -= 1 
                min_permitted_count -= 1
            reports_to_return = prereturn_multiread[0][
                                    :reports_to_retain_count - left_count - 1
                                ] + \
                random.sample(prereturn_multiread[0][
                    reports_to_retain_count - left_count - 1:
                    reports_to_retain_count - left_count - 1
                    + min_permitted_count
                ], left_count + 1)
    else:
        # Report all
        reports_to_return = prereturn_multiread[0]
    '''Change MAPQ to 255 to make results consistent; TODO: CHECK unique.h in
    BT2 code to compute MAPQ consistent with that for genome alignments.'''
    reports_to_return = [alignment[:4] + ('255',) + alignment[5:] 
                            for alignment in reports_to_return]
    if prereturn_multiread_count == 1:
        NH_field = 'NH:i:' + str(len(reports_to_return))
        return ([(alignment + (NH_field,) if 'NH:i:' not in
                  [alignment[-1][:5], alignment[-2][:5]]
                  else alignment) for alignment in reports_to_return],)
    # prereturn_multiread's length is 2, meaning there are ties
    NH_field = 'NH:i:'+ str(len(reports_to_return)
                            + len(prereturn_multiread[1]))
    return ([alignment + (NH_field,) for alignment in reports_to_return],
            [alignment + (NH_field,) for alignment in prereturn_multiread[1]])

def parsed_md(md):
    """ Divides an MD string up by boundaries between ^, letters, and numbers

        md: an MD string (example: 33A^CC).

        Return value: MD string split by boundaries described above.
    """
    md_to_parse = []
    md_group = [md[0]]
    for i, char in enumerate(md):
        if i == 0: continue
        if (re.match('[A-Za-z]', char) is not None) \
            != (re.match('[A-Za-z]', md[i-1]) is not None) or \
            (re.match('[0-9]', char) is not None) \
            != (re.match('[0-9]', md[i-1]) is not None):
            if md_group:
                md_to_parse.append(''.join(md_group))
            md_group = [char]
        else:
            md_group.append(char)
    if md_group:
        md_to_parse.append(''.join(md_group))
    return [char for char in md_to_parse if char != '0']

def reference_from_seq(cigar, seq, reference_index, rname, pos):
    """ Gets appropriate stretch of reference sequence.

        A BowtieIndexReference object is needed to extract the correct
        soft-clipped bases from the reference. Cigar is used to account
        for indels.

        cigar: CIGAR string; used to extract initial soft clip
        seq: read sequence
        reference_index: object of class BowtieIndexReference
        rname: RNAME
        pos: 1-based POS of alignment

        Return value: tuple start pos, reference sequence
    """
    cigar = re.split(r'([MINDS])', cigar)[:-1]
    del_count = sum([int(cigar[i-1]) for i in xrange(len(cigar))
                        if cigar[i] == 'D'])
    insert_count = sum([int(cigar[i-1]) for i in xrange(len(cigar))
                        if cigar[i] == 'I'])
    if cigar[1] == 'S':
        preclip = int(cigar[0])
    else:
        preclip = 0
    base_count = len(seq) - insert_count + del_count
    new_pos = pos - preclip
    if new_pos < 1:
        new_pos = 1
    reference_seq = reference_index.get_stretch(rname, new_pos - 1,
                                                base_count)
    # Clip Ns from ends of reference sequence; workaround for BT2 bug
    i = 0
    while reference_seq[i] == 'N':
        i += 1
    j = -1
    while reference_seq[j] == 'N':
        j -= 1
    if j == -1:
        return (new_pos + i, reference_seq[i:])
    return (new_pos + i, reference_seq[i:j+1])

def indels_introns_and_exons(cigar, md, pos, seq):
    """ Computes indels, introns, and exons from CIGAR, MD string,
        and POS of a given alignment.

        cigar: CIGAR string
        md: MD:Z string
        pos: position of first aligned base
        seq: read sequence

        Return value: tuple (insertions, deletions, introns, exons). Insertions
            is a list of tuples (last genomic position before insertion, 
                                 string of inserted bases). Deletions
            is a list of tuples (first genomic position of deletion,
                                 string of deleted bases). Introns is a list
            of tuples (intron start position (inclusive),
                       intron end position (exclusive),
                       left_diplacement, right_displacement). Exons is a list
            of tuples (exon start position (inclusive),
                       exon end position (exclusive)).
    """
    insertions, deletions, introns, exons = [], [], [], []
    cigar = re.split(r'([MINDS])', cigar)[:-1]
    md = parsed_md(md)
    seq_size = len(seq)
    cigar_chars, cigar_sizes = [], []
    cigar_index, md_index, seq_index = 0, 0, 0
    max_cigar_index = len(cigar)
    while cigar_index != max_cigar_index:
        if cigar[cigar_index] == 0:
            cigar_index += 2
            continue
        if cigar[cigar_index+1] == 'M':
            aligned_base_cap = int(cigar[cigar_index])
            aligned_bases = 0
            while True:
                try:
                    aligned_bases += int(md[md_index])
                    if aligned_bases <= aligned_base_cap:
                        md_index += 1
                except ValueError:
                    # Not an int, but should not have reached a deletion
                    assert md[md_index] != '^', '\n'.join(
                                                ['cigar and md:',
                                                 ''.join(cigar), ''.join(md)]
                                            )
                    if aligned_bases + len(md[md_index]) > aligned_base_cap:
                        md[md_index] = md[md_index][
                                            :aligned_base_cap-aligned_bases
                                        ]
                        aligned_bases = aligned_base_cap
                    else:
                        aligned_bases += len(md[md_index])
                        md_index += 1
                if aligned_bases > aligned_base_cap:
                    md[md_index] = aligned_bases - aligned_base_cap
                    break
                elif aligned_bases == aligned_base_cap:
                    break
            # Add exon
            exons.append((pos, pos + aligned_base_cap))
            pos += aligned_base_cap
            seq_index += aligned_base_cap
        elif cigar[cigar_index+1] == 'N':
            skip_increment = int(cigar[cigar_index])
            # Add intron
            introns.append((pos, pos + skip_increment,
                            seq_index, seq_size - seq_index))
            # Skip region of reference
            pos += skip_increment
        elif cigar[cigar_index+1] == 'I':
            # Insertion
            insert_size = int(cigar[cigar_index])
            insertions.append(
                    (pos - 1, seq[seq_index:seq_index+insert_size])
                )
            seq_index += insert_size
        elif cigar[cigar_index+1] == 'D':
            assert md[md_index] == '^', '\n'.join(
                                                ['cigar and md:',
                                                 ''.join(cigar), ''.join(md)]
                                            )
            # Deletion
            delete_size = int(cigar[cigar_index])
            md_delete_size = len(md[md_index+1])
            assert md_delete_size >= delete_size
            deletions.append((pos, md[md_index+1][:delete_size]))
            if md_delete_size > delete_size:
                # Deletion contains an intron
                md[md_index+1] = md[md_index+1][delete_size:]
            else:
                md_index += 2
            # Skip deleted part of reference
            pos += delete_size
        else:
            # Soft clip
            assert cigar[cigar_index+1] == 'S'
            # Advance seq_index
            seq_index += int(cigar[cigar_index])
        cigar_index += 2
    '''Merge exonic chunks/deletions; insertions/introns could have chopped
    them up.'''
    new_exons = []
    last_exon = exons[0]
    for exon in exons[1:]:
        if exon[0] == last_exon[1]:
            # Merge ECs
            last_exon = (last_exon[0], exon[1])
        else:
            # Push last exon to new exon list
            new_exons.append(last_exon)
            last_exon = exon
    new_exons.append(last_exon)
    return insertions, deletions, introns, new_exons

class AlignmentPrinter(object):
    """ Encapsulates methods for printing alignment information. """

    def __init__(self, manifest_object, reference_index,
                 output_stream=sys.stdout, bin_size=5000, exon_ivals=False,
                 exon_diffs=True):
        """
            manifest_object: object of type LabelsAndIndices; see manifest.py
            reference_index: object of type BowtieIndexReference; see bowtie.py
            output_stream: where to print output
            exon_ivals: True iff exon_ivals should be output
            exon_diffs: True iff exon_diffs should be output
        """
        self.manifest_object = manifest_object
        self.reference_index = reference_index
        self.exon_ivals = exon_ivals
        self.exon_diffs = exon_diffs
        self.bin_size = bin_size
        self.output_stream = output_stream

    def print_unmapped_read(self, qname, seq, qual):
        """ Prints an unmapped read from a qname, qual, and seq.

            First field is "sam", as for reported alignments from the function
            print_alignment_data().

             NOTE THAT UNPAIRED READ FLAG YT:Z:UU IS HARD-CODED; MUST UPDATE IF
             PAIRED-END ALIGNMENT IS EVER PERFORMED.

            qname: full qname in the format <original_qname> + '\x1d' + 
                <short hash of original_qname + seq + sample label> + '\x1d' +
                <sample_label>
            seq: seq to print 
            qual: qual to print

            Return value: 1, the number of lines output
        """
        print >>self.output_stream, (
                'sam\t%s\t%s\t%012d\t%s\t4\t0\t*\t*\t0\t0\t%s\t%s\tYT:Z:UU'
            ) % (self.manifest_object.label_to_index[
                        qname.rpartition('\x1d')[2]
                    ], self.reference_index.rname_to_string['*'],
                    0, qname.partition('\x1d')[0], seq, qual)
        return 1

    def print_alignment_data(self, multiread_reports_and_ties):
        """ Prints almost-SAM alignments, introns/indels, and exonic coverage.

            Descriptions of output:

            Alignments
            
            (sam_intron_ties) output only for ties in alignment score (which
                are within some tie margin as decided by multiread_to_report);
            this is the first element of multiread_reports_and_ties
            tab-delimited output tuple columns:
            Standard SAM output except fields are in different order -- and the
            first four fields include sample/intron information. If an
            alignment overlaps k introns, k lines are output. The order of the
            fields is as follows.
            1. The character 'N' so the line can be matched up
                with intron bed lines
            2. Sample index
            3. Number string representing RNAME; see BowtieIndexReference class
                in bowtie_index for conversion information
            4. Intron start position
            5. Intron end position
            6. '+' or '-' indicating which strand is sense strand
            7. '-' to ensure that the line follows all intron lines
            8. POS
            9. QNAME
            10. FLAG
            11. MAPQ
            12. CIGAR
            13. RNEXT
            14. PNEXT
            15. TLEN
            16. SEQ
            17. QUAL
            ... + optional fields

            (sam_clip_ties) output only for ties in alignment score when
            no introns are overlapped -- these alignments are almost invariably
            soft-clipped
            [SAME AS SAM FIELDS; see SAM format specification]

            (sam) output only for alignments to be reported (first element
            of tuple multiread_reports_and_ties)
            score +/- tie_margin); tab-delimited output tuple columns:
            Standard SAM output except fields are in different order to
            faciliate partitioning by sample/RNAME and coordinate sorting. The
            order of the fields is as follows.
            1. Sample index
            2. Number string representing RNAME; see BowtieIndexReference class
                in bowtie_index for conversion information
            3. POS
            4. QNAME
            5. FLAG
            6. MAPQ
            7. CIGAR
            8. RNEXT
            9. PNEXT
            10. TLEN
            11. SEQ
            12. QUAL
            ... + optional fields

            Exonic chunks (aka ECs; two formats, any or both of which may be
                emitted -- only if primary alignment is present among
                alignments to be reported):

            Exonic chunks in interval format (exon_ival);
            tab-delimited output tuple columns:
            1. Reference name (RNAME in SAM format) + ';' + bin number
            2. Sample label
            3. EC start (inclusive) on forward strand
            4. EC end (exclusive) on forward strand

            Exonic chunks in diff format (exon_diff) -- only if primary
                alignment is present among alignments to be reported;
            tab-delimited output tuple columns:
            1. Reference name (RNAME in SAM format) + ';' + 
                max(EC start, bin start) (inclusive) on forward strand IFF diff
                is positive and EC end (exclusive) on forward strand IFF diff
                is negative
            2. Bin number
            3. Sample label
            4. +1 or -1.

            Introns (intron_bed) / insertions/deletions (indel_bed);
            tab-delimited output tuple columns:
            1. 'I', 'D', or 'N' for insertion, deletion, or intron line
            2. Sample index
            3. Number string representing RNAME
            4. Start position (Last base before insertion, first base of
                                deletion, or first base of intron)
            5. End position (Last base before insertion, last base of deletion
                (exclusive), or last base of intron (exclusive))
            6. '+' or '-' indicating which strand is the sense strand for
                introns, inserted sequence for insertions, or deleted sequence
                for deletions
            ----Next fields are for introns only; they are '\x1c' for indels---
            7. Number of nucleotides between 5' end of intron and 5' end of
                read from which it was inferred, ASSUMING THE SENSE STRAND IS
                THE FORWARD STRAND. That is, if the sense strand is the reverse
                strand, this is the distance between the 3' end of the read and
                the 3' end of the intron.
            8. Number of nucleotides between 3' end of intron and 3' end of
                read from which it was inferred, ASSUMING THE SENSE STRAND IS
                THE FORWARD STRAND.
            -------------------------------------------------------------------
            9. Number of instances of intron, insertion, or deletion in sample;
                this is always +1 before bed_pre combiner/reducer

            multiread_reports_and_ties: either:
                1) 2-tuple whose second element is a list of "tied" alignments
                    and whose first element is a list of "resolved" alignments
                    that are ready to be output as secondary SAM lines; tied
                    alignments will be resolved by determining primary in a
                    subsequent step. No alignment is primary yet.
                2) Tuple whose sole element is a list of resolved alignments,
                    where one alignment is a primary.
                Every QNAME takes the form
                <original_qname> + '\x1d' + 
                <short hash of original_qname + seq + sample label> + '\x1d' +
                <sample_label>
            manifest_object: object of type LabelsAndIndices; see manifest.py
            reference_index: object of type BowtieIndexReference; see bowtie.py
            output_stream: where to print output
            exon_ivals: True iff exon_ivals should be output
            exon_diffs: True iff exon_diffs should be output

            Return value: output line count
        """
        output_line_count = 0
        try:
            primary_flag = int(multiread_reports_and_ties[0][0][1])
        except IndexError:
            # No alignments to report
            pass
        else:
            sample_index = self.manifest_object.label_to_index[
                    multiread_reports_and_ties[0][0][0].rpartition('\x1d')[2]
                ]
            if not (primary_flag & 256):
                '''First alignment to report is a primary, so output exons,
                introns, and indels.'''
                alignment = multiread_reports_and_ties[0][0]
                cigar = alignment[5]
                rname = alignment[2]
                pos = int(alignment[3])
                seq = alignment[9]
                md = [field for field in alignment
                            if field[:5] == 'MD:Z:'][0][5:]
                insertions, deletions, introns, exons \
                    = indels_introns_and_exons(cigar, md, pos, seq)
                # Output indels
                for insert_pos, insert_seq in insertions:
                    print >>self.output_stream, (
                           ('indel_bed\tI\t%s\t%s\t%012d\t%012d\t%s'
                            '\t\x1c\t\x1c\t1')
                            % (sample_index, self.reference_index.\
                                rname_to_string[rname],
                                insert_pos, insert_pos,
                                insert_seq)
                        )
                    output_line_count += 1
                for del_pos, del_seq in deletions:
                    print >>self.output_stream, (
                           ('indel_bed\tD\t%s\t%s\t%012d\t%012d\t%s'
                            '\t\x1c\t\x1c\t1')
                            % (sample_index, self.reference_index.\
                                rname_to_string[rname],
                                del_pos, del_pos + len(del_seq),
                                del_seq)
                        )
                    output_line_count += 1
                # Output exonic chunks
                if self.exon_ivals:
                    for exon_pos, exon_end_pos in exons:
                        partitions = partition.partition(
                                rname, exon_pos, exon_end_pos, self.bin_size
                            )
                        for partition_id, _, _ in partitions:
                            print >>self.output_stream, \
                                'exon_ival\t%s\t%012d\t' \
                                '%012d\t%s' \
                                % (partition_id,
                                    exon_pos, exon_end_pos, 
                                    sample_index)
                            output_line_count += 1
                if self.exon_diffs:
                    for exon_pos, exon_end_pos in exons:
                        partitions = partition.partition(
                                rname, exon_pos, exon_end_pos, self.bin_size
                            )
                        for (partition_id, partition_start, 
                                partition_end) in partitions:
                            assert exon_pos < partition_end
                            # Print increment at interval start
                            print >>self.output_stream, \
                                'exon_diff\t%s\t%s\t%012d\t1' \
                                % (partition_id,
                                    sample_index,
                                    max(partition_start, exon_pos))
                            output_line_count += 1
                            assert exon_end_pos > partition_start
                            if exon_end_pos <= partition_end:
                                '''Print decrement at interval end 
                                iff exon ends before partition
                                ends.'''
                                print >>self.output_stream, \
                                    'exon_diff\t%s\t%s\t' \
                                    '%012d\t-1' \
                                    % (partition_id, 
                                        sample_index,
                                        exon_end_pos)
                                output_line_count += 1
                try:
                    reverse_strand_string = [field for field in alignment
                                    if field[:5] == 'XS:A:'][0][5:]
                    # Output introns
                    for (intron_pos, intron_end_pos,
                            left_displacement, right_displacement) \
                        in introns:
                        print >>self.output_stream, (
                                ('intron_bed\tN\t%s\t%s\t%012d\t%012d\t%s\t'
                                 '%d\t%d\t1')
                                 % (sample_index, 
                                    self.reference_index.\
                                    rname_to_string[rname],
                                    intron_pos, intron_end_pos,
                                    reverse_strand_string,
                                    left_displacement,
                                    right_displacement)
                            )
                        output_line_count += 1
                except IndexError:
                    # No introns
                    pass
            # Write SAM output
            for alignment in multiread_reports_and_ties[0]:
                print >>self.output_stream, 'sam\t' \
                        + '\t'.join(
                            (sample_index,
                                self.reference_index.rname_to_string[
                                        alignment[2]
                                    ], '%012d' % int(alignment[3]),
                                alignment[0].partition('\x1d')[0],
                                alignment[1]) + alignment[4:]
                        )
        try:
            ties_to_print = multiread_reports_and_ties[1]
        except IndexError:
            # No ties
            pass
        else:
            for alignment in ties_to_print:
                qname = alignment[0]
                flag = alignment[1]
                cigar = alignment[5]
                rname = alignment[2]
                pos = int(alignment[3])
                seq = alignment[9]
                md = [field for field in alignment
                            if field[:5] == 'MD:Z:'][0][5:]
                insertions, deletions, introns, exons \
                    = indels_introns_and_exons(cigar, md, pos, seq)
                try:
                    sense = [field[5:] for field in alignment
                            if field[:5] == 'XS:A:'][0]
                except IndexError:
                    pass
                if introns:
                    for intron in introns:
                        print >>self.output_stream, (
                                        ('sam_intron_ties\tN\t%s\t%s\t'
                                         '%012d\t%012d\t%s\t_'
                                         '\t%012d\t%s\t%s\t') % (
                                self.manifest_object.label_to_index[
                                    qname.rpartition('\x1d')[2]
                                ], self.reference_index.rname_to_string[rname],
                                intron[0], intron[1], sense, pos,
                                qname, flag)) + '\t'.join(alignment[4:])
                else:
                    print >>self.output_stream, '\t'.join(('sam_clip_ties',) \
                                                            + alignment)
        return output_line_count

if __name__ == '__main__':
    import unittest

    class TestIndelsIntronsAndExons(unittest.TestCase):
        """ Tests indels_introns_and_exons(); needs no fixture 

            Some examples are ripped from:
            http://onetipperday.blogspot.com/2012/07/
            deeply-understanding-sam-tags.html; others are from actual
            SAM output of a dmel simulation
        """
        def test_read_1(self):
            """ Fails if example doesn't give expected indels/introns/exons."""
            self.assertEquals(([], [(18909816, 'GG')], [], 
                               [(18909796, 18909816), (18909818, 18909827)]),
                         indels_introns_and_exons(
                                '20M2D9M', '20^GG7A1', 18909796,
                                'TAGCCTCTGTCAGCACTCCTGAGTTCAGA')
                    )

        def test_read_2(self):
            """ Fails if example doesn't give expected indels/introns/exons."""
            self.assertEquals(([], [(73888560, 'GG')], [],
                               [(73888540, 73888560), (73888562, 73888571)]),
                         indels_introns_and_exons(
                                '20M2D9M', '20^GG8C0', 73888540,
                                'TAGCCTCTGTCAGCACTCCTGAGTTCAGA')
                    )

        def test_read_3(self):
            """ Fails if example doesn't give expected indels/introns/exons."""
            self.assertEquals(([(20620369, 'CA')], [(20620365, 'GT')],
                               [(20620167, 20620318, 20, 56)],
                               [(20620147, 20620167), (20620318, 20620365),
                                (20620367, 20620374)]),
                         indels_introns_and_exons(
                                '20M151N47M2D3M2I4M', '67^GT3T2C0', 20620147,
                                'CCGCACCCGTACTGCTACAGATTTCCATCATCGCCACCCGCGGGC'
                                'ATTCTGAAAAAGAGCGACGAAGAAGCAACCT')
                    )

        def test_read_4(self):
            """ Fails if example doesn't give expected indels/introns/exons."""
            self.assertEquals(([(20620155, 'CT')], [],
                               [(20620219, 20620289, 74, 2)],
                               [(20620147, 20620219), (20620289, 20620291)]),
                         indels_introns_and_exons(
                                '9M2I63M70N2M', '1A2C1A0G1G1C1C0C1G2A54',
                                 20620147,
                                'TTCTNCCTGCTTGTATGACCGTGTTGGGCGTGAGTGGCTTGTCCC'
                                'TCAAGTAGAGACCATAGCGAGATGGGTACCT')
                    )

    unittest.main()
