#!/usr/bin/env python
"""
cigar_parse.py
Part of Rail-RNA

Includes a function that outputs indels, introns, and exons from a genome
position, CIGAR string, and MD string and a function that inserts introns
in a CIGAR string.
"""
import re
import random

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
            of SAM representing an alignment.
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
        if not tokens[-1]:
            # No introns can be found
            pos = offset + int(tokens[1])
            rname = tokens[0]
            '''Use second field in each element of new_multiread to store which
            items should be tested to find whether two alignments are
            "identical".'''
            new_multiread.append(
                    ([qname.partition('\x1d')[0], flag | 256, rname,
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
        # Count number of samples in which intron combo was initially detected
        new_multiread.append(
                    ([alignment[0].partition('\x1d')[0], flag | 256,
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
    deduped_multiread = [new_multiread[0][0]]
    for i in xrange(1, len(new_multiread)):
        if new_multiread[i][1] == new_multiread[i-1][1]:
            continue
        deduped_multiread.append(new_multiread[i][0])
    deduped_multiread = [(alignment,
                             int([field for field in alignment
                                    if str(field)[:5] == 'AS:i:'][0][5:]))
                            for alignment in deduped_multiread]
    deduped_multiread.sort(key=lambda alignment: alignment[1])
    max_score = max([alignment[1] for alignment in deduped_multiread])
    ties = [alignment[0] for alignment in deduped_multiread
                            if alignment[1] == max_score]
    tie_count = len(ties)
    if tie_count > 1:
        random.seed(qname + seq + qual)
        to_primary = random.randint(0, tie_count - 1)
    else:
        to_primary = 0
    ties[to_primary][1] &= ~256
    multiread_to_return = [ties[to_primary]] + [ties[i] 
                                                for i in xrange(len(ties))
                                                if i != to_primary]
    multiread_to_return += [alignment[0] for alignment in deduped_multiread
                                if alignment[1] != max_score]
    alignment_count = len(multiread_to_return)
    for i in xrange(alignment_count):
        multiread_to_return[i][1] = str(multiread_to_return[i][1])
    NH_field = 'NH:i:%d' % alignment_count
    return [alignment + [NH_field] for alignment in multiread_to_return]

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
