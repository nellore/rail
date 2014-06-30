#!/usr/bin/env python
"""
cigar_parse.py
Part of Rail-RNA

Includes a function that outputs indels, introns, and exons from a genome
position, CIGAR string, and MD string.
"""
import re

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

def reference_from_seq(cigar, md, seq):
    """ Recovers reference sequence from read sequence and MD string.

        md: an MD string (example: 33A^CC)
        cigar: CIGAR string; used to extract initial soft clip

        Return value: reference sequence
    """
    md = parsed_md(md)
    cigar = re.split(r'([MINDS])', cigar)[:-1]
    if cigar[1] == 'S':
        read_index = int(cigar[0])
    else:
        read_index = 0
    md_index = 0
    md_size = len(md)
    reference_seq = []
    while md_index != md_size:
        try:
            bases_to_cover = int(md[md_index])
            reference_seq.append(seq[read_index:read_index+bases_to_cover])
            read_index += bases_to_cover
            md_index += 1
        except ValueError:
            if md[md_index] == '^':
                # Deletion from reference
                reference_seq.append(md[md_index+1])
                md_index += 2
            else:
                # Substitution
                reference_seq.append(md[md_index])
                md_index += 1
                read_index += 1
    return ''.join(reference_seq)

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
