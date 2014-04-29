""" Parse GTF files, combine gene annotation information with DNA sequences
    from FASTA files, and extract substrings from transcripts. """

__author__ = "Ben Langmead"

from collections import defaultdict
import string
import re

_revcomp_trans = string.maketrans("ACGT", "TGCA")


def revcomp(s):
    return s[::-1].translate(_revcomp_trans)


class GeneAnnotation(object):
    """
    Encapsulates gene annotation information, telling where each
    component of each gene lies with respect to the reference genome.
    """

    class Exon(object):
        """
        Offsets are 0-based and intervals are right-open.  All
        intervals are with respect to the Watson refrence strand.
        """

        def __init__(self, ref_id, st, en, fw):
            self.ref_id = ref_id
            self.st = st  # start position w/r/t reference Watson
            self.en = en  # end position w/r/t reference Watson
            self.fw = fw
            self.donor_after = self.acceptor_before = None

        def set_donor_after(self, donor):
            assert len(donor) == 2
            self.donor_after = donor

        def set_acceptor_before(self, acceptor):
            assert len(acceptor) == 2
            self.acceptor_before = acceptor

        def __len__(self):
            return self.en - self.st

        def __lt__(self, other):
            if self.ref_id != other.ref_id:
                return self.ref_id < other.ref_id
            elif self.fw != other.fw:
                return self.fw
            elif self.st != other.st:
                return self.st < other.st
            else:
                return self.en < other.en

        def __repr__(self):
            return "%s%s:[%d-%d)" % (self.ref_id, '+' if self.fw else '-', self.st, self.en)

    class Transcript(object):
        def __init__(self, exons, cdss=None):
            """ Create new transcript from list of exons.  TODO: also
                keep track of coding sequences.
            """
            assert len(exons) > 0
            assert all([exons[i].fw == exons[0].fw for i in xrange(len(exons))])
            assert all([exons[i].ref_id == exons[0].ref_id for i in xrange(len(exons))])
            self.fw = exons[0].fw
            self.exons = sorted(exons)  # exons
            self.ref_id = self.exons[0].ref_id
            self.ref_st, self.ref_en = self.exons[0].st, self.exons[-1].en

            assert self.ref_en > self.ref_st
            self.cdss = cdss  # coding sequences
            if self.cdss is not None:
                self.cdss.sort()
            self.nuc_length = sum(map(len, self.exons))
            self.seq = None  # nucleotide sequence

        def __iter__(self):
            """ Return iterator over exons in transcript """
            return iter(self.exons)

        def exon_length(self):
            """ Return the number of exons in transcript """
            return len(self.exons)

        def nucleotide_length(self):
            """ Return number of nucleotides covered by transcript """
            return self.nuc_length

        def unspliced_substring_from_5p(self, offset, length, ref=None):
            """ Return a substring of this transcript.  The substring
                begins at 0-based 'offset' from 5' end.  Offset can be
                negative, in which case the substring starts before
                the 5' edge.
            """
            if offset >= 0 and self.seq is not None:
                if offset + length <= self.nuc_length:
                    return self.seq[offset:offset+length]
                else:
                    # add some poly-A tail
                    return self.seq[offset:] + ('A' * (offset + length - self.nuc_length))
            else:
                # Either it starts upstream of the TSS or we haven't parsed
                # this transcript's sequence from the FASTA.  Either way, we
                # have to consult the FASTA here.
                assert ref is not None
                if self.seq is None:
                    self.set_sequence(ref)
                if self.fw:
                    seq_before = ref.get(self.ref_id, self.ref_st + offset, -offset)
                else:
                    seq_before = revcomp(ref.get(self.ref_id, self.ref_en, -offset))
                if length < -offset:
                    return seq_before[:length]
                else:
                    return seq_before + self.seq[:length+offset]

        def spliced_substring_from_5p(self, offset, length, ref=None):
            """ Return substring along with some information about
                introns spanned.  Intron spanning information TODO: should it be w/r/t substring or w/r/t 5' end of transcript?
            """
            seq = self.unspliced_substring_from_5p(offset, length, ref)
            intron_offsets, ref_ivals = [], []
            exon_begin_5p, exon_end_5p = 0, 0
            if self.fw:
                ref_5p = self.exons[0].st + offset
            else:
                ref_5p = self.exons[-1].en - offset
            for i, exon in enumerate(self.exons if self.fw else self.exons[::-1]):
                exon_begin_5p = exon_end_5p
                exon_end_5p = exon_begin_5p + len(exon)
                exon_olap = 0
                if i > 0 and offset < exon_begin_5p < offset + length:
                    intron_offsets.append(exon_begin_5p - offset)
                    diff = (self.exons[i].st - self.exons[i-1].en)
                    intron_ref_st = ref_5p - (0 if self.fw else diff)
                    intron_ref_en = ref_5p + (diff if self.fw else 0)
                    ref_ivals.append((False, intron_ref_st, intron_ref_en))
                if (exon_begin_5p <= offset < exon_end_5p) or\
                        (exon_begin_5p <= (offset + length) < exon_end_5p) or\
                        (exon_begin_5p >= offset and exon_end_5p <= offset + length):
                    olap_begin, olap_end = max(exon_begin_5p, offset) - exon_begin_5p,\
                                           min(exon_end_5p, offset + length) - exon_begin_5p
                    exon_olap = olap_end - olap_begin
                    # TODO: handle substrings that start upstream of TSS or end in poly-A
                    if self.fw:
                        ref_ivals.append((True, olap_begin + exon.st, olap_end + exon.st))
                    else:
                        ref_ivals.append((True, exon.en - olap_end, exon.en - olap_begin))
                diff = exon_olap
                if i == 0 and offset < 0:
                    diff -= offset
                if i > 0:
                    diff += (self.exons[i].st - self.exons[i-1].en)
                if self.fw:
                    ref_5p += diff
                else:
                    ref_5p -= diff
            return seq, intron_offsets, ref_ivals if self.fw else ref_ivals[::-1]

        def set_sequence(self, ref):
            """ Given reference genome sequence, fill in information
                about the sequence of this transcript, each intron's
                donor and acceptor sequences, etc.
            """
            # Set transcript sequence and donor/acceptor sequences
            seqs = []
            if self.fw:
                for i, ex in enumerate(self.exons):
                    seqs.append(ref.get(ex.ref_id, ex.st, ex.en - ex.st).upper())
                    if i != 0:
                        ex.set_acceptor_before(ref.get(ex.ref_id, ex.st - 2, 2).upper())
                    if i != len(self.exons) - 1:
                        ex.set_donor_after(ref.get(ex.ref_id, ex.en, 2).upper())
                self.seq = ''.join(seqs)
            else:
                for i, ex in enumerate(self.exons):
                    seqs.append(revcomp(ref.get(ex.ref_id, ex.st, ex.en - ex.st).upper()))
                    if i != len(self.exons) - 1:
                        ex.set_acceptor_before(revcomp(ref.get(ex.ref_id, ex.en, 2).upper()))
                    if i != 0:
                        ex.set_donor_after(revcomp(ref.get(ex.ref_id, ex.st - 2, 2).upper()))
                self.seq = ''.join(seqs[::-1])
                # Set start and stop codon sequences
            if self.cdss is not None:
                pass  # not implemented
            assert len(self.seq) == self.nuc_length

    def __init__(self, genes):
        """ Initialize from a data structure: gene_id -> transcript_id -> Transcript
        """
        self.genes = genes
        self.sequence_set = False
        self.xscripts = {}
        self.exons = []
        for gene_id in genes.iterkeys():
            for xscript_id in genes[gene_id].iterkeys():
                assert xscript_id not in self.xscripts
                self.xscripts[xscript_id] = (gene_id, genes[gene_id][xscript_id])
                self.exons.extend(genes[gene_id][xscript_id])

    def summarize(self):
        xscripts_per_gene = defaultdict(int)
        for gene_id, xscripts in self.genes.iteritems():
            xscripts_per_gene[len(xscripts)] += 1
        print "Num genes: %d" % len(self.genes)
        print "Num transcripts: %d" % len(self.xscripts)
        print "Histogram of transcripts per gene:"
        for i, num in sorted(xscripts_per_gene.iteritems()):
            print "  %d: %d" % (i, num)
        exons_per_xscript = defaultdict(int)
        exon_lengths = defaultdict(int)
        for xscript_tup in self.xscripts.itervalues():
            _, xscript = xscript_tup
            for exon in xscript:
                exon_lengths[len(exon)] += 1
            exons_per_xscript[xscript.exon_length()] += 1
        print "Exons per transcript:"
        for i, num in sorted(exons_per_xscript.iteritems()):
            print "  %d: %d" % (i, num)
            # Donors and acceptors
        if self.sequence_set:
            donors, acceptors = defaultdict(int), defaultdict(int)
            for exon in self.exons:
                if exon.donor_after is not None:
                    donors[exon.donor_after] += 1
                if exon.acceptor_before is not None:
                    acceptors[exon.acceptor_before] += 1
            print "Donors:"
            for dmotif, dnum in sorted(donors.iteritems()):
                print "  %s: %d" % (dmotif, dnum)
            print "Acceptors:"
            for amotif, anum in sorted(acceptors.iteritems()):
                print "  %s: %d" % (amotif, anum)

    @classmethod
    def from_gtf_file(cls, fns):
        gene_id_re = re.compile("gene_id \"([^\"]+)\"")
        xscript_id_re = re.compile("transcript_id \"([^\"]+)\"")
        genes = defaultdict(lambda: defaultdict(list))
        for fn in fns:
            with open(fn) as fh:
                for ln in fh:
                    if ln[0] == '#':
                        continue  # ignore comments
                    ln = ln.rstrip()
                    if len(ln) == 0:
                        continue  # ignore blank lines
                    toks = ln.split('\t')
                    assert len(toks) == 9
                    if toks[2] == 'exon':
                        st, en = int(toks[3]), int(toks[4])
                        assert en >= st
                        gene_id_ma = gene_id_re.search(toks[8])
                        gene_id = gene_id_ma.group(1)
                        xscript_id_ma = xscript_id_re.search(toks[8])
                        xscript_id = xscript_id_ma.group(1)
                        strand = toks[6]
                        assert strand in '+-'
                        genes[gene_id][xscript_id].append(GeneAnnotation.Exon(toks[0], st - 1, en, strand == '+'))
        for gk in genes.iterkeys():
            for xk in genes[gk].iterkeys():
                genes[gk][xk] = GeneAnnotation.Transcript(genes[gk][xk])
        return GeneAnnotation(genes)

    def set_sequence(self, ref):
        genes = self.genes
        for gene_id in genes.iterkeys():
            for xscript_id in genes[gene_id].iterkeys():
                genes[gene_id][xscript_id].set_sequence(ref)
        self.sequence_set = True

if __name__ == '__main__':

    import unittest

    class ReferenceString(object):

        def __init__(self, s):
            self.s = s

        def get(self, _, start, ln):
            return self.s[start:start+ln]

    class TestExon(unittest.TestCase):

        def test_basics(self):
            ex = GeneAnnotation.Exon('chr1', 10, 100, True)
            self.assertEqual(len(ex), 90)
            self.assertLess(ex, GeneAnnotation.Exon('chr2', 10, 100, True))
            self.assertLess(ex, GeneAnnotation.Exon('chr1', 10, 100, False))
            self.assertLess(ex, GeneAnnotation.Exon('chr1', 11, 100, True))
            self.assertLess(ex, GeneAnnotation.Exon('chr1', 10, 101, True))

    class TestTranscript(unittest.TestCase):

        def test_basic_members_fw(self):
            exons = [GeneAnnotation.Exon('chr1', 10, 100, True),
                     GeneAnnotation.Exon('chr1', 110, 200, True),
                     GeneAnnotation.Exon('chr1', 210, 400, True)]
            xt = GeneAnnotation.Transcript(exons)
            self.assertEqual(3, xt.exon_length())
            self.assertEqual(180 + 190, xt.nucleotide_length())
            self.assertEqual(3, len([x for x in xt]))

        def test_basic_members_rc(self):
            exons = [GeneAnnotation.Exon('chr1', 210, 400, False),
                     GeneAnnotation.Exon('chr1', 110, 200, False),
                     GeneAnnotation.Exon('chr1', 10, 100, False)]
            xt = GeneAnnotation.Transcript(exons)
            self.assertEqual(3, xt.exon_length())
            self.assertEqual(180 + 190, xt.nucleotide_length())
            self.assertEqual(3, len([x for x in xt]))

        def test_set_sequence_fw(self):
            ref = ReferenceString(('A' * 25) + ('C' * 20) + ('G' * 30))
            exons = [GeneAnnotation.Exon('chr1', 10, 20, True),
                     GeneAnnotation.Exon('chr1', 30, 40, True),
                     GeneAnnotation.Exon('chr1', 50, 60, True)]
            xt = GeneAnnotation.Transcript(exons)
            xt.set_sequence(ref)
            self.assertEquals('AAAAAAAAAACCCCCCCCCCGGGGGGGGGG', xt.seq)
            self.assertEquals(None, xt.exons[0].acceptor_before)
            self.assertEquals('AA', xt.exons[0].donor_after)
            self.assertEquals('CC', xt.exons[1].acceptor_before)
            self.assertEquals('CC', xt.exons[1].donor_after)
            self.assertEquals('GG', xt.exons[2].acceptor_before)
            self.assertEquals(None, xt.exons[2].donor_after)

        def test_set_sequence_rc(self):
            ref = ReferenceString(('A' * 25) + ('C' * 20) + ('G' * 30))
            exons = [GeneAnnotation.Exon('chr1', 50, 60, False),
                     GeneAnnotation.Exon('chr1', 30, 40, False),
                     GeneAnnotation.Exon('chr1', 10, 20, False)]
            xt = GeneAnnotation.Transcript(exons)
            xt.set_sequence(ref)
            self.assertEquals('CCCCCCCCCCGGGGGGGGGGTTTTTTTTTT', xt.seq)
            # Note: exons sorted along ref Watson in Transcript constructor
            self.assertEquals(None, xt.exons[2].acceptor_before)
            self.assertEquals('CC', xt.exons[2].donor_after)
            self.assertEquals('GG', xt.exons[1].acceptor_before)
            self.assertEquals('GG', xt.exons[1].donor_after)
            self.assertEquals('TT', xt.exons[0].acceptor_before)
            self.assertEquals(None, xt.exons[0].donor_after)

        def test_unspliced_substring_from_5p_fw(self):
            ref = ReferenceString(('A' * 25) + ('C' * 20) + ('G' * 30))
            exons = [GeneAnnotation.Exon('chr1', 10, 20, True),
                     GeneAnnotation.Exon('chr1', 30, 40, True),
                     GeneAnnotation.Exon('chr1', 50, 60, True)]
            for i in xrange(2):
                xt = GeneAnnotation.Transcript(exons)
                if i == 1:
                    xt.set_sequence(ref)
                # whole thing
                self.assertEquals('AAAAAAAAAACCCCCCCCCCGGGGGGGGGG', xt.unspliced_substring_from_5p(0, 30, ref))
                # negative offset
                self.assertEquals('AAAAAAAAAAACCCCCCCCCCGGGGGGGGG', xt.unspliced_substring_from_5p(-1, 30, ref))
                self.assertEquals('AAAAAAAAAAAAAAACCCCCCCCCCGGGGG', xt.unspliced_substring_from_5p(-5, 30, ref))
                # falls off and acquires poly-A tail
                self.assertEquals('AAAAAAAAAACCCCCCCCCCGGGGGGGGGGAAAAA', xt.unspliced_substring_from_5p(0, 35, ref))
                # from the middle
                self.assertEquals('ACCCCCCCCCCG', xt.unspliced_substring_from_5p(9, 12, ref))

        def test_unspliced_substring_from_5p_rc(self):
            ref = ReferenceString(('A' * 25) + ('C' * 20) + ('G' * 30))
            exons = [GeneAnnotation.Exon('chr1', 50, 60, False),
                     GeneAnnotation.Exon('chr1', 30, 40, False),
                     GeneAnnotation.Exon('chr1', 10, 20, False)]
            for i in xrange(2):
                xt = GeneAnnotation.Transcript(exons)
                if i == 1:
                    xt.set_sequence(ref)
                # whole thing
                self.assertEquals('CCCCCCCCCCGGGGGGGGGGTTTTTTTTTT', xt.unspliced_substring_from_5p(0, 30, ref))
                # negative offset
                self.assertEquals('CCCCCCCCCCCGGGGGGGGGGTTTTTTTTT', xt.unspliced_substring_from_5p(-1, 30, ref))
                self.assertEquals('CCCCCCCCCCCCCCCGGGGGGGGGGTTTTT', xt.unspliced_substring_from_5p(-5, 30, ref))
                # falls off and acquires poly-A tail
                self.assertEquals('CCCCCCCCCCGGGGGGGGGGTTTTTTTTTTAAAAA', xt.unspliced_substring_from_5p(0, 35, ref))
                # from the middle
                self.assertEquals('CGGGGGGGGGGT', xt.unspliced_substring_from_5p(9, 12, ref))

        def test_spliced_substring_from_5p_fw(self):
            ref = ReferenceString(('A' * 25) + ('C' * 20) + ('G' * 30))
            exons = [GeneAnnotation.Exon('chr1', 10, 20, True),
                     GeneAnnotation.Exon('chr1', 30, 40, True),
                     GeneAnnotation.Exon('chr1', 50, 60, True)]
            for i in xrange(2):
                xt = GeneAnnotation.Transcript(exons)
                if i == 1:
                    xt.set_sequence(ref)
                # whole thing
                seq, intron_offsets, ref_ivals = xt.spliced_substring_from_5p(0, 30, ref)
                self.assertEquals([10, 20], intron_offsets)
                self.assertEquals([(True, 10, 20), (False, 20, 30), (True, 30, 40), (False, 40, 50), (True, 50, 60)],
                                  ref_ivals)
                self.assertEquals('AAAAAAAAAACCCCCCCCCCGGGGGGGGGG', seq)
                # negative offset
                seq, intron_offsets, ref_ivals = xt.spliced_substring_from_5p(-1, 30, ref)
                self.assertEquals([11, 21], intron_offsets)
                self.assertEquals([(True, 10, 20), (False, 20, 30), (True, 30, 40), (False, 40, 50), (True, 50, 59)],
                                  ref_ivals)
                self.assertEquals('AAAAAAAAAAACCCCCCCCCCGGGGGGGGG', seq)
                seq, intron_offsets, ref_ivals = xt.spliced_substring_from_5p(-5, 30, ref)
                self.assertEquals([15, 25], intron_offsets)
                self.assertEquals([(True, 10, 20), (False, 20, 30), (True, 30, 40), (False, 40, 50), (True, 50, 55)],
                                  ref_ivals)
                self.assertEquals('AAAAAAAAAAAAAAACCCCCCCCCCGGGGG', seq)
                # falls off and acquires poly-A tail
                seq, intron_offsets, ref_ivals = xt.spliced_substring_from_5p(0, 35, ref)
                self.assertEquals([10, 20], intron_offsets)
                self.assertEquals([(True, 10, 20), (False, 20, 30), (True, 30, 40), (False, 40, 50), (True, 50, 60)],
                                  ref_ivals)
                self.assertEquals('AAAAAAAAAACCCCCCCCCCGGGGGGGGGGAAAAA', seq)
                # from the middle
                seq, intron_offsets, ref_ivals = xt.spliced_substring_from_5p(9, 12, ref)
                self.assertEquals([1, 11], intron_offsets)
                self.assertEquals([(True, 19, 20), (False, 20, 30), (True, 30, 40), (False, 40, 50), (True, 50, 51)],
                                  ref_ivals)
                self.assertEquals('ACCCCCCCCCCG', seq)

        def test_spliced_substring_from_5p_rc(self):
            ref = ReferenceString(('A' * 25) + ('C' * 20) + ('G' * 30))
            exons = [GeneAnnotation.Exon('chr1', 50, 60, False),
                     GeneAnnotation.Exon('chr1', 30, 40, False),
                     GeneAnnotation.Exon('chr1', 10, 20, False)]
            for i in xrange(2):
                xt = GeneAnnotation.Transcript(exons)
                if i == 1:
                    xt.set_sequence(ref)
                # whole thing
                seq, intron_offsets, ref_ivals = xt.spliced_substring_from_5p(0, 30, ref)
                self.assertEquals([10, 20], intron_offsets)
                self.assertEquals([(True, 10, 20), (False, 20, 30), (True, 30, 40), (False, 40, 50), (True, 50, 60)],
                                  ref_ivals)
                self.assertEquals('CCCCCCCCCCGGGGGGGGGGTTTTTTTTTT', seq)
                # negative offset
                seq, intron_offsets, ref_ivals = xt.spliced_substring_from_5p(-1, 30, ref)
                self.assertEquals([11, 21], intron_offsets)
                self.assertEquals([(True, 11, 20), (False, 20, 30), (True, 30, 40), (False, 40, 50), (True, 50, 60)],
                                  ref_ivals)
                self.assertEquals('CCCCCCCCCCCGGGGGGGGGGTTTTTTTTT', seq)
                seq, intron_offsets, ref_ivals = xt.spliced_substring_from_5p(-5, 30, ref)
                self.assertEquals([15, 25], intron_offsets)
                self.assertEquals([(True, 15, 20), (False, 20, 30), (True, 30, 40), (False, 40, 50), (True, 50, 60)],
                                  ref_ivals)
                self.assertEquals('CCCCCCCCCCCCCCCGGGGGGGGGGTTTTT', seq)
                # falls off and acquires poly-A tail
                seq, intron_offsets, ref_ivals = xt.spliced_substring_from_5p(0, 35, ref)
                self.assertEquals([10, 20], intron_offsets)
                self.assertEquals([(True, 10, 20), (False, 20, 30), (True, 30, 40), (False, 40, 50), (True, 50, 60)],
                                  ref_ivals)
                self.assertEquals('CCCCCCCCCCGGGGGGGGGGTTTTTTTTTTAAAAA', seq)
                # from the middle
                seq, intron_offsets, ref_ivals = xt.spliced_substring_from_5p(9, 12, ref)
                self.assertEquals([1, 11], intron_offsets)
                self.assertEquals([(True, 19, 20), (False, 20, 30), (True, 30, 40), (False, 40, 50), (True, 50, 51)],
                                  ref_ivals)
                self.assertEquals('CGGGGGGGGGGT', seq)

    unittest.main()