"""
pick.py

Routines for aggregating all spliced-alignment information within a cluster
and applying some classifier in order to pick one junction.
"""

import motifs
from collections import defaultdict
import numpy
import os
import site

base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for subdir in [ 'fasta' ]: site.addsitedir(os.path.join(base_path, subdir))

import fasta

def scanExactSpliceMotifs(lfseq, rtseq, lfmotif, rtmotif):
    """ Look for exact motif occurrences, with one end in the st window
        and the other in the en window """
    stmap, enmap = defaultdict(list), defaultdict(list)
    sawi = set()
    for seq, mp, mmap in [(lfseq, stmap, lfmotif), (rtseq, enmap, rtmotif)]:
        for i in xrange(len(seq)-1):
            s = seq[i:i+2]
            if s in mmap:
                for _, motifi, pri in mmap[s]:
                    mp[(motifi, pri)].append(i)
                    sawi.add((motifi, pri))
    for motifi, pri in sawi:
        if (motifi, pri) in stmap and (motifi, pri) in enmap:
            for sti in stmap[(motifi, pri)]:
                for eni in enmap[(motifi, pri)]:
                    yield pri, motifi, sti, eni

class SiteSelector(object):
    """ Encapsulates info needed to pick a splice site from a cluster """
    
    def __init__(self, fastafh, mfs=motifs.sa):
        """ Set up FASTA file handle and motifs """
        self.motifs = mfs
        self.lmotifmap, self.rmotifmap = defaultdict(list), defaultdict(list)
        self.motifs = mfs
        i = 0
        for pri, x, y in mfs:
            self.lmotifmap[x].append((y, i, pri))
            self.rmotifmap[y].append((x, i, pri))
            i += 1
        self.fastafh = fastafh
    
    def handleCluster(self, refid, sts, ens):
        """ Given a collection of genome intervals (intronic chunks), which we
            assume correspond to the same splice junction, guess where the
            junction is.  We use information about (a) splice site motifs and
            (b) interval end-points.  We could also use (c) other motifs, like
            branch site, landing site, (d) reads that aligned partially and
            which we can try to extend over the putative site.  Note that the
            sts and ens lists need not be parallel nor need they have the same
            # elts.  The ens offsets are exclusive. """
        # Move interval ends left by 2 so that instead of being just after the
        # motif, they're at the beginning of the motif
        ens = map(lambda x: x - 2, ens)
        minst, maxst = min(sts), max(sts)
        minen, maxen = min(ens), max(ens)
        minst, minen = max(minst - 2, 0), max(minen - 2, 0)
        maxst, maxen = maxst + 2, maxen + 2
        assert maxst >= minst and maxen >= minen, "%d >= %d, %d >= %d" % (maxst, minst, maxen, minen)
        stseq = self.fastafh.fetch_sequence(refid, minst+1, maxst+1)
        enseq = self.fastafh.fetch_sequence(refid, minen+1, maxen+1)
        sites = [ p for p in scanExactSpliceMotifs(stseq, enseq, self.lmotifmap, self.rmotifmap) ]
        # Calculate a mean and standard deviation for the intronic chunk
        # starts and ends
        stmn, stsd = numpy.mean(sts), max(numpy.std(sts), 1e-6)
        enmn, ensd = numpy.mean(ens), max(numpy.std(ens), 1e-6)
        # Assign score to each candidate site, equal to sum of absolute value
        # of Z-scores for both ends
        scoredSites = []
        for pri, motifi, st, en in sites:
            st += minst; en += minen
            stzsc = abs((st - stmn) / stsd)
            enzsc = abs((en - enmn) / ensd)
            z = stzsc + enzsc
            scoredSites.append((pri, z, self.motifs[motifi], st, en))
        scoredSites.sort()
        return scoredSites

if __name__ == '__main__':
    import unittest
    
    class TestPick(unittest.TestCase):
        
        def test1(self):
            sites = [ site for site in scanExactSpliceMotifs("GT", "AG", {"GT": [("AG", 0, 0)]}, {"AG" : [("GT", 0, 0)]}) ]
            self.assertEqual(1, len(sites))
            self.assertEqual((0, 0, 0, 0), sites[0])
        
        def test2(self):
            sites = [ site for site in scanExactSpliceMotifs("AG", "GT", {"GT": [("AG", 0, 0)]}, {"AG" : [("GT", 0, 0)]}) ]
            self.assertEqual(0, len(sites))
        
        def test3(self):
            sites = [ site for site in scanExactSpliceMotifs("GGGGGTTTTTT", "AAAAAAGGGGGG", {"GT": [("AG", 0, 0)]}, {"AG" : [("GT", 0, 0)]}) ]
            self.assertEqual(1, len(sites))
            self.assertEqual((0, 0, 4, 5), sites[0])
        
        def test4(self):
            sites = [ site for site in scanExactSpliceMotifs("GTGT", "AGAG", {"GT": [("AG", 0, 0)]}, {"AG" : [("GT", 0, 0)]}) ]
            self.assertEqual(4, len(sites))
            self.assertEqual((0, 0, 0, 0), sites[0])
            self.assertEqual((0, 0, 0, 2), sites[1])
            self.assertEqual((0, 0, 2, 0), sites[2])
            self.assertEqual((0, 0, 2, 2), sites[3])
        
        def test5(self):
            sites = [ site for site in scanExactSpliceMotifs("GTCT", "AGAC", {"GT": [("AG", 0, 0)], "CT" : [("AC", 1, 1)]}, {"AG" : [("GT", 0, 0)], "AC" : [("CT", 1, 1)]}) ]
            self.assertEqual(2, len(sites))
            self.assertEqual((0, 0, 0, 0), sites[0])
            self.assertEqual((1, 1, 2, 2), sites[1])
        
        def test6(self):
            # GTs at offsets 10 and 20, AGs at offsets 30 and 40
            fafn, _ = fasta.writeIndexedFasta(\
                "ref1", "GGGGGGGGGGGTGGGGGGGGGTAAAAAAAAAGAAAAAAAAAG")
            #                      ^10       ^20       ^30       ^40
            fafh = fasta.fasta(fafn)
            sel = SiteSelector(fafh, motifs.sa)
            sts = [ 10, 10, 10, 10 ]
            ens = [ 32, 32, 32, 32 ]
            ssite = sel.handleCluster("ref1", sts, ens)
            self.assertEqual(1, len(ssite))
            pri, _, motif, st, en = ssite[0]
            self.assertEqual(0, pri)
            self.assertEqual((0, 'GT', 'AG'), motif)
            self.assertEqual(10, st)
            self.assertEqual(30, en)
        
        def test7(self):
            # GTs at offsets 10 and 20, AGs at offsets 30 and 40
            fafn, _ = fasta.writeIndexedFasta(\
                "ref1", "GGGGGGGGGGGTGGGGGGGGGTAAAAAAAAAGAAAAAAAAAG")
            #                      ^10       ^20       ^30       ^40
            fafh = fasta.fasta(fafn)
            sel = SiteSelector(fafh, motifs.sa)
            sts = [ 0, 10, 15, 20 ]
            ens = [ 27, 29, 29, 31, 31, 32, 32, 33, 34, 35, 37, 39, 42 ]
            ssite = sel.handleCluster("ref1", sts, ens)
            self.assertEqual(4, len(ssite))
            pri, _, motif, st, en = ssite[0]
            self.assertEqual(0, pri)
            self.assertEqual((0, 'GT', 'AG'), motif)
            self.assertEqual(10, st)
            self.assertEqual(30, en)
            pri, _, motif, st, en = ssite[1]
            self.assertEqual(0,  pri)
            self.assertEqual((0, 'GT', 'AG'), motif)
            self.assertEqual(20, st)
            self.assertEqual(30, en)
            pri, _, motif, st, en = ssite[2]
            self.assertEqual(0,  pri)
            self.assertEqual((0, 'GT', 'AG'), motif)
            self.assertEqual(10, st)
            self.assertEqual(40, en)
            pri, _, motif, st, en = ssite[3]
            self.assertEqual(0,  pri)
            self.assertEqual((0, 'GT', 'AG'), motif)
            self.assertEqual(20, st)
            self.assertEqual(40, en)
        
        def test8(self):
            # GTs at offsets 10 and 20, AGs at offsets 30 and 40
            fafn, _ = fasta.writeIndexedFasta(\
                "ref1", "GGGGGGGGGGGTGGGGGGGGGTAAAAAAAAAGAGGG")
            #                      ^10       ^20       ^30
            #                                            ^32
            fafh = fasta.fasta(fafn)
            sel = SiteSelector(fafh, motifs.sa)
            sts = [ 10, 10, 10, 10, 10 ]
            ens = [ 27, 28, 29, 30, 30, 30, 31, 31, 31, 31, 31, 31, 31, 31, 31, 32, 32, 32, 32, 32, 32, 33, 34, 35 ]
            ssite = sel.handleCluster("ref1", sts, ens)
            self.assertEqual(2, len(ssite))
            pri, _, motif, st, en = ssite[0]
            self.assertEqual(0, pri)
            self.assertEqual((0, 'GT', 'AG'), motif)
            self.assertEqual(10, st)
            self.assertEqual(30, en)
            pri, _, _, st, _ = ssite[1]
            self.assertEqual(0, pri)
            self.assertEqual(10, st)
            # Even though 32 is nearby, 30 should win because it's just to the
            # left of the intronic interval end-point mean at 31
    
    unittest.main()
