"""
Parser/container for flux simulator .bed files
"""
import os
import site
from collections import defaultdict

base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
site.addsitedir(os.path.join(base_path, "struct"))
site.addsitedir(os.path.join(base_path, "annotation"))
site.addsitedir(os.path.join(base_path, "fasta"))

import search
import gtf
import fasta

import unittest

class read:
    def __init__(self,chrom,chrom_st,chrom_end,name,score,strand,thick_start,thick_end,itemRGB,blockCount,blockSizes,blockStarts):
        self.refid = chrom
        self.st1 = chrom_st
        self.en1 = chrom_end
        self.name = name
        self.score = score
        self.strand = strand
        self.thick_start = thick_start
        self.thick_end = thick_end
        self.itemRGB = itemRGB
        self.blockCount = blockCount
        self.blockSizes = blockSizes
        self.blockStarts = blockStarts
    def __str__(self):
        return "%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(self.refid,self.st1,self.en1,self.name,self.score,self.strand,self.thick_start,self.thick_end,self.itemRGB,self.blockCount,self.blockSizes,self.blockStarts)

class bed(object):
    def __init__(self,fname,xscripts):
        self.fname = fname
        self.junctions = defaultdict( list )
        self.findAnnotatedJunctions(xscripts)
    def getSites(self):
        sites = []
        for juncs in self.junctions.itervalues():#remove duplicates
            sites+=juncs
        return sites

    def findAnnotatedJunctions(self,xscripts):
        """
        Finds all splice sites that are spanned by a read
        """

        xscripts = {x.xscript_id: x for x in xscripts}
        with open(self.fname,'r') as fh:
            for ln in fh:
                ln = ln.rstrip()
                toks = ln.split("\t")
                assert len(toks)==12

                chrom,chrom_st,chrom_end,name,score,strand,thick_start,thick_end,itemRGB,blockCount,blockSizes,blockStarts = toks
                
                chrom_st,chrom_end = int(chrom_st),int(chrom_end)
                xid = name.split(":")[2]
                sites = xscripts[xid].getOverlapSites(chrom_st,chrom_end)
                self.junctions[chrom] = self.junctions[chrom] + sites
        for key in self.junctions.iterkeys():#remove duplicates
            self.junctions[key] = list( set(self.junctions[key]))
            self.junctions[key].sort(key=lambda tup:tup[0])
    """
    Note: the guess site must be in the following tuple format
    (st,en,refid,'')
    """
    def find(self,refid,guess_site):
        assert self.junctions[refid]>0
        return search.find_tuple(self.junctions[refid],guess_site)

def createTestFasta(fname,refid,refseq):
    fastaH = open(fname,'w')
    fastaIdx = open(fname+".fai",'w')
    fastaH.write(">%s\n%s\n"%(refid,refseq))
    fastaIdx.write("%s\t%d\t%d\t%d\t%d\n"%(refid,len(refseq),len(refid)+2,len(refseq),len(refseq)+1))
    fastaH.close()
    fastaIdx.close()

def createTestGTF(fname,annots):
    gtfH = open(fname,'w')
    gtfH.write(annots)
    gtfH.close()

def createTestBed(fname,frags):
    libH = open(fname,'w')
    libH.write(frags)
    libH.close()

class TestBedFunctions1(unittest.TestCase):

    def setUp(self):
        """       [AACTGTGAT CAAGGATGT] GTCGCTTGTG AAACGAGCAG [TGGATCCG] GCGAAGCACA TTGGCAC[AA GATCGCGC]"""
        """       CAACTGTGAT CAAGGATGTC TTCGCTTGTG AAACGAGCAG CTGGATCCGC CTGAAGCACA TTGGCACTAA GATCGCGCA"""
        #                  ^10        ^20        ^30        ^40        ^50        ^60        ^70       ^80
        refseq="""CAACTGTGATCAAGGATGTCTTCGCTTGTGAAACGAGCGTCTGGATCCGCCTGAAGCACATTGGCACTAAGATCGCGCA"""
        annots="\n".join( map( str, [gtf.Annot(refid='chr2R',st0=11,en0=20,orient='-',feature='exon',
                                               score='.',frame='.',attrs='CG17528\tNM_143888'),
                                     gtf.Annot(refid='chr2R',st0=51,en0=60,orient='-',feature='exon',
                                               score='.',frame='.',attrs='CG17528\tNM_143888')]))
        reads = []
        #reads.append( read(5,15,"chr2L:106903-114432W:NM_205880:30:2387:281:549/2"))
        reads.append(read('chr2R',5,15,'chrR:21471627-21476294C:NM_143888:369:3135:2411:2786/2','0','+','.','.','0,0,0','1','76',0))
        read_str = '\n'.join( map( str, reads))
        self.fasta = "test.fa"
        self.faidx = "test.fa.fai"
        self.gtf   = "test.gtf"
        self.bed   = "test.bed"
        createTestFasta(self.fasta,"chr2R",refseq)
        createTestGTF(self.gtf,annots)
        createTestGTF(self.gtf,annots)
        createTestBed(self.bed,read_str)

    def tearDown(self):
        os.remove(self.fasta)
        os.remove(self.faidx)
        os.remove(self.gtf)
        os.remove(self.bed)

    def testFind1(self):
        annots   = gtf.parseGTF([self.gtf])
        fastadb  = gtf.parseFASTA([self.fasta])
        xscripts = gtf.assembleTranscripts(annots,fastadb)
        testbed  = bed(self.bed,xscripts)
        site = testbed.find("chr2R",(18,19,"chr2R",'NM_143888'))
        self.assertEquals(site,(20,21,"chr2R",'NM_143888'))

class TestBedFunctions2(unittest.TestCase):
    def setUp(self):
        """
        AACTGTGATC 10AAGGATGTCT 20TCGCTTGTGA 30AACGAACGTC 40TGGATCCGCC 50TGAAGCATAT 60TGGCAATAAG 70ATCGCGCACT 80CCGTCCCCTA 90TATCAGACCA 100GTAAGGATCT 110GGAAACTCAT 120AAATGCCAGA 130GATAATAGCG 140TCGAATAGCG 150GCTCCTGTTG 160GTTGTCAGGC 170GCTACAAACG 180GCGGAAAACC 190GCATAGAAGT 200ATATACAAAA 210TGATTCCAGC 220GGCCCAAACG 230TCAATCTTTA 240AGAATATATA 250AGTATATTAA 260TTTTTAAGAA 270AGCTATTATT 280TACTATTACC 290TTTAGCCCAT 300ATCCGACTTC 310TAACAATATC 320TCTGGTGCTA 330CATATGTGGG 340AGTCCCGCAG 350ACAGCATATA 360GAAGATCGTT 370TACCTCACAT 380GCTAATCCAA 390AATCAGCAAG 400CTTTAGTTCA 410AGAACGTTTC 420CATGCTCATC 430CAATTTCACC 440TGCAAATAGT 450TCCAATGATT 460AGTCCTATTA 470TTGCTCAGTA 480TTTTTCATTA 490GTCATGTCAT 500ACTAGAAGAT 510TCTCAGGTTT 520TATATCTCGA 530TGCACAATGC 540CCATTGAATG 550CAGATAAGTC 560ATGGCCGCAC 570CCAAATGTCT 580AATCATAATG 590CGCGACTGGT 600TTTCTGAGAA 610CCTCGTTACC 620TGAGTTATAG 630CGTCGAATAA 640ATCACCACCT 650GTTACAAATT 660TAACAACAAA 670TTTAAGCAGA 680AAGAAGCTGT 690GAATGTATTC 700ACATTCAACT 710TACCACTTAC 720ATATTCCAGT 730ACAAGATACA 740TATTTGTATT 750TTGGTCTACA 760CTCAAAATAA 770GCGAAATTAT 780ATGCGGATGA 790TTTAATTTTT 800TCATGACGCG 810AACTTCCGCG 820TCTATGTAGT 830GTTCCTTGCC 840CTTGCATTTG 850TTCTTGTCTA 860TTATTTTTAG 870GGCGTAGGAA 880TGACCAGTTT 890GGCGGTGCTT 900AATCTTAAAC 910ACAATGGCAA 920AGTTGCCGTC 930GCCAATTATT 940CTTCCCAGAG 950AGTAAGTGTT 960ACGAATATTC 970GATGGCAATT 980CATTGATCTC 990CATACCAGTG 1000
        """
        refseq="""AACTGTGATCAAGGATGTCTTCGCTTGTGAAACGAACGTCTGGATCCGCCTGAAGCATATTGGCAATAAGATCGCGCACTCCGTCCCCTATATCAGACCAGTAAGGATCTGGAAACTCATAAATGCCAGAGATAATAGCGTCGAATAGCGGCTCCTGTTGGTTGTCAGGCGCTACAAACGGCGGAAAACCGCATAGAAGTATATACAAAATGATTCCAGCGGCCCAAACGTCAATCTTTAAGAATATATAAGTATATTAATTTTTAAGAAAGCTATTATTTACTATTACCTTTAGCCCATATCCGACTTCTAACAATATCTCTGGTGCTACATATGTGGGAGTCCCGCAGACAGCATATAGAAGATCGTTTACCTCACATGCTAATCCAAAATCAGCAAGCTTTAGTTCAAGAACGTTTCCATGCTCATCCAATTTCACCTGCAAATAGTTCCAATGATTAGTCCTATTATTGCTCAGTATTTTTCATTAGTCATGTCATACTAGAAGATTCTCAGGTTTTATATCTCGATGCACAATGCCCATTGAATGCAGATAAGTCATGGCCGCACCCAAATGTCTAATCATAATGCGCGACTGGTTTTCTGAGAACCTCGTTACCTGAGTTATAGCGTCGAATAAATCACCACCTGTTACAAATTTAACAACAAATTTAAGCAGAAAGAAGCTGTGAATGTATTCACATTCAACTTACCACTTACATATTCCAGTACAAGATACATATTTGTATTTTGGTCTACACTCAAAATAAGCGAAATTATATGCGGATGATTTAATTTTTTCATGACGCGAACTTCCGCGTCTATGTAGTGTTCCTTGCCCTTGCATTTGTTCTTGTCTATTATTTTTAGGGCGTAGGAATGACCAGTTTGGCGGTGCTTAATCTTAAACACAATGGCAAAGTTGCCGTCGCCAATTATTCTTCCCAGAGAGTAAGTGTTACGAATATTCGATGGCAATTCATTGATCTCCATACCAGTG"""

        annots="\n".join( map( str, [gtf.Annot(refid='chr2R',st0=11,en0=20,orient='-',feature='exon',
                                               score='.',frame='.',attrs='CG17528\tNM_001042999'),
                                     gtf.Annot(refid='chr2R',st0=31,en0=40,orient='-',feature='exon',
                                               score='.',frame='.',attrs='CG17528\tNM_001042999'),
                                     gtf.Annot(refid='chr2R',st0=51,en0=60,orient='-',feature='exon',
                                               score='.',frame='.',attrs='CG17528\tNM_001042999'),
                                     gtf.Annot(refid='chr2R',st0=71,en0=80,orient='-',feature='exon',
                                               score='.',frame='.',attrs='CG17528\tNM_001042999')
                                               ]))
        reads = []
        reads.append(read('chr2R',5,15,'chrR:21471627-21476294C:NM_001042999:369:3135:2411:2786/2','0','+','.','.','0,0,0','1','76',0))
        reads.append(read('chr2R',15,25,'chrR:21471627-21476294C:NM_001042999:369:3135:2411:2786/2','0','+','.','.','0,0,0','1','76',0))
        reads.append(read('chr2R',25,35,'chrR:21471627-21476294C:NM_001042999:369:3135:2411:2786/2','0','+','.','.','0,0,0','1','76',0))
        read_str = ''.join( map( str, reads))
        self.fasta = "test.fa"
        self.faidx = "test.fa.fai"
        self.gtf   = "test.gtf"
        self.bed   = "test.bed"
        createTestFasta(self.fasta,"chr2R",refseq)
        createTestGTF(self.gtf,annots)
        createTestGTF(self.gtf,annots)
        createTestBed(self.bed,read_str)
    def tearDown(self):
        os.remove(self.fasta)
        os.remove(self.faidx)
        os.remove(self.gtf)
        os.remove(self.bed)

    def testFind1(self):
        annots   = gtf.parseGTF([self.gtf])
        fastadb  = gtf.parseFASTA([self.fasta])
        xscripts = gtf.assembleTranscripts(annots,fastadb)
        testbed  = bed(self.bed,xscripts)
        site = testbed.find("chr2R",(18,19,"chr2R",'NM_001042999'))
        self.assertEquals(site,(20,21,"chr2R",'NM_001042999'))
        site = testbed.find("chr2R",(27,28,"chr2R",'NM_001042999'))
        self.assertEquals(site,(28,29,"chr2R",'NM_001042999'))
        site = testbed.find("chr2R",(40,41,"chr2R",'NM_001042999'))
        self.assertEquals(site,(40,41,"chr2R",'NM_001042999'))
        site = testbed.find("chr2R",(47,48,"chr2R",'NM_001042999'))
        self.assertEquals(site,(48,49,"chr2R",'NM_001042999'))


if __name__=="__main__":
    unittest.main()
