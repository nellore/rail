"""
display.py

Displays information from validation.py
"""

import re
import os
import site
import argparse
import sys
import math
import pickle
import bisect
import string
from collections import Counter
from collections import defaultdict
base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

import unittest

site.addsitedir(os.path.join(base_path, "annotation"))
site.addsitedir(os.path.join(base_path, "struct"))
site.addsitedir(os.path.join(base_path, "fasta"))
site.addsitedir(os.path.join(base_path, "statsmath"))

import gtf
import search
import fasta
import window

#Formats the sequence into a more readable format
def format_seq(S):
    return "".join(list(S))
#Formats lists into a more readable format
def format_list(L):
    return " ".join(["%.1f" % i for i in L])
    #Discretizes the floats into ints
    #return "".join(["%d" % int((i+5)*10) for i in L])



#No flanking site required
def printLeftSite(seqid,exon,site,annot_site,win_radius,display_st,display_end,fnh):
    """
    display_window                 (----------------------)
    exon                       |--------------|
    intron                                    |--------------------------->
    ref_start |----------- ... ------------------------------------------->
    """
    #convert to display window coordinate frame for exon information
    win_length = 2*win_radius
    display_st,display_end = site[0] - win_radius, site[0] + win_radius
    display_exon_st   = display_st - exon.st0+1
    display_exon_end  = exon.en0 - exon.st0
    exon_seq = exon.seq[display_exon_st:display_exon_end]
    #print >> sys.stderr,"Whole exon",exon.seq
    eLen = len(exon_seq)

    #Stay in reference coordinates for intron
    remaining = win_length - eLen
    display_intron_st = display_st+eLen+1
    display_intron_end= display_intron_st + remaining
    intron_seq = fnh.fetch_sequence(seqid,display_intron_st+1,display_intron_end+1)
    site_seq = intron_seq[0:2]
    swin_radius = win_radius/4
    #Get sliding window scores
    #_, norm_score, win_score, _ = window.slide_left(seqid, [flank_end], site_seq, fnh, swin_radius)

    site_st = site[0] - display_st - 1

    print "Region   ","%s:%d-%d"%(site[2],display_st,display_end)
    print "Site pos ","%s:%d-%d"%(site[2],site[0],site[1])
    print "Annotated","%s:%d-%d"%(annot_site[2],annot_site[0],annot_site[1])
    print "Exon     ",format_seq(exon_seq)
    print "Intron   ",format_seq(" "*eLen + intron_seq)
    print "Site     ",format_seq(" "*site_st+"**")
    #print "Normals  ",format_seq(" "*(site_st-swin_radius)+format_list(norm_score))
    #print "Slides   ",format_seq(" "*(site_st-swin_radius)+format_list(win_score))

#No flanking sequence required
def printRightSite(seqid,exon,site,annot_site,win_radius,display_st,display_end,fnh):
    """
    display_window                 (----------------------)
    exon                                      |--------------|
    intron        <---------------------------|
    ref_start     <--------------------------------------------------...---|
    """

    win_length = 2*win_radius
    #display_st,display_end = site[0] - win_radius, site[0] + win_radius
    display_intron_st  = display_st
    display_intron_end = exon.st0-1
    intron_seq = fnh.fetch_sequence(seqid,display_intron_st+1,display_intron_end+1)
    iLen = len(intron_seq)
    remaining = win_length - iLen
    display_exon_st  = 0
    display_exon_end = remaining
    exon_seq = exon.seq[display_exon_st:display_exon_end]

    site_st  = site[0] - display_st
    print "Region   ","%s:%d-%d"%(site[2],display_st,display_end)
    print "Site pos ","%s:%d-%d"%(site[2],site[0],site[1])
    print "Annotated","%s:%d-%d"%(annot_site[2],annot_site[0],annot_site[1])
    print "Exon    ",format_seq(" "*iLen +exon_seq)
    print "Intron  ",format_seq(intron_seq)
    print "Site    ",format_seq(" "*site_st+"**")


def printLeftSeq(seqid,flank,exon,site,annot_site,win_radius,display_st,display_end,fnh):
    """
    display_window                 (----------------------)
    flank                              |------|
    exon                       |--------------|
    intron                                    |--------------------------->
    ref_start |----------- ... ------------------------------------------->
    """
    printLeftSite(seqid,exon,site,annot_site,win_radius,display_st,display_end,fnh)
    flank_seq,flank_end = flank  #Note
    fLen = len(flank_seq)
    flank_st = (flank_end-fLen) - display_st - 1  #starting position of flanking sequence
    print "Flanks   ",format_seq(" "*flank_st + flank_seq),"\n"

def printRightSeq(seqid,flank,exon,site,annot_site,win_radius,display_st,display_end,fnh):
    """
    display_window                 (----------------------)
    flank                                     |------|
    exon                                      |--------------|
    intron        <---------------------------|
    ref_start     <--------------------------------------------------...---|
    """
    flank_seq,flank_st = flank  #Note
    printRightSite(seqid,exon,site,annot_site,win_radius,display_st,display_end,fnh)
    #Place flank_st and site_st in display coordinate frame
    flank_st = (flank_st) -  display_st
    print "Flanks  ",format_seq(" "*flank_st + flank_seq),"\n"

def printShortExon(display_st,display_end,xscript,site,fnh):
    """
    display_window                 (----------------------)
    exon                                 |---------|
    intron        <----------------------|         |-----------------...--->
    ref_start     <--------------------------------------------------...--->
    """
    short_exon=""
    #First find the exon
    for exon in xscript.exons:
        if display_st<exon.st0 and display_end>exon.en0:
                short_exon = exon
    assert short_exon!=""
    seqid = xscript.seqid
    left_in_start, left_in_end = display_st, short_exon.st0 - 1    #Left intron
    right_in_start, right_in_end = short_exon.en0, display_end     #Right intron
    display_exon_st = short_exon.st0-display_st
    display_exon_end = short_exon.en0-display_st
    left_in_seq = fnh.fetch_sequence(seqid,left_in_start+1,left_in_end+1)
    right_in_seq = fnh.fetch_sequence(seqid,right_in_start+1,right_in_end+1)
    exon_seq = exon.seq[display_exon_st:display_exon_end]
    i1Len,eLen,i2Len = len(left_in_seq), len(exon_seq), len(right_in_seq)
    display_i1_st, display_i2_st = 0, right_in_start-display_st
    eLen = len(exon_seq)
    site_st  = site[0] - display_st
    print "Region   ","%s:%d-%d"%(site[2],display_st,display_end)
    print "Site pos ","%s:%d-%d"%(site[2],site[0],site[1])
    print "Annotated","%s:%d-%d"%(site[2],site[0],site[1])
    print "Exon    ",format_seq(" "*display_exon_st +exon_seq)
    print "Intron  ",format_seq(left_in_seq+" "*eLen+right_in_seq)
    print "Site    ",format_seq(" "*site_st+"**")
    return


"""
Given a dictionary of flanking sequences and a dictionary of xscripts
For each false positive, display
1) Flanking sequences
2) Annotations
3) Sites
"""
def falsePositiveDisplay(flankDict,xscriptDict,site,annotDict,fnh):
    #print "annotations",annot_sites
    _,_,seqid,_ = site
    close = search.find_tuple(annotDict[seqid],site)
    #xscript_id = close[3]
    #x = xscript[xscript_id]
    annot_site = close
    key = (close[0],close[1],close[2])
    flanks = flankDict[key]
    xscript_id = close[3]
    xscript = xscriptDict[xscript_id]

    #Display everything around the splice site by a 50 bp radius
    win_radius = 50
    pos1,pos2 = site[0],site[1]  #positions of the estimated splice site
    #Get indexes of displayed exons.  Note that one of them should be -1
    display_st,display_end = site[0] - win_radius, site[0] + win_radius
    exon_li, exon_ri = xscript.getExon(display_st), xscript.getExon(display_end)
    print "flanks",flanks
    for flank in flanks:
        #flank_seq,flank_st = flank  #Note
        if exon_li!=-1:
            exon = xscript.exons[exon_li]
            printLeftSeq(xscript.seqid,flank,exon,site,annot_site,win_radius,display_st,display_end,fnh)
        if exon_ri!=-1:
            exon = xscript.exons[exon_ri]
            printRightSeq(xscript.seqid,flank,exon,site,annot_site,win_radius,display_st,display_end,fnh)


"""
Prints the flanking sequences, the annotated region and the site
"""
def falseNegativeDisplay(flankDict,xscriptDict,site,fnh):
    #print "annotations",annot_sites
    _,_,seqID,xscriptID = site
    key = (site[0],site[1],site[2])
    if key not in flankDict:
        flanks = []
    else:
        flanks = flankDict[key]
    xscript = xscriptDict[xscriptID]

    #Display everything around the splice site by a 50 bp radius
    win_radius = 50
    pos1,pos2 = site[0],site[1]  #positions of the estimated splice site
    #Get indexes of displayed exons.  Note that one of them should be -1 unless one of the exons are really short
    display_st,display_end = site[0] - win_radius, site[0] + win_radius
    exon_li, exon_ri = xscript.getExon(display_st), xscript.getExon(display_end)
    if exon_li==-1 and exon_ri==-1:  #most likely a really short exon
        printShortExon(display_st,display_end,xscript,site,fnh)
        return

    print "Exon indexs",exon_li,exon_ri
    if len(flanks)==0:
        if exon_li!=-1:
            exon = xscript.exons[exon_li]
            printLeftSite(xscript.seqid,exon,site,site,win_radius,display_st,display_end,fnh)
        if exon_ri!=-1:
            exon = xscript.exons[exon_ri]
            printRightSite(xscript.seqid,exon,site,site,win_radius,display_st,display_end,fnh)
    else:
        for flank in flanks:
            #flank_seq,flank_st = flank  #Note
            if exon_li!=-1:
                exon = xscript.exons[exon_li]
                printLeftSeq(xscript.seqid,flank,exon,site,site,win_radius,display_st,display_end,fnh)
            if exon_ri!=-1:
                exon = xscript.exons[exon_ri]
                printRightSeq(xscript.seqid,flank,exon,site,site,win_radius,display_st,display_end,fnh)

pattern = re.compile("(\S+):(\d+)-(\d+)") #parses chromosome name and positions

def incorrect(fps,fns,flankDict,xscriptDict,annotDict,region,fnh):

    if region!="":
        seqid,st,end = pattern.findall(region)[0]
        for fp in fps:
            if sseqid==seqid and spos1>=st and spos2<=end:
                falsePositiveDisplay(flankDict,xscriptDict,fp,annotDict,fnh)
        for fn in fns:
            if sseqid==seqid and spos1>=st and spos2<=end:
                falseNegativeDisplay(flankDict,xscriptDict,fn,fnh)
    else:
        for fp in fps:
            falsePositiveDisplay(flankDict,xscriptDict,fp,annotDict,fnh)
        for fn in fns:
            print "False negative",fn
            falseNegativeDisplay(flankDict,xscriptDict,fn,fnh)


_revcomp_trans = string.maketrans("ACGT", "TGCA")
def revcomp(s):
        return s[::-1].translate(_revcomp_trans)

def readableFormat(s):
    return "\t".join([s[i:i+10] + " " + str(i+10) for i in range(0,len(s),10)])

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

class TestDisplayFunctions(unittest.TestCase):
    def setUp(self):
        """
        AACTGTGATC 10AAGGATGTCT 20TCGCTTGTGA 30AACGAACGTC 40TGGATCCGCC 50TGAAGCATAT 60TGGCAATAAG 70ATCGCGCACT 80CCGTCCCCTA 90TATCAGACCA 100GTAAGGATCT 110GGAAACTCAT 120AAATGCCAGA 130GATAATAGCG 140TCGAATAGCG 150GCTCCTGTTG 160GTTGTCAGGC 170GCTACAAACG 180GCGGAAAACC 190GCATAGAAGT 200ATATACAAAA 210TGATTCCAGC 220GGCCCAAACG 230TCAATCTTTA 240AGAATATATA 250AGTATATTAA 260TTTTTAAGAA 270AGCTATTATT 280TACTATTACC 290TTTAGCCCAT 300ATCCGACTTC 310TAACAATATC 320TCTGGTGCTA 330CATATGTGGG 340AGTCCCGCAG 350ACAGCATATA 360GAAGATCGTT 370TACCTCACAT 380GCTAATCCAA 390AATCAGCAAG 400CTTTAGTTCA 410AGAACGTTTC 420CATGCTCATC 430CAATTTCACC 440TGCAAATAGT 450TCCAATGATT 460AGTCCTATTA 470TTGCTCAGTA 480TTTTTCATTA 490GTCATGTCAT 500ACTAGAAGAT 510TCTCAGGTTT 520TATATCTCGA 530TGCACAATGC 540CCATTGAATG 550CAGATAAGTC 560ATGGCCGCAC 570CCAAATGTCT 580AATCATAATG 590CGCGACTGGT 600TTTCTGAGAA 610CCTCGTTACC 620TGAGTTATAG 630CGTCGAATAA 640ATCACCACCT 650GTTACAAATT 660TAACAACAAA 670TTTAAGCAGA 680AAGAAGCTGT 690GAATGTATTC 700ACATTCAACT 710TACCACTTAC 720ATATTCCAGT 730ACAAGATACA 740TATTTGTATT 750TTGGTCTACA 760CTCAAAATAA 770GCGAAATTAT 780ATGCGGATGA 790TTTAATTTTT 800TCATGACGCG 810AACTTCCGCG 820TCTATGTAGT 830GTTCCTTGCC 840CTTGCATTTG 850TTCTTGTCTA 860TTATTTTTAG 870GGCGTAGGAA 880TGACCAGTTT 890GGCGGTGCTT 900AATCTTAAAC 910ACAATGGCAA 920AGTTGCCGTC 930GCCAATTATT 940CTTCCCAGAG 950AGTAAGTGTT 960ACGAATATTC 970GATGGCAATT 980CATTGATCTC 990CATACCAGTG 1000
        """
        refseq="""AACTGTGATCAAGGATGTCTTCGCTTGTGAAACGAACGTCTGGATCCGCCTGAAGCATATTGGCAATAAGATCGCGCACTCCGTCCCCTATATCAGACCAGTAAGGATCTGGAAACTCATAAATGCCAGAGATAATAGCGTCGAATAGCGGCTCCTGTTGGTTGTCAGGCGCTACAAACGGCGGAAAACCGCATAGAAGTATATACAAAATGATTCCAGCGGCCCAAACGTCAATCTTTAAGAATATATAAGTATATTAATTTTTAAGAAAGCTATTATTTACTATTACCTTTAGCCCATATCCGACTTCTAACAATATCTCTGGTGCTACATATGTGGGAGTCCCGCAGACAGCATATAGAAGATCGTTTACCTCACATGCTAATCCAAAATCAGCAAGCTTTAGTTCAAGAACGTTTCCATGCTCATCCAATTTCACCTGCAAATAGTTCCAATGATTAGTCCTATTATTGCTCAGTATTTTTCATTAGTCATGTCATACTAGAAGATTCTCAGGTTTTATATCTCGATGCACAATGCCCATTGAATGCAGATAAGTCATGGCCGCACCCAAATGTCTAATCATAATGCGCGACTGGTTTTCTGAGAACCTCGTTACCTGAGTTATAGCGTCGAATAAATCACCACCTGTTACAAATTTAACAACAAATTTAAGCAGAAAGAAGCTGTGAATGTATTCACATTCAACTTACCACTTACATATTCCAGTACAAGATACATATTTGTATTTTGGTCTACACTCAAAATAAGCGAAATTATATGCGGATGATTTAATTTTTTCATGACGCGAACTTCCGCGTCTATGTAGTGTTCCTTGCCCTTGCATTTGTTCTTGTCTATTATTTTTAGGGCGTAGGAATGACCAGTTTGGCGGTGCTTAATCTTAAACACAATGGCAAAGTTGCCGTCGCCAATTATTCTTCCCAGAGAGTAAGTGTTACGAATATTCGATGGCAATTCATTGATCTCCATACCAGTG"""
        annots="""chr2R\tunknown\texon\t290\t439\t.\t-\t.\tgene_id "CG17528"; gene_name "CG17528"; p_id "P21588"; transcript_id "NM_001042999"; tss_id "TSS13109";\nchr2R\tunknown\texon\t503\t648\t.\t-\t.\tgene_id "CG17528"; gene_name "CG17528"; p_id "P21588"; transcript_id "NM_001042999"; tss_id "TSS13109";\n"""
        self.fasta = "test.fa"
        self.faidx = "test.fa.fai"
        self.gtf   = "test.gtf"
        createTestFasta(self.fasta,"chr2R",refseq)
        createTestGTF(self.gtf,annots)
        #print >> sys.stderr,"refseq\n",readableFormat(refseq)

    def testPrint1(self):
        print "Test print 1"
        #site = (440,441,'chr2R',"NM_001042999")
        site = (439,440,'chr2R',"NM_001042999")
        # print >> sys.stderr,'fasta file',self.fasta
        # print >> sys.stderr,'gtf file  ',self.gtf
        annots = gtf.parseGTF([self.gtf])
        fastadb = gtf.parseFASTA([self.fasta])
        xscripts = gtf.assembleTranscripts(annots,fastadb)
        fnh = fasta.fasta(self.fasta)

        # print >> sys.stderr,"Fasta db      ",fastadb
        # print >> sys.stderr,"Annotations   \n",annots[0],'\n',annots[1]
        # print >> sys.stderr,"Transcripts",str(xscripts)
        #print "annotations",xscripts

        flanks = [("CCAATTTCAC",439),("CCAATTTCAC",439)]
        key = (site[0],site[1],site[2])

        flankDict = { key: flanks  }
        annotDict = {x.seqid: x.getSites() for x in xscripts}
        xscriptDict = {x.xscript_id: x for x in xscripts}
        falsePositiveDisplay(flankDict,xscriptDict,site,annotDict,fnh)
        print


    def testPrint2(self):
        print "Test print 2"
        site = (500,501,'chr2R',"NM_001042999")
        flanks = [("TAGAAGATTC",502),("TAGAAGATTC",502)]
        annots = gtf.parseGTF([self.gtf])
        fastadb = gtf.parseFASTA([self.fasta])
        xscripts = gtf.assembleTranscripts(annots,fastadb)
        fnh = fasta.fasta(self.fasta)
        annot_sites = {x.seqid: x.getSites() for x in xscripts}

        key = (site[0],site[1],site[2])

        flankDict = { key: flanks  }
        annotDict = {x.seqid: x.getSites() for x in xscripts}
        xscriptDict = {x.xscript_id: x for x in xscripts}
        falsePositiveDisplay(flankDict,xscriptDict,site,annotDict,fnh)
        print

    def testFalseNegative1(self):
        print "Test false negative 1"
        site = (500,501,'chr2R',"NM_001042999")
        flanks = []
        annots = gtf.parseGTF([self.gtf])
        fastadb = gtf.parseFASTA([self.fasta])
        xscripts = gtf.assembleTranscripts(annots,fastadb)
        fnh = fasta.fasta(self.fasta)
        flankDict = {}
        xscriptDict = {x.xscript_id: x for x in xscripts}
        falseNegativeDisplay(flankDict,xscriptDict,site,fnh)
        print
    def testShortExon(self):
        site = (500,501,'chr2R',"NM_001042999")
        #exon bounds (503,648)
        display_st,display_end = 495,650
        annots = gtf.parseGTF([self.gtf])
        fastadb = gtf.parseFASTA([self.fasta])
        xscripts = gtf.assembleTranscripts(annots,fastadb)
        fnh = fasta.fasta(self.fasta)
        xscript = xscripts[0]
        printShortExon(display_st,display_end,xscript,site,fnh)

    ###Note:  Need to test for the case where exon.st0 > display window.  aka, exon is way too short

    def tearDown(self):
        os.remove(self.fasta)
        os.remove(self.faidx)
        os.remove(self.gtf)

if __name__=="__main__":
    unittest.main()
