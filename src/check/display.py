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
#def formatSeq():

def printLeftSeq(seqid,flank,exon,site,win_radius,display_st,display_end,fnh):
    """
    display_window                 (----------------------)
    flank                              |------|
    exon                       |--------------|
    intron                                    |--------------------------->
    ref_start |----------- ... ------------------------------------------->
    """
    #convert to display window coordinate frame for exon information
    flank_seq,flank_end = flank  #Note
    win_length = 2*win_radius
    display_st,display_end = site[0] - win_radius, site[0] + win_radius
    display_exon_st   = display_st - exon.st0
    display_exon_end  = exon.en0 - exon.st0 -1
    exon_seq = exon.seq[display_exon_st:display_exon_end]
    #print >> sys.stderr,"Whole exon",exon.seq
    eLen,fLen = len(exon_seq), len(flank_seq)
    #Stay in reference coordinates for intron
    remaining = win_length - eLen
    display_intron_st = display_st+eLen+1
    display_intron_end= display_intron_st + remaining
    intron_seq = fnh.fetch_sequence(seqid,display_intron_st+1,display_intron_end+1)

    flank_st = (flank_end-fLen) - display_st - 1  #starting position of flanking sequence
    site_st = site[0] - display_st - 1

    print "Region: ",display_st,'-',display_end
    print "Site    ",site
    print "Flanks  "," "*flank_st + flank_seq
    print "Exon    ",exon_seq
    print "Intron  "," "*eLen + intron_seq
    print "Site    "," "*site_st+"**\n"

def printRightSeq(seqid,flank,exon,site,win_radius,display_st,display_end,fnh):
    """
    display_window                 (----------------------)
    flank                                     |------|
    exon                                      |--------------|
    intron        <---------------------------|
    ref_start     <--------------------------------------------------...---|
    """
    flank_seq,flank_st = flank  #Note

    win_length = 2*win_radius
    #display_st,display_end = site[0] - win_radius, site[0] + win_radius
    display_intron_st  = display_st
    display_intron_end = exon.st0
    intron_seq = fnh.fetch_sequence(seqid,display_intron_st+1,display_intron_end+1)
    iLen = len(intron_seq)
    remaining = win_length - iLen
    display_exon_st  = 0
    display_exon_end = remaining
    exon_seq = exon.seq[display_exon_st:display_exon_end]
    #Before, flank_st actually pointed to the end of the splice_site
    flank_st = (flank_st) -  display_st
    site_st  = site[0] - display_st

    print "Region: ",display_st,'-',display_end
    print "Site    ",site
    print "Flanks  "," "*flank_st + flank_seq
    print "Exon    "," "*iLen +exon_seq
    print "Intron  ",intron_seq
    print "Site    "," "*site_st+"**"

"""
Prints the flanking sequences, the annotated region and the site
"""
def printSeqs(flanks,xscript,site,fnh):
    #Display everything around the splice site by a 50 bp radius
    win_radius = 50
    pos1,pos2 = site[0],site[1]  #positions of the estimated splice site
    #Get indexes of displayed exons.  Note that one of them should be -1
    display_st,display_end = site[0] - win_radius, site[0] + win_radius
    exon_li, exon_ri = xscript.getExon(display_st), xscript.getExon(display_end)
    #print >> sys.stderr,"xscript seq",xscript.seq
    #print >> sys.stderr,"exon 1",xscript.exons[0].seq

    #flank_seq = flanks[0]

    for flank in flanks:
        #flank_seq,flank_st = flank  #Note
        if exon_li!=-1:
            exon = xscript.exons[exon_li]
            printLeftSeq(xscript.seqid,flank,exon,site,win_radius,display_st,display_end,fnh)
        if exon_ri!=-1:
            exon = xscript.exons[exon_ri]
            printRightSeq(xscript.seqid,flank,exon,site,win_radius,display_st,display_end,fnh)


pattern = re.compile("(\S+):(\d+)-(\d+)") #parses chromosome name and positions

"""
Given a dictionary of flanking sequences and a dictionary of xscripts
For each false positive, display
1) Flanking sequences
2) Annotations
3) Sites
"""

def falsePositives(fp,flanks,xscript,annot_sites,fnh):
    _,_,seqid,_ = fp
    close = search.find_tuple(annot_sites[seqid],fp)
    #xscript_id = close[3]
    #x = xscript[xscript_id]
    key = (close[0],close[1],close[2])
    flank_seqs = flanks[key]
    printSeqs(flank_seqs,xscript,fp,fnh)

#TODO: Add false negative printing
def incorrect(fps,fns,flanks,xscripts,annot_sites,region,fnh):
    print fps
    if region!="":
        seqid,st,end = pattern.findall(region)[0]
        for s in fps:
            spos1,spos2,sseqid,_ = s
            if sseqid==seqid and spos1>=st and spos2<=end:
                falsePositives(s,flanks,xscripts,annot_sites)
    else:
        for s in fps:
            _,_,sseqid,_ = s
            falsePositives(s,flanks,xscripts[sseqid],annot_sites,fnh)


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
        refseq="""CAACTGTGATCAAGGATGTCTTCGCTTGTGAAACGAACGTCTGGATCCGCCTGAAGCATATTGGCAATAAGATCGCGCACTCCGTCCCCTATATCAGACCAGTAAGGATCTGGAAACTCATAAATGCCAGAGATAATAGCGTCGAATAGCGGCTCCTGTTGGTTGTCAGGCGCTACAAACGGCGGAAAACCGCATAGAAGTATATACAAAATGATTCCAGCGGCCCAAACGTCAATCTTTAAGAATATATAAGTATATTAATTTTTAAGAAAGCTATTATTTACTATTACCTTTAGCCCATATCCGACTTCTAACAATATCTCTGGTGCTACATATGTGGGAGTCCCGCAGACAGCATATAGAAGATCGTTTACCTCACATGCTAATCCAAAATCAGCAAGCTTTAGTTCAAGAACGTTTCCATGCTCATCCAATTTCACCTGCAAATAGTTCCAATGATTAGTCCTATTATTGCTCAGTATTTTTCATTAGTCATGTCATACTAGAAGATTCTCAGGTTTTATATCTCGATGCACAATGCCCATTGAATGCAGATAAGTCATGGCCGCACCCAAATGTCTAATCATAATGCGCGACTGGTTTTCTGAGAACCTCGTTACCTGAGTTATAGCGTCGAATAAATCACCACCTGTTACAAATTTAACAACAAATTTAAGCAGAAAGAAGCTGTGAATGTATTCACATTCAACTTACCACTTACATATTCCAGTACAAGATACATATTTGTATTTTGGTCTACACTCAAAATAAGCGAAATTATATGCGGATGATTTAATTTTTTCATGACGCGAACTTCCGCGTCTATGTAGTGTTCCTTGCCCTTGCATTTGTTCTTGTCTATTATTTTTAGGGCGTAGGAATGACCAGTTTGGCGGTGCTTAATCTTAAACACAATGGCAAAGTTGCCGTCGCCAATTATTCTTCCCAGAGAGTAAGTGTTACGAATATTCGATGGCAATTCATTGATCTCCATACCAGTG"""
        annots="""chr2R\tunknown\texon\t290\t439\t.\t-\t.\tgene_id "CG17528"; gene_name "CG17528"; p_id "P21588"; transcript_id "NM_001042999"; tss_id "TSS13109";\nchr2R\tunknown\texon\t503\t648\t.\t-\t.\tgene_id "CG17528"; gene_name "CG17528"; p_id "P21588"; transcript_id "NM_001042999"; tss_id "TSS13109";\n"""
        self.fasta = "test.fa"
        self.faidx = "test.fa.fai"
        self.gtf   = "test.gtf"
        createTestFasta(self.fasta,"chr2R",refseq)
        createTestGTF(self.gtf,annots)
        #print >> sys.stderr,"refseq\n",readableFormat(refseq)

    def testPrint1(self):
        site = (440,441,'chr2R',"NM_001042999")
        flanks = [("CCAATTTCAC",440),("CCAATTTCAC",440)]
        # print >> sys.stderr,'fasta file',self.fasta
        # print >> sys.stderr,'gtf file  ',self.gtf
        annots = gtf.parseGTF([self.gtf])
        fastadb = gtf.parseFASTA([self.fasta])
        xscripts = gtf.assembleTranscripts(annots,fastadb)
        fnh = fasta.fasta(self.fasta)

        # print >> sys.stderr,"Fasta db      ",fastadb
        # print >> sys.stderr,"Annotations   \n",annots[0],'\n',annots[1]
        # print >> sys.stderr,"Transcripts",str(xscripts)

        printSeqs(flanks,xscripts[0],site,fnh)


    def testPrint2(self):
        site = (501,502,'chr2R',"NM_001042999")

        flanks = [("TAGAAGATTC",503),("TAGAAGATTC",503)]
        annots = gtf.parseGTF([self.gtf])
        fastadb = gtf.parseFASTA([self.fasta])
        xscripts = gtf.assembleTranscripts(annots,fastadb)
        fnh = fasta.fasta(self.fasta)
        printSeqs(flanks,xscripts[0],site,fnh)


    def tearDown(self):
        os.remove(self.fasta)
        os.remove(self.faidx)
        os.remove(self.gtf)

if __name__=="__main__":
    unittest.main()