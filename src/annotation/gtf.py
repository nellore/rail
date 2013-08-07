"""
gtf.py

Parses gtf files
"""

import string
import argparse
import re
from operator import itemgetter
from collections import defaultdict
import random

#parser = argparse.ArgumentParser(description='Parse a GTF file.')
def addArgs(parser):
    parser.add_argument(\
        '--fasta', metavar='path', type=str, nargs='+',
        help='FASTA file(s) containing reference genome sequences')
    parser.add_argument(\
        '--gtf', metavar='path', type=str, nargs='+', required=True,
        help='GTF file(s) containing gene annotations')
    parser.add_argument(\
        '--introns-out', metavar='PATH', type=str, required=False,
        help='File to write intron information to')
    parser.add_argument(\
        '--exons-out', metavar='PATH', type=str, required=False,
        help='File to write exon information to')

    #args = parser.parse_args()

gene_id_re = re.compile("gene_id \"([^\"]+)\"")
xscript_id_re = re.compile("transcript_id \"([^\"]+)\"")
"""
Stores objects such as exons and introns
"""
class Annot(object):
    def __init__(self, refid, st0, en0, orient, feature, score, frame, attrs):
        self.refid = refid
        self.st0 = st0 # 0-based
        self.en0 = en0 # 0-based, exclusive
        self.orient = orient
        self.feature = feature
        self.score = score
        self.frame = frame
        self.attrs = attrs

    def __len__(self):
        return self.en0-self.st0

    def __str__(self):
        return "%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s"%(self.refid, self.st0, self.en0, self.orient, self.feature, self.score, self.frame,self.attrs)
"""
Stores a chain of annotation objects
"""
class Transcript(object):
    """ Right now I ignore all associations with TSSs and promoters """
    def __init__(self, exons, start=None, stop=None, tssId=None):
        self.exons = exons # list of exon Annots
        for exon in exons:
            assert exon.feature == "exon"
        self.start = start # start coden Annot
        self.stop = stop
        self.tssId = tssId 
        self.exons.sort(key = lambda x: x.st0)
        self.orient=self.exons[0].orient
        self.st0 = self.exons[0].st0
        self.en0 = self.exons[-1].en0
        self.seqid = self.exons[0].refid
        self.gene_id,self.xscript_id = self.exons[0].attrs.split('\t')
       
    #Pass in a dictionary of fastaseqs
    def buildSeq(self,fastaseqs):
        seqs = []
        for E in self.exons:
            seqid = E.refid
            start = E.st0
            end = E.en0
            if seqid not in fastaseqs:
                continue

            seqs.append(fastaseqs[seqid][start:end].upper())
        self.seq = "".join(seqs)

    """
    Assumption: All splice sites are 2bp long
    returns a list of start positions of splice sites in the transcripts
    """
    def getSites(self):
        sites = []
        for i in range(1,len(self.exons)):
            site5 = self.exons[i-1].en0
            site3 = self.exons[i].st0
            sites.append(site5)
            sites.append(site5+1)
            sites.append(site3-1)
            sites.append(site3-2)
        return sites
    """
    Get all splice sites in terms of transcript sequence
    """
    def getXcriptSites(self):
        sites = []
        total = 0
        for i in range(1,len(self.exons)):
            total+=len(self.exons[i-1])
            sites.append(total)
            sites.append(total+1)
        return sites

    def incorporateVariants(self,mm_rate,indel_rate,var_handle):
        lseq = list(self.seq)
        for i in range(0,len(lseq)):
            r = random.random()
            if r<mm_rate:
                c = "ACGT"[random.randint(0,3)]
                while c==lseq[i]: #Can't be equal
                    lseq[i] = c
                    c = "ACGT"[random.randint(0,3)]
                var_handle.write("%s\tsnp%d"%(self.gene_id,i))
            else:
                r = random.random()
                if r<indel_rate:
                    indel = "ID"[random.randint(0,1)]
                    if indel=="I":
                        ins = "ACGT"[random.randint(0,3)]
                        seq1 = lseq[:i]
                        seq2 = lseq[i:]
                        lseq = seq1+[ins]+seq2
                        var_handle.write("%s\tinsert%d"%(self.gene_id,i))
                    else:
                        seq1 = lseq[:i]
                        seq2 = lseq[i+1:]
                        lseq = seq1+seq2
                        var_handle.write("%s\tdelete%d"%(self.gene_id,i))
        self.seq = "".join(lseq)

    def __str__(self):
        lns = []
        lns.append("%d\t%s\t%d\n"%(self.st0,self.seqid,self.en0))
        for e in self.exons:
            lns.append("%s\n"%(str(e))) 
        lns.append("%s\n"%(str(self.seq)))
        return "".join(lns)
"""
Constructs a dictionary of fasta seqs
"""
def parseFASTA(fns):
    seqs = dict()
    lns = []
    seqid = ""
    numseqs = 0
    for fn in fns:
        with open(fn, 'r') as fh:
            for ln in fh:
                if ln[0] == '>':
                    if seqid!="":
                        seqs[seqid] = ''.join(lns)
                    seqid = ln.split(" ")[0][1:].rstrip()
                    lns = []
                    numseqs+=1
                    continue
                else:
                    lns.append(ln.rstrip())
    return seqs


"""
Returns a list of exons
"""
def parseGTF(fns):
    exons = list()
    for gtf in fns:
        with open(gtf, 'r') as gtfFh:
            last_gene_id = None
            for line in gtfFh:
                gene_id, xscript_id = None, None
                refid, source, feature, start1, end1, score, orient, frame, attr = string.split(line, '\t')
                if feature == "CDS":
                    continue
                if feature == "start_codon":
                    continue
                if feature == "stop_codon":
                    continue
                assert feature == "exon"
                start1, end1 = int(start1), int(end1)
                assert end1 >= start1
                ma = gene_id_re.search(attr)
                if ma is None:
                    raise RuntimeError("No gene_id!:\n" + line)
                gene_id = ma.group(1)
                ma = xscript_id_re.search(attr)
                if ma is None:
                    raise RuntimeError("No transcript_id!:\n" + line)
                xscript_id = ma.group(1)
                if gene_id != last_gene_id:
                    last_gene_id = gene_id
                attrs = "%s\t%s"%(gene_id,xscript_id)
                annot = Annot(refid, start1-1, end1,orient,feature,score,frame,attrs)
                exons.append(annot)
    return exons

"""                                                                                                                        
Input: sorted exons and a dictionary of fasta seqs                                                                         
Chain all exons into transcript via transcript id                                                                           
"""
def assembleTranscripts(exons,fastaseqs):
    exonlist = defaultdict(list)
    #Bin all exons by transcript id
    for E in exons:
        gene_id,xscript_id = E.attrs.split("\t")
        exonlist[xscript_id].append(E)
    xscripts = list()
    for tr_id,exons in exonlist.iteritems():
        x = Transcript(exons)
        x.buildSeq(fastaseqs)
        xscripts.append(x)
    return xscripts
