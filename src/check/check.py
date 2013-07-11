"""
check.py

Simple checking script that checks splice output


"""

import pickle
import os
import site
import re
import sys
import argparse
from operator import itemgetter
from collections import defaultdict

base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
print base_path
site.addsitedir(os.path.join(base_path, "annotation"))

import gtf
parser = argparse.ArgumentParser(description=\
                                     'Splice simulator tester')

parser.add_argument(\
    '--align-out', metavar='path', type=str, required=False,
    help='File containing spliced alignments')

parser.add_argument(\
    '--output-prefix', metavar='path', type=str, required=False,
    help='Prefix for output read files')

parser.add_argument(\
    '--weight-file', metavar='path', type=str, required=False,
    help='Contains the probability weights in a pickle file')

parser.add_argument(\
    '--xscripts-file', metavar='path', type=str, required=False,
    help='Contains transcript information in a pickle file')

parser.add_argument(\
    '--refseq', metavar='path', type=str, required=False,
    help='Contains the reference sequence')

parser.add_argument(\
    '--site-file', metavar='path', type=str, required=False,
    help='Contains the splice junctions')

args = parser.parse_args()


def parseSites(fn):
    sites = []
    with open(fn,'r') as fh:
        for ln in fh:
            toks = ln.rstrip().split("\t")
            seq,st,en,lab,freq = toks[0],int(toks[1]),int(toks[2]),toks[3],int(toks[4])
            sites.append((seq,st,en,lab,freq))
    return sites
    
"""
Create a genome dictionary indexed by position to visualize the genome

First, create a genome wide visualization of all splice sites
start of transcript: \
end of transcript:   /
exon start: ^
exon end: v
splice sites: < and >

exon
<
>
intron
>
<
exon
"""
def genomeSiteDict(xscripts,sites):
    genome = dict()
    for x in xscripts:
        xstart = "012%d"%(x.st0)
        xend = "012%d"%(x.en0)
        for E in x.exons:
            start = "012%d"%(E.st0)
            end = "012%d"%(E.en0)
            genome[start] = "v"
            genome[end] = "^"
            
        genome[xstart] = "\\" #Note, this will overwrite any changes
        genome[xend] = "/"
    for s in sites:
        seq,st,en,lab,freq = s
        start = "%012d"%(st)
        start1 = "%012d"%(st+1)
        end = "%012d"%(en)
        end1 = "%012d"%(en+1)
        if start in genome:
            genome[start] = "|"
        else:
            genome[start] = "<"
        if start1 in genome:
            genome[start1] = "|"
        else:
            genome[start1] = ">"
        if end in genome:
            genome[end] = "|"
        else:
            genome[end] = ">"
        if end1 in genome:
            genome[end1] = "|"
        else:
            genome[end1] = "<"
    return genome

"""
Levels out symbols
"""
def preappendSym(genome,st,en):
    cst,cen = 0,0
    if st in genome:
        cst = len(genome[st])
    if en in genome:
        cen = len(genome[en])
    
    if cst == cen:
        return genome
    elif cst > cen:
        while cst > cen:
            genome[en].append("|")
            cst = len(genome[st])
            cen = len(genome[en])
    else:
        while cst < cen:
            genome[st].append("|")
            cst = len(genome[st])
            cen = len(genome[en])
    return genome
    

"""
exons = ^ ... v
introns = u ... n
"""
def genomeContigDict(contigs):
    genome = defaultdict(list)
    for c in contigs:
        ctype,pt,st,en,refid,lab = c
        start = "%012d"%(st)
        end = "%012d"%(en)
        genome = preappendSym(genome,start,end)
        if ctype == "exon":
            genome[st].append("^")
            genome[en].append("v")
        else:
            genome[st].append("n")
            genome[en].append("u")
    return genome

def compareSites():
    xscripts = pickle.load(open(args.xscripts_file,'rb'))
    #weights = pickle.load(open(args.weight_file,'rb'))
    sites = parseSites(args.site_file)
    genome = genomeSiteDict(xscripts,sites)
    gsyms = [(int(pos),sym) for pos,sym in genome.iteritems()]
    gsyms.sort(key=itemgetter(0))
    for g in gsyms:
        pos,sym = g
        print "%012d"%pos,sym

def parseAlignments(fn):
    contigs = []
    with open(fn,'r') as fh:
        for ln in fh:
            toks = ln.split("\t")
            assert len(toks)==6
            ctype,pt,st,en,refid,lab = toks[0],toks[1],int(toks[2]),int(toks[3]),toks[4],toks[5]
            contigs.append((ctype,pt,st,en,refid,lab))
    return contigs

def compareAlignments(fn):
    contigs = parseAlignments(fn)
    genome = genomeContigDict(contigs)
    gsyms = [(int(pos),sym) for pos,sym in genome.iteritems()]
    gsyms.sort(key=itemgetter(0))
    for g in gsyms:
        pos,syms = g
        print "%012d"%pos,"".join(syms)

if __name__=="__main__":    
    print args.align_out
    compareAlignments(args.align_out)
    
