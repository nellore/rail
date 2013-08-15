"""
sim_splice.py

Simulates differiential gene expression using annotated genes
"""
import string
import os
import site
import argparse
import sys
import random
import math
import re
import bisect
import pickle
from operator import itemgetter
from collections import defaultdict

base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
site.addsitedir(os.path.join(base_path, "annotation"))
site.addsitedir(os.path.join(base_path, "fasta"))
site.addsitedir(os.path.join(base_path, "util"))

import gtf
import counter

parser = argparse.ArgumentParser(description=\
'Counts all of the splice sites in a transcriptome')

gtf.addArgs(parser)
args = parser.parse_args()

def getSiteSeq(site,fastaseq):
    #print >> sys.stderr,"fastaseq",fastaseq[  site[0]:site[1]+2 ].upper()
    return fastaseq[ site[0]:site[1]+1 ].upper()

def getSiteDistribution(xscripts,fastadb):

    cnts = counter.Counter()
    for x in xscripts:
        site_pairs = x.getSitePairs()
        for sp in site_pairs:
            left_site,right_site = sp
            if x.seqid in fastadb:
                fastaseq = fastadb[x.seqid]
                leftseq=getSiteSeq(left_site,fastaseq)
                rightseq=getSiteSeq(right_site,fastaseq)
                cnts[(leftseq,rightseq)]+=1
    print "\n".join(["%s-%s\t%d"%(k[0],k[1],v) for k,v in cnts.items()])

if __name__=="__main__":
    annots = gtf.parseGTF(args.gtf)
    fastadb = gtf.parseFASTA(args.fasta)
    xscripts = gtf.assembleTranscripts(annots,fastadb)
    getSiteDistribution(xscripts,fastadb)

