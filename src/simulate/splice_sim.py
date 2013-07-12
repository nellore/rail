"""
annotate_sim.py

Simulates differiential gene expression using annotated genes
"""

import os
import site
import argparse
import sys
import random
import math
import re
import bisect
from operator import itemgetter
from collections import defaultdict

base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
site.addsitedir(os.path.join(base_path, "annotation"))

import gtf

parser = argparse.ArgumentParser(description=\
                                     'Transcript simulator')
parser.add_argument(\
    '--output-prefix', metavar='path', type=str, required=False,
    help='Prefix for output read files')
parser.add_argument(\
    '--read-len', metavar='int', action='store', type=int, default=100,
    help='Read length to simulate')
parser.add_argument(\
    '--num-replicates', metavar='int', action='store', type=int, default=8,
    help='Number of replicates per group')
parser.add_argument(\
    '--seed', metavar='int', action='store', type=int, default=874,
    help='Pseudo-random seed')
parser.add_argument(\
    '--num-nucs', metavar='int', action='store', type=float, default=1e7,
    help='Number of total nucleotides of reads to generate')
parser.add_argument(\
    '--num-xscripts', metavar='int', action='store', type=float, default=1e7,
    help='Number of transcripts to simulate')

gtf.addArgs(parser)
args = parser.parse_args()


class WeightedRandomGenerator(object):

    def __init__(self, weights):
        self.totals = []
        running_total = 0
        for w in weights:
            running_total += w
            self.totals.append(running_total)

    def next(self):
        rnd = random.random() * self.totals[-1]
        return bisect.bisect_right(self.totals, rnd)

    def __call__(self):
        return self.next()

#Make weighted random generator for transcriptome
def makeWeights(xscripts):
    if args.num_xscripts==0:
        num_xscripts= len(xscripts)
    else:
        num_xscripts = args.num_xscripts
    weights = [1.0] * len(xscripts)

    num = 0
    for i in range(0,len(weights)):
        if num<num_xscripts:
            weights[i] = random.random()
        else:
            weights[i] = 0
        num+=1
    return WeightedRandomGenerator(weights),weights

def simulate(xscripts,readlen,targetNucs,fastaseqs):
    gen,weights = makeWeights(xscripts)
    n = 0
    seqs = []
    while n<targetNucs:
        #Pick a transcript at weighted random
        i = gen.next()
        x = xscripts[i]
        start,end,seqid = 0,len(x.seq)-readlen,x.seqid
        if end<start or len(x.seq)<readlen:
            continue
        i = random.randint(start,end)
        read = x.seq[i:i+readlen]
        seqs.append(read)
        n+=readlen
    return seqs,weights

#Just prints out one file
def writeReads(seqs, fnPre,manifestFn):
    """ Only unpaired for now """
    fn = "%s.seqs.tab6"%(fnPre)
    with open(manifestFn,'w') as manFh:
        manFh.write("%s\t0\tsplice-0-0\n"%(fn))
    with open(fn,'w') as fh:
        for i in xrange(0,len(seqs)):
            seq = seqs[i]
            if len(seq)==0:
                continue
            qual = "I" * len(seq)
            nm = "r_n%d;LB:splice-0-0" % (i)
            fh.write("%s\t%s\t%s\n" % (nm, seq, qual))
        
def writeWeights(weights,xscripts,fout):
    with open(fout,'w') as fh:
        for i in range(0,len(weights)):
            fh.write("weight:%s\n%s\n"%(weights[i],str(xscripts[i])))
    

if __name__=="__main__":
    annots = gtf.parseGTF(args.gtf)
    fastadb = gtf.parseFASTA(args.fasta)
    xscripts = gtf.assembleTranscripts(annots,fastadb)
    seqs,weights = simulate(xscripts,args.read_len,args.num_nucs,args.fasta)
    writeReads(seqs,args.output_prefix,args.output_prefix+".manifest")
    writeWeights(weights,xscripts,args.output_prefix+".weights")
    
