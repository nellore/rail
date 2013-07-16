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
import pickle
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
    '--num-replicates', metavar='int', action='store', type=int, default=2,
    help='Number of replicates per group')
parser.add_argument(\
    '--seed', metavar='int', action='store', type=int, default=874,
    help='Pseudo-random seed')
parser.add_argument(\
    '--num-nucs', metavar='int', action='store', type=float, default=1e7,
    help='Number of total nucleotides of reads to generate')
parser.add_argument(\
    '--num-xscripts', metavar='int', action='store', type=float, default=2,
    help='Number of transcripts to simulate')
parser.add_argument(\
    '--stranded', metavar='int', action='store', type=int, default=False,
    help='Indicates if boths strands need to be simulated')

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
            num+=1
        else:
            weights[i] = 0
    return WeightedRandomGenerator(weights),weights
"""
Randomly gets two transcripts: one from each strand
"""
def makeStrandedWeights(xscripts):
    num_xscripts = args.num_xscripts
    weights = [1.0] * len(xscripts)
    last_strand = "+"
    num = 0
    for i in range(0,len(weights)):
        if num<num_xscripts and xscripts[i].orient!=last_strand:
            weights[i] = random.random()
            print xscripts[i].xscript_id
            last_strand = xscripts[i].orient
            num+=1
        else:
            weights[i] = 0
    return WeightedRandomGenerator(weights),weights

def simulate(xscripts,readlen,targetNucs,fastaseqs):
    if args.stranded:
        gen,weights = makeStrandedWeights(xscripts)
    else:
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

def replicateize(seqs1, seqs2, nreps):
    """ Take all the sequence reads for groups 1 and 2 and split them into into
        a bunch of replicates with varying average coverage. """
    weights1 = [1.0] * nreps
    weights2 = [1.0] * nreps
    for i in xrange(0, nreps):
        weights1[i] += (random.random() - 0.5)
        weights2[i] += (random.random() - 0.5)
    # Each replicate has a weight, so now we assign sequences to replicates
    seqs1rep = [ [] for _ in xrange(0, nreps) ]
    seqs2rep = [ [] for _ in xrange(0, nreps) ]
    gen1 = WeightedRandomGenerator(weights1)
    gen2 = WeightedRandomGenerator(weights2)
    for seq in seqs1:
        i = gen1.next()
        seqs1rep[i].append(seq)
    for seq in seqs2:
        i = gen2.next()
        seqs2rep[i].append(seq)
    return seqs1rep, seqs2rep


#Just prints out one file
def writeReads(seqs1rep,seqs2rep,fnPre,manifestFn):
    """ Only unpaired for now """
    seqsrep = [ seqs1rep, seqs2rep ]
    fn = "%s.seqs.tab6"%(fnPre)
    with open(manifestFn, 'w') as manFh:
        for group in xrange(0, 2):
            sr = seqsrep[group]
            for rep in xrange(0, len(sr)):
                fn = "%s.group%d.rep%d.tab6" % (fnPre, group, rep)
                manFh.write("%s\t0\tsplice-%d-%d\n" % (fn, group, rep))
                with open(fn, 'w') as fh:
                    for i in xrange(0, len(sr[rep])):
                        seq = sr[rep][i]
                        qual = "I" * len(seq)
                        nm = "r_n%d;LB:splice-%d-%d" % (i, group, rep)
                        fh.write("%s\t%s\t%s\n" % (nm, seq, qual))

    
if __name__=="__main__":
    annots = gtf.parseGTF(args.gtf)
    fastadb = gtf.parseFASTA(args.fasta)
    xscripts = gtf.assembleTranscripts(annots,fastadb)
    seqs,weights = simulate(xscripts,args.read_len,args.num_nucs,args.fasta)
    seqs1,seqs2 = replicateize(seqs,seqs,args.num_replicates)
    writeReads(seqs1,seqs2,args.output_prefix,args.output_prefix+".manifest")
    #pickle.dump(weights,open(args.output_prefix+".weights",'wb'))
    #pickle.dump(xscripts,open(args.output_prefix+".xscripts",'wb'))
    
