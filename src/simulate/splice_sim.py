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

import gtf
import chrsizes

parser = argparse.ArgumentParser(description=\
                                     'Transcript simulator')
parser.add_argument(\
    '--output-prefix', metavar='path', type=str, required=False,
    help='Prefix for output read files')
parser.add_argument(\
    '--output-annotations', metavar='path', type=str, required=False,
    help='Prefix for simulated transcripts')
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
parser.add_argument(\
    '--alternative-spliced', metavar='int', action='store', type=int, default=False,
    help='Indicates if alternatively spliced transcripts should be simulated')
parser.add_argument(\
    '--paired-end', metavar='int', action='store', type=int, default=False,
    help='Indicates if the reads are paired end')
parser.add_argument(\
    '--readmm_rate', metavar='float', action='store', type=float, default=0.01,
    help='The rate of mismatches')
parser.add_argument(\
    '--snp_rate', metavar='float', action='store', type=float, default=0.001,
    help='The rate of snps')
parser.add_argument(\
    '--indel_rate', metavar='float', action='store', type=float, default=0.0002,
    help='The rate of inserts or deletions')
parser.add_argument(\
    '--variants_file', metavar='PATH', action='store', type=str, default="variants.txt",
    help='Stores a list of variants for each transcript')
parser.add_argument(\
    '--chrsizes', type=str, required=True,
    help='The sizes of each chromosome')


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
def makeWeights(xscripts,seq_sizes,annots_handle):
    if args.num_xscripts==0:
        num_xscripts= len(xscripts)
    else:
        num_xscripts = args.num_xscripts
    weights = [1.0] * len(xscripts)
    num = 0
    for i in range(0,len(weights)):
        if num<num_xscripts and xscripts[i].seqid in seq_sizes: #Can't have no existent chromosome
            weights[i] = random.random()
            if num==0:
                annots_handle.write(xscripts[i].xscript_id)
            else:
                annots_handle.write("\n"+xscripts[i].xscript_id)
            num+=1
        else:
            weights[i] = 0
    return WeightedRandomGenerator(weights),weights

"""
Randomly picks transcripts such that half of the transcripts come from one strand and the other half come from the other strand
"""
def makeStrandedWeights(xscripts,seq_sizes,annots_handle):
    num_xscripts = args.num_xscripts
    weights = [1.0] * len(xscripts)
    last_strand = "+"
    num = 0
    for i in range(0,len(weights)):
        if num<num_xscripts and xscripts[i].orient!=last_strand and xscripts[i].seqid in seq_sizes:  #Can't have no existent chromosome
            weights[i] = random.random()
            last_strand = xscripts[i].orient
            if num==0:
                annots_handle.write(xscripts[i].xscript_id)
            else:
                annots_handle.write("\n"+xscripts[i].xscript_id)
            num+=1
        else:
            weights[i] = 0
    return WeightedRandomGenerator(weights),weights

total_reads = 0
total_mismatches = 0
"""
Incorporates sequencing error into reads
"""
def sequencingError(read,mm_rate):
    global total_reads
    global total_mismatches
    length = len(read)
    lread = list(read)
    for i in range(1,len(lread)-1):
        r = random.random()
        if r<mm_rate:
            lread[i] = "ACGT"[random.randint(0,3)]
    nread = "".join(lread)
    if nread!=read:
        total_mismatches+=1
    total_reads+=1
    return nread

"""
Incorporants variants into transcripts with nonzero weights
"""
def incorporateVariants(weights,xscripts,mm_rate,indel_rate,var_handle):
    for i in range(0,len(weights)):
        if weights[i]!=0:
            oldseq = xscripts[i].seq
            xscripts[i].incorporateVariants(mm_rate,indel_rate,var_handle)
            if xscripts[i].seq!=oldseq:
                print >> sys.stderr, "No variants!"
    return xscripts

_revcomp_trans = string.maketrans("ACGT", "TGCA")
def revcomp(s):
    return s[::-1].translate(_revcomp_trans)

def overlapping_sites(xscript,read_st,read_end):
    sites = xscript.getSites()         #in the genome coordinate frame
    xsites = xscript.getXcriptSites()  #in the transcript coordinate frame
    # sites.sort(key=lambda tup:tup[1])
    # sites.sort(key=lambda tup:tup[0])
    # xsites.sort(key=lambda tup:tup[1])
    # xsites.sort(key=lambda tup:tup[0])
    overlaps = []
    for i in range(0,len(xsites)):
        if read_st<=xsites[i][0] and read_end>=xsites[i][1]:
            overlaps.append(sites[2*i])
            overlaps.append(sites[2*i+1])
    return overlaps

def simulateSingle(xscript,sites,readlen):
    start,end,seqid = 0,len(xscript.seq)-readlen,xscript.seqid
    if end<start or len(xscript.seq)<readlen:
        print >> sys.stderr,"end<start",(end<start)
        print >> sys.stderr,"len(x.seq)<readlen",(len(xscript.seq)<readlen)
        return None,None
        #continue
    i = random.randint(start,end)
    if xscript.orient=="+":
        read = xscript.seq[i:i+readlen]
    else:
        read = revcomp(xscript.seq[i:i+readlen])
    overlaps = overlapping_sites(xscript,i,i+readlen)
    sites = sites.union(overlaps)
    read = sequencingError(read,args.readmm_rate)
    return read,sites

def simulatePairedEnd(xscript,sites,readlen):
    start,end,seqid = 0,len(xscript.seq)-readlen,xscript.seqid
    if end<start or len(xscript.seq)<readlen:
        #print >> sys.stderr, xscript.seq
        #print >> sys.stderr,"end<start",(end<start)
        #print >> sys.stderr,"len(x.seq)<readlen",(len(xscript.seq)<readlen)
        return None,None
        #continue
    i = random.randint(start,end)
    j = random.randint(start,end)
    while j < i-readlen or j>i+readlen: #get other mate
        j = random.randint(start,end)

    if xscript.orient=="+":
        mate1 = xscript.seq[i:i+readlen]
        mate2 = xscript.seq[j:j+readlen]
    else:
        mate1 = revcomp(xscript.seq[i:i+readlen])
        mate2 = revcomp(xscript.seq[j:j+readlen])

    overlaps1 = overlapping_sites(xscript,i,i+readlen)
    overlaps2 = overlapping_sites(xscript,j,j+readlen)
    sites = sites.union(overlaps1)
    sites = sites.union(overlaps2)
    mate1 = sequencingError(mate1,args.readmm_rate)
    mate2 = sequencingError(mate2,args.readmm_rate)

    return (mate1,mate2),sites

def simulate(xscripts,readlen,targetNucs,fastaseqs,var_handle,seq_sizes,annots_handle):
    if args.stranded:
        gen,weights = makeStrandedWeights(xscripts,seq_sizes,annots_handle)
    else:
        gen,weights = makeWeights(xscripts,seq_sizes,annots_handle)
    n = 0
    seqs = []
    incorporateVariants(weights,xscripts,args.snp_rate,args.indel_rate,var_handle)
    sim_xscripts = set()
    sites = set()
    while n<targetNucs:
        #Pick a transcript at weighted random
        i = gen.next()
        x = xscripts[i]
        if x not in sim_xscripts:
            sim_xscripts.add(x)
        #print x.gene_id,x.xscript_id,x.seqid
        if args.paired_end:
            tmp_reads,tmp_sites = simulatePairedEnd(x,sites,readlen)
            if tmp_reads==None and tmp_sites==None:
                continue
            else:
                reads,sites = tmp_reads,tmp_sites
            n+=readlen+readlen
            seqs.append(reads) #Appends a pair of reads
        else:
            tmp_read,tmp_sites = simulateSingle(x,sites,readlen)
            if tmp_reads==None and tmp_sites==None:
                continue
            else:
                reads,sites = tmp_reads,tmp_sites
            n+=readlen
            seqs.append(read) #Appends just one read

    return seqs,weights,list(sim_xscripts),sites

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

def writeSingleReads(seqs1rep,seqs2rep,fnPre,manifestFn):
    """ Only unpaired for now """
    seqsrep = [ seqs1rep, seqs2rep ]
    fn = "%s.seqs.tab"%(fnPre)
    with open(manifestFn, 'w') as manFh:
        for group in xrange(0, 2):
            sr = seqsrep[group]
            for rep in xrange(0, len(sr)):
                fn = "%s.group%d.rep%d.tab" % (fnPre, group, rep)
                manFh.write("%s\t0\tsplice-%d-%d\n" % (fn, group, rep))
                with open(fn, 'w') as fh:
                    for i in xrange(0, len(sr[rep])):
                        seq = sr[rep][i]
                        qual = "I" * len(seq)
                        nm = "r_n%d;LB:splice-%d-%d" % (i, group, rep)
                        fh.write("%s\t%s\t%s\n" % (nm, seq, qual))
"""
Uses the 5 tok format in align.py
"""
def writePairedEndReads(seqs1rep,seqs2rep,fnPre,manifestFn):
    """ Only unpaired for now """
    seqsrep = [ seqs1rep, seqs2rep ]
    fn = "%s.seqs.tab"%(fnPre)
    with open(manifestFn, 'w') as manFh:
        for group in xrange(0, 2):
            sr = seqsrep[group]
            for rep in xrange(0, len(sr)):
                fn = "%s.group%d.rep%d.tab" % (fnPre, group, rep)
                manFh.write("%s\t0\tsplice-%d-%d\n" % (fn, group, rep))
                with open(fn, 'w') as fh:
                    for i in xrange(0, len(sr[rep])):
                        pair = sr[rep][i]
                        mate1,mate2 = pair
                        qual = "I" * len(mate1)
                        nm = "r_n%d;LB:splice-%d-%d" % (i, group, rep)
                        fh.write("%s\t%s\t%s\t%s\t%s\n" % (nm, mate1, qual, mate2, qual))

"""
Get all transcripts that exhibit alternative splicing
"""
def test_alternativeSplicing(xscripts):
    genes = defaultdict(list)
    axscripts = []
    for x in xscripts:
        genes[x.gene_id].append(x)
    for gene_ids,isoforms in genes.iteritems():
        if len(isoforms)>1:
            axscripts+=isoforms
            return axscripts
    #return axscripts

if __name__=="__main__":
    print >> sys.stderr,"Mismatch Rate",args.readmm_rate
    print >> sys.stderr,"SNP Rate",args.snp_rate
    print >> sys.stderr,"Indel Rate",args.indel_rate
    seq_sizes = chrsizes.getSizes(args.chrsizes)
    annots_handle = open(args.output_annotations,'w')
    var_handle = open(args.variants_file,'w')
    annots = gtf.parseGTF(args.gtf)
    fastadb = gtf.parseFASTA(args.fasta)
    xscripts = gtf.assembleTranscripts(annots,fastadb)
    if args.alternative_spliced:
        xscripts = test_alternativeSplicing(xscripts)
    seqs,weights,xscripts,sites = simulate(xscripts,args.read_len,args.num_nucs,args.fasta,var_handle,seq_sizes,annots_handle)
    seqs1,seqs2 = replicateize(seqs,seqs,args.num_replicates)
    if args.paired_end:
        writePairedEndReads(seqs1,seqs2,args.output_prefix,args.output_prefix+".manifest")
    else:
        writeSingleReads(seqs1,seqs2,args.output_prefix,args.output_prefix+".manifest")
    #This stores the list in pickle files for serialization
    #pickle.dump(weights,open(args.output_prefix+".weights",'wb'))
    pickle.dump(xscripts,open(args.output_prefix+".xscripts",'wb'))
    #sites = list(sites)
    real_sites,sim_sites = xscripts[0].getSites(),list(sites)
    real_sites.sort()
    sim_sites.sort()
    real_sites = "\t".join(map(str,real_sites))
    sim_sites = "\t".join(map(str,sim_sites))
    #print >> sys.stderr,"Transcripts    ",real_sites
    #print >> sys.stderr,"Simulated sites",sim_sites

    pickle.dump(sites,open(args.output_prefix+".sites",'wb'))

    #site_handle = open(args.output_prefix+".sites",'w')
    #sites = map( str, sites)
    #site_handle.write("\t".join(sites))

    print >> sys.stderr,"Total number of reads",total_reads
    print >> sys.stderr,"Total number of mismatched reads",total_mismatches

    #print >> sys.stderr,"Sites",sites
