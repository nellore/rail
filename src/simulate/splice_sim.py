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
import chrsizes
import counter

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
    '--chrsizes', type=str, required=False,
    help='The sizes of each chromosome')
parser.add_argument(\
    '--test', action='store_const', const=True, default=False,
    help='Run unit tests')
parser.add_argument(\
    '--profile', action='store_const', const=True, default=False,
    help='Profile simulation generation')

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

#Creates an error model using a binomial random number generator based off of read length=N
def errorPMF(N,mm_rate):
    return WeightedRandomGenerator([ (1-mm_rate)**(N-i) * mm_rate**i for i in range(0,N)])

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
def sequencingError(read,errModel):
    global total_reads
    global total_mismatches
    length = len(read)
    lread = list(read)
    #num = random.expovariate(1.0/mm_rate)
    #print >> sys.stderr,"expo num",num,"mean",(1.0/mm_rate)
    errs = errModel.next()
    if errs==0:
        total_reads+=1
        return "".join(lread)

    for i in range(0,errs):
        r = random.randint(1,len(lread)-1)
        bases = ["A","C","G","T"]
        bases.remove(lread[r])
        lread[r] = bases[random.randint(0,2)]
        
    nread = "".join(lread)
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
    st,en = read_st+xscript.st0,read_end+xscript.st0
    for i in range(0,len(xsites)):
        if read_st<=xsites[i][0] and read_end>=xsites[i][1]:
            overlaps.append(sites[2*i])
            overlaps.append(sites[2*i+1])
        diff = sites[2*i+1][1] - sites[2*i][0] + 1
        if read_st>=xsites[i][0]:
            st+=diff
        if read_end>=xsites[i][0]:
            en+=diff
    return overlaps,st,en


def simulateSingle(xscript,readlen,errModel):
    start,end,seqid = 0,len(xscript.seq)-readlen,xscript.seqid
    if end<start or len(xscript.seq)<readlen:
        return None,None,None
    i = random.randint(start,end)
    if xscript.orient=="+":
        read = xscript.seq[i:i+readlen]
    else:
        read = revcomp(xscript.seq[i:i+readlen])
    sites,st,end = overlapping_sites(xscript,i,i+readlen)
    read = sequencingError(read,errModel)
    return read,sites,(st,end)

def simulatePairedEnd(xscript,readlen,errModel):
    start,end,seqid = 0,len(xscript.seq)-readlen,xscript.seqid
    if end<start or len(xscript.seq)<readlen:
        #print >> sys.stderr, xscript.seq
        #print >> sys.stderr,"end<start",(end<start)
        #print >> sys.stderr,"len(x.seq)<readlen",(len(xscript.seq)<readlen)
        return None,None,None,None
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

    overlaps1,st1,end1 = overlapping_sites(xscript,i,i+readlen)
    overlaps2,st2,end2 = overlapping_sites(xscript,j,j+readlen)
    #sites = sites.union(overlaps1)
    #sites = sites.union(overlaps2)
    sites = set(overlaps1)
    sites |= set(overlaps2)
    mate1 = sequencingError(mate1,errModel)
    mate2 = sequencingError(mate2,errModel)

    return (mate1,mate2),sites,(st1,end1),(st2,end2)

def simulate(xscripts,readlen,targetNucs,fastaseqs,var_handle,seq_sizes,annots_handle):
    if args.stranded:
        gen,weights = makeStrandedWeights(xscripts,seq_sizes,annots_handle)
    else:
        gen,weights = makeWeights(xscripts,seq_sizes,annots_handle)
    n = 0
    seqs = []
    cov_sts, cov_ends = counter.Counter(), counter.Counter()
    incorporateVariants(weights,xscripts,args.snp_rate,args.indel_rate,var_handle)
    sim_xscripts = set()
    sites = set()
    errModel = errorPMF(readlen,args.readmm_rate)
    while n<targetNucs:
        #Pick a transcript at weighted random
        i = gen.next()
        x = xscripts[i]
        if x not in sim_xscripts:
            sim_xscripts.add(x)
        #print x.gene_id,x.xscript_id,x.seqid
        if args.paired_end:
            tmp_reads,tmp_sites,bounds1,bounds2,overlaps= simulatePairedEnd(x,readlen,errModel)
            if tmp_reads==None and tmp_sites==None:
                continue
            else:
                reads = tmp_reads
                sites|= set(tmp_sites)
            n+=readlen+readlen
            cov_sts[ bounds1[0] ]+=1
            cov_sts[ bounds2[0] ]+=1
            cov_ends[ bounds1[1] ]+=1
            cov_ends[ bounds2[1] ]+=1
            seqs.append(reads) #Appends a pair of reads
        else:
            tmp_reads,tmp_sites,bounds = simulateSingle(x,readlen,errModel)
            if tmp_reads==None and tmp_sites==None:
                continue
            else:
                reads = tmp_reads
                sites|= set(tmp_sites)
            n+=readlen
            cov_sts[ bounds[0] ]+=1
            cov_ends[ bounds[1] ]+=1
            seqs.append(reads) #Appends just one read

    return seqs,weights,list(sim_xscripts),sites,cov_sts,cov_ends

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

def go():
    print >> sys.stderr,"Mismatch Rate",args.readmm_rate
    print >> sys.stderr,"SNP Rate",args.snp_rate
    print >> sys.stderr,"Indel Rate",args.indel_rate
    seq_sizes = chrsizes.getSizes(args.chrsizes)
    annots_handle = open(args.output_annotations,'w')
    var_handle = open(args.output_prefix+".variants",'w')
    annots = gtf.parseGTF(args.gtf)
    fastadb = gtf.parseFASTA(args.fasta)
    xscripts = gtf.assembleTranscripts(annots,fastadb)
    if args.alternative_spliced:
        xscripts = test_alternativeSplicing(xscripts)
    seqs,weights,xscripts,sites,cov_sts,cov_ends = simulate(xscripts,args.read_len,args.num_nucs,args.fasta,var_handle,seq_sizes,annots_handle)
    seqs1,seqs2 = replicateize(seqs,seqs,args.num_replicates)
    if args.paired_end:
        writePairedEndReads(seqs1,seqs2,args.output_prefix,args.output_prefix+".manifest")
    else:
        writeSingleReads(seqs1,seqs2,args.output_prefix,args.output_prefix+".manifest")
    #This stores the list in pickle files for serialization
    pickle.dump(xscripts,open(args.output_prefix+".xscripts",'wb'))
    ###BIG NOTE:  This pickles a tuple of counter objects
    pickle.dump( (cov_sts,cov_ends) ,open(args.output_prefix+".cov",'wb'))
    real_sites,sim_sites = xscripts[0].getSites(),list(sites)
    real_sites.sort()
    sim_sites.sort()
    real_sites = "\t".join(map(str,real_sites))
    sim_sites = "\t".join(map(str,sim_sites))
    
    pickle.dump(sites,open(args.output_prefix+".sites",'wb'))
    print >> sys.stderr,"Total number of reads",total_reads
    print >> sys.stderr,"Total number of mismatched reads",total_mismatches

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

def readableFormat(s):
    return "\t".join([s[i:i+10] + " " + str(i+10) for i in range(0,len(s),10)])

if __name__=="__main__":
    if not args.test and not args.profile:
        go() 
    elif args.profile:
        import cProfile
        cProfile.run('go()')
    else:
        del sys.argv[1:]
        import unittest
        class TestSimulationFunctions(unittest.TestCase):
            def setUp(self):
                """       CAACTGTGAT (CAAGGATGTC) TTCGCTTGTG (AAACGAACGT) CTGGATCCGC (CTGAAGCATA) TTGGCAATAA GATCGCGCA"""
                refseq="""CAACTGTGATCAAGGATGTCTTCGCTTGTGAAACGAACGTCTGGATCCGCCTGAAGCATATTGGCAATAAGATCGCGCA"""
                annots="""chr2R\tunknown\texon\t11\t20\t.\t-\t.\tgene_id "CG17528"; gene_name "CG17528"; p_id "P21588"; transcript_id "NM_001042999"; tss_id "TSS13109";\nchr2R\tunknown\texon\t31\t40\t.\t-\t.\tgene_id "CG17528"; gene_name "CG17528"; p_id "P21588"; transcript_id "NM_001042999"; tss_id "TSS13109";\nchr2R\tunknown\texon\t51\t60\t.\t-\t.\tgene_id "CG17528"; gene_name "CG17528"; p_id "P21588"; transcript_id "NM_001042999"; tss_id "TSS13109";\n"""
                self.fasta = "test.fa"
                self.faidx = "test.fa.fai"
                self.gtf   = "test.gtf"
                createTestFasta(self.fasta,"chr2R",refseq)
                createTestGTF(self.gtf,annots)
            def tearDown(self):
                os.remove(self.fasta)
                os.remove(self.faidx)
                os.remove(self.gtf)
            def test_overlap1(self):
                annots = gtf.parseGTF([self.gtf])
                fastadb = gtf.parseFASTA([self.fasta])
                xscripts = gtf.assembleTranscripts(annots,fastadb)
                print readableFormat(fastadb["chr2R"])
                read_st,read_end = 5,15
                _,st,end = overlapping_sites(xscripts[0],read_st,read_end)
                print >> sys.stderr,"read",read_st,read_end
                print >> sys.stderr,"ref",st,end
                self.assertEquals(15,st)
                self.assertEquals(35,end)
            def test_overlap2(self):
                annots = gtf.parseGTF([self.gtf])
                fastadb = gtf.parseFASTA([self.fasta])
                xscripts = gtf.assembleTranscripts(annots,fastadb)
                print readableFormat(fastadb["chr2R"])
                read_st,read_end = 5,25
                _,st,end = overlapping_sites(xscripts[0],read_st,read_end)
                print >> sys.stderr,"read",read_st,read_end
                print >> sys.stderr,"ref",st,end
                self.assertEquals(15,st)
                self.assertEquals(55,end)
            def test_errorPMF(self):
                N = 1000000
                n = 100
                rate = 0.01
                cnts = counter.Counter()
                model = errorPMF(n,rate)
                for i in range(0,N):
                    cnts[model.next()]+=1
                print >> sys.stderr,"Histogram",cnts
                self.assertGreater(cnts[0],cnts[1]) 
                self.assertGreater(cnts[0],cnts[2])
                self.assertGreater(cnts[1],cnts[2])

        unittest.main()
