"""
sim_splice.py

Simulates differiential gene expression using annotated genes.

Outputs
=======

The simulator outputs several kinds of objects:

1. Simulated reads
2. Manifest file
3. File describing sequence variants added
4. Transcripts sampled
5. Coverage
6. List of sites traversed

Strandedness
============

A few notes about strands and strandedness of RNA-seq protocols.  For a given
gene, the DNA strand that contains the actual translated codon sequences is
called the "sense" strand, and its complement is the "anti-sense" strand.  For
a given gene, the sense strand might be either the Watson or the Crick strand.
It varies from gene to gene.

An RNA-seq protocol (lab procedures used to generate the RNA-seq data) is
either "stranded" or not.  In a stranded protocol, the reads always (or almost
always, since no protocol is 100% efficient) comes from the sense strand.  In a
non-stranded protocol, the read is equally likely to come from either strand
(either sense or anti-sense or, equivalently, either Watson or Crick).  Most
RNA-seq datasets were generated with non-stranded protocols.

See also:

Levin JZ, Yassour M, Adiconis X, Nusbaum C, Thompson DA, Friedman N, Gnirke A,
Regev A. Comprehensive comparative analysis of strand-specific RNA sequencing
methods. Nat Methods. 2010 Sep;7(9):709-15.

The --stranded option tells us to sample reads only from the sense strand.
When --stranded is not specified, we sample reads from both the sense and
anti-sense strands approximately equally

Paired-end
==========

A brief note about paired-end sequencing.  When doing paired-end sequencing,
we're really sequencing in from both ends of a larger fragment, i.e.:

 GGTCATCTGCACGTAT
 >>>>>>>>>>>>>>>>
 GGTCATCTGCACGTATCGTACGTACGTTCACTACTCTGGACTACGTACGAG
                                  <<<<<<<<<<<<<<<<<<
                                  CTCTGGACTACGTACGAG

Mate #1 on top, mate #2 on bottom.  Two important things to note.  First: the
fragment itself might have come from the Watson or Crick stand of the genome.
Second: mate #2 is reverse-complemented with respect to mate #1.

TODO:
- Output a bed file with the simulated reads w/ splices.  Stick
  "name=junctions" at the top so that viewers know to show connecting arcs,
  per: http://www.broadinstitute.org/igv/splice_junctions
- Fix issue in test_alternativeSplicing (see my comment)
- Clarify purpose of st and en in overlapping_sites (see my comment)

"""

import string
import os
import site
import argparse
import sys
import random
import bisect
import pickle
from collections import defaultdict

base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
site.addsitedir(os.path.join(base_path, "annotation"))
site.addsitedir(os.path.join(base_path, "fasta"))
site.addsitedir(os.path.join(base_path, "util"))

import gtf
import chrsizes
import counter
import fasta

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
    '--paired-end', action='store_const', const=True, default=False,
    help='Generate paired-end reads')
parser.add_argument(\
    '--fragment-mean', metavar='float', action='store', type=float, default=300.0,
    help='Mean of (gaussian) fragment distribution')
parser.add_argument(\
    '--fragment-sd', metavar='float', action='store', type=float, default=15.0,
    help='Standard deviation of (gaussian) fragment distribution')
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
    '--canonical-sites', action='store_const', const=True, default=False,
    help='Only simulates reads spanning canonical sites')
parser.add_argument(\
    '--noncanonical-sites', action='store_const', const=True, default=False,
    help='Only simulates reads spanning noncanonical sites')
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

def factorial(x):
    return 1 if x==0 else x*factorial(x-1)

# n choose k
# n = read length, k = number of sequencing errors
# binomial coefficient
b_coeff = [ factorial(args.read_len)/(factorial(k)*factorial(args.read_len-k)) for k in range(0,args.read_len+1)]
# Creates an error model using a binomial random number generator based off of read length=N
def errorPMF(N,mm_rate):
    return WeightedRandomGenerator([ b_coeff[i] * (1-mm_rate)**(N-i) * mm_rate**i for i in range(0,N)])

"""
Picks transcripts arbitrarily (i.e. starting from the beginning of the
xscripts list) such that half of the transcripts come from the Watson strand
and the other half from Crick.
"""
def makeWeights(xscripts, seq_sizes, annots_handle):
    num_xscripts = args.num_xscripts
    weights = [1.0] * len(xscripts)
    last_strand = "+"
    num, i = 0, 0
    # We're picking the first num_xscripts transcripts, but we're constrained
    # to alternate between Watson and Crick being the sense strand.
    while i < len(weights) and num < num_xscripts:
        if xscripts[i].orient != last_strand and xscripts[i].seqid in seq_sizes:
            weights[i] = random.random()
            last_strand = xscripts[i].orient
            annots_handle.write(xscripts[i].xscript_id)
            annots_handle.write('\n')
            num += 1
        i += 1
    return WeightedRandomGenerator(weights), weights

total_reads = 0
total_mismatches = 0

"""
Incorporates sequencing error into reads
"""
_otherNucs = { 'A' : ['C', 'G', 'T'],
               'C' : ['A', 'G', 'T'],
               'G' : ['A', 'C', 'T'],
               'T' : ['A', 'C', 'G'],
               'N' : ['N'] }

def sequencingError(read,errModel):
    global total_reads
    global total_mismatches

    length = len(read)
    lread = list(read)
    #num = random.expovariate(1.0/mm_rate)
    #print >> sys.stderr,"expo num",num,"mean",(1.0/mm_rate)
    errs = errModel.next()
    if errs == 0:
        total_reads += 1
        return "".join(lread)

    for _ in xrange(errs):
        r = random.randint(0, length - 1)
        c = lread[r]
        if c not in _otherNucs: c = 'N'
        lread[r] = random.choice(_otherNucs[c])

    total_mismatches += 1
    total_reads += 1
    return "".join(lread)

"""
Incorporates variants into transcripts with nonzero weights
"""
def incorporateVariants(weights,xscripts,mm_rate,indel_rate,var_handle):
    if mm_rate==0 and indel_rate==0:
        return xscripts
    for i in xrange(len(weights)):
        if weights[i] != 0:
            oldseq = xscripts[i].seq
            xscripts[i].incorporateVariants(mm_rate,indel_rate,var_handle)
            if xscripts[i].seq != oldseq:
                print >> sys.stderr, "No variants!"
    return xscripts

_revcomp_trans = string.maketrans("ACGT", "TGCA")
def revcomp(s):
    return s[::-1].translate(_revcomp_trans)

def overlapCanonical(xscript,read_st,read_end,can_sites):
    """Checks to see if the read spans a canonical splice site"""
    st, en = read_st + xscript.st0, read_end + xscript.st0
    sites = xscript.getSitePairs()                          #in the genome coordinate frame
    #can_sites = set(xscript.getCanonicalSites(fastaHandle)) #canonical sites in the genome coordinate frame
    xsites = xscript.getXcriptSites()                       #in the transcript coordinate frame
    for i in xrange(len(xsites)):
        if read_st <= xsites[i][0] and read_end >= xsites[i][1]:
            if sites[i] in can_sites:
                return True
    return False

def overlapNonCanonical(xscript,read_st,read_end,ncan_sites):
    """Checks to see if the read spans a canonical splice site"""
    st, en = read_st + xscript.st0, read_end + xscript.st0
    sites = xscript.getSitePairs()                  #in the genome coordinate frame
    #ncan_sites = set(xscript.getNonCanonicalSites(fastaHandle)) #noncanonical sites in the genome coordinate frame
    xsites = xscript.getXcriptSites()               #in the transcript coordinate frame
    for i in xrange(len(xsites)):
        if read_st <= xsites[i][0] and read_end >= xsites[i][1] and sites[i] in ncan_sites:
            return True
    return False

def overlapSpecificSites(xscript,read_st,read_end,siteType, fastaHandle):
    """Returns true if the read overlaps a set of either canonical or noncanonical sites as specified by siteType"""
    if siteType=="canonical":
        can_sites = set(xscript.getCanonicalSites(fastaHandle)) #canonical sites in the genome coordinate frame
        return overlapCanonical(xscript,read_st,read_end,can_sites)
    else:
        ncan_sites = set(xscript.getNonCanonicalSites(fastaHandle)) #canonical sites in the genome coordinate frame
        return overlapNonCanonical(xscript,read_st,read_end,ncan_sites)

def overlapping_sites(xscript, read_st, read_end):
    """ Given a transcript and an interval on the transcript, return a list of
        splice sites that are spanned by the interval """
    sites = xscript.getSites()         #in the genome coordinate frame
    xsites = xscript.getXcriptSites()  #in the transcript coordinate frame
    # sites is twice as long as xsites??
    # sites is twice as long as xsite for the following reason
    """
    =======|^^ ^^|==========   reference genome
    -------^     ^----------   transcript
    Splice sites aren't present in the transcript, so only two break points are necessary.
    For completeness sake, the splice sites in the reference genome are included, making sites = 2*xsites
    """

    # sites.sort(key=lambda tup:tup[1])
    # sites.sort(key=lambda tup:tup[0])
    # xsites.sort(key=lambda tup:tup[1])
    # xsites.sort(key=lambda tup:tup[0])
    overlaps = []
    st, en = read_st + xscript.st0, read_end + xscript.st0
    for i in xrange(len(xsites)):
        if read_st <= xsites[i][0] and read_end >= xsites[i][1]:
            overlaps.append(sites[2*i])
            overlaps.append(sites[2*i+1])
        diff = sites[2*i+1][1] - sites[2*i][0] + 1
        if read_st >= xsites[i][0]:
            st += diff
        if read_end >= xsites[i][0]:
            en += diff
    # BTL: what are st and en?
    # Jamie: They are the adjusted starts and ends of the reads
    # It is adjusted because the read positions in transcriptome coordinates
    # are being converted to genome coordinates
    return overlaps, st, en

def simulateSpanningSingle(xscript, readlen, errModel, siteType, fastaHandle):
    """Simulate a single unpaired read from a given transcript that overlaps a either a canonical site ir a noncanonical site"""
    start, end = 0, len(xscript.seq) - readlen
    if end < start or len(xscript.seq) < readlen:
        return None,None,None
    if (( len(xscript.getCanonicalSites(fastaHandle))==0 and args.canonical_sites ) or
        ( len(xscript.getNonCanonicalSites(fastaHandle))==0 and args.noncanonical_sites )):
        return None,None,None

    i = random.randint(start, end)
    while not overlapSpecificSites(xscript,start,end,siteType, fastaHandle):
        i = random.randint(start, end)

    # Read is always taken from sense strand?
    read = xscript.seq[i:i+readlen]

    if args.stranded:
        if xscript.orient == "-":
            read = revcomp(read)
    elif random.random() < 0.5:
        read = revcomp(read)
    sites, st, end = overlapping_sites(xscript, i, i + readlen)
    read = sequencingError(read, errModel)
    return read, sites, (st, end)


def simulateSingle(xscript, readlen, errModel):
    """ Simulate a single unpaired read from the given transcript """
    start, end = 0, len(xscript.seq) - readlen
    if end < start or len(xscript.seq) < readlen:
        return None,None,None
    i = random.randint(start, end)
    # Read is always taken from sense strand?
    read = xscript.seq[i:i+readlen]
    if args.stranded:
        if xscript.orient == "-":
            read = revcomp(read)
    elif random.random() < 0.5:
        read = revcomp(read)

    sites, st, end = overlapping_sites(xscript, i, i + readlen)
    read = sequencingError(read, errModel)
    return read, sites, (st, end)

def simulatePairedEnd(xscript, readlen, errModel, fraglenGen=lambda: None):
    fraglen = fraglenGen()
    if fraglen is None:
        fraglen = readlen * 2.5
    if len(xscript.seq) < fraglen:
        return None, None, None, None
    start, end = 0, len(xscript.seq) - fraglen

    frag_i = random.randint(start, end)
    i, j = frag_i, frag_i + fraglen - readlen
    mate1 = xscript.seq[i : i + readlen]
    mate2 = revcomp(xscript.seq[j : j + readlen])
    if (args.stranded and xscript.orient == "-") or \
       (not args.stranded and random.random() < 0.5):
        # If mate1 and mate2 are different lengths, i and j must be adjusted
        mate1, mate2 = revcomp(mate2), revcomp(mate1)

    overlaps1, st1, end1 = overlapping_sites(xscript, i, i + readlen)
    overlaps2, st2, end2 = overlapping_sites(xscript, j, j + readlen)
    #sites = sites.union(overlaps1)
    #sites = sites.union(overlaps2)
    sites = set(overlaps1)
    sites |= set(overlaps2)
    mate1 = sequencingError(mate1,errModel)
    mate2 = sequencingError(mate2,errModel)

    return (mate1,mate2),sites,(st1,end1),(st2,end2)

def simulate(xscripts,readlen,targetNucs,fastaseqs,var_handle,seq_sizes,annots_handle):
    #TODO: Think about using fasta dictionary instead ...
    fastaHandle = fasta.fasta(args.fasta[0])#Just look at first fasta file for now

    #
    # Step 1: Assign weights to isoforms
    #
    gen, weights = makeWeights(xscripts,seq_sizes,annots_handle)
    n = 0
    seqs = []
    cov_sts, cov_ends = counter.Counter(), counter.Counter()

    #
    # Step 2: Incorporate sequence variants
    #
    incorporateVariants(weights,xscripts,args.snp_rate,args.indel_rate,var_handle)

    #
    # Step 3: Generate sequence reads
    #
    sim_xscripts = set()
    sites = set()
    errModel = errorPMF(readlen,args.readmm_rate)
    while n < targetNucs:
        # Pick a transcript at weighted random
        i = gen.next()
        x = xscripts[i]
        if x not in sim_xscripts:
            sim_xscripts.add(x)
        # print x.gene_id,x.xscript_id,x.seqid
        if args.canonical_sites or args.noncanonical_sites:
            siteType = "canonical" if args.canonical_sites else "noncanonical"
            tmp_reads,tmp_sites,bounds = simulateSpanningSingle(x,readlen,errModel,siteType,fastaHandle)
            if tmp_reads is None and tmp_sites is None:
                continue
            else:
                reads = tmp_reads
                sites |= set(tmp_sites)
        elif args.paired_end:
            def fraglenGen():
                return int(random.gauss(args.fragment_mean, args.fragment_sd))
            tmp_mates,tmp_sites,bounds1,bounds2 = simulatePairedEnd(x, readlen, errModel, fraglenGen)
            if tmp_mates is None and tmp_sites is None:
                continue
            else:
                reads = tmp_mates
                sites |= set(tmp_sites)
            n += readlen + readlen
            cov_sts[ bounds1[0] ]+=1
            cov_sts[ bounds2[0] ]+=1
            cov_ends[ bounds1[1] ]+=1
            cov_ends[ bounds2[1] ]+=1
            seqs.append(tmp_mates) #Appends a pair of reads
        else:
            tmp_reads,tmp_sites,bounds = simulateSingle(x,readlen,errModel)
            if tmp_reads is None and tmp_sites is None:
                continue
            else:
                reads = tmp_reads
                sites |= set(tmp_sites)
            n+=readlen
            # What's up with cov_sts and cov_ends??
            # Jamie: cov_sts and cov_ends keeps a count of start and end positions
            cov_sts[ bounds[0] ]+=1
            cov_ends[ bounds[1] ]+=1
            seqs.append(reads) #Appends just one read

    return seqs, weights, list(sim_xscripts), sites, cov_sts, cov_ends

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
Randomly gets the first set of transcript isoforms from the same transcript that exhibit alternative splicing
"""
def test_alternativeSplicing(xscripts):
    genes = defaultdict(list)
    axscripts = []
    for x in xscripts:
        genes[x.gene_id].append(x)
    # BTL: this doesn't look right.  Return axscripts with only one element?
    # Jamie: Yeah, it isn't right.  I addressed it in the above comment
    for _, isoforms in genes.iteritems():
        if len(isoforms) > 1:
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
    seqs, _, xscripts, sites, cov_sts, cov_ends = simulate(xscripts,args.read_len,args.num_nucs,args.fasta,var_handle,seq_sizes,annots_handle)
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
        class TestSimulationFunctions1(unittest.TestCase):
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
                for _ in xrange(N): cnts[model.next()] += 1
                print >> sys.stderr,"Histogram",cnts
                self.assertGreater(cnts[1],cnts[0])
                self.assertGreater(cnts[0],cnts[2])
                self.assertGreater(cnts[1],cnts[2])

        class TestSimulationFunctions2(unittest.TestCase):
            def setUp(self):
                """           ^^                        ^^      ^^        """
                """       CAACTGATGT CTTCGCTTGT GAAACGAACC ATATTGATCGC GCA"""
                """       0          1          2          3          4          5          6          7         8        """
                """       [AAC]GTGAT CAAG[ATGTC TTCGCTTGTG AAACGAA]GT CTGGATCCGC CTGAAG[ATA TT]GCAATAA GATCGCGCAG [TCGCGC]"""
                """       CAACTGTGAT CAAGGATGTC TTCGCTTGTG AAACGAACGT CTGGATCCGC CTGAAGCATA TTGGCAATAA GATCGCGCAG ATCGCGCA"""
                refseq="""CAACTGTGATCAAGGATGTCTTCGCTTGTGAAACGAACGTCTGGATCCGCCTGAAGCATATTGGCAATAAGATCGCGCAGATCGCGCA"""
                annots="""chr2R\tunknown\texon\t1\t5\t.\t-\t.\tgene_id "CG17528"; gene_name "CG17528"; p_id "P21588"; transcript_id "NM_001042999"; tss_id "TSS13109";\nchr2R\tunknown\texon\t15\t38\t.\t-\t.\tgene_id "CG17528"; gene_name "CG17528"; p_id "P21588"; transcript_id "NM_001042999"; tss_id "TSS13109";\nchr2R\tunknown\texon\t57\t63\t.\t-\t.\tgene_id "CG17528"; gene_name "CG17528"; p_id "P21588"; transcript_id "NM_001042999"; tss_id "TSS13109";\nchr2R\tunknown\texon\t80\t88\t.\t-\t.\tgene_id "CG17528"; gene_name "CG17528"; p_id "P21588"; transcript_id "NM_001042999"; tss_id "TSS13109";\n"""
                self.fasta = "test.fa"
                self.faidx = "test.fa.fai"
                self.gtf   = "test.gtf"
                createTestFasta(self.fasta,"chr2R",refseq)
                createTestGTF(self.gtf,annots)

            def tearDown(self):
                os.remove(self.fasta)
                os.remove(self.faidx)
                os.remove(self.gtf)

            def test_canonical_overlap1(self):
                annots = gtf.parseGTF([self.gtf])
                fastadb = gtf.parseFASTA([self.fasta])
                xscripts = gtf.assembleTranscripts(annots,fastadb)
                print readableFormat(fastadb["chr2R"])
                read_st,read_end = 1,10
                fastaHandle = fasta.fasta(self.fasta)
                didOverlap = overlapSpecificSites(xscripts[0],read_st,read_end,"canonical",fastaHandle)
                print >> sys.stderr,"read",read_st,read_end,didOverlap
                self.assertEquals(didOverlap, True)

                read_st,read_end = 22,32
                didOverlap = overlapSpecificSites(xscripts[0],read_st,read_end,"canonical",fastaHandle)
                print >> sys.stderr,"read",read_st,read_end,didOverlap
                self.assertEquals(didOverlap, True)

                read_st,read_end = 32,42
                didOverlap = overlapSpecificSites(xscripts[0],read_st,read_end,"canonical",fastaHandle)
                print >> sys.stderr,"read",read_st,read_end,didOverlap
                self.assertEquals(didOverlap, False)

            def test_noncanonical_overlap1(self):
                annots = gtf.parseGTF([self.gtf])
                fastadb = gtf.parseFASTA([self.fasta])
                xscripts = gtf.assembleTranscripts(annots,fastadb)
                print readableFormat(fastadb["chr2R"])
                read_st,read_end = 1,10
                fastaHandle = fasta.fasta(self.fasta)
                didOverlap = overlapSpecificSites(xscripts[0],read_st,read_end,"noncanonical",fastaHandle)
                print >> sys.stderr,"read",read_st,read_end,didOverlap
                self.assertEquals(didOverlap, False)

                read_st,read_end = 22,32
                didOverlap = overlapSpecificSites(xscripts[0],read_st,read_end,"noncanonical",fastaHandle)
                print >> sys.stderr,"read",read_st,read_end,didOverlap
                self.assertEquals(didOverlap, False)

                read_st,read_end = 32,42
                didOverlap = overlapSpecificSites(xscripts[0],read_st,read_end,"noncanonical",fastaHandle)
                print >> sys.stderr,"read",read_st,read_end,didOverlap
                self.assertEquals(didOverlap, True)



            def test_errorPMF(self):
                N = 1000000
                n = 100
                rate = 0.01
                cnts = counter.Counter()
                model = errorPMF(n,rate)
                for _ in xrange(N): cnts[model.next()] += 1
                print >> sys.stderr,"Histogram",cnts
                self.assertGreater(cnts[1],cnts[0])
                self.assertGreater(cnts[0],cnts[2])
                self.assertGreater(cnts[1],cnts[2])

        unittest.main()
