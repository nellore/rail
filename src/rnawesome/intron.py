import math
import numpy as np
"""
Tab-delimited input tuple columns:
1. Partition ID for partition overlapped by interval (also includes strand information)
2. Interval start
3. Interval end (exclusive)
4. Reference ID
5. Sample label
6. Readlet Sequence before 5' site
7. Readlet Sequence after 5' site
8. Readlet Sequence before 3' site
9. Readlet Sequence after 3' site

Tab-delimited splice site output tuple columns:
1. Reference ID
2. 5' start
3. 3' start
4. Sample label
5. Read frequency (number of times sample read overlapped junction)

Tab-delimited cooccurence output tuple columns:
1. rdid
2. refID
3. left_site
4. right_site

TODO:
1) Fix the sliding window skewness problem
2) Test against a mismatch rate of 0.01% for small dataset
3) Test against genetic variation for small dataset
4) Test against huge data set using the validation problem for verification
5) Run against human data set
"""
import os
import sys
import argparse
import site
import time
import re
import string
from collections import defaultdict
timeSt = time.clock()

base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
site.addsitedir(os.path.join(base_path, "fasta"))
site.addsitedir(os.path.join(base_path, "read"))
site.addsitedir(os.path.join(base_path, "alignment"))
site.addsitedir(os.path.join(base_path, "statsmath"))
site.addsitedir(os.path.join(base_path, "util"))

import chrsizes
import fasta
import readlet
import needlemanWunsch
import histogram
import counter
import window

parser = argparse.ArgumentParser(description=\
                                     'Reports splice junction information')
parser.add_argument(\
    '--test', action='store_const', const=True, default=False, help='Run unit tests')
parser.add_argument(\
    '--refseq', type=str, required=False,
    help='The fasta sequence of the reference genome. The fasta index of the reference genome is also required')
parser.add_argument(\
    '--radius', type=int, required=False,default=10,
    help='The radius of the clustering algorithm and the radius of the sliding window algorithm.')
parser.add_argument(\
    '--scores-file', type=str, required=False,default="",
    help='The location of the scores output file')


readlet.addArgs(parser)

args = parser.parse_args()

ninp = 0                   # # lines input so far
nout = 0                   # # lines output so far

"""
Conducts radial clustering.
All introns that have start and end positions within 2 readletIvals within each other are binned together
"""
def cluster(ivals):
    points = []
    bins = defaultdict(list)
    rIval = args.radius
    p = ivals[0]
    points.append(p)
    #key = "%s,%s"%(p[0],p[1])
    bins[p].append(ivals[0])
    notFound = True
    for i in range(1,len(ivals)):
        for j in range(0,len(points)): #Check all the neighborhood of all points
            p = points[j]
            ival = ivals[i]
            if (ival[0]>(p[0]-rIval) and ival[0]<(p[0]+rIval) and
                ival[1]>(p[1]-rIval) and ival[1]<(p[1]+rIval) ):
                #key = "%s,%s"%(p[0],p[1])
                #bins[key].append(ivals[i])
                bins[p].append(ivals[i])
                notFound = False
                break
        if notFound:
            p = ivals[i]
            points.append(p)
            #key = "%s,%s"%(p[0],p[1])
            bins[p].append(ivals[i])
        notFound = True
    return bins
"""
Clusters by start and intron length instead of start and end
"""
def diagonal_cluster(ivals):
    points = []
    bins = defaultdict(list)
    rIval = args.radius
    p = ivals[0]
    points.append(p)
    #key = "%s,%s"%(p[0],p[1])
    bins[p].append(ivals[0])
    notFound = True
    for i in range(1,len(ivals)):
        for j in range(0,len(points)): #Check all the neighborhood of all points
            ival = ivals[i]
            p = points[j]
            st1,end1 = ival[0],ival[1]
            st2,end2 = p[0],p[1]
            if ( abs(st1-st2)<rIval and (end1-st1)==(end2-st2) ):
                bins[p].append(ivals[i])
                notFound = False
                break
        if notFound:
            p = ivals[i]
            points.append(p)
            bins[p].append(ivals[i])
        notFound = True
    return bins


"""
Just a fancier way to print out lists
"""
def format_list(L):
    return "\t".join(["%.2f"%i for i in L])
def format_seq(L):
    return "\t".join( map( str,L))


"""
Note: site is formatted as follows: XX-XX (e.g. GT-AG)
Returns the 5' and 3' splice sites within multiple intervals
n is length of the histogram window, which is specified by user
"""
def sliding_window(refID, sts,ens, site, n, fastaF):
    #n,r = 2*args.radius, args.radius

    in_start, in_end = min(sts),max(ens)
    toks = site.split("-")
    assert len(toks)==2
    site5p,site3p = toks[0],toks[1]

    hist5 = histogram.hist_score(sts,in_start,"5",2*n+1)
    hist3 = histogram.hist_score(ens,in_end,"3",2*n+1)
    mean5,std5 = hist5.index(max(hist5))+1,   histogram.stddev(hist5)/2 #offset bias correction for 5' end
    mean3,std3 = hist3.index(max(hist3)),   histogram.stddev(hist3)/2
    #Create a normal distributed scoring scheme based off of candidates
    cost,win_length = -3,2*n+1
    h5 = histogram.normal_score(win_length,mean5,std5)
    h3 = histogram.normal_score(win_length,mean3,std3)
    """Remember that fasta index is base 1 indexing"""
    seq5 = fastaF.fetch_sequence(refID,in_start-n,in_start+n).upper()
    seq3 = fastaF.fetch_sequence(refID,in_end-n,in_end+n).upper()

    #score5 = window.score(seq5,site5p,h5,cost)
    #score3 = window.score(seq3,site3p,h3,cost)

    score5 = window.match(seq5,site5p,h5,cost)
    score3 = window.match(seq3,site3p,h3,cost)

    # print >> sys.stderr,"Site",site
    # print >> sys.stderr,"Region",in_start-n,in_start+n
    # print >> sys.stderr,"Seq left   \t",format_seq(seq5)
    # print >> sys.stderr,"Hist left  \t",format_list(h5)
    # print >> sys.stderr,"Score left \t",format_list(score5)
    # print >> sys.stderr,"Seq right  \t",format_seq(seq3)
    # print >> sys.stderr,"Hist right \t",format_list(h3)
    # print >> sys.stderr,"Score right\t",format_list(score3)
    #Find candidates in sliding window scores
    maxwin_5,score_5 = window.findSite(score5,"5")
    maxwin_3,score_3 = window.findSite(score3,"3")
    #Convert candidates into reference genome coordinates
    junc5, junc3 = maxwin_5+in_start-n-1, maxwin_3+(in_end-n-1)
    return junc5,score_5,junc3,score_3  #returned transformed coordinates of junction sites


cigar_pattern = re.compile(r"(\d+)(\S)")

_revcomp_trans = string.maketrans("ACGT", "TGCA")
def revcomp(s):
    return s[::-1].translate(_revcomp_trans)

"""
Generates a distribution and finds the most likely occurrance
"""
def findMode(sites):
    hist = counter.Counter()
    for s in sites:
        hist[s]+=1
    return (hist.most_common(1)[0])[0]

"""
This calculates the corrected position of the site given the cigar alignment
Note that read_site must be in terms of read coordinates"""
def cigar_correct(read_site,cigar,site5,site3):
    left_site,right_site = site5,site3
    align = cigar_pattern.findall(cigar)
    i,j = 0,0  #pointers in Needleman Wunsch trace back matrix.  i = ref, j = read
    for k in range(0,len(align)):
        cnt,char = int(align[k][0]),align[k][1]
        if char=="M" or char=="R":
            i,j = i+cnt,j+cnt
        elif char=="D":
            if read_site>j and read_site<j+cnt:
                diff = j+cnt-read_site
                left_site,right_site = left_site+diff,right_site-diff
            j+=cnt
        elif char=="I":
            if read_site>i and read_site<i+cnt:
                diff = i+cnt-read_site
                left_site,right_site = left_site-diff,right_site+diff
            i+=cnt
    #left_site,right_site indicate starting positions of splice site
    return left_site,right_site

"""
Applies the Needleman-Wunsch algorithm to provide a list of candiadtes
"""
def nw_correct(refID,site5,site3,introns,strand,fastaF):
    sites5,sites3 = [],[]
    M = needlemanWunsch.lcsCost()
    for intr in introns:
        in_st,in_en,lab,rdseq5_flank,rdseq3_flank,_ = intr
        rdseq5_over,rdseq3_over = rdseq3_flank,rdseq5_flank
        rdseq5 = rdseq5_flank+rdseq5_over
        rdseq3 = rdseq3_flank+rdseq3_over
        n = len(rdseq5_flank)
        overlap = len(rdseq5_over)
        st,en = site5-n,site3+n
        refseq5_flank = fastaF.fetch_sequence(refID,st-1,site5-1).upper()
        refseq3_flank = fastaF.fetch_sequence(refID,site3+3,en+3).upper()
        refseq5_over = refseq3_flank[:overlap]
        refseq3_over = refseq5_flank[-overlap:]
        refseq5 = refseq5_flank+refseq5_over
        refseq3 = refseq3_flank+refseq3_over
        _,cigar5 = needlemanWunsch.needlemanWunschXcript(refseq5,rdseq5,M)
        nsite5_1,nsite3_1 = cigar_correct(len(rdseq5_flank),cigar5,site5,site3)
        sites5.append(nsite5_1)
        sites3.append(nsite3_1+2) #Adjust for 3' end offset bias in nw alignment
    del M
    return sites5,sites3


def findHistogramLen(refID,sts,ends,radius,fastaF):
    min_st, max_st = min(sts), max(sts)
    min_end, max_end = min(ends), max(ends)
    stR = max_st-min_st
    endR = max_end-max_end

    n = max( 2*stR, 2*endR, 2*radius )
    #print >> sys.stderr,stR,endR,radius

    lengths = chrsizes.getSizes(fastaF.fasta_file+".fai")
    #Check bounds to make sure that sliding window won't over-extend genome coordinates
    if min_st-n<0 or max_end+n>lengths[refID]:
        n = 2*args.radius

    return n

"""Weighs different canonical and non-canonical sites and weighs them
Note that fw_site,rev_site and weight are zipped up in each site
"""
def findBestSite(refID,sts,ens,sites,introns,strand,fastaF):
    bs5,bs3 = 0,0  #Best sites
    bscore = 0     #Best score
    bseq = ""      #Best seq
    for s in sites:
        seq = s[0] if strand=='+' else s[1]
        w = s[2] #weight
        N = findHistogramLen(refID,sts,ens,args.radius,fastaF)
        site5,s5,site3,s3 = sliding_window(refID,sts,ens,seq,N,fastaF)
        #sites5,sites3   = nw_correct(refID,site5,site3,introns,strand,fastaF)
        #nsts,nens = sites5+list(sts),sites3+list(ens)
        #site5,s5,site3,s3 = sliding_window(refID,sites5,sites3,seq,fastaF) #Retrain using nw
        #N = findHistogramLen(refID,nsts,nens,args.radius,fastaF)
        #site5,s5,site3,s3 = sliding_window(refID,nsts,nens,seq,N,fastaF) #Retrain using nw
        if (s5+s3)*w > bscore:
            bscore = s5+s3
            bs5,bs3 = site5,site3
            bseq = seq
    return bs5,bs3,bseq,bscore

def known_noncanonical(refID,st,en):

    radius = 100
    if ((refID=='chr2R' and ( abs(st-14644850)<=radius or abs(en-14645050)<=radius)) or #False negatives
        (refID=='chr2R' and ( abs(st-2642301)<=radius or abs(en-2642501)<=radius)) or
        (refID=='chr2R' and ( abs(st-2652141)<=radius or abs(en-2652341)<=radius)) or
        (refID=='chr2R' and ( abs(st-2652938)<=radius or abs(en-2653138)<=radius)) or
        (refID=='chr3L' and ( abs(st-13433423)<=radius or abs(en-13433623)<=radius)) or
        (refID=='chr2L' and ( abs(st-20796030)<=radius or abs(en-20796031)<=radius)) or #False positives
        (refID=='chr2L' and ( abs(st-7709286)<=radius or abs(en-7709287)<=radius)) or
        (refID=='chr2R' and ( abs(st-1589725)<=radius or abs(en-1589726)<=radius)) or
        (refID=='chr2R' and ( abs(st-2652219)<=radius or abs(en-2652220)<=radius)) or
        (refID=='chr2R' and ( abs(st-2653058)<=radius or abs(en-2653059)<=radius)) or
        (refID=='chr3L' and ( abs(st-23087777)<=radius or abs(en-23087778)<=radius)) or
        (refID=='chr3L' and ( abs(st-3178206)<=radius or abs(en-3178207)<=radius)) or
        (refID=='chr3R' and ( abs(st-22302865)<=radius or abs(en-22302866)<=radius)) or
        (refID=='chr3R' and ( abs(st-25637180)<=radius or abs(en-25637181)<=radius))):
        print >> sys.stderr,"Noncanonical",st,en
        return True
    #print >> sys.stderr,"Canonical",st,en
    return False


"""
Finds canonical sites (e.g GT-AG sites)
"""
def getJunctionSites(pt,refID,bins,fastaF):
    global nout
    strand = pt[-1]
    sites5, sites3 = [],[]
    for coords,introns in bins.iteritems():
        samples = counter.Counter()
        coOccurences = defaultdict( list )
        splice_site = "GT-AG" if strand=="+" else "CT-AC"  #only consider canonical sites
        sts,ens,labs,_,_,rdids = zip(*introns)
        N = findHistogramLen(refID,sts,ens,args.radius,fastaF)

        site5,s5,site3,s3 = sliding_window(refID,sts,ens,splice_site,N,fastaF)
        #sites5,sites3   = nw_correct(refID,site5,site3,introns,strand,fastaF)
        #site5,s5,site3,s3 = sliding_window(refID,sites5,sites3,splice_site,N,fastaF) #Retrain using nw
        # threshold = 1.0

        #if s5<threshold or s3<threshold:
        splice_sites = [("GC-AC","CT-GC",1.0),
                        ("AT-AC","GT-AT",1.0)]
        nsite5,nsite3,ncseq,nc = findBestSite(refID,sts,ens,splice_sites,introns,strand,fastaF)
        site_chr = "N" if known_noncanonical(refID,nsite5,nsite3) else "C"
        if args.scores_file!="" and nc>0 and (s5+s3)>0:
            handle = open(args.scores_file,'a')
            handle.write("%lf\t%lf\t%s\n"%( (s5+s3),nc,site_chr) )

        for intr in introns:
            _,_,lab,_,_,rdid = intr
            coOccurences[rdid].append( (site5,site3) )
            samples[lab]+=1

        #Output for bed sites
        for sam,counts in samples.items():
            print "site\t%s\t%012d\t%d\t%s\t%d"%(refID,site5,site3,sam,counts)
            nout+=1

        #Output for co-occurences
        for rdid,sites in coOccurences.items():
            if len(sites)>1:
                for s in sites:
                    left_site,right_site = s
                    print "cooccurence\t%s\t%s\t%d\t%d"%(rdid,refID,left_site,right_site)

def go():

    global ninp
    starts = []  #Contains starting positions of introns
    ends = []    #Contains ending positions of introns
    labs = []    #Sample labels of introns
    rdids = []   #The read ids
    seq5_flanks,seq3_flanks = [],[]
    last_pt = "\t"
    fnh = fasta.fasta(args.refseq)
    last_ref = "\t"

    for ln in sys.stdin:
        # Parse next read
        ln = ln.rstrip()
        toks = ln.split('\t')
        assert len(toks)>=8
        pt, st, en, refid, lab, seq5_flank, seq3_flank, rdid = toks[0], int(toks[1]), int(toks[2]), toks[3], toks[4], toks[5], toks[6], toks[7]
        if last_pt=='\t':
            last_pt, last_ref = pt, refid
        elif last_pt!=pt:
            intron_ivals = zip(starts,ends,labs,seq5_flanks,seq3_flanks,rdids)
            #Cluster all introns with similar start and end positions
            #bins = cluster(intron_ivals)
            bins = diagonal_cluster(intron_ivals)
            #Apply sliding windows to find splice junction locations
            getJunctionSites(last_pt,last_ref,bins,fnh)
            starts,ends,labs = [],[],[]
            seq5_flanks,seq3_flanks = [],[]

        starts.append(st)
        ends.append(en)
        labs.append(lab)
        seq5_flanks.append(seq5_flank)
        seq3_flanks.append(seq3_flank)
        rdids.append(rdid)
        last_pt,last_ref = pt,refid
        ninp+=1

    if last_pt!='\t':
        #Handle last partition
        intron_ivals = zip(starts,ends,labs,seq5_flanks,seq3_flanks,rdids)
        #Cluster all introns with similar start and end positions
        #bins = cluster(intron_ivals)
        bins = diagonal_cluster(intron_ivals)
        #Apply sliding windows to find splice junction locations
        getJunctionSites(last_pt,last_ref,bins,fnh)

    # Done
    timeEn = time.clock()
    print >>sys.stderr, "DONE with intron.py; in/out = %d/%d; time=%0.3f secs" % (ninp, nout, timeEn-timeSt)

def createTestFasta(fname,refid,refseq):
    fastaH = open(fname,'w')
    fastaIdx = open(fname+".fai",'w')
    fastaH.write(">%s\n%s\n"%(refid,refseq))
    fastaIdx.write("%s\t%d\t%d\t%d\t%d\n"%(refid,len(refseq),len(refid)+2,len(refseq),len(refseq)+1))
    fastaH.close()
    fastaIdx.close()


if not args.test:
    if args.scores_file!="":
        handle = open(args.scores_file,'w')
    go()
else:
    del sys.argv[1:]
    import unittest


    class TestIntronFunctions1(unittest.TestCase):
        def setUp(self):
            self.refseq="""TCGATGTCGATGGGTCCAAGCTGCTCAAATATCCCGCTGCCGGTGGATGCAACACCGGGTCCCCCTTGCAGCCAGATTACCAGCGGCCTCTCTATGAAATGAGATACATTGGCAGTGGTGTATAGAAGCCAGTAAAAGAGGTGAGCGCCCTTCCGAACTTCCACATAGTCCCATTCCTGTACTCCAGGTCCCAGACCAACACGTCCTGCAACGAAATAACTAAGACTTTTGGAGTATTTCTCTCAAACATCGAAACTTATAAATGACCCCATTTAGTAGATTTTAATTAACCTCAATATGGCAACCACAACTACGCCATTTTTTTCACTTTGGTAACCATACCACATTTATGTCTCAGAAAACGTACACACCTTGCACGCAGATCAGTGATAAAAAAAAAC"""

            self.rightseqs = ["TCCTTGCACG","TGCACGCAGA","TGCACGCAGA","TCCTTGCACG","TTGCACGCAG","CCTTGCACGC","TTGCACGCAG","TCCTTGCACG","TGCACGCAGA","CTTGCACGCA","TTGCACGCAG","CTTGCACGCA","CCTTGCACGC","TTGCACGCAG","TCCTTGCACG","CTTGCACGCA","CTTGCACGCA","TTGCACGCAG","CCTTGCACGC","TTGCACGCAG","CCTTGCACGC","TGCACGCAGA","CTTGCACGCA","TGCACGCAGA","TTGCACGCAG","TTGCACGCAG","TCCTTGCACG","TCCTTGCACG","TGCACGCAGA","CTTGCACGCA","CCTTGCACGC","CTTGCACGCA","CTTGCACGCA","TGCACGCAGA","TGCACGCAGA","CTTGCACGCA","TTGCACGCAG","CCTTGCACGC","CCTTGCACGC","TCCTTGCACG","TGCACGCAGA","TCCTTGCACG","TCCTTGCACG","CCTTGCACGC","CTTGCACGCA","TGCACGCAGA","CCTTGCACGC","TTGCACGCAG","TGCACGCAGA","TGCACGCAGA","TTGCACGCAG","TCCTTGCACG","CTTGCACGCA","CTTGCACGCA","CCTTGCACGC","CCTTGCACGC","TCCTTGCACG","TGCACGCAGA","TTGCACGCAG","TCCTTGCACG","TCCTTGCACG","TGCACGCAGA","CCTTGCACGC","TTGCACGCAG","CTTGCACGCA","TCCTTGCACG","CTTGCACGCA","TTGCACGCAG","CCTTGCACGC","CCTTGCACGC","TGCACGCAGA","TCCTTGCACG","TTGCACGCAG","CTTGCACGCA","TCCTTGCACG","CCTTGCACGC","TTGCACGCAG","TTGCACGCAG","CCTTGCACGC","TGCACGCAGA","TTGCACGCAG","TTGCACGCAG","CCTTGCACGC","TCCTTGCACG","TGCACGCAGA","CTTGCACGCA","CCTTGCACGC","CCTTGCACGC","TCCTTGCACG","TTGCACGCAG","TGCACGCAGA","TCCTTGCACG","TGCACGCAGA","TGCACGCAGA","TCCTTGCACG","TTGCACGCAG","CTTGCACGCA","CCTTGCACGC","CTTGCACGCA","TGCACGCAGA","CTTGCACGCA","TTGCACGCAG","TGCACGCAGA","TGCACGCAGA","CTTGCACGCA","CCTTGCACGC","TCCTTGCACG","TCCTTGCACG","CCTTGCACGC","TGCACGCAGA","CTTGCACGCA","CTTGCACGCA","TGCACGCAGA","TGCACGCAGA","CTTGCACGCA","CTTGCACGCA"]

            self.leftseqs  = ["GACCAACACG","AACACGTCCT","AACACGTCCT","GACCAACACG","CAACACGTCC","ACCAACACGT","CAACACGTCC","GACCAACACG","AACACGTCCT","CCAACACGTC","CAACACGTCC","CCAACACGTC","ACCAACACGT","CAACACGTCC","GACCAACACG","CCAACACGTC","CCAACACGTC","CAACACGTCC","ACCAACACGT","CAACACGTCC","ACCAACACGT","AACACGTCCT","CCAACACGTC","AACACGTCCT","CAACACGTCC","CAACACGTCC","GACCAACACG","GACCAACACG","AACACGTCCT","CCAACACGTC","ACCAACACGT","CCAACACGTC","CCAACACGTC","AACACGTCCT","AACACGTCCT","CCAACACGTC","CAACACGTCC","ACCAACACGT","ACCAACACGT","GACCAACACG","AACACGTCCT","GACCAACACG","GACCAACACG","ACCAACACGT","CCAACACGTC","AACACGTCCT","ACCAACACGT","CAACACGTCC","AACACGTCCT","AACACGTCCT","CAACACGTCC","GACCAACACG","CCAACACGTC","CCAACACGTC","ACCAACACGT","ACCAACACGT","GACCAACACG","AACACGTCCT","CAACACGTCC","GACCAACACG","GACCAACACG","AACACGTCCT","ACCAACACGT","CAACACGTCC","CCAACACGTC","GACCAACACG","CCAACACGTC","CAACACGTCC","ACCAACACGT","ACCAACACGT","AACACGTCCT","GACCAACACG","CAACACGTCC","CCAACACGTC","GACCAACACG","ACCAACACGT","CAACACGTCC","CAACACGTCC","ACCAACACGT","AACACGTCCT","CAACACGTCC","CAACACGTCC","ACCAACACGT","GACCAACACG","AACACGTCCT","CCAACACGTC","ACCAACACGT","ACCAACACGT","GACCAACACG","CAACACGTCC","AACACGTCCT","GACCAACACG","AACACGTCCT","AACACGTCCT","GACCAACACG","CAACACGTCC","CCAACACGTC","ACCAACACGT","CCAACACGTC","AACACGTCCT","CCAACACGTC","CAACACGTCC","AACACGTCCT","AACACGTCCT","CCAACACGTC","ACCAACACGT","GACCAACACG","GACCAACACG","ACCAACACGT","AACACGTCCT","CCAACACGTC","CCAACACGTC","AACACGTCCT","AACACGTCCT","CCAACACGTC","CCAACACGTC"]
            self.fasta = "test.fa"
            self.faidx = "test.fa.fai"
            createTestFasta(self.fasta,"test",self.refseq)

        def tearDown(self):
            os.remove(self.fasta)
            os.remove(self.faidx)

        def test_nw_correct1(self):
            print >> sys.stderr,"NW test 1"


            n = len(self.leftseqs)
            sts,ends,labs,rdids = [205]*n, [371]*n, ["test_labs"]*n, map( str, range(0,n))
            fnh = fasta.fasta("test.fa")

            refID, splice_site, strand= "test","CT-AC","-"
            left_site,_,right_site,_ = sliding_window(refID,sts,ends,splice_site,args.radius,fnh)
            # print "left site",left_site,205
            # print "right site",right_site,369
            self.assertEquals( left_site, 205)
            self.assertEquals( right_site, 369)
            print >> sys.stderr,"Sliding window test passed !"
            introns = zip(sts,ends,labs,self.leftseqs,self.rightseqs,rdids)
            sites5,sites3   = nw_correct(refID,left_site,right_site,introns,strand,fnh)
            left_site,_,right_site,_ = sliding_window(refID,sites5,sites3,splice_site,args.radius,fnh)

            # print "left site ",left_site,205
            # print "left histogram ",sites5
            # print "right site",right_site,369
            # print "right histogram",sites3

            self.assertEquals( left_site, 205)
            self.assertEquals( right_site, 369)
            print >> sys.stderr,"Needleman Wunsch test passed ! \n"

    class TestIntronFunctions2(unittest.TestCase):
        def setUp(self):
            self.leftseqs =["CGTGTGCACG","CGTGTGCACG","CGTGTGCACG","CGTGTGCACG","CGTGTGCACG","CGTGTGCACG","CGTGTGCACG","CGTGTGCACG","CGTGTGCACG","CGTGTGCACG","CGTGTGCACG","CGTGTGCACG","CGTGTGCACG","CGTGTGCACG","CGTGTGCACG","CGTGTGCACG","CGTGTGCACG","CGTGTGCACG","CGTGTGCACG","CGTGTGCACG","CGTGTGCACG","CGTGTGCACG","CGTGTGCACG","CGTGTGCACG","CGTGTGCACG","CGTGTGCACG","CGTGTGCACG","CGTGTGCACG","CGTGTGCACG","CGTGTGCACG","CGTGTGCACG","CGTGTGCACG","CGTGTGCACG","CGTGTGCACG","CGTGTGCACG","CGTGTGCACG","CGTGTGCACG","CGTGTGCACG","CGTGTGCACG","CGTGTGCACG","CGTGTGCACG","CGTGTGCACG","CGTGTGCACG","CGTGTGCACG","CGTGTGCACG","CGTGTGCACG","CGTGTGCACG","CGTGTGCACG","CGTGTGCACG","CGTGTGCACG"]

            self.rightseqs=["AAACGCCCAG","AAACGCCCAG","AAACGCCCAG","AAACGCCCAG","AAACGCCCAG","AAACGCCCAG","AAACGCCCAG","AAACGCCCAG","AAACGCCCAG","AAACGCCCAG","AAACGCCCAG","AAACGCCCAG","AAACGCCCAG","AAACGCCCAG","AAACGCCCAG","AAACGCCCAG","AAACGCCCAG","AAACGCCCAG","AAACGCCCAG","AAACGCCCAG","AAACGCCCAG","AAACGCCCAG","AAACGCCCAG","AAACGCCCAG","AAACGCCCAG","AAACGCCCAG","AAACGCCCAG","AAACGCCCAG","AAACGCCCAG","AAACGCCCAG","AAACGCCCAG","AAACGCCCAG","AAACGCCCAG","AAACGCCCAG","AAACGCCCAG","AAACGCCCAG","AAACGCCCAG","AAACGCCCAG","AAACGCCCAG","AAACGCCCAG","AAACGCCCAG","AAACGCCCAG","AAACGCCCAG","AAACGCCCAG","AAACGCCCAG","AAACGCCCAG","AAACGCCCAG","AAACGCCCAG","AAACGCCCAG","AAACGCCCAG"]

            self.refseq = "CGACGACACCGACGACGCCAAAGTTGCCACAGGAAACGGAAATCTGAGCGTGTGCACGTGTGTGTGTGCGCGCACATGGCGTTCATATTTATTTATTTCTTTTTCGGTACAGGAAACGCCCAGCAGGATTAAGAATGGAGTAGTCTTGTGACCATCGGGAACTTTTCGGGGGACAGCCATAAGTGTCAAGACTTAAAGCTG"
            self.fasta = "test.fa"
            self.faidx = "test.fa.fai"
            createTestFasta(self.fasta,"test",self.refseq)
        def tearDown(self):
            os.remove(self.fasta)
            os.remove(self.faidx)

        def test_nw_correct2(self):
            print >> sys.stderr,"NW test 2"

            n = len(self.leftseqs)
            sts,ends,labs,rdids = [57]*n, [112]*n,["test_labs"]*n, map( str, range(0,n))

            fnh = fasta.fasta("test.fa")
            refID, splice_site,strand= "test","GT-AG",'+'
            left_site,_,right_site,_ = sliding_window(refID,sts,ends,splice_site,args.radius,fnh)
            # print "left site",left_site,57
            # print "right site",right_site,111
            self.assertEquals( left_site,57)
            self.assertEquals( right_site,110)
            print >> sys.stderr,"Sliding window test passed !"
            introns = zip(sts,ends,labs,self.leftseqs,self.rightseqs,rdids)
            sites5,sites3   = nw_correct(refID,left_site,right_site,introns,strand,fnh)
            left_site,_,right_site,_ = sliding_window(refID,sites5,sites3,splice_site,args.radius,fnh)

            # print "left site ",left_site,57
            # print "left histogram ",sites5
            # print "right site",right_site,110
            # print "right histogram",sites3
            # print "left  seq",self.refseq[:left_site]
            # print "right seq",self.refseq[right_site:]
            self.assertEquals( left_site,57)
            self.assertEquals( right_site,110)
            print >> sys.stderr,"Needleman Wunsch test passed ! \n"

    class TestIntronFunctions3(unittest.TestCase):
        def setUp(self):
            self.leftseqs =["GCGTGTGCAC","CGTGTGCACG"]
            self.rightseqs=["GAAACGCCCA","AAACGCCCAG"]
            "                                                    CG TGTGCACG                                                             AAACGCC CAG"
            "                                                   GCG TGTGCAC                                                             GAAACGCC CA"
            "CGACGACACC GACGACGCCA AAGTTGCCAC AGGAAACGGA AATCTGAGCG TGTGCACGTG TGTGTGTGCG CGCACATGGC GTTCATATTT ATTTATTTCT TTTTCGGTAC AGGAAACGCC CAGCAGGATT AAGAATGGAG TAGTCTTGTG ACCATCGGGA ACTTTTCGGG GGACAGCCAT AAGTGTCAAG ACTTAAAGCT G"
            self.refseq = "CGACGACACCGACGACGCCAAAGTTGCCACAGGAAACGGAAATCTGAGCGTGTGCACGTGTGTGTGTGCGCGCACATGGCGTTCATATTTATTTATTTCTTTTTCGGTACAGGAAACGCCCAGCAGGATTAAGAATGGAGTAGTCTTGTGACCATCGGGAACTTTTCGGGGGACAGCCATAAGTGTCAAGACTTAAAGCTG"
            self.fasta = "test.fa"
            self.faidx = "test.fa.fai"
            createTestFasta(self.fasta,"test",self.refseq)

        def tearDown(self):
            os.remove(self.fasta)
            os.remove(self.faidx)

        def test_noncanonical1(self):
            print >> sys.stderr,"Non canonical 1"

            n = len(self.leftseqs)
            sts,ends,labs,rdids = [57]*n, [112]*n,["test_labs"]*n, map( str, range(0,n))

            splice_sites = [("GT-AG","CT-AC",1.0),
                            ("GC-AC","CT-GC",0.0),
                            ("AT-AC","GT-AT",0.0)]


            fnh = fasta.fasta("test.fa")
            refID,strand= "test",'+'
            introns = zip(sts,ends,labs,self.leftseqs,self.rightseqs,rdids)
            left_site,right_site,seq,_ = findBestSite(refID,sts,ends,splice_sites,introns,strand,fnh)
            self.assertEquals( left_site,57)
            self.assertEquals( right_site,110)
            self.assertEquals( seq,"GT-AG")

    class TestIntronFunctions4(unittest.TestCase):
        def setUp(self):
            self.leftseqs =["GCGTGTGCAC","CGTGTGCACC"]
            self.rightseqs=["GAAACGCCCA","AAACGCCCAG"]
            "                                                    CG TGTGCACG                                                             AAACGCC CAG"
            "                                                   GCG TGTGCAC                                                             GAAACGCC CA "
            "CGACGACACC GACGACGCCA AAGTTGCCAC AGGAAACGGA AATCTGAGCG TGTGCACCTG TGTGTGTGCG CGCACATGGC GTTCATATTT ATTTATTTCT TTTTCGGTAC ACGAAACGCC CAGCAGGATT AAGAATGGAG TAGTCTTGTG ACCATCGGGA ACTTTTCGGG GGACAGCCAT AAGTGTCAAG ACTTAAAGCT G"
            self.refseq = "CGACGACACCGACGACGCCAAAGTTGCCACAGGAAACGGAAATCTGAGCGTGTGCACCTGTGTGTGTGCGCGCACATGGCGTTCATATTTATTTATTTCTTTTTCGGTACACGAAACGCCCAGCAGGATTAAGAATGGAGTAGTCTTGTGACCATCGGGAACTTTTCGGGGGACAGCCATAAGTGTCAAGACTTAAAGCTG"
            self.fasta = "test.fa"
            self.faidx = "test.fa.fai"
            createTestFasta(self.fasta,"test",self.refseq)

        def tearDown(self):
            os.remove(self.fasta)
            os.remove(self.faidx)

        def test_noncanonical2(self):
            print >> sys.stderr,"Non canonical 2"

            n = len(self.leftseqs)
            sts,ends,labs,rdids = [57]*n, [112]*n,["test_labs"]*n, map( str, range(0,n))

            splice_sites = [("GT-AG","CT-AC",1.0),
                            ("GC-AC","CT-GC",0.0),
                            ("AT-AC","GT-AT",0.0)]

            fnh = fasta.fasta("test.fa")
            refID,strand= "test",'-'
            introns = zip(sts,ends,labs,self.leftseqs,self.rightseqs,rdids)
            left_site,right_site,seq,_ = findBestSite(refID,sts,ends,splice_sites,introns,strand,fnh)
            self.assertEquals( left_site,57)
            self.assertEquals( right_site,110)
            self.assertEquals( seq,"CT-AC")


    class TestIntronFunctions5(unittest.TestCase):
        def setUp(self):
            self.leftseqs =["GCGTGTGCAC","CGTGTGCACG"]
            self.rightseqs=["GAAACGCCCA","AAACGCCCAG"]
            "                                                    CG TGTGCACG                                                             AAACGCC CAG"
            "                                                   GCG TGTGCAC                                                             GAAACGCC CA"
            "CGACGACACC GACGACGCCA AAGTTGCCAC AGGAAACGGA AATCTGAGCG CGCGCACGCG CGCGCGCGCG CGCACATGGC GTTCATATTT ATTTATTTCT TTTTCGGTAC ACGAAACGCC CAGCAGGATT AAGAATGGAG TAGTCTTGTG ACCATCGGGA ACTTTTCGGG GGACAGCCAT AAGTGTCAAG ACTTAAAGCT G"

            self.refseq = "CGACGACACCGACGACGCCAAAGTTGCCACAGGAAACGGAAACCCGAGCGCGCGCACGCGCGCGCGCGCGCGCACATGGCGTTCATATTTATTTATTTCTTTTTCGGTACACGAAACGCCCAGCAGGATTAAGAATGGAGTAGTCTTGTGACCATCGGGAACTTTTCGGGGGACAGCCATAAGTGTCAAGACTTAAAGCTG"

            self.fasta = "test.fa"
            self.faidx = "test.fa.fai"
            createTestFasta(self.fasta,"test",self.refseq)
        def tearDown(self):
            os.remove(self.fasta)
            os.remove(self.faidx)
        def test_noncanonical3(self):
            print >> sys.stderr,"Non canonical 3"
            n = len(self.leftseqs)
            sts,ends,labs,rdids = [57]*n, [112]*n,["test_labs"]*n, map( str, range(0,n))
            splice_sites = [("GT-AG","CT-AC",1.0),
                            ("GC-AC","CT-GC",1.0),
                            ("AT-AC","GT-AT",1.0)]
            fnh = fasta.fasta("test.fa")
            refID,strand= "test",'+'
            introns = zip(sts,ends,labs,self.leftseqs,self.rightseqs,rdids)
            left_site,right_site,seq,_ = findBestSite(refID,sts,ends,splice_sites,introns,strand,fnh)
            self.assertEquals( left_site,57)
            self.assertEquals( right_site,110)
            self.assertEquals( seq,"GC-AC")


    class TestIntronFunctions6(unittest.TestCase):
        def setUp(self):
            self.leftseqs =["GCGTGTGCAC","CGTGTGCACG"]
            self.rightseqs=["GAAACGCCCA","AAACGCCCAG"]
            "                                                    CG TGTGCACG                                                             AAACGCC CAG"
            "                                                   GCG TGTGCAC                                                             GAAACGCC CA"
            "CGACGACACC GACGACGCCA AAGTTGCCAC AGGAAACGGA AATCTGAGCG TGTGCACGTG TGTGTGTGCG CGCACATGGC GTTCATATTT ATTTATTTCT TTTTCGGTAC AGGAAACGCC CAGCAGGATT AAGAATGGAG TAGTCTTGTG ACCATCGGGA ACTTTTCGGG GGACAGCCAT AAGTGTCAAG ACTTAAAGCT G"
            self.refseq = "CGACGACACCGACGACGCCAAAGTTGCCACAGGAAACGGAAATCTGAGCGTGTGCACGTGTGTGTGTGCGCGCACATGGCGTTCATATTTATTTATTTCTTTTTCGGTACAGGAAACGCCCAGCAGGATTAAGAATGGAGTAGTCTTGTGACCATCGGGAACTTTTCGGGGGACAGCCATAAGTGTCAAGACTTAAAGCTG"
            self.fasta = "test.fa"
            self.faidx = "test.fa.fai"
            createTestFasta(self.fasta,"test",self.refseq)

        def tearDown(self):
            os.remove(self.fasta)
            os.remove(self.faidx)

        def test_noncanonical4(self):
            print >> sys.stderr,"Non canonical 4"

            n = len(self.leftseqs)
            sts,ends,labs,rdids = [57]*n, [112]*n,["test_labs"]*n, map( str, range(0,n))

            splice_sites = [("GT-AG","CT-AC",1.0),
                            ("GC-AC","CT-GC",1.0),
                            ("AT-AC","GT-AT",1.0)]


            fnh = fasta.fasta("test.fa")
            refID,strand= "test",'+'
            introns = zip(sts,ends,labs,self.leftseqs,self.rightseqs,rdids)
            left_site,right_site,seq,_ = findBestSite(refID,sts,ends,splice_sites,introns,strand,fnh)
            self.assertEquals( left_site,57)
            self.assertEquals( right_site,110)
            self.assertEquals( seq,"GT-AG")

    unittest.main()

