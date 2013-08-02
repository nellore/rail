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

Tab-delimited output tuple columns:
1. Reference ID
2. 5' start
3. 3' start
4. Sample label
5. Read frequency (number of times sample read overlapped junction)


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

import fasta
import readlet
import sw
import nw
import histogram
import counter

parser = argparse.ArgumentParser(description=\
                                     'Reports splice junction information')

parser.add_argument(\
    '--refseq', type=str, required=False,
    help='The fasta sequence of the reference genome. The fasta index of the reference genome is also required')
parser.add_argument(\
    '--radius', type=int, required=False,default=10,
    help='The radius of the clustering algorithm and the radius of the sliding window algorithm.')


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
    key = "%s,%s"%(p[0],p[1])
    bins[key].append(ivals[0])
    notFound = True
    for i in range(1,len(ivals)):
        for j in range(0,len(points)): #Check all the neighborhood of all points
            p = points[j]
            ival = ivals[i]
            if (ival[0]>(p[0]-rIval) and ival[0]<(p[0]+rIval) and
                ival[1]>(p[1]-rIval) and ival[1]<(p[1]+rIval) ):
                key = "%s,%s"%(p[0],p[1])
                bins[key].append(ivals[i])
                notFound = False
                break
        if notFound:
            p = ivals[i]
            points.append(p)
            key = "%s,%s"%(p[0],p[1])
            bins[key].append(ivals[i])
        notFound = True
    return bins


"""
Scores a set of windows based off of splice site
"""
def score(seq, site, hist):
    wsize = len(site) # window size
    nwins = len(seq)-wsize+1
    wins = [0]*nwins

    for i in range(0,nwins):
        for j in range(0,len(site)):
            s = 1 if site[j]==seq[i+j] else -3
            wins[i]+=s*hist[i+j]
    return wins


"""
Returns the site by finding the maximum in the scores
To break ties it uses the direction.
If direction=="5", that means its a 5' end and it will return the score closest to the 5' end (aka. left)
The vice versa happens with direction=="3"

Note that this just returns offsets wrt to window frame
"""
def findSite(scores,direction):
    count = -1 if direction=="5" else 1
    i = len(scores)-1 if direction=="5" else 0
    m, ind = -1, -1
    while i>=0 and i<len(scores):
        if m < scores[i]:
            ind = i
            m = scores[i]
        i+=count
    return ind,scores[ind]

"""
Just a fancier way to print out lists
"""
def format_list(L):
    s = ""
    for i in L:
        s+="%.3f "%i
    return s

"""
Note: site is formatted as follows: XX-XX (e.g. GT-AG)
Returns the 5' and 3' splice sites within multiple intervals
"""
def sliding_window(refID, sts,ens, site, fastaF):
    n,r = 2*args.radius, args.radius
    in_start, in_end = min(sts),max(ens)
    toks = site.split("-")
    assert len(toks)==2
    site5p,site3p = toks[0],toks[1]
    hist5, hist3 = histogram.hist_score(sts,in_start,"5",2*n+1), histogram.hist_score(sts,in_start,"3",2*n+1)
    #mean5,std5 = histogram.average(hist5),2*histogram.stddev(hist5)
    #mean3,std3 = histogram.average(hist3),2*histogram.stddev(hist3)
    mean5,std5 = hist5.index(max(hist5)),2*histogram.stddev(hist5)
    mean3,std3 = hist3.index(max(hist3)),2*histogram.stddev(hist3)
    #Create a normal distributed scoring scheme based off of candidates
    h5,h3 = histogram.normal_score(2*n+1,mean5,std5), histogram.normal_score(2*n+1,mean5,std5)
    """Remember that fasta index is base 1 indexing"""
    seq5 = fastaF.fetch_sequence(refID,in_start-n,in_start+n)
    seq3 = fastaF.fetch_sequence(refID,in_end-n,in_end+n)
    score5,score3 = score(seq5,site5p,h5),score(seq3,site3p,h3)
    j5,s5 = findSite(score5,"5")
    j3,s3 = findSite(score3,"3")
    # print >> sys.stderr,"Seq 5",seq5
    # print >> sys.stderr,"Seq 3",seq3
    # print >> sys.stderr,"Histogram 5\t",format_list(hist5)
    # print >> sys.stderr,"Histogram 3\t",format_list(hist3)
    # print >> sys.stderr,"Score 5 \t",format_list(h5)
    # print >> sys.stderr,"Score 3 \t",format_list(h3)
    # print >> sys.stderr, "Mean 5\t",histogram.average(hist5),"Mode 5",hist5.index(max(hist5))
    # print >> sys.stderr, "Mean 3\t",histogram.average(hist3),"Mode 3",hist5.index(max(hist3))
    return j5+in_start-n-1,s5,j3+(in_end-n-1),s3  #returned transformed coordinates of junction sites


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
Note that read_site must be in terms of read coordinates
"""
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
Applies the Needleman-Wunsch algorithm to provide a list of candiates
"""
def nw_correct(refID,site5,site3,introns,strand,fastaF):
    sites5,sites3 = [],[]
    for intr in introns:
        in_st,in_en,lab,rdseq5_flank,rdseq5_over,rdseq3_flank,rdseq3_over = intr
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
        _,cigar5 = nw.c_needlemanWunschXcript(refseq5,rdseq5,nw.lcsCost)
        nsite5_1,nsite3_1 = cigar_correct(len(rdseq5_flank),cigar5,site5,site3)
        sites5.append(nsite5_1)
        sites3.append(nsite3_1)
    return sites5,sites3

"""
Finds canonical sites (e.g GT-AG sites)
"""
def getJunctionSites(pt,refID,bins,fastaF):
    global nout
    strand = pt[-1]
    samples = counter.Counter()
    sites5, sites3 = [],[]
    for coords,introns in bins.iteritems():
        splice_site = "GT-AG" if strand=="+" else "CT-AC"  #only consider canonical sites
        sts,ens,labs,_,_,_,_ = zip(*introns)
        site5,_,site3,_ = sliding_window(refID,sts,ens,splice_site,fastaF)
        sites5,sites3 = nw_correct(refID,site5,site3,introns,strand,fastaF)
        site5,_,site3,_ = sliding_window(refID,sites5,sites3,splice_site,fastaF) #Retrain using nw
        for intr in introns:
            lab = intr[2]
            samples[lab]+=1
        for sam,counts in samples.items():
            print "%s\t%012d\t%d\t%s\t%d"%(refID,site5,site3,sam,counts)
            nout+=1

starts = []  #Contains starting positions of introns
ends = []    #Contains ending positions of introns
labs = []    #Sample labels of introns
seq5_flanks,seq3_flanks,seq5_overs,seq3_overs = [],[],[],[]
last_pt = "\t"
fnh = fasta.fasta(args.refseq)
last_ref = "\t"

for ln in sys.stdin:
    # Parse next read
    ln = ln.rstrip()
    toks = ln.split('\t')
    assert len(toks)>=9
    pt, st, en, refid, lab, seq5_flank, seq5_over, seq3_flank, seq3_over = toks[0], int(toks[1]), int(toks[2]), toks[3], toks[4], toks[5], toks[6], toks[7], toks[8]
    if last_pt=='\t':
        last_pt, last_ref = pt, refid
    elif last_pt!=pt:
        intron_ivals = zip(starts,ends,labs,seq5_flanks,seq5_overs,seq3_flanks,seq3_overs)
        #Cluster all introns with similar start and end positions
        bins = cluster(intron_ivals)
        #Apply sliding windows to find splice junction locations
        getJunctionSites(last_pt,last_ref,bins,fnh)
        starts,ends,labs = [],[],[]
        seq5_flanks,seq3_flanks,seq5_overs,seq3_overs = [],[],[],[]

    starts.append(st)
    ends.append(en)
    labs.append(lab)
    seq5_flanks.append(seq5_flank)
    seq3_flanks.append(seq3_flank)
    seq5_overs.append(seq5_over)
    seq3_overs.append(seq3_over)
    last_pt,last_ref = pt,refid
    ninp+=1

if last_pt!='\t':
    #Handle last partition
    intron_ivals = zip(starts,ends,labs,seq5_flanks,seq5_overs,seq3_flanks,seq3_overs)
    #Cluster all introns with similar start and end positions
    bins = cluster(intron_ivals)
    #Apply sliding windows to find splice junction locations
    getJunctionSites(last_pt,last_ref,bins,fnh)

# Done
timeEn = time.clock()
print >>sys.stderr, "DONE with intron.py; in/out = %d/%d; time=%0.3f secs" % (ninp, nout, timeEn-timeSt)
