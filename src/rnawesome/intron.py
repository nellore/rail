import math
import numpy as np
"""                               
Tab-delimited input tuple columns:  
1. Partition ID for partition overlapped by interval 
2. Interval start  
3. Interval end (exclusive) 
4. Reference ID    
5. Sample label     

Tab-delimited output tuple columns:           
1. Reference ID
2. 5' start
3. 3' start
4. Sample label
5. Read frequency (number of times sample read overlapped junction)

Questions to consider
1) Should the splice junction sites be scored according to a histogram?
2) How much can be assume about intron coverage?  This will break if there aren't enough spanning introns
"""
import os
import sys
import argparse
import site
import time
from collections import defaultdict
from collections import Counter
timeSt = time.clock()

base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
site.addsitedir(os.path.join(base_path, "fasta"))
site.addsitedir(os.path.join(base_path, "read"))
site.addsitedir(os.path.join(base_path, "alignment"))

import fasta
import readlet
import sw

parser = argparse.ArgumentParser(description=\
                                     'Reports splice junction information')

parser.add_argument(\
    '--refseq', type=str, required=False,
    help='The fasta sequence of the reference genome. The fasta index of the reference genome is also required')
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
    rIval = args.readletIval
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
            s = 1 if site[j]==seq[i+j] else -2
            wins[i]+=s*hist[i+j]
    return wins            

"""
Scores based off of a histogram of boundary locations
"""
def hist_score(coords,offset,endtype,N):
    #offset is placed in the middle of the histogram
    hist = [1.0]*N #pseudo counts
    n = N/2
    for c in coords:
        if abs(offset-c)>N:
            print>>sys.error,"Out of bounds coordinate"
            continue
        ind = (c-offset)+n if endtype=="5" else (N-(offset-c)-1)-n
        hist[ind]+=1
    total = sum(hist)
    hist = [h/float(total) for h in hist]
    return hist

"""
Assigns scores based off of how close it is to the end based off of a p=2 series
Sites closer to the ends get higher scores
"""
def exp_score(endtype,N):
    if endtype=="5":
        hist = [.5**n for n in range(0,N)]
    else: 
        hist = [.5**(N-n+1) for n in range(0,N)]
    return hist

"""
Assigns scores based off of a harmonic series (aka. logarithmic score)
"""
def harmonic_score(endtype,N):
    if endtype=="5":
        hist = [1.0/(n+1) for n in range(0,N)]
    else: 
        hist = [1.0/(N-n) for n in range(0,N)]
    return hist
"""
Assigns scores based off of a normal distribution with mean=center of window and std=window_size/2
"""
def normal(x,m,s):
    return (1.0/(math.sqrt(2*math.pi*s*s)))*math.exp(-(x-m)*(x-m)/(2*s*s))

"""
Score using a normal distribution s.t. positions in the center will be weighted higher
"""
def normal_score(N):
    m = N/2
    s = m
    hist = [normal(i,m,s) for i in range(0,N)]
    return hist
    
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
def print_list(L):
    s = ""
    for i in L:
        s+="%.2f "%i
    return s
"""
Note: site is formatted as follows: XX-XX (e.g. GT-AG)
Returns the 5' and 3' splice sites within multiple intervals
"""
def sliding_window(refID, ivals, site, fastaF):
    n,r = 2*args.readletIval, args.readletIval
    sts,ens,labs,seqs = zip(*ivals)
    in_start, in_end = min(sts),max(ens)
    toks = site.split("-")
    assert len(toks)==2
    site5p,site3p = toks[0],toks[1]
    h5,h3 = normal_score(2*n+1),normal_score(2*n+1)
    """Remember that fasta index is base 1 indexing"""
    seq5 = fastaF.fetch_sequence(refID,in_start-n,in_start+n)
    seq3 = fastaF.fetch_sequence(refID,in_end-n,in_end+n)
    score5,score3 = score(seq5,site5p,h5),score(seq3,site3p,h3)
    j5,s5 = findSite(score5,"5")
    j3,s3 = findSite(score3,"3")
    return j5+in_start-n-1,s5,j3+(in_end-n-1),s3  #returned transformed coordinates of junction sites

"""
Applies the Smith-Waterman algorithm to correct the initial splice site estimate
"""
def sw_correct(refID,site5,site3,introns):
    starts,ends,labs,seqs = zip(*introns)
    for rdseq in seqs:
        n = len(s)/2
        st = site5-n
        en = site3+n
        seq5 = fastaF.fetch_sequence(refID,st,site5)
        seq3 = fastaF.fetch_sequence(refID,site3,en)
        score5,cigar5 = sw.smithWatermanXcript(seq5,rdseq)
        score3,cigar3 = sw.smithWatermanXcript(seq5,rdseq)
        
        
"""
Finds canonical sites (e.g GT-AG sites)
"""
def getJunctionSites(refID,bins,fastaF):
    global nout
    samples = Counter()
    sites5, sites3 = [],[]
    for coords,introns in bins.iteritems():
        site5,score5,site3,score3 = sliding_window(refID,introns,"GT-AG",fastaF)
        site5r,score5r,site3r,score3r = sliding_window(refID,introns,"CT-AC",fastaF) #Reverse complement
        if score5r>score5 and score3r>score3:
            site5,site3 = site5r,site3r

        
        for intr in introns:
            lab = intr[2]
            samples[lab]+=1
        for sam,counts in samples.items():
            print "%s\t%012d\t%d\t%s\t%d"%(refID,site5,site3,sam,counts)
            nout+=1

starts = []  #Contains starting positions of introns
ends = []    #Contains ending positions of introns
labs = []    #Sample labels of introns
seqs = []
last_pt = "\t"
fnh = fasta.fasta(args.refseq)
last_ref = "\t"

for ln in sys.stdin:
    # Parse next read
    ln = ln.rstrip()
    toks = ln.split('\t')
    assert len(toks)>=5
    pt, st, en, refid, lab, seq = toks[0], int(toks[1]), int(toks[2]), toks[3], toks[4], toks[5]
    if last_pt=='\t':
        last_pt, last_ref = pt, refid
    elif last_pt!=pt:
        intron_ivals = zip(starts,ends,labs,seqs)
        #Cluster all introns with similar start and end positions   
        bins = cluster(intron_ivals)
        #Apply sliding windows to find splice junction locations
        getJunctionSites(last_ref,bins,fnh)

        starts,ends,labs = [],[],[]
        
    starts.append(st)
    ends.append(en)
    labs.append(lab)
    seqs.append(seq)
    last_pt,last_ref = pt,refid
    ninp+=1

#Handle last partition
intron_ivals = zip(starts,ends,labs)
bins = cluster(intron_ivals)
getJunctionSites(last_ref,bins,fnh)

# Done                                                                                                                      
timeEn = time.clock()
print >>sys.stderr, "DONE with intron.py; in/out = %d/%d; time=%0.3f secs" % (ninp, nout, timeEn-timeSt)
