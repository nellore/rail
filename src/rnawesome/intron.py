import math

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
"""
import os
import sys
import argparse
import site
import time
from collections import defaultdict
timeSt = time.clock()

base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
site.addsitedir(os.path.join(base_path, "fasta"))
site.addsitedir(os.path.join(base_path, "read"))

import fasta
import readlet

parser = argparse.ArgumentParser(description=\
                                     'Reports splice junction information')

# parser.add_argument(\
#     '--readletIval', type=int, required=False,
#     help='If readlets are desired, interval between readlet starts')

parser.add_argument(\
    '--refseq', type=str, required=False,
    help='The fasta sequence of the reference genome. The fasta index of the reference genome is also required')
readlet.addArgs(parser)
args = parser.parse_args()

"""
Conducts radial clustering
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
def score(seq, site):
    wsize = len(site) # window size
    nwins = len(seq)-wsize
    wins = [0]*nwins
    for i in range(0,nwins):
        for j in range(0,len(site)):
            s = 1 if site[j]==seq[i] else -1
            wins[i]+=s
    return wins

"""
Returns the site by finding the maximum in the scores
To break ties it uses the direction.  
If direction=="5", that means its a 5' end and it will return the score closest to the 5' end (aka. left)
The vice versa happens with direction=="3"
"""
def findSite(scores,direction):
    count = -1 if direction=="5" else 1
    i = len(scores)-1 if direction=="5" else 0
    m, ind = -1, -1
    while i>=0 and i<len(scores):
        if m > scores[i]:
            ind = i
        i+=count        
    return ind

"""
Note: site is formatted as follows: XX-XX (e.g. GT-AG)
Returns the 5' and 3' splice sites within a multiple interval
"""
def sliding_window(refID, ivals, site,fastaF):
    print ival
    assert len(ival)==2
    start,end = ival[0],ival[1]
    toks = site.split("-")
    assert len(toks)==2
    site5p,site5p = toks[0],toks[1]
    seq = fastaF.fetch_sequence(refID,start,end)
    prime5 = seq[:2*args.readletIval]  #5' intron end
    prime3 = seq[-2*args.readletIval:] #3' intron end
    scores5,scores3 = score(prime5,site5p), score(prime3,site3p)
    return findSite(scores5,"5"),findSite(scores5,"3")


#Get counts of splice junctions per sample
def getJunctionSites(refID,bins,fastaF):
    sites5, sites3 = [],[]
    for coords,introns in bins.iteritems():
        print coords,introns
        site5,site3 = sliding_window(refID,introns,"GT-AG",fastaF)
        sites5.append(site5)
        sites3.append(site3)
    return sites5,sites3



starts = []  #Contains starting positions of introns
ends = []    #Contains ending positions of introns
labs = []    #Sample labels of introns
last_pt = "\t"
fnh = fasta.fasta(args.refseq)
last_ref = "\t"

for ln in sys.stdin:
    # Parse next read                                                                                                       
    ln = ln.rstrip()
    toks = ln.split('\t')
    assert len(toks)==5
    pt, st, en, refid, lab = toks[0], int(toks[1]), int(toks[2]), toks[3], toks[4]
    if last_pt=='\t':
        last_pt, last_ref = pt, refid
    elif last_pt!=pt:
        intron_ivals = zip(starts,ends,labs)
        #Cluster all introns with similar start and end positions   
        bins = cluster(intron_ivals)
        #Apply sliding windows to find splice junction locations
        sites5,sites3 = getJunctionSites(last_ref,bins,fnh)
        print "Site 5",sites5
        print "Site 3",sites3

        starts,ends,labs = [],[],[]

    starts.append(st)
    ends.append(en)
    labs.append(lab)
    last_pt,last_ref = pt,refid
        
intron_ivals = zip(starts,ends,labs)
bins = cluster(intron_ivals)
