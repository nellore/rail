"""
Tab-delimited input tuple columns:
1. Partition ID for partition overlapped by interval (also includes strand information)
2. Interval start
3. Interval end (exclusive)
4. Sample label
5. Readlet Sequence before 5' site
6. Readlet Sequence after 5' site
7. Readlet Sequence before 3' site
8. Readlet Sequence after 3' site

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
"""

import os
import sys
import argparse
import site
import time
import re
import string
from collections import defaultdict
timeSt = time.time()

base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for subdir in [ 'fasta', 'cluster', 'read', 'alignment', 'statsmath', 'util', 'interval' ]:
    site.addsitedir(os.path.join(base_path, subdir))

import partition
import chrsizes
import fasta
import readlet
import needlemanWunsch
import histogram
import counter
import window
from strip import Strip, cluster

parser = argparse.ArgumentParser(description=\
                                     'Reports splice junction information')
parser.add_argument(\
    '--test', action='store_const', const=True, default=False, help='Run unit tests')
parser.add_argument(\
    '--refseq', type=str, required=False,
    help='The fasta sequence of the reference genome. The fasta index of the reference genome is also required')
parser.add_argument(\
    '--cluster-radius', type=int, required=False, default=10,
    help='The radius of the clustering algorithm and the radius of the sliding window algorithm.')
parser.add_argument(\
    '--sliding-window-radius', type=int, required=False, default=10,
    help='The amount to search on either side of the mode in sliding-window algorithm.')
parser.add_argument(\
    '--intron-partition-overlap', type=int, required=False, default=20,
    help='Amount that partitions overlap their left and right neighbors by.')
parser.add_argument(\
    '--scores-file', type=str, required=False,default="",
    help='The location of the scores output file')

readlet.addArgs(parser)
partition.addArgs(parser)

args = parser.parse_args()

ninp, nout = 0, 0 # # lines input/output so far

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
        n = 2 * radius

    return n

"""Weighs different canonical and non-canonical sites and weighs them
Note that fw_site,rev_site and weight are zipped up in each site
"""
def findBestSite(refID,sts,ens,sites,introns,strand,radius,fastaF):
    bs5,bs3 = 0,0  #Best sites
    bscore = 0     #Best score
    bseq = ""      #Best seq
    for s in sites:
        seq = s[0] if strand=='+' else s[1]
        w = s[2] #weight
        N = findHistogramLen(refID, sts, ens, radius, fastaF)
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

"""
Finds canonical sites (e.g GT-AG sites)

bins: (seed st, seed en) -> (st, en, lab, rdid)
"""
def getJunctionSites(refID, strand, bins, fastaF, radius):
    global nout
    for introns in bins.itervalues():
        samples = counter.Counter()
        coOccurences = defaultdict( list )
        splice_site = "GT-AG" if strand=="+" else "CT-AC"  #only consider canonical sites
        sts, ens, lab, rdid = zip(*introns)
        N = findHistogramLen(refID,sts,ens,radius,fastaF)

        site5,s5,site3,s3 = sliding_window(refID,sts,ens,splice_site,N,fastaF)
        #sites5,sites3   = nw_correct(refID,site5,site3,introns,strand,fastaF)
        #site5,s5,site3,s3 = sliding_window(refID,sites5,sites3,splice_site,N,fastaF) #Retrain using nw
        # threshold = 1.0

        #if s5<threshold or s3<threshold:
        splice_sites = [("GC-AC","CT-GC",1.0),
                        ("AT-AC","GT-AT",1.0)]
        _, _, _,nc = findBestSite(refID, sts, ens, splice_sites, introns, strand, radius, fastaF)
        
        # I don't get why findBestSite can't return "N" or "C"
        site_chr = "C"
        if args.scores_file!="" and nc>0 and (s5+s3)>0:
            handle = open(args.scores_file,'a')
            handle.write("%lf\t%lf\t%s\n"%( (s5+s3),nc,site_chr) )
        
        for intr in introns:
            _, _, lab, rdid = intr
            coOccurences[rdid].append( (site5,site3) )
            samples[lab]+=1
        
        # Output for bed sites
        for sam,counts in samples.items():
            print "site\t%s\t%012d\t%d\t%s\t%d"%(refID,site5,site3,sam,counts)
            nout+=1
        
        # Output for co-occurences
        for rdid,sites in coOccurences.items():
            if len(sites)>1:
                for s in sites:
                    left_site,right_site = s
                    print "cooccurence\t%s\t%s\t%d\t%d"%(rdid,refID,left_site,right_site)

def go():

    global ninp
    
    starts, ends, labs, rdids = [], [], [], []
    seq5_flanks, seq3_flanks = [], []
    last_pt, last_refid = '\t', '\t'
    last_pt_st, last_pt_en = None, None
    fnh = fasta.fasta(args.refseq)
    binsz = partition.binSize(args)
    mat = {}
    
    def binIntervals(clustering):
        bins = {}
        for st, en, lab, rdid in zip(starts, ends, labs, rdids):
            assert (st, en) in clustering.cmap
            idx = clustering.cmap[(st, en)]
            if idx in bins:
                bins[idx].append((st, en, lab, rdid))
            else:
                bins[idx] = [(st, en, lab, rdid)]
        return bins
    
    def handlePartition():
        # Step 1. Make a strip
        strip = Strip.fromHistogram(mat)
        
        # Step 2. Cluster splice sites within strip
        clustering = cluster(strip, N=args.cluster_radius)
        clustering.limitTo(last_pt_st, last_pt_en)
        
        # Step 3. Build a bins object
        bins = binIntervals(clustering)
        
        # Step 4: Apply sliding windows to find splice junction locations
        getJunctionSites(last_refid, last_pt[-1], bins, fnh, args.sliding_window_radius)
    
    for ln in sys.stdin:
        # Parse next read
        toks = ln.rstrip().split('\t')
        assert len(toks) >= 7
        pt, st, en, lab, seq5_flank, seq3_flank, rdid = \
            toks[0], int(toks[1]), int(toks[2]), toks[3], toks[4], toks[5], toks[6]
        refid, pt_st, pt_en = partition.parse(pt[:-1], binsz)
        assert st >= pt_st and st < pt_en, "Intron start %d not in partition [%d, %d), partition id=%s" % (st, pt_st, pt_en, pt)
        
        if last_pt != pt and last_pt != '\t':
            handlePartition()
            starts, ends, labs, rdids = [], [], [], []
            seq5_flanks, seq3_flanks = [], []
            mat = {}
        
        # Compile histogram of how many times a start, end 
        mat[(st, en)] = mat.get((st, en), 0) + 1
        
        starts.append(st)
        ends.append(en)
        labs.append(lab)
        seq5_flanks.append(seq5_flank)
        seq3_flanks.append(seq3_flank)
        rdids.append(rdid)
        last_pt, last_refid = pt, refid
        last_pt_st, last_pt_en = pt_st, pt_en
        ninp += 1
    
    if last_pt != '\t': handlePartition()
    
    timeEn = time.time()
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

