"""
validate.py

Reads in the bed file containing the estimated splice sites and the pickle file containing the transcripts and provides the following statistics
1.  Number of exactly correct splice junctions
2.  Number of splice junctions within 1 radius (specified by user)
3.  Number of splice junctions completely off
4.  Plot of the distribution of error
"""
import re
import os
import site
import argparse
import sys
import math
import pickle
import bisect
from collections import Counter
from collections import defaultdict
base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
site.addsitedir(os.path.join(base_path, "annotation"))
site.addsitedir(os.path.join(base_path, "struct"))
site.addsitedir(os.path.join(base_path, "fasta"))
site.addsitedir(os.path.join(base_path, "statsmath"))

import gtf
import search
import fasta
import window
import display

parser = argparse.ArgumentParser(description=\
                                     'Splice junction validator')
parser.add_argument(\
    '--xscripts-file', metavar='path', type=str, required=True,
    help='Path of the transcripts pickle file')
parser.add_argument(\
    '--sites-file', metavar='path', type=str, required=True,
    help='Path of the annotated splice sites pickle file')
parser.add_argument(\
    '--bed-file', metavar='path', type=str, required=True,
    help='Path of the estimated splice sites bed file')
parser.add_argument(\
    '--radius', type=int, required=False,default=10,
    help='The radius of tolerance for identifying splice site neighborhoods')
parser.add_argument(\
    '--refseq', type=str, required=True,
    help='The reference sequence')
parser.add_argument(\
                    '--false-positives',action='store_const', const=True, default=False, help='Indicates if false positives are printed')
parser.add_argument(\
                    '--false-negatives',action='store_const', const=True, default=False, help='Indicates if false positives are printed')
parser.add_argument(\
                    '--flank-seqs', type=str,required=True,help='The flanking sequences surrounding the intron')
parser.add_argument(\
                    '--region',type=str,required=False,default="",help='The coordinates of the sites to be displayed (e.g. chrX:1-100)')

args = parser.parse_args()


"""
Given a dictionary of sites, obtains reference sequence and returns a Counter object
"""
def siteDistribution(sites,fh):
    seq_hist = Counter()
    for k,s in sites.iteritems():
        for i in range(0,len(sites[k]),2):
            print >> sys.stderr,len(sites[k]),i,s[i],s[i+1]
            st,en = s[i],s[i+1]
            refseq = fh.fetch_sequence(k,st+1,en+1).upper()
            seq_hist[refseq]+=1
    return seq_hist
"""
Finds all annotated splice sites based off of annotated transcripts
and returns a list of splice sites binned by reference sequence id
"""
def annotated_sites(xscripts):
    sites = defaultdict(list)
    for x in xscripts:
        sites[x.seqid]=sites[x.seqid]+x.getSites()
    for k,v in sites.iteritems():
        sites[k] = list(set(sites[k]))
        sites[k].sort(key=lambda tup:tup[1])
        sites[k].sort(key=lambda tup:tup[0])
    return sites

"""
Bins all splice sites in bed file and bins them by reference sequence id
"""
def readBedSites(bedfile):
    sites = defaultdict(list)
    i = 0
    with open(bedfile,'r') as fh:
        for ln in fh:
            if i%2==1:
                sites[seqid].append( (st,en,seqid,"") )
            else:
                line = ln.rstrip()
                toks = line.split("\t")
                seqid,st,en = toks[0],int(toks[1]),int(toks[2])
            i+=1
    for k,v in sites.iteritems():
        sites[k] = list( set(sites[k])) #Remove redundancies
        sites[k].sort(key=lambda tup:tup[1])
        sites[k].sort(key=lambda tup:tup[0])
    return sites

#Removes all duplicates
def unique_sites(sites):
    total_sites = set()
    for k,v in sites.iteritems():
        total_sites = total_sites.union( set(v) )
    return set(total_sites)
"""
Compares the simulated sites and the output from the pipeline
"""
def compare(bed_sites,annot_sites,radius):
    correct = 0
    nearby  = 0
    incorrect = 0
    total_sites = unique_sites(annot_sites)
    missed_sites = unique_sites(annot_sites)
    found_sites = set()
    close_sites = set()
    false_sites = set()
    total = len(missed_sites)
    for k,v in bed_sites.iteritems():
        for guess in v:
            if len(annot_sites[k])==0:
                continue
            exact = search.find_tuple(annot_sites[k],guess)
            if (guess[0],guess[1],guess[2]) == (exact[0],exact[1],exact[2]):
                correct+=1
                found_sites.add(exact)
                missed_sites.discard(exact)
            elif abs(guess[0]-exact[0])<=radius:
                if exact not in close_sites:
                    close_sites.add(exact)
                    missed_sites.discard(exact)
                else:
                    false_sites.add(guess)
                    incorrect+=1
            else:
                false_sites.add(guess)
                incorrect+=1
    incorrect_sites = found_sites.intersection(close_sites)
    nearby = len(close_sites.difference(found_sites))
    incorrect+=len(incorrect_sites)
    return found_sites,close_sites,false_sites,missed_sites,total_sites



#def readOverlappedSite(fname):
#return set( map( int,open(fname,'r').readline().rstrip().split("\t")))

# """
# hist_obj is keyed by site and contains a base distribution of A,G,C,T,N,...
# """
# def flankingHist(site,flank_seq,hist_obj):
#     key = (site[0],site[1],site[2]) #Don't want xscript_id for key
#     for i in range(0,len(flank_seq)):

"""
Takes all flanking sequences, finds the closest annotated site and bins them
"""
def binFlanks(sites,flanks_file):
    flanks_dict = deflault_dict(list)
    with open(flanks_file,'r') as fnh:
        for ln in fnh:
            ln = ln.rstrip()
            toks = ln.split("\t")
            assert len(toks)==8
            st,end,seqid,flank_left,flank_right = int(toks[2]), int(toks[3]), toks[4], toks[6], toks[7]
            #Search for nearby sites
            guess5,guess3 = (st,st+1,seqid,""), (end-1,end,seqid,"")
            site5 = search.find_tuple(sites,guess5)
            site3 = search.find_tuple(sites,guess3)
            #Create key and bin flanking sequences
            key5 = (site5[0],site5[1],site5[2])
            key3 = (site3[0],site3[1],site3[2])
            flanks[key5].append(flank_left)
            flanks[key3].append(flank_right)
    return flanks


"""
Prints out false positives and false negatives specified by a region
as well as flanking sequences
"""
def display(fps,fns,flanks_dict,xscripts,annot_sites,region):
    #First convert xscripts into dictionary
    xscriptDict = {x.seqid: x for x in xscripts}  #Only available for >=python2.7

    # xscriptDict = dict()
    # for x in xscripts:
    #     xscriptDict[x.seqid] = x


    return



if __name__=="__main__":
    #sites = readOverlappedSites(args.site_file)
    xscripts = pickle.load(open(args.xscripts_file,'rb'))
    sites = pickle.load(open(args.sites_file,'rb'))
    bed_sites = readBedSites(args.bed_file)
    annot_sites = annotated_sites(xscripts)

    #print annot_sites
    found_sites,close_sites,false_sites,missed_sites,total_sites = compare(bed_sites,annot_sites,args.radius)
    fastaH = fasta.fasta(args.refseq)

    intersect_sites = list(missed_sites.intersection(sites))
    missed_sites    = list(missed_sites)
    total_sites     = list(total_sites)
    sites           = list(sites)
    found_sites     = list(found_sites)
    false_sites     = list(false_sites)
    close_sites     = list(close_sites)

    #Sort all of the lists with respect to start position and chromosome
    total_sites.sort(key=lambda tup:tup[1])
    total_sites.sort(key=lambda tup:tup[0])
    sites.sort(key=lambda tup:tup[1])
    sites.sort(key=lambda tup:tup[0])
    found_sites.sort(key=lambda tup:tup[1])
    found_sites.sort(key=lambda tup:tup[0])
    missed_sites.sort(key=lambda tup:tup[1])
    missed_sites.sort(key=lambda tup:tup[0])
    false_sites.sort(key=lambda tup:tup[1])
    false_sites.sort(key=lambda tup:tup[0])
    intersect_sites.sort(key=lambda tup:tup[1])
    intersect_sites.sort(key=lambda tup:tup[0])
    close_sites.sort(key=lambda tup:tup[1])
    close_sites.sort(key=lambda tup:tup[0])

    print >>sys.stderr, "Close Sites ",close_sites
    print >>sys.stderr, "Missed Sites",missed_sites
    print >>sys.stderr, "False Sites ",false_sites
    print >>sys.stderr, "Intersect   ",intersect_sites

    missed = len(missed_sites)
    total = len(total_sites)
    sim_total = len(sites)
    correct = len(found_sites)
    nearby = len(close_sites)
    incorrect = len(false_sites)

    print "Total annot sites   \t",total
    print "Num sim sites       \t",sim_total
    print "Correct             \t",correct
    print "Nearby              \t",nearby
    print "False positives     \t",incorrect
    print "False negatives     \t",missed
    #print "Bed site stats      \t",bed_site_stats
    #print "Annotated site stats\t",annot_site_stats


