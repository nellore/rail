"""
validate.py

Reads in the bed file containing the estimated splice sites and the pickle file containing the transcripts and provides the following statistics
1.  Number of exactly correct splice junctions
2.  Number of splice junctions within 1 radius (specified by user)
3.  Number of splice junctions completely off
4.  Plot of the distribution of error
"""
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

import gtf
import search
import fasta

parser = argparse.ArgumentParser(description=\
                                     'Splice junction validator')
parser.add_argument(\
    '--xscripts-file', metavar='path', type=str, required=True,
    help='Path of the transcripts pickle file')
parser.add_argument(\
    '--sites-file', metavar='path', type=str, required=True,
    help='Path of the splice sites pickle file')
parser.add_argument(\
    '--bed-file', metavar='path', type=str, required=True,
    help='Path of the transcripts pickle file')
parser.add_argument(\
    '--radius', type=int, required=False,default=10,
    help='The radius of tolerance for identifying splice site neighborhoods')
parser.add_argument(\
    '--refseq', type=str, required=True,
    help='The reference sequence')
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
        sites[k].sort()
    return sites

"""
Bins all splice sites in bed file and bins them by reference sequence id
"""
def readBedSites(bedfile):
    sites = defaultdict(list)
    with open(bedfile,'r') as fh:
        for ln in fh:
            line = ln.rstrip()
            toks = line.split("\t")
            seq,st = toks[0],int(toks[1])
            sites[seq].append(st)
    for k,v in sites.iteritems():
        sites[k] = list(set(sites[k]))
        sites[k].sort()
    return sites

def union_sites(sites):
    total_sites = set()
    for k,v in sites.iteritems():
        total_sites = total_sites.union(set(v))
    return set(total_sites)

def compare(bed_sites,annot_sites,radius):
    correct = 0
    nearby  = 0
    incorrect = 0
    total_sites = union_sites(annot_sites)
    missed_sites = union_sites(annot_sites)
    found_sites = set()
    close_sites = set()
    false_sites = set()
    total = len(missed_sites)
    for k,v in bed_sites.iteritems():
        for guess in v:
            if len(annot_sites[k])==0:
                continue
            exact = search.find(annot_sites[k],guess)
            if guess==exact:
                #print "Correct","Guess",guess,"Exact",exact
                correct+=1
                found_sites.add(exact)
                missed_sites.discard(exact)
            elif abs(guess-exact)<=radius:
                #print "Nearby","Guess",guess,"Exact",exact
                if exact not in close_sites:
                    close_sites.add(exact)
                    missed_sites.discard(exact)
                else:
                    false_sites.add(guess)
                    incorrect+=1
            else:
                #print "Incorrect","Guess",guess,"Exact",exact
                false_sites.add(guess)
                incorrect+=1
    incorrect_sites = found_sites.intersection(close_sites)
    nearby = len(close_sites.difference(found_sites))
    incorrect+=len(incorrect_sites)
    return found_sites,close_sites,false_sites,missed_sites,total_sites #since we looking at 2x sites

if __name__=="__main__":
    xscripts = pickle.load(open(args.xscripts_file,'rb'))
    sites = pickle.load(open(args.sites_file,'rb'))
    bed_sites = readBedSites(args.bed_file)
    annot_sites = annotated_sites(xscripts)

    #print annot_sites
    found_sites,close_sites,false_sites,missed_sites,total_sites = compare(bed_sites,annot_sites,args.radius)
    fastaH = fasta.fasta(args.refseq)

    intersect_sites = list(missed_sites.intersection(sites))
    missed_sites = list(missed_sites)
    total_sites = list(total_sites)
    sites = list(sites)
    found_sites = list(found_sites)
    false_sites = list(false_sites)
    close_sites = list(close_sites)

    total_sites.sort()
    sites.sort()
    found_sites.sort()
    missed_sites.sort()
    false_sites.sort()
    intersect_sites.sort()
    close_sites.sort()
    
    print >>sys.stderr, "Annot Sites ",total_sites
    print >>sys.stderr, "Sim Sites   ",sites
    print >>sys.stderr, "Found Sites ",found_sites
    print >>sys.stderr, "Close Sites ",close_sites
    print >>sys.stderr, "Missed Sites",missed_sites
    print >>sys.stderr, "False Sites ",false_sites
    print >>sys.stderr, "Intersect   ",intersect_sites

    missed = len(missed_sites)/2
    #bed_site_stats = siteDistribution(bed_sites,fastaH)
    #annot_site_stats = siteDistribution(annot_sites,fastaH)
    total = len(total_sites)/2
    sim_total = len(sites)/2
    correct = len(found_sites)/2
    nearby = len(close_sites)/2
    incorrect = len(false_sites)/2

    print "Total annot sites   \t",total
    print "Num sim sites       \t",sim_total
    print "Correct             \t",correct
    print "Nearby              \t",nearby
    print "False positives     \t",incorrect
    print "False negatives     \t",missed
    #print "Bed site stats      \t",bed_site_stats
    #print "Annotated site stats\t",annot_site_stats


