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

base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
site.addsitedir(os.path.join(base_path, "annotation"))
site.addsitedir(os.path.join(base_path, "struct"))

import gtf
import search

parser = argparse.ArgumentParser(description=\
                                     'Splice junction validator')
parser.add_argument(\
    '--xscripts-file', metavar='path', type=str, required=True,
    help='Path of the transcripts pickle file')
parser.add_argument(\
    '--bed-file', metavar='path', type=str, required=True,
    help='Path of the transcripts pickle file')
parser.add_argument(\
    '--radius', type=int, required=False,default=10,
    help='The radius of tolerance for identifying splice site neighborhoods')
args = parser.parse_args()

"""
Finds all annotated splice sites based off of annotated transcripts
"""
def annotated_sites(xscripts):
    sites = []
    for x in xscripts:
        sites+=x.getSites()
    sites = list(set(sites))
    sites.sort()
    return sites

def readBedSites(bedfile):
    sites = []
    with open(bedfile,'r') as fh:
        for ln in fh:
            line = ln.rstrip()
            toks = line.split("\t")
            sites.append(int(toks[1]))
    return sites


def compare(bed_sites,annot_sites,radius):
    correct = 0
    nearby  = 0
    incorrect = 0
    for guess in bed_sites:
        exact = search.find(annot_sites,guess)
        if guess==exact:
            #print "Correct","Guess",guess,"Exact",exact
            correct+=1
        elif abs(guess-exact)<=radius:
            #print "Nearby","Guess",guess,"Exact",exact
            nearby+=1
        else:
            #print "Incorrect","Guess",guess,"Exact",exact
            incorrect+=1
    return correct/2,nearby/2,incorrect/2 #since we looking at 2x sites

if __name__=="__main__":
    xscripts = pickle.load(open(args.xscripts_file,'rb'))
    bed_sites = readBedSites(args.bed_file)
    annot_sites = annotated_sites(xscripts)
    #print annot_sites
    correct,nearby,incorrect = compare(bed_sites,annot_sites,args.radius)
    total=(correct+nearby+incorrect)
    print "Correct  \t",correct
    print "Nearby   \t",nearby
    print "Incorrect\t",incorrect
    print "Total    \t",total


