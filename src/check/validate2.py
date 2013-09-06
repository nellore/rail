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
import copy
from collections import Counter
from collections import defaultdict
base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
site.addsitedir(os.path.join(base_path, "annotation"))
site.addsitedir(os.path.join(base_path, "struct"))
site.addsitedir(os.path.join(base_path, "fasta"))
site.addsitedir(os.path.join(base_path, "statsmath"))
site.addsitedir(os.path.join(base_path, "util"))
site.addsitedir(os.path.join(base_path, "simulation"))

import gtf
import search
import fasta
import window
import display
import counter
import library

parser = argparse.ArgumentParser(description=\
                                     'Splice junction validator')
parser.add_argument(\
    '--bed-file', metavar='path', type=str, required=False, default=""
    help='Path of the estimated splice sites bed file')
parser.add_argument(\
    '--radius', type=int, required=False,default=10,
    help='The radius of tolerance for identifying splice site neighborhoods')
parser.add_argument(\
    '--window-radius', type=int, required=False,default=50,
    help='The radius of display window')
parser.add_argument(\
    '--refseq', type=str, required=False, default=""
    help='The reference sequence')
parser.add_argument(\
    '--gtf', type=str, required=False, default=""
    help='The annotated transcripts file')
parser.add_argument(\
    '--lib-file', type=str, required=False, default=""
    help='The library file containing all of the correct positions of the fragments')
parser.add_argument(\
    '--out-dir', type=str, required=False, default=""
    help='The output directory of the false positive and false negative regions')
parser.add_argument(\
    '--profile', action='store_const', const=True, default=False,
    help='Profile simulation generation')
parser.add_argument(\
    '--test', action='store_const', const=True, default=False,
    help='Run unittests')

display.addArgs(parser)

args = parser.parse_args()

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

"""
Compares the simulated sites and the output from the pipeline
"""
#TODO: Modify for flux simulator
def compare(bed_sites,annot_sites,radius):
    correct = 0
    nearby  = 0
    incorrect = 0
    total_sites = unique_sites(annot_sites)
    missed_sites = unique_sites(annot_sites) #false negatives
    found_sites = set() #correct sites
    false_sites = set() #false positives
    total = len(missed_sites)
    for k,v in bed_sites.iteritems():
        for guess in v:
            if len(annot_sites[k])==0:
                continue
            exact = search.find_tuple(annot_sites[k],guess)
            if (guess[0],guess[1],guess[2]) == (exact[0],exact[1],exact[2]):
                found_sites.add(exact)
                missed_sites.discard(guess)
            else:
                #NOTE: This contains both false positives and false negatives
                false_sites.add(guess)

    return found_sites,close_sites,false_sites,missed_sites,total_sites

def go():

    #Step 1: Isolate all read in .lib file that span splice junctions.  These are the annotated splice junctions
    annots = gtf.parseGTF(args.gtf)
    lib_frags = library.library(args.lib_file,annots)
    bed_sites = readBedSites(args.bed_file)
    #Step 2: Compare detected splice junctions to annotated splice junctions

    #Step 3: Output two files: false positive regions and false negative regions


if __name__=="__main__":
    if args.profile:
        import cProfile
        cProfile.run('go()')
    else:
        go_flux()
