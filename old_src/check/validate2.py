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
import itertools
from collections import Counter
from collections import defaultdict
base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
site.addsitedir(os.path.join(base_path, "annotation"))
site.addsitedir(os.path.join(base_path, "struct"))
site.addsitedir(os.path.join(base_path, "fasta"))
site.addsitedir(os.path.join(base_path, "statsmath"))
site.addsitedir(os.path.join(base_path, "util"))
site.addsitedir(os.path.join(base_path, "simulate"))

import gtf
import search
import fasta
import window
import display
import counter

parser = argparse.ArgumentParser(description=\
                                     'Splice junction validator')
parser.add_argument(\
    '--bed-file', metavar='path', type=str, required=False, default="",
    help='Path of the estimated splice sites bed file')
parser.add_argument(\
    '--radius', type=int, required=False,default=10,
    help='The radius of tolerance for identifying splice site neighborhoods')
parser.add_argument(\
    '--window-radius', type=int, required=False,default=50,
    help='The radius of display window')
parser.add_argument(\
    '--refseq', type=str, required=True, nargs='+', default="",
    help='FASTA file(s) containing reference genome sequences')
parser.add_argument(\
    '--gtf', metavar='path', type=str, nargs='+', required=False,
    help='GTF file(s) containing gene annotations')
parser.add_argument(\
    '--lib-file', type=str, required=False, default="",
    help='The library file containing all of the correct positions of the fragments')
parser.add_argument(\
    '--actual-sites', type=str, required=False, default="",
    help='The bed file containing all of the correct positions of the reads')
parser.add_argument(\
    '--out-dir', type=str, required=False, default="",
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

def writeBedSites(bedfile,sites):
    handle = open(bedfile,'w')
    for site in sites:
        handle.write("%s\t%d\t%d\t1\n"%(site[2],site[0],site[1]+1))
    handle.close()

"""
Compares the simulated sites and the output from the pipeline
"""
def compare(bed_sites,annot_sites):
    correct = 0
    nearby  = 0
    incorrect = 0
    a = [v for v in annot_sites.itervalues()]
    missed_sites = set(itertools.chain.from_iterable(a))    #false negatives
    found_sites  = set()                                    #correct sites
    false_sites  = set()                                    #false positives
    for refid,guesses in bed_sites.iteritems():
        for guess in guesses:
            if len(annot_sites[refid])==0:
                continue
            exact = search.find_tuple(annot_sites[refid],guess)
            #print "Exact",exact,"Guess",guess
            if (guess[0],guess[1],guess[2]) == (exact[0],exact[1],exact[2]):
                found_sites.add(exact)
                missed_sites.discard(exact)
            else:
                false_sites.add(guess)
                missed_sites.discard(exact)

    return found_sites,false_sites,missed_sites

def go():

    #Step 1: Isolate all read in .red file that span splice junctions.  These are the annotated splice junctions
    annots = gtf.parseGTF(args.gtf)
    fastadb  = gtf.parseFASTA(args.refseq)
    xscripts = gtf.assembleTranscripts(annots,fastadb)
    annot_sites = readBedSites(args.actual_sites)
    bed_sites = readBedSites(args.bed_file)
    #Step 2: Compare detected splice junctions to annotated splice junctions
    found_sites,false_sites,missed_sites = compare(bed_sites,annot_sites)
    #Step 3: Sort sites
    found_sites  = list(found_sites)
    false_sites  = list(false_sites)
    missed_sites = list(missed_sites)
    annot_sites  = list( set( itertools.chain.from_iterable([v for v in annot_sites.itervalues()])))
    annot_sites.sort(key=lambda tup:tup[0])
    found_sites.sort(key=lambda tup:tup[0])
    false_sites.sort(key=lambda tup:tup[0])
    missed_sites.sort(key=lambda tup:tup[0])
    #Step 4: Output two files: false positive regions and false negative regions
    writeBedSites("%s/false_positives.bed"%(args.out_dir),false_sites)
    writeBedSites("%s/false_negatives.bed"%(args.out_dir),missed_sites)
    writeBedSites("%s/annotated.bed"%(args.out_dir),annot_sites)
    print "Correct sites  ", len(found_sites)
    print "False positives", len(false_sites)
    print "False negatives", len(missed_sites)

if __name__=="__main__":
    if args.profile:
        import cProfile
        cProfile.run('go()')
    else:
        go()
