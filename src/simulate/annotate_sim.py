"""
annotate_sim.py

Simulates differiential gene expression using annotated genes
"""

import os
import site
import argparse
import sys
import random
import math
from operator import itemgetter
from collections import defaultdict
base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
site.addsitedir(os.path.join(base_path, "fasta"))

import fasta
import re

"""
Input: sorted contigs and a dictionary of fasta seqs
Chain all exons into transcript via transcript id
"""
def assembleTranscripts(exons,fastaseqs):
    xscripts = defaultdict(list)
    for E in exons:
        tr_id = E[9]
        xscripts[tr_id].append(E)
        
    for tr_id,exons in xscripts.iteritems():
        exons.sort(key=itemgetter(3))
        for E in exons:
            seqid = E[0]
            start = E[3]
            end = E[4]
            seq = fastaseqs[seqid][start:end]
    return transcripts

def parseFASTA(fns):
    seqs = dict()
    lns = []
    seqid = ""
    numseqs = 0
    for fn in fns:
        if ln[0] == '>':
            seqid = ln[0].split(' ')[0][1:]
            if numseqs>0:
                seqs[seqid] = ''.join(lns)
            lns = []
            numseqs+=1
            continue
        else:
            lns.append(ln.rstrip())
    return seqs

reg = re.compile(r"(\S+)\s(\S+)\s(\S+)\s(\d+)\s(\d+)\s(\S+)\s(\S+)\s(\S+)\sgene_id\s\S(\S+)\S;\sgene_name\s\S\S+\S;\s*\S*\s*\S*;*\stranscript_id\s\S(\S+)\S;")

def parseGTF(fn):
    contigs = []
    with open(fn,'r') as fh:
        for ln in fh:
            l = ln.rstrip()
            seqname,source,feature,start,end,score,strand,frame,gene_id,transcript_id = reg.findall(ln)[0]
            if feature=="exon":
                contigs.append((seqname,source,feature,start,end,score,strand,frame,gene_id,transcript_id))
    return contigs

    
if __name__ =="__main__":
    parser = argparse.ArgumentParser(description=\
                                         'Transcript simulator')
    parser.add_argument(\
        '--fasta', metavar='path', type=str, required=True,
        help='FASTA file(s) containing reference genome sequences')
    parser.add_argument(\
        '--annotations', metavar='path', type=str, required=True,
        help='GTF file containing all fo the gene annotations')
    parser.add_argument(\
        '--output-prefix', metavar='path', type=str, required=False,
        help='Prefix for output read files')
    parser.add_argument(\
        '--read-len', metavar='int', action='store', type=int, default=100,
        help='Read length to simulate')
    parser.add_argument(\
        '--num-replicates', metavar='int', action='store', type=int, default=8,
        help='Number of replicates per group')
    parser.add_argument(\
        '--seed', metavar='int', action='store', type=int, default=874,
        help='Pseudo-random seed')
    parser.add_argument(\
        '--num-nucs', metavar='int', action='store', type=float, default=1e7,
        help='Number of total nucleotides of reads to generate')
    args = parser.parse_args()

    
    contigs = parseGTF(args.annotations)
    # contigs.sort(key=itemgetter(3))
    # contigs.sort(key=itemgetter(9))
    # contigs.sort(key=itemgetter(0))
    print contigs
    
