#!/usr/bin/env python

"""
align.py

Alignment script, usually to be used in Hadoop apps.

Wraps Bowtie.  Has features for (1) optionally extracting readlets,
(2) optionally truncating reads or omitting mates, (3) hashing the read
name or partitioning the genome position.
"""

import os
import sys
import argparse
import site

base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
site.addsitedir(os.path.join(base_path, "interval"))

import interval

parser = argparse.ArgumentParser(description=\
    'Take readlets aligned using Bowtie and combine into spliced alignments.')

args = parser.parse_args()

def handleRead(rdid, ivals):
    ''' Given all the intervals where readlets from a given read
        aligned, compose a spliced alignment record. '''
    for k in ivals.iterkeys():
        for iv in iter(ivals[k]):
            # No reason to print rdid
            print "%s\t%d\t%d" % (k, iv.start, iv.end)

last_rdid = "\t"
nlines = 0

ivals = dict()

for ln in sys.stdin:
    ln = ln.rstrip()
    toks = ln.split('\t')
    rdid, fw, refid, refoff, seq, qual = toks[0], toks[1], toks[2], int(toks[3]), toks[4], toks[5]
    if rdid != last_rdid:
        handleRead(rdid, ivals)
        ivals = dict()
        last_rdid = rdid
    if refid not in ivals:
        ivals[refid] = interval.FlatIntervals()
    ivals[refid].add(interval.Interval(refoff, refoff+len(seq)))
    nlines += 1

handleRead(rdid, ivals)

# Done
print >>sys.stderr, "DONE with splice.py; processed %d lines" % nlines
