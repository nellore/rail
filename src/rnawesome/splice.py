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
site.addsitedir(os.path.join(base_path, "sample"))

import interval
import partition
import sample

parser = argparse.ArgumentParser(description=\
    'Take readlets aligned using Bowtie and combine into spliced alignments.')

partition.addArgs(parser)

args = parser.parse_args()

binsz = partition.binSize(args)

def handleRead(rdid, ivals):
    ''' Given all the intervals where readlets from a given read
        aligned, compose a spliced alignment record. '''
    for k in ivals.iterkeys():
        for iv in iter(ivals[k]):
            st, en = iv.start, iv.end
            # Keep stringing rdid along because it contains the label string
            # Add a partition id that combines the ref id and some function of
            # the offsets
            for pt in iter(partition.partition(k, st, en, binsz)):
                print "%s\t%s\t%d\t%d\t%s" % (pt, k, st, en, sample.parseLab(rdid))

last_rdid = "\t"
nlines = 0
rdid = None
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

if rdid is not None:
    handleRead(rdid, ivals)

# Done
print >>sys.stderr, "DONE with splice.py; processed %d lines" % nlines
