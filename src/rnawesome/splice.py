#!/usr/bin/env python

"""
splice.py

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
site.addsitedir(os.path.join(base_path, "manifest"))

import interval
import partition
import sample
import manifest

parser = argparse.ArgumentParser(description=\
    'Take readlets aligned using Bowtie and combine into spliced alignments.')

partition.addArgs(parser)
manifest.addArgs(parser)

args = parser.parse_args()

binsz = partition.binSize(args)

nout = 0

def handleRead(rdid, ivals):
    ''' Given all the intervals where readlets from a given read
        aligned, compose a spliced alignment record. '''
    global nout
    for k in ivals.iterkeys():
        for iv in iter(ivals[k]):
            st, en = iv.start, iv.end
            # Keep stringing rdid along because it contains the label string
            # Add a partition id that combines the ref id and some function of
            # the offsets
            for pt in iter(partition.partition(k, st, en, binsz)):
                print "%s\t%d\t%d\t%s\t%s" % (pt, st, en, k, sample.parseLab(rdid))
                nout += 1

last_rdid = "\t"
ninp = 0
rdid = None
ivals = dict()

for ln in sys.stdin:
    # Parse next read
    ln = ln.rstrip()
    toks = ln.split('\t')
    rdid, fw, refid, refoff, seq, qual = toks[0], toks[1], toks[2], int(toks[3]), toks[4], toks[5]
    # If we just moved on to readlets for the next read...
    if rdid != last_rdid:
        if last_rdid != "\t":
            # Process the intervals accumulated for the previous read
            handleRead(last_rdid, ivals)
        # Start from scratch for next read
        ivals = dict()
        last_rdid = rdid
    # Add this interval to the flattened interval collection for current read
    if refid not in ivals:
        ivals[refid] = interval.FlatIntervals()
    ivals[refid].add(interval.Interval(refoff, refoff+len(seq)))
    ninp += 1

# Process the intervals accumulated for the final read
if last_rdid != "\t":
    handleRead(last_rdid, ivals)

# Done
print >>sys.stderr, "DONE with splice.py; in/out = %d/%d" % (ninp, nout)
