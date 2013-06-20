#!/usr/bin/env python

"""
splice.py
(after align.py, before merge.py)

Given alignments from align.py, merge all intervals covered by a readlet into a
minimal set of intervals covered by at least one readlet.  For output tuples,
we generate a partition id based on the partition size given by --ntasks and
--genomeLen.

Tab-delimited input tuple columns:
1. Read name
2. Orientation +/-
3. Reference ID
4. 0-based reference offset
5. Nucleotide sequence
6. Quality sequence

Binning/sorting prior to this step:
1. Binned by read name

Tab-delimited output tuple columns:
1. Partition ID for partition overlapped by interval
2. Interval start
3. Interval end (exclusive)
4. Reference ID
5. Sample label
"""

import os
import sys
import argparse
import site
import time
timeSt = time.clock()

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
ndigits = 12 #maxinum number of digits of start position number 

def handleRead(rdid, ivals):
    ''' Given all the intervals where readlets from a given read
        aligned, compose a spliced alignment record. '''
    global nout
    global ndigits
    for k in ivals.iterkeys():
        for iv in iter(ivals[k]):
            st, en = iv.start, iv.end
            # Keep stringing rdid along because it contains the label string
            # Add a partition id that combines the ref id and some function of
            # the offsets
            for pt in iter(partition.partition(k, st, en, binsz)):
                start = str(st)
                assert len(start)<=ndigits
                start = (ndigits-len(start))*"0"+start
                print "%s\t%s\t%d\t%s\t%s" % (pt, start, en, k, sample.parseLab(rdid))
                #print "%s\t%d\t%d\t%s\t%s" % (pt, st, en, k, sample.parseLab(rdid))
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
timeEn = time.clock()
print >>sys.stderr, "DONE with splice.py; in/out = %d/%d; time=%0.3f secs" % (ninp, nout, timeEn-timeSt)
