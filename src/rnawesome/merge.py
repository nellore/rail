#!/usr/bin/env python

"""
merge.py
(after splice.py, before walk_prenorm.py)

If there are a lot of intervals that are identical given a particular label &
partition, merge them down into a single interval with counter >1.

Tab-delimited input tuple columns:
 1. Partition ID for partition overlapped by interval
 2. Interval start
 3. Interval end (exclusive)
 4. Reference ID
 5. Sample label

Binning/sorting prior to this step:
 1. Binned by partition
 2. Bins sorted by Interval start

Tab-delimited output tuple columns (only column 5 is new):
 1. Partition ID for partition overlapped by interval
 2. Interval start
 3. Interval end (exclusive)
 4. Reference ID
 5. Interval count
 6. Sample label
"""

import sys
import argparse
import time
timeSt = time.clock()

parser = argparse.ArgumentParser(description=\
    'Merge reads that cover the same interval within a given partition.')
args = parser.parse_args()

ninp, nout = 0, 0 # # lines input/output so far
last_pt = "\t"    # last partition id
last_refid = "\t" # last ref id
last_st = -1      # last start pos
last_en = -1      # last end pos
cnts = dict()     # per-label counts for a given rdid
ndigits = 12      #maxinum number of digits of start position number 

def flushCnts(pt, refid, st, en):
    global nout
    for k, v in cnts.iteritems():
        # k is group label, v is # times the interval occurred
        start = str(st)
        assert len(start)<=ndigits
        print "%s\t%012d\t%d\t%s\t%d\t%s" % (pt, st, en, refid, v, k)
        nout += 1

for ln in sys.stdin:
    ln = ln.rstrip()
    toks = ln.split('\t')
    assert len(toks) == 5, "Bad input:\n" + ln
    pt, st, en, refid, lab = toks
    st, en = int(st), int(en)
    assert pt != last_pt or st >= last_st
    if pt == last_pt and st == last_st and en == last_en:
        # In a run of same-partition, same-end-point reads.  The
        # previous may or may not have been for the same label.
        cnts[lab] = cnts.get(lab, 0) + 1
    else:
        # Flush previous dict
        if last_st >= 0:
            flushCnts(last_pt, last_refid, last_st, last_en)
        cnts = {lab:1}
    last_pt, last_refid, last_st, last_en = pt, refid, st, en
    ninp += 1

if last_st >= 0:
    flushCnts(last_pt, last_refid, last_st, last_en)

# Done
timeEn = time.clock()
print >>sys.stderr, "DONE with merge.py; in/out = %d/%d; time=%0.3f secs" % (ninp, nout, timeEn-timeSt)
