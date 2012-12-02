#!/usr/bin/env python

"""
merge.py

If there are a lot of label, interval pairs that are identical, merge
them down into a single pair with a larger counter.
"""

import sys
import argparse

parser = argparse.ArgumentParser(description=\
    'Merge reads that cover the same interval within a given partition.')
args = parser.parse_args()

ninp, nout = 0, 0 # # lines input/output so far
last_pt = "\t"    # last partition id
last_refid = "\t" # last ref id
last_st = -1      # last start pos
last_en = -1      # last end pos
cnts = dict()     # per-label counts for a given rdid

def flushCnts(pt, refid, st, en):
    global nout
    for k, v in cnts.iteritems():
        print "%s\t%d\t%d\t%s\t%d\t%s" % (pt, st, en, refid, v, k)
        nout += 1

for ln in sys.stdin:
    ln = ln.rstrip()
    toks = ln.split('\t')
    assert len(toks) == 5
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
print >>sys.stderr, "DONE with merge.py; in/out = %d/%d" % (ninp, nout)
