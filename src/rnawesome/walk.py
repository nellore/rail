'''
walk.py

Walk along a partition of the genome and compile all the spiced
alignments into a series of coverage vectors.  In a way, we're
realizing the positions x samples coverage matrix one row at a time.

How do we know the extents of the partition?

TODO:
- Need to have all the labels known
- Optionally do run length encoding
- Optionally omit all-0 lines
'''

import os
import sys
import site
import argparse

base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
site.addsitedir(os.path.join(base_path, "interval"))
site.addsitedir(os.path.join(base_path, "struct"))
site.addsitedir(os.path.join(base_path, "manifest"))

import circular
import partition
import manifest

parser = argparse.ArgumentParser(description=\
    'Take spliced alignments binned and sorted into partitions and '
    'construct coverage vectors for each sample.')

partition.addArgs(parser)
manifest.addArgs(parser)

args = parser.parse_args()

binsz = partition.binSize(args)

nlines = 0                 # # lines seen so far
last_pt = "\t"             # id of previous partition
last_st = -1               # start offset of previous fragment
part_st, part_en = -1, -1  # start/end offsets of current partition
ends = dict()              # maps sample names to circular buffers
cov = dict()               # current coverage in each sample
omitAllZero = True         # whether to omit rows that are all 0s

labs = manifest.labels(args)
ls = sorted(labs)

def handleInterval(last_st, st):
    # Wind all the buffers forward to just before this read's starting
    # position
    print >>sys.stderr, "  handling [%d, %d)" % (last_st, st)
    while last_st > -1 and st > last_st:
        if last_st >= part_st and last_st < part_en:
            elts = []
            # For all labels
            tot = 0
            for l in ls:
                if l in ends:
                    ends[l].advanceTo(last_st)
                    # Take into account reads ending at this position
                    nen = ends[l].get(last_st)
                    assert nen <= cov[l], "%d reads ended at %d but coverage was only %d" % (nen, last_st, cov[l])
                    cov[l] -= ends[l].get(last_st)
                    tot += cov[l]
                    elts.append(str(cov[l]))
                else:
                    elts.append("0")
            if not omitAllZero or tot > 0:
                print ("%d\t" % last_st) + "\t".join(elts)
        last_st += 1

def finishPartition(last_st, part_st, part_en):
    global ends, cov
    print >>sys.stderr, "Finished partition [%d, %d), last read start %d" % (part_st, part_en, last_st)
    handleInterval(last_st, part_en)
    ends = dict()
    cov = dict()

for ln in sys.stdin:
    ln = ln.rstrip()
    toks = ln.split('\t')
    assert len(toks) == 6
    pt, st, en, refid, weight, lab = toks
    st, en, weight = int(st), int(en), int(weight)
    if pt != last_pt:
        # We moved on to a new partition
        last_part_st, last_part_en = part_st, part_en
        if last_part_st >= 0:
            finishPartition(last_st, last_part_st, last_part_en)
        refid2, part_st, part_en = partition.parse(pt, binsz)
        print >>sys.stderr, "Started partition [%d, %d); first read: [%d, %d)" % (part_st, part_en, st, en)
        assert refid == refid2
        assert part_en > part_st
        last_st = -1 
    assert part_st >= 0
    assert part_en > part_st
    if lab not in ends:
        # We only use the count buffer to record where reads *end* so
        # we don't need to extend past either end of the partition
        ends[lab] = circular.CircularCountBuffer(part_st, binsz)
        assert len(ends[lab]) == binsz
        cov[lab] = 0
    if en >= part_st and en < part_en:
        ends[lab].add(weight, en) # increment a counter at the position where the read ends
    # Wind all the buffers forward to just before this read's starting
    # position
    if last_st > -1 and st > last_st:
        handleInterval(last_st, st)
    elif last_st == -1 and st > part_st:
        handleInterval(part_st, st)
    cov[lab] += weight
    last_pt, last_st = pt, st
    nlines += 1

if part_st > -1:
    finishPartition(last_st, part_st, part_en)

# Done
print >>sys.stderr, "DONE with walk.py; processed %d lines" % nlines
