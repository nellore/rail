'''
walk.py

Walk along a partition of the genome and compile all the spiced
alignments into a series of coverage vectors.  In a way, we're
realizing the positions x samples coverage matrix one row at a time.
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

parser.add_argument(\
    '--columnize', action='store_const', const=True,
    help='Emit per-sample, per-position counts, instead of per-sample vectors')

partition.addArgs(parser)
manifest.addArgs(parser)

args = parser.parse_args()

binsz = partition.binSize(args)

ninp = 0                   # # lines input so far
nout = 0                   # # lines output so far
last_pt = "\t"             # id of previous partition
last_st = -1               # start offset of previous fragment
part_st, part_en = -1, -1  # start/end offsets of current partition
ends = dict()              # maps sample names to circular buffers
cov = dict()               # current coverage in each sample
omitAllZero = True         # whether to omit rows that are all 0s

# if true, each line is a per-sample per-position count. otherwise,
# each is per-position count vector for all samples.  Note that if
# this is true, a downstream aggregator will have to re-aggregate by
# genome window in rebuild the per-position count vector. 
columnize = args.columnize

# Get the set of all labels by parsing the manifest file, given on the
# filesystem or in the Hadoop file cache
labs = manifest.labels(args)
ls = sorted(labs)

def handleInterval(last_st, st):
    # Wind all the buffers forward to just before this read's starting
    # position
    global nout
    while last_st > -1 and st > last_st:
        if last_st >= part_st and last_st < part_en:
            # For all labels
            tot = 0
            if columnize:
                for l in ls:
                    if l in ends:
                        ends[l].advanceTo(last_st)
                        # Take into account reads ending at this position
                        nen = ends[l].get(last_st)
                        assert nen <= cov[l], "%d reads ended at %d but coverage was only %d" % (nen, last_st, cov[l])
                        cov[l] -= ends[l].get(last_st)
                        tot += cov[l]
                    if l in cov and cov[l] > 0:
                        # Print the position then all the coverage numbers
                        print "%s\t%d\t%d" % (l, last_st, cov[l])
                        nout += 1
            else:
                elts = []
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
                # 'elts' now contains a list of stringized coverage numbers
                if not omitAllZero or tot > 0:
                    # Print the position then all the coverage numbers
                    print ("%d\t" % last_st) + "\t".join(elts)
                    nout += 1
        last_st += 1

def finishPartition(last_st, part_st, part_en, verbose=False):
    global ends, cov
    if verbose:
        print >>sys.stderr, "Finished partition [%d, %d), last read start %d" % (part_st, part_en, last_st)
    handleInterval(last_st, part_en)
    ends = dict()
    cov = dict()

verbose = False

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
        if verbose:
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
    ninp += 1

if part_st > -1:
    finishPartition(last_st, part_st, part_en)

# Done
print >>sys.stderr, "DONE with walk.py; in/out = %d/%d" % (ninp, nout)
