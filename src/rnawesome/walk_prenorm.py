'''
walk_prenorm.py
(after merge.py, before normalize.py)

Walk along a partition of the genome.  Given all the merged spliced-alignment
intervals, collapse them down into a series of per-position coverage vectors.
Then print those vectors on a per-group basis.  When printing output, we don't
have to distinguish between positions or reference IDs because all that matters
to normalization is the group label and the count.

Todo: Right now, the entire partition name is being used for the chromosome name

Tab-delimited input tuple columns:
 1. Partition ID for partition overlapped by interval
 2. Interval start
 3. Interval end (exclusive)
 4. Reference ID
 5. Interval count
 6. Sample label

Binning/sorting prior to this step:
 1. Binned by partition
 2. Bins sorted by Interval start

Tab-delimited output tuple columns:
 1. Sample label
 2. Chromosome name
 3. Genome Position
 4. Count (at some position)

'''

import os
import sys
import site
import argparse
import time
timeSt = time.clock()

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

ninp = 0                   # # lines input so far
nout = 0                   # # lines output so far
last_pt = "\t"             # id of previous partition
last_st = -1               # start offset of previous fragment
part_st, part_en = -1, -1  # start/end offsets of current partition
ends = dict()              # maps sample names to circular buffers
cov = dict()               # current coverage in each sample
omitAllZero = True         # whether to omit rows that are all 0s

# Get the set of all labels by parsing the manifest file, given on the
# filesystem or in the Hadoop file cache
labs = manifest.labels(args)
ls = sorted(labs)    

def finish(partId, sampleId, covtup):
    st, cov = covtup
    nout = 0
    chrid = partId[:partId.rfind(';')]
    for i in xrange(0, len(cov)):
        if cov[i] > 0:
            print "%s\t%s\t%012d\t%d" % (sampleId, chrid, st+i, cov[i])
            nout += 1
    return nout

def finishAll(pt, cov):
    nout = 0
    for samp, covbuf in cov.iteritems():
        covst, covl = covbuf.finalize()
        if len(covl) > 0:
            nout += finish(pt, samp, (covst, covl))
    return nout

verbose = False
maxlen = 100

for ln in sys.stdin:
    ln = ln.rstrip()
    toks = ln.split('\t')
    assert len(toks) == 6
    pt, st, en, refid, weight, lab = toks  
    st, en, weight = int(st), int(en), int(weight)
    maxlen = max(en - st, maxlen)
    if pt != last_pt:
        # We moved on to a new partition
        last_part_st, last_part_en = part_st, part_en
        nout += finishAll(last_pt, cov)
        cov = {}
        refid2, part_st, part_en = partition.parse(pt, binsz)
        if verbose:
            print >>sys.stderr, "Started partition [%d, %d); first read: [%d, %d)" % (part_st, part_en, st, en)
        assert refid == refid2
        assert part_en > part_st
        last_st = -1 
    assert part_st >= 0
    assert part_en > part_st
    if lab not in cov:
        cov[lab] = circular.CircularCoverageBuffer(part_st, part_en, maxlen)
    covst, covl = cov[lab].add(st, en, weight)
    if len(covl) > 0:
        nout += finish(pt, lab, (covst, covl))
    last_pt, last_st = pt, st
    ninp += 1

if part_st > -1:
    nout += finishAll(last_pt, cov)

# Done
timeEn = time.clock()
print >>sys.stderr, "DONE with walk_prenorm.py; in/out = %d/%d; time=%0.3f secs" % (ninp, nout, timeEn-timeSt)
