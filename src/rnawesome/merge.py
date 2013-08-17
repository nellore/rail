#!/usr/bin/env python

"""
merge.py
(after splice.py, before walk_prenorm.py)

If there are a lot of intervals that are identical given a particular label &
partition, merge them down into a single interval with counter >1.

Also, we re-partition them into smaller partitions (cleverly) in an effort to
make partitions with approximately equal interval membership.  This gives the
following walk_prenorm and walk_fit steps much better load balance.

Tab-delimited input tuple columns:
 1. Partition ID for partition overlapped by interval
 2. Interval start
 3. Interval end (exclusive)
 4. Sample label

Binning/sorting prior to this step:
 1. Binned by partition
 2. Bins sorted by Interval start

Tab-delimited output tuple columns (only column 5 is new):
 1. New Partition ID
 2. Interval start
 3. Interval end (exclusive)
 4. Interval count
 5. Sample label
"""

import os
import sys
import argparse
import time
import site

base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
site.addsitedir(os.path.join(base_path, "interval"))

import partition

parser = argparse.ArgumentParser(description=\
    'Merge reads that cover the same interval within a given partition.')

parser.add_argument(\
    '--partition-stats', action='store_const', const=True, help='Output statistics about bin sizes, time taken per bin, number of bins per reducer')
parser.add_argument(\
    '--profile', action='store_const', const=True, default=False, help='Profile the code')
parser.add_argument(\
    '--merge', action='store_const', const=True, default=False, help='Merge intervals with same start/end/sample')

parser.add_argument(\
    '--repartition', action='store_const', const=True, default=False, help='Repartition with finer-grain keys')
parser.add_argument(\
    '--target-cost', type=int, default=10000,
    help='For repartitioning: When deciding where to draw bin boundaries for output from this step, shoot for this per-bin cost')
parser.add_argument(\
    '--min-bin-size', type=int, default=750, help='For repartitioning: Minimum permissible bin size')
parser.add_argument(\
    '--nuc-cost', type=float, default=0.02,
    help='For repartitioning: Marginal cost of including one additional nucleotide in a bin')
parser.add_argument(\
    '--ival-cost', type=float, default=1.0,
    help='For repartitioning: Marginal cost of including one additional interval in a bin')

partition.addArgs(parser)
args = parser.parse_args()

prefix = "o\t" if args.partition_stats else ""

def reportRange(pt, st, en, amt, lab):
    # k is group label, v is # times the interval occurred
    assert len(str(st)) <= 12
    # we have the option of attaching a new partition label here, perhaps
    # to load balance the tuples better
    print "%s%s\t%012d\t%d\t%d\t%s" % (prefix, pt, st, en, amt, lab)
    return 1

class Writer(object):
    
    def __init__(self):
        pass
    
    def newBin(self, pt, pt_st, pt_en):
        self.pt = pt
    
    def next(self, st, en, amt, lab):
        return reportRange(self.pt, st, en, amt, lab)

class RepartitionWriter(object):
    """ This writer repartitions the genome in a way that attempts to equalize
        the 'cost' of the partitions.  'cost' is defined as ax + by, where
        a = # nucs, x = per-nuc cost, b = # intervals, y = per-interval cost.
    """
    
    def __init__(self, targetCost=10000, ivalCost=1, nucCost=0.02, minBin=750, costThresh=0.9):
        self.minBin = minBin
        self.ivalCost = ivalCost
        self.nucCost = nucCost
        self.costThresh = costThresh * targetCost
    
    def newBin(self, pt, pt_st, pt_en):
        self.pt = pt
        self.pt_st, self.pt_en = pt_st, pt_en
        self.lRepartSt, self.rRepartSt = pt_st, None
        self.lCostSoFar, self.rCostSoFar = 0, None
        self.max_en = -1
        self.bin = 0
        self.lRepartId, self.rRepartId = self.pt + '_0', None
        
        # home for records that are being buffered up until we know the end boundary for the bin
        self.buf = []
    
    def next(self, st, en, amt, lab):
        
        self.max_en = max(self.max_en, en)
        olapsLeft, olapsRight = True, False
        
        if self.rRepartSt is not None:
            if en >  self.rRepartSt: olapsRight = True
            if st >= self.rRepartSt: olapsLeft = False
        
        if olapsLeft: self.lCostSoFar += self.ivalCost
        if olapsRight: self.rCostSoFar += self.ivalCost
        
        if not olapsLeft:
            assert olapsRight
            self.lRepartSt,  self.rRepartSt  = self.rRepartSt,  None
            self.lCostSoFar, self.rCostSoFar = self.rCostSoFar, None
            self.lRepartId,  self.rRepartId  = self.rRepartId,  None
            self.bin += 1
            olapsLeft, olapsRight = True, False
        
        # Check if it's time to create a new bin
        sz = en - self.lRepartSt
        if self.rRepartSt is None and \
           en >= self.max_en and \
           en <= self.pt_en - 100 and \
           sz >= self.minBin and \
           self.lCostSoFar + sz * self.nucCost > self.costThresh:
            self.rRepartSt = self.max_en
            self.rCostSoFar = 0
            self.rRepartId = '_'.join([self.pt, str(self.bin+1)])
            assert not olapsRight
        
        # Now we can emit this interval
        n = 0
        if olapsLeft:
            if self.rRepartSt is not None:
                reportRange(self.lRepartId, st, en, amt, lab)
            else:
                self.buf.append((st, en, amt, lab))
            n += 1
        if olapsRight:
            reportRange(self.rRepartId, st, en, amt, lab)
            n += 1
        return n

def flushCnts(pt, st, en, cnts, wer):
    return sum([wer.next(st, en, v, k) for k, v in cnts.iteritems()])

def go():
    ninp, nout = 0, 0 # # lines input/output so far
    last_pt = "\t"    # last partition id
    last_st = -1      # last start pos
    last_en = -1      # last end pos
    cnts = dict()     # per-label counts for a given rdid
    nbin = 1
    timeSt = time.time()
    timeBinSt = timeSt
    nInBin = 0
    binsz = partition.binSize(args)
    
    wer = Writer()
    if args.repartition:
        wer = RepartitionWriter(targetCost=args.target_cost,
                                ivalCost=args.ival_cost,
                                nucCost=args.nuc_cost,
                                minBin=args.min_bin_size,
                                costThresh=0.9)
    
    for ln in sys.stdin:
        toks = ln.rstrip().split('\t')
        assert len(toks) == 4, "Bad input:\n" + ln
        pt, st, en, lab = toks[0], int(toks[1]), int(toks[2]), toks[3]
        assert pt != last_pt or st >= last_st
        
        if pt != last_pt:
            # Just transitioned to a new bin
            _, part_st, part_en = partition.parse(pt, binsz)
            wer.newBin(pt, part_st, part_en)
            if last_pt != '\t' and args.partition_stats:
                # Report some summary statistics about previous bin
                timeBin = time.time() - timeBinSt
                print '\t'.join(['partstats', str(nInBin), str(timeBin)])
                timeBinSt = time.time()
                nInBin = 0
            nbin += 1
        nInBin += 1
        
        if args.merge:
            # Possibly merge with previous range
            if pt == last_pt and st == last_st and en == last_en:
                # In a run of same-partition, same-end-point reads.  The
                # previous may or may not have been for the same label.
                cnts[lab] = cnts.get(lab, 0) + 1
            else:
                # Flush previous dict
                if last_st >= 0:
                    nout += flushCnts(last_pt, last_st, last_en, cnts, wer)
                cnts = {lab:1}
        else:
            # Output the range, without attempting to merge it with any other
            nout += wer.next(st, en, 1, lab)
        
        last_pt, last_st, last_en = pt, st, en
        ninp += 1
    
    if last_pt != '\t' and args.partition_stats:
        # Report some summary statistics about final bin
        timeBin = time.time() - timeBinSt
        print '\t'.join(['partstats', str(nInBin), str(timeBin)])
    
    if last_st >= 0:
        nout += flushCnts(last_pt, last_st, last_en, cnts, wer)
    
    timeEn = time.time()
    
    if args.partition_stats:
        # Report some summary statistics about this reduce task
        print '\t'.join(['reducerstats', str(nbin), str(ninp), str(nout), str(timeEn - timeSt)])
    
    print >>sys.stderr, "DONE with merge.py; in/out = %d/%d; time=%0.3f secs" % (ninp, nout, timeEn-timeSt)

if args.profile:
    import cProfile
    cProfile.run('go()')
else: go()
