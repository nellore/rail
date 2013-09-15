'''
normalize_pre.py
(after align.py --differential, before normalize.py)

Take the count differentials output by align.py and compile per-sample
coverage information.  We go bin-by-bin.

Tab-delimited input tuple columns:
 1. Partition ID for partition overlapped by interval
 2. Sample label
 3. Position
 4. Differential (-1 or +1, unless they've been aggregated)

Binning/sorting prior to this step:
 1. Binned by partition (field 1)
 2. Binned by sample label (field 2)
 3. Sorted by position (field 3)

Tab-delimited output tuple columns:
 1. Sample label
 2. Reference ID
 3. Position
 4. Coverage

'''

import os
import sys
import site
import argparse
import time

base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
site.addsitedir(os.path.join(base_path, "interval"))

import partition

parser = argparse.ArgumentParser(description=\
    'Take coverage differentials and compile coverage information for bucket.')

partition.addArgs(parser)

parser.add_argument(\
    '--partition-stats', action='store_const', const=True, help='Output statistics about bin sizes, time taken per bin, number of bins per reducer')

args = parser.parse_args()

def go():
    
    timeSt = time.time()
    binsz = partition.binSize(args)
    ninp, nout = 0, 0
    last_pt = '\t'
    prefix = 'o\t' if args.partition_stats else '';
    nbin, ninbin, timebin = 0, 0, timeSt
    last_lab, last_refid, last_off = None, None, None
    
    def output(sampleId, refid, off, cov):
        print "%s%s\t%s\t%012d\t%d" % (prefix, sampleId, refid, off, cov)
        return 1
    
    cnt = 0 # Current count
    for ln in sys.stdin:
        ninp += 1
        toks = ln.rstrip().split('\t')
        if toks[0] == 'DUMMY':
            continue
        assert len(toks) == 4, "Bad input line:\n" + ln
        pt, lab, off, diff = toks[0], toks[1], int(toks[2]), int(toks[3])
        refid, _, _ = partition.parse(pt, binsz)
        newChunk = pt != last_pt or lab != last_lab
        if newChunk: cnt = 0
        newOff = newChunk or off != last_off
        if newOff and last_off is not None:
            nout += output(last_lab, last_refid, last_off, cnt)
        cnt += diff
        
        if pt != last_pt:
            # We moved on to a new partition
            if last_pt != '\t' and args.partition_stats:
                timeDiff = time.time() - timebin
                print '\t'.join(['partstats', str(ninbin), str(timeDiff)])
                timebin = time.time()
                ninbin = 0
            nbin += 1
        ninbin += 1
        
        last_lab, last_refid, last_off, last_pt = lab, refid, off, pt
    
    if last_pt != '\t' and args.partition_stats:
        timeDiff = time.time() - timebin
        print '\t'.join(['partstats', str(ninbin), str(timeDiff)])
    
    if last_pt != '\t':
        assert last_off is not None
        nout += output(last_lab, last_refid, last_off, cnt)
    
    timeEn = time.time()
    
    if args.partition_stats:
        print '\t'.join(['reducerstats', str(nbin), str(ninp), str(nout), str(timeEn - timeSt)])
    
    print >>sys.stderr, "DONE with normalize_pre.py; in/out = %d/%d; time=%0.3f secs" % (ninp, nout, timeEn-timeSt)

go()
