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

TODO:
 - Do input tuples really need the refid (4th) field?  Can't we recover from
   partition ID?

'''

import os
import sys
import site
import argparse
import time

base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
site.addsitedir(os.path.join(base_path, "interval"))
site.addsitedir(os.path.join(base_path, "struct"))

import circular
import partition

parser = argparse.ArgumentParser(description=\
    'Take spliced alignments binned and sorted into partitions and '
    'construct coverage vectors for each sample.')

partition.addArgs(parser)

parser.add_argument(\
    '--profile', action='store_const', const=True, help='Profile the code')
parser.add_argument(\
    '--partition-stats', action='store_const', const=True, help='Output statistics about bin sizes, time taken per bin, number of bins per reducer')

args = parser.parse_args()

def go():
    
    timeSt = time.clock()
    
    binsz = partition.binSize(args)
    
    ninp = 0                   # # lines input so far
    nout = 0                   # # lines output so far
    last_pt = "\t"             # id of previous partition
    part_st, part_en = -1, -1  # start/end offsets of current partition
    cov = dict()               # current coverage in each sample
    prefix = "o\t" if args.partition_stats else "";
    
    def finish(partId, sampleId, covtup):
        st, cov = covtup
        nout = 0
        chrid = partId[:partId.rfind(';')]
        for i in xrange(0, len(cov)):
            if cov[i] > 0:
                print "%s%s\t%s\t%012d\t%d" % (prefix, sampleId, chrid, st+i, cov[i])
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
    nInBin = 0
    nbin = 0
    timeBinSt = 0
    
    for ln in sys.stdin:
        ln = ln.rstrip()
        toks = ln.split('\t')
        assert len(toks) == 6
        pt, st, en, _, weight, lab = toks  
        st, en, weight = int(st), int(en), int(weight)
        _, part_st, part_en = partition.parse(pt, binsz)
        assert en > st
        assert en > part_st
        assert st < part_en
        maxlen = max(en - st, maxlen)
        if pt != last_pt:
            # We moved on to a new partition
            nout += finishAll(last_pt, cov)
            cov = {}
            if verbose:
                print >>sys.stderr, "Started partition [%d, %d); first read: [%d, %d)" % (part_st, part_en, st, en)
            assert part_en > part_st
            if last_pt != '\t' and args.partition_stats:
                timeBin = time.clock() - timeBinSt
                print '\t'.join(['partstats', str(nInBin), str(timeBin)])
                timeBinSt = time.clock()
                nInBin = 0
            nbin += 1
        nInBin += 1
        assert part_st >= 0
        assert part_en > part_st
        if lab not in cov:
            cov[lab] = circular.CircularCoverageBuffer(part_st, part_en, maxlen)
        covst, covl = cov[lab].add(st, en, weight)
        if len(covl) > 0:
            nout += finish(pt, lab, (covst, covl))
        last_pt = pt
        ninp += 1
    
    if last_pt != '\t' and args.partition_stats:
        timeBin = time.clock() - timeBinSt
        print '\t'.join(['partstats', str(nInBin), str(timeBin)])
    
    if part_st > -1:
        nout += finishAll(last_pt, cov)
    
    timeEn = time.clock()
    
    if args.partition_stats:
        print '\t'.join(['reducerstats', str(nbin), str(ninp), str(nout), str(timeEn - timeSt)])
    
    print >>sys.stderr, "DONE with walk_prenorm.py; in/out = %d/%d; time=%0.3f secs" % (ninp, nout, timeEn-timeSt)

if args.profile:
    import cProfile
    cProfile.run('go()')
else:
    go()
