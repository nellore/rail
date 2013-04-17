'''
normalize.py
(after walk_prenorm.py, before normalize_post.py)

Given position counts for a given sample, calculates some relevant per-sample
summary statistic, such as median, mean, total, or upper quartile.  Summary
statistics are passed along to normalize_post.py, which prints them all to a
file.

Tab-delimited input tuple columns:
 1. Sample label
 2. Count (at some position)

Binning/sorting prior to this step:
 1. Binned by sample label

Tab-delimited output tuple columns:
 1. Sample label
 2. Normalization factor

'''

import sys
import argparse
import time
timeSt = time.clock()

parser = argparse.ArgumentParser(description=\
    'Takes per-position counts for given samples and calculates a summary '
    'statistic such as median or 75% percentile.')

parser.add_argument(\
    '--percentile', metavar='FRACTION', type=float, required=False, default=0.75,
    help='For a given sample, the per-position percentile to extract as the normalization factor')

args = parser.parse_args()

ninp = 0                   # # lines input so far
nout = 0                   # # lines output so far
last_samp = "\t"           # id of previous partition
cov = dict()               # coverage histogram
totcov = 0                 # total coverage
totnonz = 0                # total positions with non-0 coverage

def percentile(cov):
    ''' Given a histogram (dictionary mapping keys to integers), return
        the upper quartile.  Same as 75th percentile.  Return a string.
        '''
    cur = 0
    for (k, v) in sorted(cov.iteritems(), reverse=True):
        cur += v
        assert cur <= totnonz
        if cur > ((1.0 - args.percentile) * totnonz):
            return str(k)
    raise RuntimeError("Should not reach this point")

for ln in sys.stdin:
    ln = ln.rstrip()
    toks = ln.split('\t')
    assert len(toks) == 2
    samp, cv = toks[0], int(toks[1])
    if samp != last_samp and last_samp != "\t":
        print "%s\t%s" % (last_samp, percentile(cov))
        nout += 1
        cov = dict()
        totcov, totnonz = 0, 0
    last_samp = samp
    ninp += 1
    cov[cv] = cov.get(cv, 0) + 1
    totcov += cv
    totnonz += 1

if last_samp != "\t":
    print "%s\t%s" % (last_samp, percentile(cov))
    nout += 1

# Done
timeEn = time.clock()
print >>sys.stderr, "DONE with normalize.py; in/out = %d/%d; time=%0.3f secs" % (ninp, nout, timeEn-timeSt)
