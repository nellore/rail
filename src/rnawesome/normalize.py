'''
normalize.py

Takes input from walk.py --columnize (each record is a sample name,
position, coverage tuple), and calculates some relevant per-sample
summary, such as median, mean, total, or upper quartile.  Summaries are
then written to an augmented manifest file on HDFS or other shared
filesystem and that file is passed to downstream stages via file
cacheing. 
'''

import sys
import argparse

parser = argparse.ArgumentParser(description=\
    'Takes per-sample bins of coverage rows (in any order) and '
    'calculates some relevant summary over all the coverage counts.')

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
    assert len(toks) == 3
    samp, pos, cv = toks[0], int(toks[1]), int(toks[2])
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
print >>sys.stderr, "DONE with normalize.py; in/out = %d/%d" % (ninp, nout)
