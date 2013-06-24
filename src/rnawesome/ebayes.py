'''
ebayes.py
(after walk_fit.py, before hmm_params.py)

Given all the unmoderated t-statistics and others results from the
linear fits, calculate the moderated t-statistics.  In general, we need
to do this both for the test statistics, and for zero or more sets of
null statistics, each corresponding to a permutation of the group
labels.  When we output the moderated t-statistics, we include a
partition id since, later in the pipeline, we will want to partition
them into genome windows.

Tab-delimited input tuple columns:
 1. Reference ID
 2. Reference offset (0-based)
 3. Mean A
 4. Degrees of freedom
 5+. Comma-delimited triples of (1) coefficient, (2) standard deviation, (3)
     sigma.  One for the data, and then N more triples for each of the N
     permutations of the data tried.

Binning/sorting prior to this step:
 (none)

Tab-delimited output tuple columns:
 1. Partition ID
 2. Reference offset (0-based)
 3+. Comma-delimited pairs of (1) moderated t-staistic, (2) log fold-change.
     One pair for the data, then N more pairs for each of N permutations.

'''

import argparse
import sys
import os
import site
import math
import numpy as np
import time
timeSt = time.clock()

base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
site.addsitedir(os.path.join(base_path, "statsmath"))
site.addsitedir(os.path.join(base_path, "interval"))

# Some key functions ported from R
from statsmath import digamma, trigamma, trigammaInverse, is_infinite

import partition

parser = argparse.ArgumentParser(description=\
    'Take all or some of the unmoderated T-statistics and moderate them .')

parser.add_argument(\
    '--sample', metavar='FRACTION', type=float,
    default=1.1, help='Randomly sample this fraction of input t-stats')
parser.add_argument(\
    '--hmm-overlap', dest='hmmolap', metavar='INT', type=int, default=1000,
    help='Number of observations into previous bin to begin')

partition.addArgs(parser)
args = parser.parse_args()

ninp, nout = 0, 0

def eBayes(tt, df):
    ''' Given list of tstats tuples, moderate them '''
    coef = np.array(map(lambda x: x[0], tt)) # get the coefficients
    stddev = np.array(map(lambda x: x[1], tt)) # get the coefficients
    sg = np.array(map(lambda x: x[2], tt)) # get the sigmas
    sg2 = sg * sg
    sg2med = np.median(sg2)
    sg2 = np.maximum(sg2, 1e-05 * sg2med)
    zg = np.log(sg2)
    dfd2 = df / 2.0
    eg = zg - digamma(dfd2) + math.log(dfd2)
    G = len(sg)
    ebar = np.mean(eg)
    d0 = 2 * trigammaInverse( np.mean( ((eg - ebar) ** 2) * (G/(G-1.0)) - trigamma(dfd2) ) )
    d0d2 = d0 / 2.0
    s02 = math.exp(ebar + digamma(d0d2) - math.log(d0d2))
    sgtilde = None
    if is_infinite(d0): sgtilde = s02
    else: sgtilde = (d0 * s02 + df * sg2) / (d0 + df)
    assert sgtilde is not None
    ttmod = coef / (np.sqrt(sgtilde) * stddev)
    logfchange = coef
    return zip(ttmod, logfchange)

# tts is a list of { lists of 3-tuples }.  One outermost list element
# per test/null permutation.  Inner list is per position.  Elements of
# each 3-tuple are (1) coefficient, (2) standard deviation, (3) sigma.
tts = []
# parallel lists of referene ids and positions
refids, poss = [], []

first = True
df = None
for ln in sys.stdin:
    ln = ln.rstrip()
    toks = ln.split('\t')
    ninp += 1
    # Token is pos(tab)row-mean(tab)degrees-of-freedom(tab)...
    # Where ... is a variable-length list of coef,stddev,sigma triples joined by tabs
    assert len(toks) >= 5
    refid, pos, mn, mydf = toks[0], int(toks[1]), float(toks[2]), int(toks[3])
    if df is None:
        df = mydf
    assert df == mydf
    for i in xrange(4, len(toks)):
        toks2 = toks[i].split(',')
        assert len(toks2) == 3
        coef, stdddev, sigma = float(toks2[0]), float(toks2[1]), float(toks2[2])
        if first: tts.append([(coef, stdddev, sigma)])
        else: tts[i-4].append((coef, stdddev, sigma))
    poss.append(pos)
    refids.append(refid)
    first = False

tts_mod = [ eBayes(ttsx, df) for ttsx in tts ]
# tts_mod is a list of { 2-tuples of { a list of floats, another list of floats } }

for i in xrange(0, len(tts_mod[0])):
    # Make output string
    ttstr_list = []
    for ttup in tts_mod:
        ttmod, logfchange = ttup[i]
        ttstr_list.append("%f,%f" % (ttmod, logfchange))
    # Place the position within a partition
    for pt in partition.partition(refids[i], poss[i], poss[i] + args.hmmolap, partition.binSize(args)):
        print(pt + ("\t%012d\t" % poss[i]) + "\t".join(ttstr_list))
        nout += 1

# Done
timeEn = time.clock()
print >>sys.stderr, "DONE with ebayes.py; in/out = %d/%d; time=%0.3f secs" % (ninp, nout, timeEn-timeSt)
