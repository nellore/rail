'''
ebayes.py

Given all the unmoderated t-statistics and others results from the
linear fits, calculate the moderated t-statistics.  In general, we need
to do this both for the test statistics, and for zero or more sets of
null statistics, each corresponding to a permutation of the group
labels.
'''

import argparse
import random
import sys
import os
import site

base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
site.addsitedir(os.path.join(base_path, "statsmath"))

import rpy2
from rpy2.robjects.packages import importr
import rpy2.robjects as robjects

import statsmath

parser = argparse.ArgumentParser(description=\
    'Take all or some of the unmoderated T-statistics and moderate them .')

parser.add_argument(\
    '--sample', metavar='FRACTION', type=float,
    default=1.1, help='Randomly sample this fraction of input t-stats')

ninp, nout = 0, 0

def moderate(t):
    ''' Moderate a single t statistic '''
    pass

def getTstats(tt):
    ''' Given list of tstats tuples, moderate them '''
    pass

# We need some special math functions to do the moderation, including digamma, trigamma and 

tts = []

first = True
for ln in sys.stdin:
    ln = ln.rstrip()
    toks = ln.split('\t')
    ninp += 1
    # Token is pos(tab)row-mean(tab)degrees-of-freedom(tab)...
    # Where ... is a variable-length list of coef,stddev,sigma triples joined by tabs
    assert len(toks) >= 4
    pos, mn, df = int(toks[0]), float(toks[1]), int(toks[2])
    for i in xrange(3, len(toks)):
        toks2 = toks[i].split(',')
        assert len(toks2) == 3
        coef, stdddev, sigma = float(toks2[0]), float(toks2[1]), float(toks2[2])
        if first:
            tts.append([(coef, stdddev, sigma)])
        else:
            tts[i-3].append((coef, stdddev, sigma))
    first = False

mod_tts = [ getTstats(ttsx) for ttsx in tts ]

# Done
print >>sys.stderr, "DONE with ebayes.py; in/out = %d/%d" % (ninp, nout)
