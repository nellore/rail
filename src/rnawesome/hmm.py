'''
hmm.py

Given the moderated t-statistics, run the HMM.
'''

import sys
import ghmm
import argparse

parser = argparse.ArgumentParser(description=\
    'Given HMM parameters and sorted bins of moderated t-statistics, '
    'run the HMM and emit state string.')

parser.add_argument(\
    '--params', metavar='PATH', type=str, required=True,
    help='Randomly sample this fraction of input t-stats')
parser.add_argument(\
    '--hmm-overlap', dest='hmmolap', metavar='INT', type=int, default=1000,
    help='Number of observations into previous bin to begin')

args = parser.parse_args()

ninp, nout = 0, 0

for ln in sys.stdin:
    ln = ln.rstrip()
    toks = ln.split('\t')
    assert len(toks) >= 4
    ninp += 1

# Done
print >>sys.stderr, "DONE with hmm.py; in/out = %d/%d" % (ninp, nout)
