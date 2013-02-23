'''
hmm.py
(after hmm_params.py)

Given the moderated t-statistics and HMM parameters, run HMM on the original
data, as well as on all the permuted datasets.

(not yet implemented)

Tab-delimited input tuple columns:
 1. Partition ID
 2. Reference ID
 3. Reference offset (0-based)
 4+. Comma-delimited pairs of (1) moderated t-staistic, (2) log fold-change.
     One pair for the data, then N more pairs for each of N permutations.

Binning/sorting prior to this step:
 1. Binned by partition ID
 2. Sorted by reference offset

Tab-delimited output tuple columns:
 1. Reference ID
 2. Reference offset (0-based)
 3+. HMM state.  One column for the original data, then N more columns for each
     of the N permutations.

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
