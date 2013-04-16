'''
hmm_params.py
(after ebayes.py, before hmm.py)

Given moderated t-statistics, fit HMM parameters.  This might be run just on
the test statistics, or (rarely) on both test and null statistics.  DERfinder's
behavior is to fit just on the test statistics.

The HMM has 3 hidden states, as discussed in section 3.2 of the DERfinder
paper.  D(l) = 0 is meant for nucleotides where there is basically no
expression in any sample.  D(l) = 1 is meant for nucleotides where there is
expression, but no differential expression.  D(l) = 2 is meant for nucleotides
where there is differential expression.  The HMM emissions are the moderated
t-statistics.  The t-statistics are modeled as being drawn from one of three
gaussians depending on the hidden state.  So that transition matrix is 3 x 3
and the emission distributes are 3 gaussians.

In this stage, we use the R code in source files R/paramHelpers.R,
R/getParams.R, and R/locfdrFit.R.

Tab-delimited input tuple columns:
 1. Partition ID
 2. Reference ID
 3. Reference offset (0-based)
 4+. Comma-delimited pairs of (1) moderated t-staistic, (2) log fold-change.
     One pair for the data, then N more pairs for each of N permutations.

Binning/sorting prior to this step:
 (none)

Tab-delimited output tuple columns:
 1. Either "test" for data, or "nullX" where X=integer for permutations
 2. Comma-delimited HMM parameters
'''

import os
import sys
import argparse

# rpy2 is the glue between Python and R
import rpy2
from rpy2.robjects.packages import importr
import rpy2.robjects as robjects

#
# Parse arguments
#

parser = argparse.ArgumentParser(description=\
    'Given moderated t-statistics, fit HMM parameters.')

parser.add_argument(\
    '--null', action='store_const', const=True,
    help='Get HMM parameters for null as well as test statistics')
parser.add_argument(\
    '--out', metavar='PATH', type=str, required=False, default=None,
    help='File to write output to')

args = parser.parse_args()

#
# Set up R functions
#

R_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "R")
R_paramHelpers = os.path.join(R_dir, 'paramHelpers.R')
R_getParams = os.path.join(R_dir, 'getParams.R')
R_locfdrFit = os.path.join(R_dir, 'locfdrFit.R')
assert os.path.exists(R_getParams)
assert os.path.exists(R_locfdrFit)

R_paramHelpers_fh = open(R_paramHelpers, 'r')
R_paramHelpers_str = R_paramHelpers_fh.read()
R_paramHelpers_fh.close()

R_getParams_fh = open(R_getParams, 'r')
R_getParams_str = R_getParams_fh.read()
R_getParams_fh.close()

R_locfdrFit_fh = open(R_locfdrFit, 'r')
R_locfdrFit_str = R_locfdrFit_fh.read()
R_locfdrFit_fh.close()

# Load the paramHelpers helper functions
robjects.r(R_paramHelpers_str)
# Load the getParams function
robjects.r(R_getParams_str)
# Load the loadfdrFit function
robjects.r(R_locfdrFit_str)

r_getParams = robjects.r['getParams']

ninp, nout = 0, 0

def getParams(ttmods):
    ''' Given moderated t statistics, calculate HMM parameters for use
        in the next step. '''
    res = r_getParams(robjects.FloatVector(ttmods))
    # These are the initial probabilities for the 4 hidden states
    stateprobs = res.rx2('stateprobs')
    # These are the means and standard deviations for the 4 emission
    # probability distributions
    mns = res.rx2('params').rx2('mean')
    sds = res.rx2('params').rx2('sd')
    return tuple(stateprobs + mns + sds)

ttmods = []
poss = []
names = []
first = True

#
# Parse moderated t-statistics into a list of lists.  The first list is
# the list of all moderated test t-statistics.  Subsequent lists are
# lists of moderated null t-statistics for the various permutations.
#
for ln in sys.stdin:
    ln = ln.rstrip()
    toks = ln.split('\t')
    assert len(toks) >= 3
    pt, pos = toks[0], int(toks[1])
    if first:
        names.append("test") # test data
    for i in xrange(2, len(toks)):
        if i > 2 and not args.null:
            break # skip null t-statistics
        toks2 = toks[i].split(',')
        assert len(toks2) == 2
        ttmod, logfchange = float(toks2[0]), float(toks2[1])
        if first:
            ttmods.append([ttmod])
            names.append("null%d" % (i-2)) # permutation data
        else:
            ttmods[i-2].append(ttmod)
    first = False
    poss.append(pos)
    ninp += 1

#
# Apply getParams to each list of moderated t-statistics.
#
params = map(getParams, ttmods)

# For each position, let tup = a tuple of tuples, where the outer
# tuples correspond to sets of test/null statistics and the inner
# tuples are the return values from getParams 
tupstrs = [ ",".join(map(str, tup)) for tup in params ]

ofh = sys.stdout
if args.out is not None:
    ofh = open(args.out, 'w')

for i in xrange(0, len(tupstrs)):
    ofh.write(names[i] + "\t" + tupstrs[i] + "\n")
    nout += 1

if args.out is not None:
    ofh.close()

# Done
print >>sys.stderr, "DONE with hmm_params.py; in/out = %d/%d" % (ninp, nout)
