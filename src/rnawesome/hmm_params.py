'''
hmm_params.py

Given moderated t-statistics, fit HMM parameters.  This might be run
just on the test statistics, or (rarely) on test and null statistics.
'''

import os
import sys
import argparse
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
R_getParams = os.path.join(R_dir, 'getParams.R')
R_locfdrFit = os.path.join(R_dir, 'locfdrFit.R')
assert os.path.exists(R_getParams)
assert os.path.exists(R_locfdrFit)

R_getParams_fh = open(R_getParams, 'r')
R_getParams_str = R_getParams_fh.read()
R_getParams_fh.close()

R_locfdrFit_fh = open(R_locfdrFit, 'r')
R_locfdrFit_str = R_locfdrFit_fh.read()
R_locfdrFit_fh.close()

# Load the getParams fucntion
robjects.r(R_getParams_str)
# Load the loadfdrFit fucntion
robjects.r(R_locfdrFit_str)

r_getParams = robjects.r['getParams']

ninp, nout = 0, 0

def getParams(ttmods):
    
    ''' Given moderated t statistics, calculate HMM parameters for use
        in the next step. '''
    #res = r_getParams(robjects.FloatVector(ttmods))
    return (0, 1, 2, 3)

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
    assert len(toks) >= 2
    pos = int(toks[0])
    if first:
        names.append("test")
    for i in xrange(1, len(toks)):
        if i > 1 and not args.null:
            break # skip null t-statistics
        toks2 = toks[i].split(',')
        assert len(toks2) == 2
        ttmod, logfchange = float(toks2[0]), float(toks2[1])
        if first:
            ttmods.append([ttmod])
            names.append("null%d" % i)
        else:
            ttmods[i-1].append(ttmod)
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
