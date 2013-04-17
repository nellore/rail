'''
hmm.py
(after hmm_params.py)

Given the moderated t-statistics and HMM parameters, run HMM on the original
data, as well as on all the permuted datasets.

(not yet implemented)

Tab-delimited input tuple columns:
 1. Partition ID
 2. Reference offset (0-based)
 3+. Comma-delimited pairs of (1) moderated t-staistic, (2) log fold-change.
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

import os
import site
import sys
import ghmm
import argparse
import numpy
import string
import time
timeSt = time.clock()

base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
site.addsitedir(os.path.join(base_path, "interval"))

import partition

parser = argparse.ArgumentParser(description=\
    'Given HMM parameters and sorted bins of moderated t-statistics, '
    'run the HMM and emit state string.')

parser.add_argument(\
    '--params', metavar='PATH', type=str, required=True,
    help='Parameter file for HMM')
parser.add_argument(\
    '--hmm-overlap', dest='hmmolap', metavar='INT', type=int, default=1000,
    help='Number of observations into previous bin to begin')

partition.addArgs(parser)
args = parser.parse_args()

ninp, nout = 0, 0

class HMM:
    def __init__(self, probs, transProb=(0.999, 1e-12)):
        """ Given: initial probabilities, a 2-value summary of the transition
            probabilities, and the parameters for the gaussians modeling the
            observed processes, calculate the complete set of parameters we
            will use to set up the HMM. """
        
        assert len(probs) == 12
        initial = probs[0:4]
        means = probs[4:8]
        sds = probs[8:12]
        stayprob, EtoDE = transProb
        EtoZero = 1.0 - stayprob - 2.0 * EtoDE
        stayrem = (1.0 - stayprob) / 3.0
        A = [[ stayprob,  stayrem,  stayrem,  stayrem],
             [ EtoZero,  stayprob,    EtoDE,    EtoDE],
             [ EtoZero,     EtoDE, stayprob,    EtoDE],
             [ EtoZero,     EtoDE,    EtoDE, stayprob]]
        
        # Emission parameters
        E = map(list, zip(means, sds))
        print >> sys.stderr, E
        
        F = ghmm.Float()  # emission domain of this model
        self.hmm = ghmm.HMMFromMatrices(\
            F,
            ghmm.GaussianDistribution(F),
            A,
            E,
            initial)
        self.A, self.E, self.I = A, E, initial
    
    def viterbi(self, x):
        return self.hmm.viterbi(ghmm.EmissionSequence(ghmm.Float(), x))
    
    def __str__(self):
        return 'A:\n%s\nE:\n%s\nI:\n%s' % (self.A, self.E, self.I)
    
    def __repr__(self):
        return str(self)

def parseParams(fn):
    """ Parse HMM parameters """
    hmms = []
    with open(fn, 'r') as fh:
        for ln in fh:
            ln = ln.rstrip()
            assert ln.count('\t') == 1
            _, params = string.split(ln, '\t')
            assert params.count(',') == 11
            # 12 params: 4 initial probs, 4 emission means, 4 emission stddevs
            params = string.split(params, ',')
            hmms.append(HMM(map(float, params)))
    return hmms

# Parse all the HMM parameters and build HMMs using ghmm
hmms = parseParams(args.params)

def writeResults(res, refid, ofh):
    for st, off in res:
        ofh.write("%s\t%d\t%s\n" % (refid, off, st))

def processPartition(tts, offs, hmm, st, en):
    """ Run Viterbi on a sequence of moderated t-statistics """
    denseTts = []
    cur = 0
    #print >> sys.stderr, (st, en, offs)
    for i in xrange(st - args.hmmolap, en):
        if cur >= len(offs) or i < offs[cur]:
            denseTts.append(0.0)
        else:
            assert i == offs[cur]
            denseTts.append(tts[cur])
            cur += 1
    assert cur == len(tts)
    results = []
    path, _ = hmm.viterbi(denseTts)
    assert len(path) == en-st + args.hmmolap, "%s" % str(path)
    cur = args.hmmolap
    for i in xrange(st, en):
        results.append((path[cur], i))
        cur += 1
    return results

tts, offs = [], []
refid, st, en = None, None, None
lastPartid = None
for ln in sys.stdin:
    ln = ln.rstrip()
    #print >>sys.stderr, ln
    toks = ln.split('\t')
    assert len(toks) >= 3
    partid, refoff = toks[0], toks[1]
    refoff = int(refoff)
    first = False
    startedPartition = lastPartid is None or partid != lastPartid
    finishedPartition = lastPartid is not None and partid != lastPartid
    if finishedPartition:
        for i in xrange(0, len(tts)):
            writeResults(processPartition(tts[i], offs[i], hmms[i], st, en), refid, sys.stdout)
        tts, offs = [], []
    if startedPartition:
        refid, st, en = partition.parse(partid, partition.binSize(args))
        first = True
    for i in xrange(2, len(toks)):
        assert toks[i].count(',') == 1
        tt, logf = string.split(toks[i], ',')
        tt, logf = float(tt), float(logf)
        ii = i - 2
        if first:
            tts.append([])
            offs.append([])
            assert len(tts) == i-1
        tts[i-2].append(tt)
        offs[i-2].append(refoff)
    lastPartid = partid
    ninp += 1

if lastPartid is not None:
    for i in xrange(0, len(tts)):
        writeResults(processPartition(tts[i], offs[i], hmms[i], st, en), refid, sys.stdout)

# Done
timeEn = time.clock()
print >>sys.stderr, "DONE with hmm.py; in/out = %d/%d; time=%0.3f secs" % (ninp, nout, timeEn-timeSt)
