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
 3. Length of run of positions with this state
 4. HMM state
 5. Dataset id (0 for test, >0 for permutations)

'''

import os
import site
import sys
#import ghmm
import argparse
import numpy
import string
import time
#import scipy.stats
import math
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
parser.add_argument(\
    '--profile', action='store_const', const=True,
    help='Profile the code')

partition.addArgs(parser)
args = parser.parse_args()

ninp, nout = 0, 0

class HMM:
    """ Encapsulates a simple HMM for finding regions of differential
        expression.  Emissions are moderated t-statistics.  States
        indicate differential expression. """
    
    def __init__(self, probs, transProb=(0.999, 1e-12)):
        """ Given: initial probabilities, a 2-value summary of the transition
            probabilities, and the parameters for the gaussians modeling the
            observed processes, calculate the complete set of parameters we
            will use to set up the HMM. """
        
        assert len(probs) == 12
        initial, means, sds = probs[:4], probs[4:8], probs[8:]
        stayprob, EtoDE = transProb
        EtoZero = 1.0 - stayprob - 2.0 * EtoDE
        stayrem = (1.0 - stayprob) / 3.0
        
        # Transition probabilities
        A = [[ stayprob,  stayrem,  stayrem,  stayrem],
             [ EtoZero,  stayprob,    EtoDE,    EtoDE],
             [ EtoZero,     EtoDE, stayprob,    EtoDE],
             [ EtoZero,     EtoDE,    EtoDE, stayprob]]
        
        # Emission probabilities
        E = map(list, zip(means, sds))
        #self.Edists = [ scipy.stats.norm(E[i][0], E[i][1]) for i in xrange(0, 4) ]
        
        #F = ghmm.Float()  # emission domain of this model
        # debugFile = open("debug.txt","w")
        # debugFile.write("Initialized emission and transition probabilities\n")
        # self.hmm = ghmm.HMMFromMatrices(
        #     F,
        #     ghmm.GaussianDistribution(F),
        #     A,
        #     E,
        #     initial)

        self.A, self.E, self.I = A, E, initial

    def Elog(self,i,x): #Emission probability of t-statistic
        #prob = self.Edists[i].pdf(x)
        m = self.E[i][0]
        s = self.E[i][1]
        prob = (1.0/(math.sqrt(2*math.pi*s*s)))*math.exp(-(x-m)*(x-m)/(2*s*s))
        if(prob==0):
            return float("-inf")
        return math.log(prob,2)

    def Ilog(self,i): #Initial probability of t-statistic
        if(self.I[i]==0):
            return float("-inf")
        return math.log(self.I[i],2)

    def Alog(self,i,j):
        if(self.A[i][j]==0):
            return float("-inf")
        return math.log(self.A[i][j],2)

    def viterbi(self, x):
        """ Return most likely sequence of states through the HMM given
            sequence of moderated t-statistics (x) """
        # debugFile = open("debug.txt","a")
        # debugFile.write("Started Viterbi algorithm\n")
        # debugFile.write(str(self.hmm)+'\n')

        #return self.hmm.viterbi(ghmm.EmissionSequence(ghmm.Float(), x))
    
        """ Given sequence of emissions, return the most probable path
        along with its log probability.  Do all calculations in log
        space to avoid underflow. """
        #x = map(self.smap.get, x) # turn emission characters into ids
        nrow, ncol = len(self.E), len(x)  #nrows = states, ncols = observations
        mat   = numpy.zeros(shape=(nrow, ncol), dtype=float) # prob
        matTb = numpy.zeros(shape=(nrow, ncol), dtype=int)   # traceback
        #print >> sys.stderr, (self.Alog, self.Elog, self.Ilog)
        # Fill in first column
        for i in xrange(0, nrow):
            mat[i, 0] = self.Elog(i, x[0]) + self.Ilog(i)
        # Fill in rest of log prob and Tb tables
        for j in xrange(1, ncol):
            for i in xrange(0, nrow):
                ep = self.Elog(i, x[j])
                mx, mxi = mat[0, j-1] + self.Alog(0, i) + ep, 0
                for i2 in xrange(1, nrow):
                    pr = mat[i2, j-1] + self.Alog(i2, i) + ep
                    if pr > mx:
                        mx, mxi = pr, i2
                mat[i, j], matTb[i, j] = mx, mxi
        # Find final state with maximal log probability
        omx, omxi = mat[0, ncol-1], 0
        for i in xrange(1, nrow):
            if mat[i, ncol-1] > omx:
                omx, omxi = mat[i, ncol-1], i
        # Traceback
        i, p = omxi, [omxi]
        for j in xrange(ncol-1, 0, -1):
            i = matTb[i, j]
            p.append(i)
        #p = map(lambda x: self.Q[x], p[::-1])
        #print "Path",p
        return p, omx # Return path and log probability

    def __str__(self):
        return 'A:\n%s\nE:\n%s\nI:\n%s' % (self.A, self.E, self.I)
    
    def __repr__(self):
        return str(self)

def parseParams(fn):
    """ Parse HMM parameters """
    hmms = []
    # debugFile = open("debug.txt","a")
    # debugFile.write("Reading HMM parameters\n")
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

def writeResults(res, refid, dataset, ofh):
    """ Given HMM path through a partition, output intervals """
    global nout
    if len(res) == 0:
        return
    sti, offi = res[0]
    for st, off in res[1:]:
        if st != sti:
            assert off <= args.genomeLen
            ofh.write("%s\t%d\t%d\t%s\t%d\n" % (refid, offi, off - offi, sti, dataset))
            sti, offi = st, off
    assert off <= args.genomeLen, "off=%d, args.genomeLen=%d" % (off, args.genomeLen)
    off = res[-1][1] + 1
    ofh.write("%s\t%d\t%d\t%s\t%d\n" % (refid, offi, off - offi, sti, dataset))
    nout += 1

def processPartition(tts, offs, hmm, st, en):
    """ Run Viterbi on a sequence of moderated t-statistics """
    denseTts = []
    cur = 0
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
    assert len(path) == en-st + args.hmmolap
    if en > args.genomeLen:
        en = args.genomeLen
    cur = args.hmmolap
    for i in xrange(st, en):
        results.append((path[cur], i))
        cur += 1
    return results

def go():
    global ninp, nout
    
    # Parse all the HMM parameters and build HMMs using ghmm
    hmms = parseParams(args.params)
    
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
        assert refoff < args.genomeLen, "refoff=%d, args.genomeLen=%d" % (refoff, args.genomeLen)
        first = False
        startedPartition = lastPartid is None or partid != lastPartid
        finishedPartition = lastPartid is not None and partid != lastPartid
        if finishedPartition:
            for i in xrange(0, len(tts)):
                writeResults(processPartition(tts[i], offs[i], hmms[i], st, en), refid, i, sys.stdout)
            tts, offs = [], []
        if startedPartition:
            refid, st, en = partition.parse(partid, partition.binSize(args))
            first = True
        for i in xrange(2, len(toks)):
            assert toks[i].count(',') == 1
            tt, logf = string.split(toks[i], ',')
            tt, logf = float(tt), float(logf)
            if first:
                tts.append([]); offs.append([])
                assert len(tts) == i-1
            tts[i-2].append(tt)
            offs[i-2].append(refoff)
        lastPartid = partid
        ninp += 1
    
    if lastPartid is not None:
        for i in xrange(0, len(tts)):
            writeResults(processPartition(tts[i], offs[i], hmms[i], st, en), refid, i, sys.stdout)

if args.profile:
    import cProfile
    cProfile.run('go()')
else:
    go()

# Done
timeEn = time.clock()
print >>sys.stderr, "DONE with hmm.py; in/out = %d/%d; time=%0.3f secs" % (ninp, nout, timeEn-timeSt)
