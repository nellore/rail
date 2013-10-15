'''
walk_fit.py
(after normalize_post.py, before ebayes.py)

Again, walk over the readlet intervals in a partition.  This time we emit
per-position *vectors* instead of per-position, per-sample elements.  Also, we
now know normalization factors.  So we have all we need to obtain the
per-position unmoderated t-statistics.

Tab-delimited input tuple columns:
 1. Partition ID for partition overlapped by interval
 2. Interval start
 3. Interval end (exclusive)
 4. Reference ID
 5. Interval count
 6. Sample label

Binning/sorting prior to this step:
 1. Binned by partition
 2. Bins sorted by Interval start

Tab-delimited output tuple columns:
 1. Reference ID
 2. Reference offset (0-based)
 3. Mean A
 4. Degrees of freedom
 5+. Comma-delimited triples of (1) coefficient, (2) standard deviation, (3)
     sigma.  One for the data, and then N more triples for each of the N
     permutations of the data tried.

TODO:
 - Investigate intermediate folder problem.  Need to find a way to regularly
   clean and replace the contents in the folder
 - Do input tuples really need the refid (4th) field?  Can't we recover from
   partition ID?

'''

import os
import sys
import site
import argparse
import math
import string
import numpy
import random
import time
timeSt = time.clock()

base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
site.addsitedir(os.path.join(base_path, "interval"))
site.addsitedir(os.path.join(base_path, "struct"))

import circular
import partition

parser = argparse.ArgumentParser(description=\
    'Take spliced alignments binned and sorted into partitions and '
    'construct coverage vectors for each sample.')

parser.add_argument(\
    '--normals', metavar='PATH', type=str, required=True,
    help='File to read normalization factors from')
parser.add_argument(\
    '--fudge-factor', metavar='PATH', type=float,
    default=32, help='Add this to counts before taking log2')
parser.add_argument(\
    '--permutations', metavar='INT', type=int,
    default=3, help='Number of permutation tests to perform')
parser.add_argument(\
    '--permutations-out', metavar='PATH', type=str, required=False,
    help='File to write permutations to')
parser.add_argument(\
    '--seed', metavar='INT', type=int,
    default=41092, help='Seed for pseudo-random number generator')
parser.add_argument(\
    '--profile', action='store_const', const=True,
    help='Profile the code')

partition.addArgs(parser)
args = parser.parse_args()
binsz = partition.binSize(args)

class LMFitter(object):
    """ Helper class for repeatedly fitting linear models to count
        vectors.  There is one LMFitter per permutation.  Thus, one
        LMFitter is associated with a (possibly permuted) group vector
        and a (possibly permuted) normalization factor vector.  We use
        a LRUCache to retrieve previously-calculated model parameters,
        since we expect certain count vectors (e.g. ones with low total
        count) to be common and this saves time. """
    
    def __init__(self, gs, ns, cacheSz = 1000):
        """ Initialize with permuted group and normalization-factor
            vectors """
        assert len(gs) == len(ns)
        self.gs, self.ns = gs, ns
        
        # Design matrix
        self.D = D = numpy.array([gs, ns, numpy.ones(len(gs))]).T
        # Element of (XTX)-1 useful for std err calculation later
        XTX = numpy.dot(D.T, D)
        XTXi = numpy.linalg.inv(XTX)
        self.XTXi00 = XTXi[0, 0]
        # Degrees of freedom
        self.df = len(gs) - 3
        
        self.useCache = True
        
        # Very simple initial cache: remember last query and answer
        self.prevY, self.prevAns = None, None
        
        # Cache is a doubly-linked list
        # Link layout:     [PREV, NEXT, KEY, RESULT]
        self.root = root = [None, None, None, None]
        self.nmiss, self.nhit = 0, 0
        self.cache = cache = {}
        last = root
        for _ in xrange(cacheSz):
            key = object()
            cache[key] = last[1] = last = [last, root, key, None]
        root[0] = last
    
    def __fit(self, y):
        """ Fit linear model for given count vector y """
        assert len(y) == len(self.gs)
        Amean = numpy.mean(y)
        fit = numpy.linalg.lstsq(self.D, y)
        coefs = fit[0]
        residuals = [ y[i] - (coefs[0] * self.gs[i] + coefs[1] * self.ns[i] + coefs[2]) for i in xrange(len(y)) ]
        rss = sum([r * r for r in residuals])
        sigma = math.sqrt(rss / self.df)
        groupStderr = math.sqrt(self.XTXi00 * (rss / self.df))
        return (Amean, self.df, coefs[0], groupStderr/sigma, sigma)
    
    def fit(self, y):
        """ Fit linear model given count vector y """
        if y == self.prevY:
            return self.prevAns # Cache hit!
        self.prevY = y
        tupy = tuple(y)
        cache = self.cache
        root = self.root
        if self.useCache:
            link = cache.get(tupy)
            if link is not None:
                # Cache hit!
                link_prev, link_next, _, ans = link
                link_prev[1] = link_next
                link_next[0] = link_prev
                last = root[0]
                last[1] = root[0] = link
                link[0] = last
                link[1] = root
                self.nhit += 1
                self.prevAns = ans
                return ans
        # Cache miss
        ans = self.__fit(y)
        if self.useCache:
            root[2] = tupy
            root[3] = ans
            oldroot = root
            root = self.root = root[1]
            root[2], oldkey = None, root[2]
            root[3], _ = None, root[3]
            del cache[oldkey]
            cache[tupy] = oldroot
            self.nmiss += 1
        self.prevAns = ans
        return ans

def go():

    ninp = 0                   # # lines input so far
    nout = 0                   # # lines output so far
    last_pt = "\t"             # id of previous partition
    part_st, part_en = -1, -1  # start/end offsets of current partition
    
    def normals(args):
        """ Reads in --normals file, places all the sample labels in 'labels'
            list and places all normalization factors in 'normals' list.
            Returns them. """
        if not os.path.isfile(args.normals):
            raise RuntimeError("No such --normals file '%s'" % args.normals)
        normals = dict()
        labels = []
        fh = open(args.normals)
        for line in fh:
            line = line.rstrip()
            toks = line.split('\t')
            print>>sys.stderr, toks
            assert len(toks) == 2
            labels.append(toks[0])
            if toks[1] != "NA":
                normals[toks[0]] = int(toks[1])
            else:
                normals[toks[0]] = 0.0
        fh.close()
        return labels, normals
    
    labs, normals = normals(args)
    ls = sorted(labs)
    
    random.seed(args.seed)
    normalsl = [ float(normals[l]) for l in ls ]
    
    def groupName(lab):
        toks = string.split(lab, '-')
        assert len(toks) >= 2
        return toks[1]
    
    def labelMapping(ls):
        idx, gidx = 0, 0
        lab2grp, grp2lab, lab2idx, idx2lab = {}, {}, {}, {}
        for l in ls:
            if l not in lab2idx:
                lab2idx[l] = idx
                idx2lab[idx] = l
                idx += 1
                gr = groupName(l)
                if gr not in lab2grp:
                    lab2grp[gr] = gidx
                    grp2lab[gidx] = gr
                    gidx += 1
        return lab2grp, grp2lab, lab2idx, idx2lab
    
    lab2grp, grp2lab, lab2idx, _ = labelMapping(ls)
    
    verbose = False
    
    idxs = [ lab2idx[x] for x in ls ]
    gs = map(lab2grp.get, [groupName(x) for x in ls])
    
    # Make the permutations
    
    idxs_permutations = []  # permutations of sample indexes
    gs_permutations = []    # permutations of groups indexes
    norm_permutations = []  # permutations of sample norm factors
    
    for _ in xrange(0, args.permutations):
        shuf = zip(idxs, gs, normalsl)
        random.shuffle(shuf)
        idxs_permutations.append(map(lambda x: x[0], shuf))
        gs_permutations.append(map(lambda x: x[1], shuf))
        norm_permutations.append(map(lambda x: x[2], shuf))
    
    # Make the objects that we'll use to fit the linear models for the
    # unpermuted (testFitter) and permuted (permFitters) coverage vectors
    testFitter = LMFitter(gs, normalsl)
    permFitters = [ LMFitter(gs, ns) for gs, ns in \
                    zip(gs_permutations, norm_permutations) ]
    
    # Write out a file describing all the permutations
    if args.permutations_out is not None:
        with open(args.permutations_out, 'w') as fh:
            fh.write(','.join(map(grp2lab.get, gs)) + '\n')
            for g in gs_permutations:
                fh.write(','.join(map(grp2lab.get, g)) + '\n')
    
    def transformCount(c):
        return math.log(c + args.fudge_factor, 2) # base-2 log
    
    def modelFit(refid, st, y, testFitter, permFitters):
        xformy = map(transformCount, y)
        # Fit the unpermuted vector
        mn, df, coef, stdev, sig = testFitter.fit(xformy) 
        fitstrs = [ "%f,%f,%f" % (coef, stdev, sig) ]
        # Fit the npermuted vectors
        for fitter in permFitters:
            _, _, coef, stdev, sig = fitter.fit(xformy)
            fitstrs.append("%f,%f,%f" % (coef, stdev, sig))
        print (("%s\t%012d\t%f\t%d\t" % (refid, st, mn, df)) + '\t'.join(fitstrs))
        return 1
    
    def analyzeWindow(refid, st, cov, testFitter, permFitters):
        tcov = zip(*cov) # transpose
        nout = 0
        for i in xrange(0, len(tcov)):
            if sum(tcov[i]) > 0:
                nout += modelFit(refid, st + i, tcov[i], testFitter, permFitters)
        return nout
    
    def finalize(refid, mcov, testFitter, permFitters):
        st, cov = mcov.finalize()
        if cov is not None:
            return analyzeWindow(refid, st, cov, testFitter, permFitters)
        return 0
    
    maxlen = 100
    nsamp = len(ls)
    mcov = None
    
    for ln in sys.stdin:
        ln = ln.rstrip()
        toks = ln.split('\t')
        assert len(toks) == 6
        pt, st, en, _, weight, lab = toks
        st, en, weight = int(st), int(en), int(weight)
        maxlen = max(en - st, maxlen)
        if pt != last_pt:
            # We moved on to a new partition
            if mcov is not None:
                last_refid, _, _ = partition.parse(last_pt, binsz)
                nout += finalize(last_refid, mcov, testFitter, permFitters)
            _, part_st, part_en = partition.parse(pt, binsz)
            mcov = circular.CircularMultiCoverageBuffer(nsamp, part_st, part_en, maxlen)
            if verbose:
                print >>sys.stderr, "Started partition [%d, %d); first read: [%d, %d)" % (part_st, part_en, st, en)
            assert part_en > part_st
        assert part_st >= 0
        assert part_en > part_st
        sampid = lab2idx[lab]
        covst, covl = mcov.add(sampid, st, en, weight)
        if covl is not None:
            refid, _, _ = partition.parse(pt, binsz)
            nout += analyzeWindow(refid, covst, covl, testFitter, permFitters)
        last_pt = pt
        ninp += 1
    
    if part_st > -1:
        refid, _, _ = partition.parse(pt, binsz)
        nout += finalize(refid, mcov, testFitter, permFitters)
    
    # Done
    timeEn = time.clock()
    print >>sys.stderr, "DONE with walk_fit.py; in/out = %d/%d; time=%0.3f secs" % (ninp, nout, timeEn-timeSt)

if args.profile:
    import cProfile
    cProfile.run('go()')
else:
    go()
