"""
simple.py

A very simple RNA-seq simulator.
"""

import sys
import random
import bisect

def makeRef(ln, stayprob, stateFn):
    stayrem = (1.0 - stayprob) / 3.0
    EtoZero = (1.0 - stayprob) / 2.0
    EtoDE = (1.0 - stayprob) / 4.0
    
    # Transition matrix
    A = [[ stayprob,  stayrem,  stayrem,  stayrem],
         [ EtoZero,  stayprob,    EtoDE,    EtoDE],
         [ EtoZero,     EtoDE, stayprob,    EtoDE],
         [ EtoZero,     EtoDE,    EtoDE, stayprob]]
    
    def transition(state):
        """ Transition from 'state' to next state in time series using
            HMM transition matrix """
        assert state >= 0 and state <= 3
        rnd = random.random()
        assert rnd >= 0.0 and rnd <= 1.0
        tot = 0.0
        for i in xrange(0, 4):
            tot += A[state][i]
            if rnd < tot:
                return i
        return random.randint(0, 3)
    
    state = 0
    coverage, path = [], []
    cov1, cov2 = 0.001, 0.001
    with open(stateFn, 'w') as fh:
        for i in xrange(0, ln):
            newstate = transition(state)
            if newstate != state:
                if newstate == 0:
                    cov1, cov2 = 0.001, 0.001
                elif newstate == 1:
                    cov1 = random.random()
                    cov2 = 1.1
                    while cov2 >= cov1:
                        cov2 = random.random()
                elif newstate == 2:
                    cov1 = random.random()
                    cov2 = cov1
                else:
                    cov2 = random.random()
                    cov1 = 1.1
                    while cov1 >= cov2:
                        cov1 = random.random()
            state = newstate
            path.append("0dED"[state])
            fh.write("0dED"[state])
            if (i+1) % 70 == 0:
                fh.write('\n')
            coverage.append((cov1, cov2))

    return path, coverage

def simulate(ref, readlen, targetNucs, stateFn):
    _, coverage = makeRef(len(ref), 0.999, stateFn)
    nucs = 0
    seqs1, seqs2 = [], []
    print >> sys.stderr, "Generating..."
    passi = 1
    while nucs < targetNucs:
        print >> sys.stderr, "  Pass %d..." % passi
        for i in xrange(0, len(coverage)):
            r1, r2 = random.random(), random.random()
            if r1 < coverage[i][0]:
                seqs1.append(ref[i:i+readlen])
                nucs += readlen
            if r2 < coverage[i][1]:
                seqs2.append(ref[i:i+readlen])
                nucs += readlen
        passi += 1
    print >> sys.stderr, "Shuffling"
    random.shuffle(seqs1)
    random.shuffle(seqs2)
    return seqs1, seqs2

class WeightedRandomGenerator(object):

    def __init__(self, weights):
        self.totals = []
        running_total = 0
        for w in weights:
            running_total += w
            self.totals.append(running_total)

    def next(self):
        rnd = random.random() * self.totals[-1]
        return bisect.bisect_right(self.totals, rnd)

    def __call__(self):
        return self.next()

def replicateize(seqs1, seqs2, nreps):
    """ Take all the sequence reads for groups 1 and 2 and split them into into
        a bunch of replicates with varying average coverage. """
    weights1 = [1.0] * nreps
    weights2 = [1.0] * nreps
    for i in xrange(0, nreps):
        weights1[i] += (random.random() - 0.5)
        weights2[i] += (random.random() - 0.5)
    # Each replicate has a weight, so now we assign sequences to replicates
    seqs1rep = [ [] for _ in xrange(0, nreps) ]
    seqs2rep = [ [] for _ in xrange(0, nreps) ]
    gen1 = WeightedRandomGenerator(weights1)
    gen2 = WeightedRandomGenerator(weights2)
    for seq in seqs1:
        i = gen1.next()
        seqs1rep[i].append(seq)
    for seq in seqs2:
        i = gen2.next()
        seqs2rep[i].append(seq)
    return seqs1rep, seqs2rep

def writeReads(seqs1rep, seqs2rep, fnPre, manifestFn):
    """ Only unpaired for now """
    seqsrep = [ seqs1rep, seqs2rep ]
    with open(manifestFn, 'w') as manFh:
        for group in xrange(0, 2):
            sr = seqsrep[group]
            for rep in xrange(0, len(sr)):
                fn = "%s.group%d.rep%d.tab5" % (fnPre, group, rep)
                manFh.write("%s\t0\tsimple-%d-%d\n" % (fn, group, rep))
                with open(fn, 'w') as fh:
                    for i in xrange(0, len(sr[rep])):
                        seq = sr[rep][i]
                        qual = "I" * len(seq)
                        nm = "r_group%d_rep%d_n%d" % (group, rep, i)
                        fh.write("%s\t%s\t%s\n" % (nm, seq, qual))

def parseFasta(fns):
    lns = []
    for fn in fns:
        with open(fn, 'r') as fh:
            for ln in fh:
                if ln[0] == '>':
                    continue
                else:
                    lns.append(ln.rstrip())
    return ''.join(lns)

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description=\
        'Simple RNA-seq simulation.')
    
    parser.add_argument(\
        '--fasta', metavar='path', type=str, nargs='+', required=True,
        help='FASTA file(s) containing reference genome sequences')
    parser.add_argument(\
        '--output-prefix', metavar='path', type=str, required=True,
        help='Prefix for output read files')
    parser.add_argument(\
        '--read-len', metavar='int', action='store', type=int, default=100,
        help='Read length to simulate')
    parser.add_argument(\
        '--num-replicates', metavar='int', action='store', type=int, default=8,
        help='Number of replicates per group')
    parser.add_argument(\
        '--num-nucs', metavar='int', action='store', type=float, default=1e7,
        help='Number of total nucleotides of reads to generate')

    args = parser.parse_args()
    
    seqs1, seqs2 = simulate(parseFasta(args.fasta), args.read_len, args.num_nucs, args.output_prefix + ".states")
    print >> sys.stderr, "  reads in group 1: %d" % len(seqs1)
    print >> sys.stderr, "  reads in group 2: %d" % len(seqs2)
    print >> sys.stderr, "Replicatizing..."
    seqs1rep, seqs2rep = replicateize(seqs1, seqs2, args.num_replicates)
    print >> sys.stderr, "  reads in group 1 replicates: " + str(map(len, seqs1rep))
    print >> sys.stderr, "  reads in group 2 replicates: " + str(map(len, seqs2rep))
    print >> sys.stderr, "Writing..."
    writeReads(seqs1rep, seqs2rep, args.output_prefix, args.output_prefix + ".manifest")
