#!/usr/bin/env python

"""
align.py

Alignment script, usually to be used in Hadoop apps.

Wraps Bowtie.  Has features for (1) optionally extracting readlets,
(2) optionally truncating reads or omitting mates, (3) hashing the read
name or partitioning the genome position.
"""

import sys
import os
import site
import argparse
import threading

base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
site.addsitedir(os.path.join(base_path, "bowtie"))
site.addsitedir(os.path.join(base_path, "read"))
site.addsitedir(os.path.join(base_path, "sample"))
site.addsitedir(os.path.join(base_path, "manifest"))

import bowtie
import readlet
import sample
import manifest

ninp = 0               # # lines input so far
nout = 0               # # lines output so far
pe = False
discardMate = None
lengths = dict()       # read lengths after truncation
rawLengths = dict()    # read legnths prior to truncation
qualCnts = dict()      # quality counts after adjustments
rawQualCnts = dict()   # quality counts before adjustments
qualAdd = None         # amt to add to qualities
truncateAmt = None     # amount to truncate reads
truncateTo = None      # amount to truncate reads

readletize = None      # if we're going to readletize, 

xformReads = qualAdd is not None or truncateAmt is not None or truncateTo is not None 

parser = argparse.ArgumentParser(description=\
    'Align reads using Bowtie, usually as the map step in a Hadoop program.')
parser.add_argument(\
    '--binReadName', metavar='NUMBINS', type=int, default=0,
    help='Emit hashed version of read name as first column')
parser.add_argument(\
    '--binPos', metavar='NUMBINS', type=int, default=0,
    help='Emit hashed version of read name as first column')

bowtie.addArgs(parser)
readlet.addArgs(parser)
manifest.addArgs(parser)

args = parser.parse_args()

def xformRead(seq, qual):
    # Possibly truncate
    newseq, newqual = "", ""
    if truncateAmt is not None:
        pass
    if truncateTo is not None:
        pass
    if qualAdd is not None:
        pass
    return newseq, newqual

def bowtieCmd(bowtieExe, bowtieIdx, bowtieArgs):
    # Check that Bowtie exists and is executable
    if not os.path.isfile(bowtieExe):
        print >>sys.stderr, "Bowtie executable '%s' does not exist" % bowtieExe
        sys.exit(1)
    if not os.access(bowtieExe, os.X_OK):
        print >>sys.stderr, "Bowtie executable '%s' exists but is not executable" % bowtieExe
        sys.exit(1)
    # Check that index is there
    for i in [".1.ebwt", ".2.ebwt", ".3.ebwt", ".4.ebwt", ".rev.1.ebwt", ".rev.2.ebwt"]:
        if not os.path.isfile(bowtieIdx + i):
            print >>sys.stderr, "Could not find index file " + (bowtieIdx + i)
            sys.exit(1)
    bowtieCmd = "%s %s --12 --mm %s -" % (bowtieExe, bowtieArgs, bowtieIdx)

bowtieOutDone = threading.Event()
bowtieErrDone = threading.Event()

def bowtieOut(st):
    ''' Process standard out (stdout) output from Bowtie '''
    global nout
    for line in st:
        line = line.rstrip()
        toks, fron = None, None
        if args.binReadName > 0:
            if toks is None:
                toks = line.split('\t')
            nm = abs(hash(toks[0]) % args.binReadName)
            fron = [str(nm)]
        if fron is not None:
            assert toks is not None
            nout += 1
            print '\t'.join(fron) + '\t' + '\t'.join(toks)
        else:
            nout += 1
            print line
    bowtieOutDone.set()

def bowtieErr(st):
    ''' Process standard error (stderr) output from Bowtie '''
    for line in st:
        print >> sys.stderr, line.rstrip()
    bowtieErrDone.set()

proc = bowtie.proc(args, outHandler=bowtieOut, errHandler=bowtieErr)

for ln in sys.stdin:
    ln = ln.rstrip()
    toks = ln.split('\t')
    ninp += 1
    pair = False
    nm, seq, qual = None, None, None
    nm1, seq1, qual1 = None, None, None
    nm2, seq2, qual2 = None, None, None
    if len(toks) == 3:
        # Unpaired read
        nm, seq, qual = toks
        sample.hasLab(nm, mustHave=True)
    elif len(toks) == 5 or len(toks) == 6:
        # Paired-end read
        if len(toks) == 5:
            # 6-token version
            nm1, seq1, qual1, seq2, qual2 = toks
            nm2 = nm1
        else:
            # 5-token version
            nm1, seq1, qual1, nm2, seq2, qual2 = toks
        sample.hasLab(nm1, mustHave=True)
        if discardMate is not None:
            if discardMate == 1:
                nm, seq, qual = nm2, seq2, qual2
            else:
                nm, seq, qual = nm1, seq1, qual1
        else:
            pair = True
    else:
        raise RuntimeError("Wrong number of tokens for line: " + ln)
    
    if pair:
        # Paired-end
        if xformReads:
            # Truncate and transform quality values
            seq1, qual1 = xformRead(seq1, qual1)
            seq2, qual2 = xformRead(seq2, qual2)
        proc.stdin.write("%s\t%s\t%s\t%s\t%s\n" % (nm1, seq1, qual1, seq2, qual2))
    else:
        # Unpaired
        if xformReads:
            # Truncate and transform quality values
            seq, qual = xformRead(seq, qual)
        if args.readletLen > 0:
            # Readletize
            for rlet in iter(readlet.readletize(args, nm, seq, qual)):
                nm_rlet, seq_rlet, qual_rlet = rlet
                proc.stdin.write("%s\t%s\t%s\n" % (nm_rlet, seq_rlet, qual_rlet))
        else:
            proc.stdin.write("%s\t%s\t%s\n" % (nm, seq, qual))

# Close and flush STDIN.
proc.stdin.close()

# Wait until the threads processing stdout and stderr are done.

# Close stdout and stderr
print >>sys.stderr, "Waiting for Bowtie stdout processing thread to finish"
bowtieOutDone.wait()
print >>sys.stderr, "Waiting for Bowtie stderr processing thread to finish"
bowtieErrDone.wait()

proc.stdout.close()
proc.stderr.close()

# Done
print >>sys.stderr, "DONE with align.py; in/out = %d/%d" % (ninp, nout)
