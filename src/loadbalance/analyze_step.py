"""
analyze_step.py

Analyze the load balancing characteristics of a pipeline step.
"""

import os
import sys
import site
import tempfile
import string

base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
site.addsitedir(os.path.join(base_path, "util"))

import url
import filemover

import argparse

parser = argparse.ArgumentParser(description='Analyzes a directory of files containing load balance information from a pipeline step.')

filemover.addArgs(parser)

parser.add_argument(\
    '--input', metavar='URL', type=str, required=True, help='Directory with partition statistics')
parser.add_argument(\
    '--temp', metavar='PATH', type=str, help='Temporary directory to put any downloaded partition statistics in')
parser.add_argument(\
    '--output', metavar='URL', type=str, required=True, help='Directory to put output files in')

args = parser.parse_args()

mover = filemover.FileMover(args=args)

def analyzePartStatsUrl(u):
    ret = []
    for ln in analyzeUrl(u):
        ln = ln.rstrip()
        if len(ln) == 0:
            continue
        nInBin, secs = string.split(ln, '\t')
        ret.append((int(nInBin), float(secs)))
    return ret

def analyzeReducerStatsUrl(u):
    ret = []
    for ln in analyzeUrl(u):
        ln = ln.rstrip()
        if len(ln) == 0:
            continue
        nbin, ninp, nout, secs = string.split(ln, '\t')
        ret.append((int(nbin), int(ninp), int(nout), float(secs)))
    return ret

def analyzeUrl(u):
    # First go grab all the files
    if u.isLocal():
        d = u.toUrl()
    else:
        print >> sys.stderr, "Getting files in %s" % u.plus('*').toNonNativeUrl()
        tmp = tempfile.mkdtemp()
        mover.get(u.plus('*'), tmp)
        d = tmp
    
    for fn in os.listdir(d):
        fn = os.path.join(d, fn)
        print >> sys.stderr, "  Analyzing file %s" % fn
        with open(fn, 'r') as fh:
            for ln in fh:
                yield ln

inpUrl = url.Url(args.input)

print >> sys.stderr, "Gathering partstats..."
ps = analyzePartStatsUrl(inpUrl.plus('partstats'))
print >> sys.stderr, "Gathering reducerstats..."
rs = analyzeReducerStatsUrl(inpUrl.plus('reducerstats'))

print >> sys.stderr, "Writing %d partstats..." % len(ps)
psFn = os.path.join(args.output, 'partstats.tsv')
with open(psFn, 'w') as fh:
    for p in ps:
        fh.write( '\t'.join(map(str, p)) )
        fh.write('\n')

print >> sys.stderr, "Writing %d reducerstats..." % len(ps)
rsFn = os.path.join(args.output, 'reducerstats.tsv')
with open(rsFn, 'w') as fh:
    for r in rs:
        fh.write( '\t'.join(map(str, r)) )
        fh.write('\n')
