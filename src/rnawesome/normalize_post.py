'''
normalize_post.py
(after normalize.py, before walk_fit.py)

Takes results from the normalize phase as a single bin and writes out a
new manifest file annotated with normalization factors.
 
Tab-delimited input tuple columns (just 1 tuple per sample label):
 1. Sample label
 2. Summary statistic

Binning/sorting prior to this step:
 (none)

Tab-delimited output tuple columns:
 (no Hadoop output - just file)

Other output:
 --out file gets a table, with a row for each sample and columns:
 1. Sample name
 2. Normalization factor
'''

import os
import sys
import site
import argparse
import time
timeSt = time.clock()
 
base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
site.addsitedir(os.path.join(base_path, "manifest"))
site.addsitedir(os.path.join(base_path, "util"))

import manifest
import url
import filemover

parser = argparse.ArgumentParser(description=\
    'Take results from normalization phase and write them to a file.')

parser.add_argument(\
    '--out', metavar='URL', type=str, required=False, default=None,
    help='URL to write output to.  Goes to stdout by default.')
parser.add_argument(\
    '--verbose', action='store_const', const=True, default=False,
    help='Prints out extra debugging statements')

manifest.addArgs(parser)
filemover.addArgs(parser)

args = parser.parse_args()

if args.verbose:
    print >> sys.stderr, " vvv Manifest file vvv"
    with open(args.manifest) as fh: print >> sys.stderr, fh.read()
    print >> sys.stderr, " ^^^ Manifest file ^^^"

# Get the set of all labels by parsing the manifest file, given on the
# filesystem or in the Hadoop file cache
labs = manifest.labels(args)
ls = sorted(labs)

ninp = 0       # # lines input so far
facts = dict() # Normalization factors

for ln in sys.stdin:
    ln = ln.rstrip()
    if len(ln) == 0: continue
    toks = ln.split('\t')
    assert len(toks) == 2
    facts[toks[0]] = int(toks[1])
    ninp += 1

ofh = sys.stdout
outFn, outUrl = None, None

if args.out is not None:
    # If --out is a local file, just write directly to that file.  Otherwise,
    # write to a temporary file that we will later upload to destination.
    outUrl = url.Url(args.out)
    if outUrl.isLocal():
        try: os.makedirs(outUrl.toUrl())
        except: pass
        outFn = os.path.join(args.out, 'normalization_factors.tsv')
    else:
        outFn = "normalize_post.temp"
    ofh = open(outFn, 'w')

for l in ls:
    if l in facts: ofh.write("%s\t%d\n" % (l, facts[l]))
    else: ofh.write("%s\tNA\n" % l)

if args.out is not None:
    ofh.close()
    if not outUrl.isLocal():
        mover = filemover.FileMover(args=args)
        mover.put(outFn, outUrl.plus("normalization_factors.tsv"))
        os.remove(outFn)

# Done
timeEn = time.clock()
print >>sys.stderr, "DONE with normalize_post.py; in = %d; time=%0.3f secs" % (ninp, timeEn-timeSt)
