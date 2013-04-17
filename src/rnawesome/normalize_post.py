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

import manifest

parser = argparse.ArgumentParser(description=\
    'Take results from normalization phase and write them to a file.')

parser.add_argument(\
    '--out', metavar='PATH', type=str, required=False, default=None,
    help='File to write output to')

manifest.addArgs(parser)
args = parser.parse_args()

# Get the set of all labels by parsing the manifest file, given on the
# filesystem or in the Hadoop file cache
labs = manifest.labels(args)
ls = sorted(labs)

ninp = 0       # # lines input so far
facts = dict() # Normalization factors

for ln in sys.stdin:
    ln = ln.rstrip()
    toks = ln.split('\t')
    assert len(toks) == 2
    facts[toks[0]] = int(toks[1])
    ninp += 1

ofh = sys.stdout
if args.out is not None:
    ofh = open(args.out, 'w')

for l in ls:
    if l in facts: ofh.write("%s\t%d\n" % (l, facts[l]))
    else: ofh.write("%s\tNA\n" % l)

if args.out is not None:
    ofh.close()    

# Done
timeEn = time.clock()
print >>sys.stderr, "DONE with normalize_post.py; in = %d; time=%0.3f secs" % (ninp, timeEn-timeSt)
