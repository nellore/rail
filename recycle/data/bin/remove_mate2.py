"""
remove_mate2.py

Remove all the mate-2 URLs from a manifest file.
"""

import sys
import string
import argparse

parser = argparse.ArgumentParser(description=\
    'Scan manifest file, summarize it, and check for obvious mistakes.')
parser.add_argument(\
    '--manifest', type=str, required=True, help='Manifest file')

args = parser.parse_args()

with open(args.manifest) as fh:
    print >> sys.stderr, 'Processing manifest file "%s"' % args.manifest
    for ln in fh:
        ln = ln.rstrip()
        if len(ln) == 0 or ln[0] == '#':
            print ln
            continue
        toks = string.split(ln, '\t')
        if len(toks) == 3:
            print ln
        else:
            print '\t'.join([toks[0], toks[1], toks[-1]])
