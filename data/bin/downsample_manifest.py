"""
downsample_manifest.py

Remove some fraction of the entries from a manifest file.
"""

import sys
import argparse
import random

parser = argparse.ArgumentParser(description=\
    'Scan manifest file, summarize it, and check for obvious mistakes.')
parser.add_argument(\
    '--manifest', type=str, required=True, help='Manifest file')
parser.add_argument(\
    '--sample', type=float, required=True, help='Fraction to sample')
parser.add_argument(\
    '--seed', type=int, default=844, help='Seed for random generator')

args = parser.parse_args()
random.seed(args.seed)

with open(args.manifest) as fh:
    print >> sys.stderr, 'Processing manifest file "%s"' % args.manifest
    for ln in fh:
        ln = ln.rstrip()
        if len(ln) == 0 or ln[0] == '#':
            print ln
            continue
        if random.random() < args.sample:
            print ln
