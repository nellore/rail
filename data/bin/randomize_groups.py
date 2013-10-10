"""
randomize_groups.py
"""

import os
import sys
import string
import time
import subprocess
import argparse
import random

parser = argparse.ArgumentParser(description=\
    'Scan manifest file and assign new, randomized group labels.')
parser.add_argument(\
    '--manifest', type=str, required=True, help='Manifest file')
parser.add_argument(\
    '--num-groups', type=int, default=2, help='Number of randomized group labels to use')
parser.add_argument(\
    '--seed', type=int, default=774, help='Seed for random generator')

args = parser.parse_args()

random.seed(args.seed)

# TODO: Try to assign the same group label to records that started with the same name
with open(args.manifest) as fh:
    print >> sys.stderr, 'Processing manifest file "%s"' % args.manifest
    for ln in fh:
        ln = ln.rstrip()
        if len(ln) == 0 or ln[0] == '#':
            print ln
            continue
        toks = string.split(ln, '\t')
        assert len(toks) == 3 or len(toks) == 5
        nametoks = string.split(toks[-1], '-')
        newname = toks[-1].replace('-', '_')
        assert len(nametoks)
        name = '-'.join(['random_%d' % random.randint(1, args.num_groups), newname, '1'])
        print '\t'.join(toks[:-1]) + '\t' + name
