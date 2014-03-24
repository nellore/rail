#!/usr/bin/env python

import sys
import os

for ln in sys.stdin:
    toks = ln.rstrip().split('\t')
    toks[0] = os.path.abspath(toks[0])
    if len(toks) == 5:
        toks[2] = os.path.abspath(toks[2])
    print '\t'.join(toks)
