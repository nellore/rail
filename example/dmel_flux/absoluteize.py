#!/usr/bin/env python

import sys
import os

for ln in sys.stdin:
    toks = ln.rstrip().split('\t')
    toks[0] = os.path.abspath(toks[0])
    print '\t'.join(toks)
