"""
bed2raw.py

Turn bed into lines with just position(tab)score, without run-length encoding.
Assumes just one chromosome.  You can optionally specify the number of
nucleotides in the chromosome via the first command-line arg.  If unspecified,
we print to the end of the rightmost bed-file interval.
"""

import sys
import string

maxpos = None
if len(sys.argv) > 1:
    maxpos = int(sys.argv[1])
first = True

for ln in sys.stdin:
    toks = string.split(ln.rstrip(), '\t')
    st, en = int(toks[1]), int(toks[2])
    if first and st > 1:
        for i in xrange(1, st):
            print '%d\t0' % (i-1)
    first = False
    for i in xrange(st, en):
        print '\t'.join([str(i-1), toks[3]])
        if i >= maxpos:
            break
    if en >= maxpos:
        break
