"""
Converts splice sites into bed format

Doesn't take into consideration sample label atm
"""

import sys

from collections import defaultdict
from operator import itemgetter
import counter

gpositions = defaultdict(counter.Counter)

for ln in sys.stdin:
    toks = ln.rstrip().split("\t")
    if len(toks) < 7:
        continue
    refid, st, en, motifl, motifr, lab, freq = toks[0], int(toks[1]), int(toks[2]),toks[3], toks[4], toks[5], int(toks[6])
    gpositions[refid][str(st)] += freq
    gpositions[refid][str(st+1)] += freq
    gpositions[refid][str(en-2)] += freq
    gpositions[refid][str(en-1)] += freq

for d,c in gpositions.iteritems():
    glist = [(int(pos),freq) for pos,freq in c.items()]
    glist.sort(key=itemgetter(0))
    for G in glist:
        pos,freq = G
        print "%s\t%d\t%d\t%d" % (d, pos, pos+1, freq)
