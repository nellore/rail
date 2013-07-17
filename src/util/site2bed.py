"""
Converts splice sites into bed format

Doesn't take into consideration sample label atm
"""
import os
import sys
import argparse
import site
import time

from collections import Counter
from collections import defaultdict
from operator import itemgetter

gpositions = defaultdict(Counter)

for ln in sys.stdin:
    ln = ln.rstrip()
    toks = ln.split("\t")
    refid,st,en,lab,freq = toks[0],int(toks[1]),int(toks[2]),toks[3],int(toks[4])
    gpositions[refid][str(st)]+=freq
    gpositions[refid][str(st+1)]+=freq
    gpositions[refid][str(en)]+=freq
    gpositions[refid][str(en+1)]+=freq

for d,c in gpositions.iteritems():
    glist = [(int(pos),freq) for pos,freq in c.items()]
    glist.sort(key=itemgetter(0))
    for G in glist:
        pos,freq = G
        print "%s\t%d\t%d\t%d"%(d,pos,pos+1,freq)
