import os
import sys
import argparse
import site
import time

#TODO: Make bedgraph output cleaner
for ln in sys.stdin:
    ln = ln.rstrip()
    toks = ln.split('\t')
    if len(toks)<11:
        continue
    seqid,junc_st,junc_end,depth,blocksz = toks[0],int(toks[1]),int(toks[2]),int(toks[4]),toks[10]
    sizes = blocksz.split(',')
    if len(sizes)>1:
        left_size,right_size = int(sizes[0]), int(sizes[1])
        left, right = junc_st+left_size, junc_end-right_size
        print "%s\t%d\t%d\t%d"%(seqid,left,left+1,depth)
        print "%s\t%d\t%d\t%d"%(seqid,left+1,left+2,depth)
        print "%s\t%d\t%d\t%d"%(seqid,right-2,right-1,depth)
        print "%s\t%d\t%d\t%d"%(seqid,right-1,right,depth)
