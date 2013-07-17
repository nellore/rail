"""
chrsizes.py 

Utility functions for handling sequence lengths from fasta index file
"""


"""
Returns a dictionary with all of the chromosome sizes given a fai file
"""
def getSizes(fname):
    sizes = dict()
    with open(fname,'r') as fh:
        for ln in fh:
            line = ln.rstrip()
            toks = line.split('\t')
            seq,length  = toks[0],int(toks[1])
            sizes[seq] = length
    return sizes
"""
Calculates the total length of all of the sequences added together given a fai file
"""
def totalLength(fname):
    sizes = getSizes(fname)
    total = 0
    for k,v in sizes.iteritems():
        total+=v
    return total
