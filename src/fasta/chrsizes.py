


"""
Returns a dictionary with all of the chromosome sizes
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

def totalLength(sizes):
    total = 0
    for k,v in sizes.iteritems():
        total+=v
    return total
