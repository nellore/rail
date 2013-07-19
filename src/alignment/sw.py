"""
sw.py Functions for Smith-Waterman which calculates local alignment value between two strings given a cost function
"""

import sys
import numpy
import re

numpy.set_printoptions(threshold=numpy.nan)

def cost(xc,yc):
    if xc=='-' or yc=='-': #gap
        return -2
    if xc==yc: #match
        return 1
    return -1  #mismatch

def smithWaterman(x,y,s):
    """ Calculate local alignment value of sequences x and y using
        dynamic programming.  Return local alignment value. """
    D = numpy.zeros((len(x)+1, len(y)+1), dtype=int)
    for i in xrange(1, len(x)+1):
        for j in xrange(1, len(y)+1):
            D[i, j] = max(D[i-1, j-1] + s(x[i-1], y[j-1]), # diagonal
                          D[i-1, j  ] + s(x[i-1], '-'),    # vertical
                          D[i  , j-1] + s('-',    y[j-1]), # horizontal
                          0)
    return D[len(x), len(y)]

cigar_pattern = re.compile(r"(\d+)(\S)")
def sumcigar(alignment):
    caligns = cigar_pattern.findall(alignment)
    total = 0
    for c in caligns:
        if c[1]=="I" or c[1]=="M":
            total+=int(c[0])
        elif c[0]=="D":
            total-=int(c[0])
    return total

"""
Returns a cigar like format of the alignment
"""
def cigar(alignment):
    count = 1
    letter = alignment[0]
    s = ""
    for i in range(1,len(alignment)):
        if letter!=alignment[i]:
            s+="%d%s"%(count,letter)
            letter = alignment[i]
            count = 1
        else:
            count+=1
    s+="%d%s"%(count,letter)
    letter = alignment[i]
    count = 1
    return s
        

def traceback(D, x, y,xi,yi,s):
    """ Trace back from bottom-right cell in edit-distance matrix D for
        strings x and y """
    #i, j = len(x), len(y)
    i,j = xi,yi
    xscript = []
    while i > 0 or j > 0:
        diag, vert, horz = -sys.maxint, -sys.maxint, -sys.maxint
        if i > 0 and j > 0:
            diag = D[i-1, j-1] + s(x[i-1], y[j-1])
        if i > 0:
            vert = D[i-1, j] + s(x[i-1], '-')
        if j > 0:
            horz = D[i, j-1] + s('-', y[j-1])

        if diag > vert and diag > horz:
            xscript.append('R' if x[i-1] != y[j-1] else 'M')
            i -= 1; j -= 1
        elif vert >= horz:
            xscript.append('I')
            i -= 1
        else:
            xscript.append('D')
            j -= 1
    return (''.join(xscript))[::-1]

def smithWatermanXcript(x, y, s):
    """ Calculate global alignment value of sequences x and y using
        dynamic programming.  Return global alignment value, optimal
        transcript. """
    D = numpy.zeros((len(x)+1, len(y)+1), dtype=int)
    for i in xrange(1, len(x)+1):
        for j in xrange(1, len(y)+1):
            D[i, j] = max(D[i-1, j-1] + s(x[i-1], y[j-1]), # diagonal
                          D[i-1, j  ] + s(x[i-1], '-'),    # vertical
                          D[i  , j-1] + s('-',    y[j-1]), # horizontal
                          0)
    
    #Find ending point
    maxI,indI = 0,len(x)
    for i in range(0,len(x)+1):
        maxI =  i if D[i,len(y)]>maxI else maxI

    return D[maxI, len(y)], cigar(traceback(D, x, y,maxI,len(y), s))
