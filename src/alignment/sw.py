"""
sw.py Functions for Smith-Waterman which calculates local alignment value between two strings given a cost function
"""

import sys
import numpy



def cost(xc,yc):
    if xc=='-' or yc=='-': #gap
        return -2 
    if xc==yc: #match
        return 2
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

def traceback(D, x, y, s):
    """ Trace back from bottom-right cell in edit-distance matrix D for
        strings x and y """
    i, j = len(x), len(y)
    xscript = []
    while i > 0 or j > 0:
        diag, vert, horz = -sys.maxint, -sys.maxint, -sys.maxint
        if i > 0 and j > 0:
            diag = D[i-1, j-1] + s(x[i-1], y[j-1])
        if i > 0:
            vert = D[i-1, j] + s(x[i-1], '-')
        if j > 0:
            horz = D[i, j-1] + s('-', y[j-1])
        if diag >= vert and diag >= horz:
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
    return D[len(x), len(y)], traceback(D, x, y, s)
