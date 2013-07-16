"""
sw.py Functions for Smith-Waterman which calculates local alignment value between two strings given a cost function
"""

import sys
import numpy


def exampleCost(xc, yc):
    """ Cost function assigning 0 to match, 2 to transition, 4 to
        transversion, and 8 to gaps """
    if xc == yc: return 0 # match
    if xc == '-' or yc == '-': return 8 # gap
    minc, maxc = min(xc, yc), max(xc, yc)
    if minc == 'A' and maxc == 'G': return 2 # transition
    elif minc == 'C' and maxc == 'T': return 2 # transition
    return 4 # transversion

def lcsCost(xc, yc):
    return -1 if xc == yc else 0


def smithWaterman(x,y,s):
    """ Calculate local alignment value of sequences x and y using
        dynamic programming.  Return local alignment value. """
    D = numpy.zeros((len(x)+1, len(y)+1), dtype=int)
    for i in xrange(1, len(x)+1):
        for j in xrange(1, len(y)+1):
            D[i, j] = min(D[i-1, j-1] + s(x[i-1], y[j-1]), # diagonal
                          D[i-1, j  ] + s(x[i-1], '-'),    # vertical
                          D[i  , j-1] + s('-',    y[j-1])) # horizontal
    return D[len(x), len(y)]

