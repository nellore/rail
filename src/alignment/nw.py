#!/usr/bin/env python

"""
nw.py: Functions for Needleman-Wunsch, which calculates global
       alignment value between two strings given a cost function that
       assigns a cost to all possible substitutions and gaps.

Author: Ben Langmead
Date: 2/26/2013
Contact: langmea@cs.jhu.edu
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

def needlemanWunsch(x, y, s):
    """ Calculate global alignment value of sequences x and y using
    dynamic programming.  Return global alignment value. """
    D = numpy.zeros((len(x)+1, len(y)+1), dtype=int)
    for j in xrange(1, len(y)+1):
        D[0, j] = j * s('-', y[j-1])
    for i in xrange(1, len(x)+1):
        D[i, 0] = i * s(x[i-1], '-')
    for i in xrange(1, len(x)+1):
        for j in xrange(1, len(y)+1):
            D[i, j] = min(D[i-1, j-1] + s(x[i-1], y[j-1]), # diagonal
                          D[i-1, j  ] + s(x[i-1], '-'),    # vertical
                          D[i  , j-1] + s('-',    y[j-1])) # horizontal
    return D[len(x), len(y)]
             
def traceback(D, x, y, s):
    """ Trace back from bottom-right cell in edit-distance matrix D for
        strings x and y """
    i, j = len(x), len(y)
    xscript = []
    while i > 0 or j > 0:
        diag, vert, horz = sys.maxint, sys.maxint, sys.maxint
        if i > 0 and j > 0:
            diag = D[i-1, j-1] + s(x[i-1], y[j-1])
        if i > 0:
            vert = D[i-1, j] + s(x[i-1], '-')
        if j > 0:
            horz = D[i, j-1] + s('-', y[j-1])
        if diag <= vert and diag <= horz:
            xscript.append('R' if x[i-1] != y[j-1] else 'M')
            i -= 1; j -= 1
        elif vert <= horz:
            xscript.append('I')
            i -= 1
        else:
            xscript.append('D')
            j -= 1
    return (''.join(xscript))[::-1]

def needlemanWunschXcript(x, y, s):
    """ Calculate global alignment value of sequences x and y using
        dynamic programming.  Return global alignment value, optimal
        transcript. """
    D = numpy.zeros((len(x)+1, len(y)+1), dtype=int)
    for j in xrange(1, len(y)+1):
        D[0, j] = j * s('-', y[j-1])
    for i in xrange(1, len(x)+1):
        D[i, 0] = i * s(x[i-1], '-')
    for i in xrange(1, len(x)+1):
        for j in xrange(1, len(y)+1):
            D[i, j] = min(D[i-1, j-1] + s(x[i-1], y[j-1]), # diagonal
                          D[i-1, j  ] + s(x[i-1], '-'),    # vertical
                          D[i  , j-1] + s('-',    y[j-1])) # horizontal
    return D[len(x), len(y)], traceback(D, x, y, s)

