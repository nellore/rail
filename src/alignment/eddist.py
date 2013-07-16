#!/usr/bin/env python

"""
eddist.py: Functions for calculating edit distance between two strings.

Author: Ben Langmead
Date: 2/14/2013
Contact: langmea@cs.jhu.edu
"""

import sys
import numpy

def edDistRecursive(x, y):
    if len(x) == 0: return len(y)
    if len(y) == 0: return len(x)
    delt = 1 if x[-1] != y[-1] else 0
    diag = edDistRecursive(x[:-1], y[:-1]) + delt 
    vert = edDistRecursive(x[:-1], y) + 1
    horz = edDistRecursive(x, y[:-1]) + 1
    return min(diag, vert, horz)

def edDistRecursiveMemo(x, y, memo=None):
    if memo is None:
        memo = {}
    if len(x) == 0:
        return len(y)
    elif len(y) == 0:
        return len(x)
    if (len(x), len(y)) in memo:
        return memo[(len(x), len(y))]
    mm = 0
    if x[-1] != y[-1]:
        mm = 1
    diag = edDistRecursiveMemo(x[:-1], y[:-1], memo) + mm 
    vert = edDistRecursiveMemo(x[:-1], y, memo) + 1
    horz = edDistRecursiveMemo(x, y[:-1], memo) + 1
    ans = min(diag, vert, horz)
    memo[(len(x), len(y))] = ans
    return ans

import numpy
def edDistDp(x, y):
    """ Calculate edit distance between sequences x and y using dynamic
        programming.  Return distance. """
    D = numpy.zeros((len(x)+1, len(y)+1), dtype=int)
    D[0, 1:] = range(1, len(y)+1)
    D[1:, 0] = range(1, len(x)+1)
    for i in xrange(1, len(x)+1):
        for j in xrange(1, len(y)+1):
            delt = 1 if x[i-1] != y[j-1] else 0
            D[i, j] = min(D[i-1, j-1]+delt, D[i-1, j]+1, D[i, j-1]+1)
    return D[len(x), len(y)]

def traceback(D, x, y):
    """ Trace back from bottom-right cell in edit-distance matrix D for
        strings x and y """
    i, j = len(x), len(y)
    xscript = []
    while i > 0 or j > 0:
        diag, vert, horz = sys.maxint, sys.maxint, sys.maxint
        delt = None
        if i > 0 and j > 0:
            delt = 0 if x[i-1] == y[j-1] else 1
            diag = D[i-1, j-1] + delt
        if i > 0:
            vert = D[i-1, j] + 1
        if j > 0:
            horz = D[i, j-1] + 1
        if diag <= vert and diag <= horz:
            xscript.append('R' if delt == 1 else 'M')
            i -= 1; j -= 1
        elif vert <= horz:
            xscript.append('I')
            i -= 1
        else:
            xscript.append('D')
            j -= 1
    return (''.join(xscript))[::-1]

import numpy
def edDistDpXcript(x, y):
    """ Calculate edit distance between sequences x and y using dynamic
        programming.  Return distance and optimal edit transcript. """
    D = numpy.zeros((len(x)+1, len(y)+1), dtype=int)
    D[0, 1:] = range(1, len(y)+1)
    D[1:, 0] = range(1, len(x)+1)
    for i in xrange(1, len(x)+1):
        for j in xrange(1, len(y)+1):
            delt = 1 if x[i-1] != y[j-1] else 0
            D[i, j] = min(D[i-1, j-1] + delt, D[i-1, j] + 1, D[i, j-1] + 1)
    return (D[len(x), len(y)], traceback(D, x, y))

def edDistMat(x, y):
    """ Calculate edit distance between sequences x and y using dynamic
        programming.  Return distance and optimal edit transcript. """
    D = numpy.zeros((len(x)+1, len(y)+1), dtype=int)
    D[0, 1:] = range(1, len(y)+1)
    D[1:, 0] = range(1, len(x)+1)
    for i in xrange(1, len(x)+1):
        for j in xrange(1, len(y)+1):
            delt = 1 if x[i-1] != y[j-1] else 0
            D[i, j] = min(D[i-1, j-1] + delt, D[i-1, j] + 1, D[i, j-1] + 1)
    return D

def edDistMatFill1(x, y):
    """ Do columns in outermost loop. """
    D = numpy.zeros((len(x)+1, len(y)+1), dtype=int)
    D[0, 1:] = range(1, len(y)+1)
    D[1:, 0] = range(1, len(x)+1)
    for j in xrange(1, len(y)+1):
        for i in xrange(1, len(x)+1):
            delt = 1 if x[i-1] != y[j-1] else 0
            D[i, j] = min(D[i-1, j-1] + delt, D[i-1, j] + 1, D[i, j-1] + 1)
    assert numpy.all(D == edDistMat(x, y))
    return D

def edDistMatFill2(x, y):
    """ Do anti-diagonals in outermost loop. """
    D = numpy.zeros((len(x)+1, len(y)+1), dtype=int)
    D[0, 1:] = range(1, len(y)+1)
    D[1:, 0] = range(1, len(x)+1)
    mn = min(len(x), len(y))
    mx = max(len(x), len(y))
    for k in xrange(0, len(x) + len(y) - 1):
        if k < mn:
            i, j = k+1, 1
            for m in xrange(0, k+1):
                delt = 1 if x[i-1] != y[j-1] else 0
                D[i, j] = min(D[i-1, j-1] + delt, D[i-1, j] + 1, D[i, j-1] + 1)
                i -= 1
                j += 1
        elif k < mx and len(x) < len(y):
            i, j = len(x), k-mn-2
            for m in xrange(0, mn):
                delt = 1 if x[i-1] != y[j-1] else 0
                D[i, j] = min(D[i-1, j-1] + delt, D[i-1, j] + 1, D[i, j-1] + 1)
                i -= 1
                j += 1
        elif k < mx:
            i, j = k+1, 1
            for m in xrange(0, mn):
                delt = 1 if x[i-1] != y[j-1] else 0
                D[i, j] = min(D[i-1, j-1] + delt, D[i-1, j] + 1, D[i, j-1] + 1)
                i -= 1
                j += 1
        else:
            for m in xrange(0, mn - (k-mx) - 1):
                pass
    print D
    assert numpy.all(D == edDistMat(x, y))
    return D

def hammingDist(x, y):
    """ Return Hamming distance between 2 same-length strings """
    assert len(x) == len(y)
    nmm = 0
    for i in xrange(0, len(x)):
        if x[i] != y[i]:
            nmm += 1
    return nmm

def edDistBounds(x, y):
    """ Return lower and upper bounds on the edit distance between two
        strings without doing complicated stuff (nothing O(mn)). """
    if x == y: return 0, 0
    absDiff = abs(len(x) - len(y))
    # if lengths differ, need at least some deletions/insertions are
    # required to make lengths the same
    lower = max(1, absDiff)
    upper = max(len(x), len(y))
    # an upper bound is the hamming distance between the shorter string
    # and a same-length prefix of the longer string, plus the number of
    # deletions/insertions needed to make them the same length
    if len(x) > len(y):
        upper = hammingDist(y, x[:len(y)]) + absDiff
    else:
        upper = hammingDist(x, y[:len(x)]) + absDiff
    return lower, upper
