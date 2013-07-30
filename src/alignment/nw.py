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
import re
import scipy.weave

def exampleCost(xc, yc):
    if type(xc)==type(1) and type(yc)==type(1):
        xc,yc = chr(xc),chr(yc)
    """ Cost function assigning 0 to match, 2 to transition, 4 to
        transversion, and 8 to gaps """
    if xc == yc: return 0 # match
    if xc == '-' or yc == '-': return 3 # gap
    minc, maxc = min(xc, yc), max(xc, yc)
    if minc == 'A' and maxc == 'G': return 1 # transition
    elif minc == 'C' and maxc == 'T': return 1 # transition
    return 2 # transversion

def lcsCost(xc, yc):
    if type(xc)==type(1) and type(yc)==type(1):
        xc,yc = chr(xc),chr(yc)
    return -1 if xc == yc else 0

def matchCost(xc,yc):
    #print >> sys.stderr,xc,yc,chr(xc),chr(yc)
    return 1 if xc==yc else -4

def inverseMatchCost(xc,yc):
    return 0 if xc==yc else 1

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


"""
An inlined C implementation of Needleman Wunsch.

Note: s MUST be a functor, otherwise very bad things will happen
Also notice that this function takes the maxinum whereas its counterpart takes the mininum
"""
def c_needlemanWunsch(px, py, s):
    assert type(px) == type("")
    assert type(py) == type("")
    D = numpy.zeros((len(px)+1, len(py)+1), dtype=int)
    code = """
    py::tuple arg(2);
    std::string x = std::string(px);
    std::string y = std::string(py);
    for(int j = 0 ; j<y.length();j++){
       arg[0] = '-'; arg[1] = y[j-1];
       D(0,j) = j * (int)s.call(arg);    }
    for(int i = 0 ; i<x.length();i++){
       arg[0] = x[i-1]; arg[1] = '-';
       D(i,0) = i * (int)s.call(arg);    }
    for(int i = 1; i<x.length()+1; i++){
       for(int j = 1; j<y.length()+1; j++){
           arg[0] = x[i-1]; arg[1] = y[j-1];
           int cc = (int)s.call(arg);
           arg[0] = x[i-1]; arg[1] = '-';
           int c_ = (int)s.call(arg);
           arg[0] = '-'; arg[1] = y[j-1];
           int _c = (int)s.call(arg);
           D(i,j) = std::max( D(i-1,j-1)+cc, 
                              std::max( D(i-1,j) + c_,
                                        D(i,j-1) + _c)); 
       }
    }
    """
    
    scipy.weave.inline(code,['D','px','py','s'],
                       type_converters=scipy.weave.converters.blitz,
                       compiler='gcc')
    return D[len(px), len(py)], D


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
                          D[i  , j-1] + s('-',    y[j-1])) # horizonta
    return D[len(x), len(y)], cigar(traceback(D, x, y, s))



"""
An inlined C implementation of Needleman Wunsch.

Note: s MUST be a functor, otherwise very bad things will happen
"""
def c_needlemanWunschXcript(px, py, s):
    assert type(px) == type("")
    assert type(py) == type("")
    D = numpy.zeros((len(px)+1, len(py)+1), dtype=int)
    code = """
    py::tuple arg(2);
    std::string x = std::string(px);
    std::string y = std::string(py);
    for(int j = 0 ; j<y.length();j++){
       arg[0] = '-'; arg[1] = y[j-1];
       D(0,j) = j * (int)s.call(arg);    }
    for(int i = 0 ; i<x.length();i++){
       arg[0] = x[i-1]; arg[1] = '-';
       D(i,0) = i * (int)s.call(arg);    }
    for(int i = 1; i<x.length()+1; i++){
       for(int j = 1; j<y.length()+1; j++){
           arg[0] = x[i-1]; arg[1] = y[j-1];
           int cc = (int)s.call(arg);
           arg[0] = x[i-1]; arg[1] = '-';
           int c_ = (int)s.call(arg);
           arg[0] = '-'; arg[1] = y[j-1];
           int _c = (int)s.call(arg);
           D(i,j) = std::min( D(i-1,j-1)+cc, 
                              std::min( D(i-1,j) + c_,
                                        D(i,j-1) + _c)); 
       }
    }
    """
    
    scipy.weave.inline(code,['D','px','py','s'],
                       type_converters=scipy.weave.converters.blitz,
                       compiler='gcc')
    return D[len(px), len(py)], cigar(traceback(D, px, py, s))
