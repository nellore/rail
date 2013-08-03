import sys
import numpy
import nw
import re

#Some commonly used substitution matrices
#IMPORTANT: All substitution matrices must be of type numpy.int32

def lcsCost():
    """ A substitution matrix where match=-1, everything else =0.  Only makes
       sense for solving LCS problems, and only if DP algorithm is taking a
       min at each cell. """
    M = numpy.zeros((6,6),numpy.int32)
    for i in range(0,6):
        M[i,i] = -1
    return M

def matchCost():
    """ Return a substitution matrix where match=+1, everything else=-1.  The
        6 symbols are 0=A, 1=C, 2=G, 3=T, 4=N, 5=- (gap). """
    M = numpy.zeros((6,6),numpy.int32)
    M.fill(-1)
    for i in range(0,6):
        M[i,i] = 1
    return M

def inverseMatchCost():
    """ Return a substitution matrix where match=0, everything else=1.  Only
        makes sense if the DP algorithm is taking a min at each cell. """
    M = numpy.zeros((6,6),numpy.int32)
    M.fill(1)
    for i in range(0,6):
        M[i,i] = 0
    return M

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

def char2index(string):
    return [ "ACGTN-".index(c) for c in string ]

def traceback(D, sx, sy, s):
    """ Trace back from bottom-right cell in edit-distance matrix D for
    strings x and y """
    x,y = char2index(sx),char2index(sy)
    i, j = len(x), len(y)
    xscript = []
    while i > 0 or j > 0:
        diag, vert, horz = sys.maxint, sys.maxint, sys.maxint
        if i > 0 and j > 0:
            diag = D[i-1, j-1] + s[x[i-1], y[j-1]]
        if i > 0:
            vert = D[i-1, j] + s[x[i-1], '-']
        if j > 0:
            horz = D[i, j-1] + s['-', y[j-1]]
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

#x,y are the alignment sequences, s is the substitution matrix
def needlemanWunsch(x,y,s):
    L = len(x)
    D = numpy.zeros((L+1,L+1),numpy.int32)
    nw.nw(D,x,y,s)
    return D[L,L]

#x,y are the alignment sequences, s is the substitution matrix
def needlemanWunschXcript(x,y,s):
    L = len(x)
    D = numpy.zeros((L+1,L+1),numpy.int32)
    nw.nw(D,x,y,s)
    return D[L,L],D
