import math
import sys

#Mean of histogram.  Array MUST be normalized (e.g. between [0,1] )
def average(p):
    result = 0
    for x in range(0,len(p)):
        result+=x*p[x]
    return result
#Second moment of histogram. Array MUST be normalized (e.g. between [0,1] )
def moment2(p):
    result = 0
    for x in range(0,len(p)):
        result+=x*x*p[x]
    return result
#Standard deviation of histogram. Array MUST be normalized (e.g. between [0,1] )
def stddev(p):
    m = average(p)
    m2 = moment2(p)
    return math.sqrt(m2-(m*m))

"""
Scores based off of a histogram of boundary locations
"""
def hist_score(coords,offset,endtype,N):
    #offset is placed in the middle of the histogram
    hist = [1.0]*N #pseudo counts
    n = N/2
    for c in coords:
        if abs(offset-c)>N:
            print >> sys.stderr,"Out of bounds coordinate"
            continue
        ind = (c-offset)+n if endtype=="5" else (N-(offset-c)-1)-n
        hist[ind]+=1
    total = sum(hist)
    hist = [h/float(total) for h in hist]
    return hist

"""
Assigns scores based off of how close it is to the end based off of a p=2 series
Sites closer to the ends get higher scores
"""
def exp_score(endtype,N):
    if endtype=="5":
        hist = [.5**n for n in range(0,N)]
    else:
        hist = [.5**(N-n+1) for n in range(0,N)]
    return hist

"""
Assigns scores based off of a harmonic series (aka. logarithmic score)
"""
def harmonic_score(endtype,N):
    if endtype=="5":
        hist = [1.0/(n+1) for n in range(0,N)]
    else:
        hist = [1.0/(N-n) for n in range(0,N)]
    return hist
"""
Assigns scores based off of a normal distribution with mean=center of window and std=window_size/2
"""
def normal(x,m,s):
    return (1.0/(math.sqrt(2*math.pi*s*s)))*math.exp(-(x-m)*(x-m)/(2*s*s))

"""
Score using a normal distribution s.t. positions in the center will be weighted higher
"""
def normal_score(N,m,s):
    #m = N/2
    #s = m
    hist = [normal(i,m,s) for i in range(0,N)]
    return hist
