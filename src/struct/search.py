"""
search.py 
Various ways to conduct a binary search on a sorted list

Copied from http://docs.python.org/2/library/bisect.html
"""
import bisect

def find(a,x):
    'Find the closest value to x'
    i = bisect.bisect_left(a,x)
    if i > 0 and i<len(a):
        return a[i-1] if abs(a[i-1]-x)< abs(a[i]-x) else a[i]
    else:
        if i==0:
            return a[i]
        elif i==len(a):
            return a[i-1]
        else:
            print "Index not found",i
            return -100000000
            
def index(a, x):
    'Locate the leftmost value exactly equal to x'
    i = bisect.bisect_left(a, x)
    if i != len(a) and a[i] == x:
        return i
    raise ValueError

def find_lt(a, x):
    'Find rightmost value less than x'
    i = bisect.bisect_left(a, x)
    if i:
        return a[i-1]
    raise ValueError

def find_ge(a, x):
    'Find leftmost item greater than or equal to x'
    i = bisect.bisect_left(a, x)
    if i != len(a):
        return a[i]
    raise ValueError

def find_le(a, x):
    'Find rightmost value less than or equal to x'
    i = bisect.bisect_right(a, x)
    if i:
        return a[i-1]
    raise ValueError

def find_gt(a, x):
    'Find leftmost value greater than x'
    i = bisect.bisect_right(a, x)
    if i != len(a):
        return a[i]
    raise ValueError

