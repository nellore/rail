'''
statsmath.py

Some functions that mirror some needed functionality in R.
'''

import math
import scipy
from scipy import special
import sys

def digamma(x):
    return scipy.special.digamma()

def trigamma(x):
    return scipy.special.polygamma(1, x)

def tetragamma(x):
    return scipy.special.polygamma(2, x)

def trigammaInverse(x, iterlimit=50):
    ''' Inverse trigamma function, along the lines of the R
        implementation below, but assumes argument is scalar. '''
    if math.isnan(x):
        return x
    if x < 0:
        return float("nan")
    if x > 1e07:
        if x == 0.0:
            return float("inf")
        else:
            return 1.0 / math.sqrt(x)
    if x < 1e-06:
        if x == 0.0:
            return float("inf")
        else:
            return 1.0 / x
    assert x >= 1e-06
    y = 0.5 + 1.0 / x
    itr = 0
    while True:
        itr += 1
        tri = trigamma(y)
        dif = tri * (1 - tri/x) / tetragamma(y)
        y += dif
        if -dif / y < 1e-08:
            break
        if itr > iterlimit:
            print >>sys.stderr, "trigammaInverse iteration limit (%s) exceeded for argument %f" % (iterlimit, x)
    return y

#trigammaInverse <- function (x) 
#{
#    if (!is.numeric(x)) 
#        stop("Non-numeric argument to mathematical function")
#    if (length(x) == 0) 
#        return(numeric(0))
#    omit <- is.na(x)
#    if (any(omit)) {
#        y <- x
#        if (any(!omit)) 
#            y[!omit] <- Recall(x[!omit])
#        return(y)
#    }
#    omit <- (x < 0)
#    if (any(omit)) {
#        y <- x
#        y[omit] <- NaN
#        warning("NaNs produced")
#        if (any(!omit)) 
#            y[!omit] <- Recall(x[!omit])
#        return(y)
#    }
#    omit <- (x > 1e+07)
#    if (any(omit)) {
#        y <- x
#        y[omit] <- 1/sqrt(x[omit])
#        if (any(!omit)) 
#            y[!omit] <- Recall(x[!omit])
#        return(y)
#    }
#    omit <- (x < 1e-06)
#    if (any(omit)) {
#        y <- x
#        y[omit] <- 1/x[omit]
#        if (any(!omit)) 
#            y[!omit] <- Recall(x[!omit])
#        return(y)
#    }
#    y <- 0.5 + 1/x
#    iter <- 0
#    repeat {
#        iter <- iter + 1
#        tri <- trigamma(y)
#        dif <- tri * (1 - tri/x)/psigamma(y, deriv = 2)
#        y <- y + dif
#        if (max(-dif/y) < 1e-08) 
#            break
#        if (iter > 50) {
#            warning("Iteration limit exceeded")
#            break
#        }
#    }
#    y
#}
