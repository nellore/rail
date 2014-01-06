#! /usr/bin/env python
import sys
import csv
import math
import re
import string
import numpy as np
import matplotlib
import pylab
import read_alignments

'''
Compares the coverage vectors from 2 alignments of a set of reads
'''

# Construct coverage arrays
print 'Initializing coverage arrays'
total_length = 0
for [a, b] in read_alignments.chr_lengths:
    total_length += b

coverage_actual = np.zeros(total_length)
coverage_predicted = np.zeros(total_length)

if len(sys.argv) != 3:
    print "Usage: ./comp_coverage.py actual predicted"
    print "  .sam and .bed files are supported"

actual = sys.argv[1]
predicted = sys.argv[2]

format1 = actual[-3:len(actual)]
format2 = predicted[-3:len(predicted)]

# Read alignments from files
print 'Reading actual alignments'
if format1 == 'sam':
    coverage_actual, hits_actual = read_alignments.readSAM(actual, coverage_actual)
elif format1 == 'bed':
    coverage_actual, hits_actual = read_alignments.readBED(actual, coverage_actual)
else:
    print 'Only .sam and .bed files are supported'
print str(hits_actual) + ' actual hits'

print 'Reading predicted alignments'
if format2 == 'sam':
    coverage_predicted, hits_predicted = read_alignments.readSAM(predicted, coverage_predicted)
elif format2 == 'bed':
    coverage_predicted, hits_predicted = read_alignments.readBED(predicted, coverage_predicted)
else:
    print 'Only .sam and .bed files are supported'
print str(hits_predicted) + ' predicted hits'


# normalize coverage vectors
#len_actual = np.linalg.norm(coverage_actual)
#len_predicted = np.linalg.norm(coverage_predicted)
#coverage_actual = coverage_actual / len_actual
#coverage_predicted = coverage_predicted / len_predicted

# Compute Euclidean distance between the coverage vectors
#dist = np.linalg.norm(np.subtract(coverage_actual, coverage_predicted))
#print 'Distance:\t' + str(dist)

# Compute the correlation between the coverage vectors
correlation = np.correlate(coverage_actual, coverage_predicted)
print 'Correlation:\t' + str(correlation)


# plot every 10th point, ignoring all (0,0) points
print 'Plotting'
xs = coverage_actual[0:len(coverage_actual):10]
ys = coverage_predicted[0:len(coverage_predicted):10]
matplotlib.pyplot.scatter(xs+ys, xs-ys, alpha=0.5)
matplotlib.pyplot.show()
print 'Done!'
