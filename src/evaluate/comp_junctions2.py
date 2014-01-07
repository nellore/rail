#! /usr/bin/env python
import sys
import csv
import math
import re
import string
import numpy as np
import read_junctions 

'''
Compares the actual and predicted junction locations from 2 files
'''

'''
Returns the chromosome and index within the chromosome of the nt 
at the given position in the 1d array
'''
def getChromosomePos(index):
    for [a,b] in read_junctions.chr_lengths:
        if index < b:
            return index, a
        else:
            index -= b
    print 'Error! Index too high!'


# Initialize empty coverage arrays
j_actual = dict()
j_predicted = dict()

if len(sys.argv) != 3 and len(sys.argv) != 4:
    print "Usage: ./comp_coverage.py actual predicted error"
    print "  .sam .bed and .gtf files are supported"

actual = sys.argv[1]
predicted = sys.argv[2]

# Optional argument error
# Allows for some error of up to 'error' bases in junction placement
error = 0
if len(sys.argv) == 4:
    error = int(sys.argv[3])

format1 = actual[-3:len(actual)]
format2 = predicted[-3:len(predicted)]

#print 'Reading actual junctions'
if format1 == 'sam':
    j_actual = read_junctions.readJunctionsSAM(actual, j_actual, 1)
elif format1 == 'bed':
    j_actual = read_junctions.readJunctionsBED(actual, j_actual, 1)
elif format1 == 'gtf':
    j_actual = read_junctions.readJunctionsGTF(actual, j_actual, 1)
else:
    print 'Only .sam, .bed, and .gtf files are supported'

#print 'Reading predicted junctions'
if format2 == 'sam':
    j_predicted = read_junctions.readJunctionsSAM(predicted, j_predicted, 1)
elif format2 == 'bed':
    j_predicted = read_junctions.readJunctionsBED(predicted, j_predicted, 1)
elif format2 == 'gtf':
    j_predicted = read_junctions.readJunctionsGTF(predicted, j_predicted, 1)
else:
    print 'Only .sam, .bed, and .gtf files are supported'

print len(j_actual)
print len(j_predicted)
tp = 0
fp = 0
fn = 0
tn = 0

jarray_actual = np.array([])
jarray_predicted = np.array([])

# Write heatmap for IGV showing tp, fp, tn, & fn bases
#   tp = 0 (white)
#   fp = 1 (red)
#   fn = -1 (blue)
#   tn = nothing (gray)
f = open('junctions.cn', 'w')
f.write('SNP\tChromosome\tPhysicalPosition\tJunctions\n')

# Count tp, fp, and fn
# First loop through predicted junctions and find 
#   matches (tp) and mismatches (fp) in the actual junctions list
for i in j_predicted:
    pos, chrom = getChromosomePos(i)
    jarray_predicted = np.append(jarray_predicted, j_predicted[i])

    found = False

    for n in xrange(-error,error+1):
        if (i+n) in j_actual:
            found = True
            tp += 1
            jarray_actual = np.append(jarray_actual, j_actual[i+n])
            f.write('SNP\t' + chrom + '\t' + str(pos) + '\t0\n')
            break

    if not found:
        fp += 1
        jarray_actual = np.append(jarray_actual, 0)
        f.write('SNP\t' + chrom + '\t' + str(pos) + '\t1\n')

# Then loop through the actual junctions list and find 
#   mismatches (fn) in the predicted list
for i in j_actual:
    found = False

    for n in xrange(-error, error+1):
        if (i+n) in j_predicted:
            found = True
            break

    if not found:
        pos, chrom = getChromosomePos(i)
        fn += 1
        jarray_actual = np.append(jarray_actual, j_actual[i])
        jarray_predicted = np.append(jarray_predicted, 0)
        
        pos, chrom = getChromosomePos(i)
        f.write('SNP\t' + chrom + '\t' + str(pos) + '\t-1\n')

f.close()        

#for (name, length) in chr_lengths:
#    tn += length
#tn = tn - tp - fp - fn

total = len(jarray_actual)

print 'True positives:\t\t' + str(tp)
print 'False positives:\t' + str(fp) 
print 'False negatives:\t' + str(fn)
#print 'True negatives:\t\t' + str(tn)
print '\n'

#print 'Sensitivity:\t\t' + str(float(tp) / (tp+fn))
#print 'Specificity:\t\t' + str(float(tn) / (tn+fp))
#print 'Accuracy:\t\t' + str(float(tp+tn) / (tp+fp+fn+tn))
print 'Precision:\t\t' + str(float(tp) / (tp+fp)) 
print 'Recall:\t\t\t' + str(float(tp) / (tp+fn))
print '\n'

#dist = np.linalg.norm(np.subtract(jarray_actual, jarray_predicted))

# Normalize junction vectors by size
len_actual = np.linalg.norm(jarray_actual)
len_predicted = np.linalg.norm(jarray_predicted)
jarray_actual = jarray_actual / len_actual
jarray_predicted = jarray_predicted / len_predicted

correlation = np.correlate(jarray_actual, jarray_predicted)
print 'Correlation:\t' + str(correlation)

print 'Done!'
