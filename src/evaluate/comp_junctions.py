#! /usr/bin/env python
import sys
import csv
import math
import re
import string
import numpy as np
import matplotlib
import pylab

# Description goes here

# Read the alignment data from a sam file
# Update the coverage arrays
def read_data(filename, form, junctions, wgt):
    print 'Reading ' + filename

    count= 0 

    with open(filename, 'r') as tsv:
        rows = []
        for line in tsv:
            row = line.strip().split('\t')

            # read SAM file
            if (form == 'sam' and len(row) >= 11):
                chr_name = row[2]
                if chr_name in chromosomes:
                    start = getIndex(int(row[3]), chr_name)

                    pattern = row[5]
                    junctionOffsets = []
                    currOffset = 0

                    origPattern = pattern

                    match = re.search("\D", pattern)
                    while match and match.start() < len(pattern)-1:
                        index = match.start()
                        currOffset += int(''.join(pattern[:index]))
                        pattern = pattern[index+1:]
                        junctionOffsets.append(currOffset)
                        match = re.search("\D", pattern)

                    #print origPattern + ' -> ' + str(junctionOffsets)
                    
                    for off in junctionOffsets:
                        i = off + start - 1
                        junctions[i] = 1
                #else:
                #    print 'Unknown chromosome ' + row[2]

            # read BED file
            elif (form == 'bed'):
                chr_name = row[0]
                if chr_name == 'dmel_mitochondrion_genome':
                    chr_name = 'M'
                if chr_name in chromosomes:
                    start = getIndex(int(row[1]), chr_name)
                    junctionLens = [int(s) for s in string.split(row[10], ',')]
                    junctionOffsets = [int(s) for s in string.split(row[11], ',')]

                    for i in xrange(len(junctionOffsets)):
                        off = junctionOffsets[i]
                        length = junctionLens[i]
                        if i > 0:
                            junctions[off+start] = 1
                        if i < (len(junctionOffsets)-1):
                            junctions[off+start+length] = 1
                #else:
                #    print 'Unknown chromosome ' + chr_name

            # read GTF file
            else:
                chr_name = row[0]
                if chr_name == 'dmel_mitochondrion_genome':
                    chr_name = 'M'

                if chr_name in chromosomes:
                    start = getIndex(int(row[3]) - 1, chr_name)
                    end = getIndex(int(row[4]), chr_name)

                    # don't count endpoints of transcript
                    if row[2] == 'transcript':
                        junctions[start] = 0
                        junctions[end] = 0

                    if row[2] == 'exon':
                        count += 1
                        if start in junctions:
                            junctions.pop(start, None)
                        else:
                            junctions[start] = 1
  
                        if end in junctions:
                            junctions.pop(end, None)
                        else:
                            junctions[end] = 1
                    #count += 1
                    #junctions[start] = 1
                    #junctions[end] = 1

    return junctions


chromosomes = ['2L', '2R', '3L', '3R', '4', 'M', 'X', '2LHet', '2RHet', '3LHet', '3RHet', 'XHet', 'YHet', 'U', 'Uextra']
chromosomes = ['chr' + el for el in chromosomes]
chr_lengths = [['2L', 23011544],
               ['2R', 21146708],
               ['3L', 24543557],
               ['3R', 27905053],
               ['4', 1351857],
               ['M', 19517],
               ['X', 22422827],
               ['2LHet', 368872],
               ['2RHet', 3288761],
               ['3LHet', 2555491],
               ['3RHet', 2517507],
               ['XHet', 204112],
               ['YHet', 347038],
               ['U', 10049037],
               ['Uextra', 29004656]]
chr_lengths = [['chr' + el[0], el[1]] for el in chr_lengths]


#chromosomes = ['2L']
#chr_lengths = [['2L', 23011544]]

# Return the position in the 1d arry of the nt at the given position in a chromosome
def getIndex(position, chromosome):
    index = 0
    for [a,b] in chr_lengths:
        if a == chromosome:
            return index + position
        else:
            index += b

def getChromosomePos(index):
    for [a,b] in chr_lengths:
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

error = 0
if len(sys.argv) == 4:
    error = int(sys.argv[3])

format1 = actual[-3:len(actual)]
format2 = predicted[-3:len(predicted)]

if (format1 != 'sam' and format1 != 'bed' and format1 != 'gtf') or (format2 != 'sam' and format2 != 'bed' and format2 != 'gtf'):
    print 'Only .sam, .bed, and .gtf files are supported'


#print 'Reading actual junctions'
j_actual = read_data(actual, format1, j_actual, 1)

#print 'Reading predicted junctions'
j_predicted = read_data(predicted, format2, j_predicted, 1)

print len(j_actual)
print len(j_predicted)
tp = 0
fp = 0
fn = 0
tn = 0

jarray_actual = np.array([])
jarray_predicted = np.array([])

f = open('junctions.cn', 'w')
f.write('SNP\tChromosome\tPhysicalPosition\tJunctions\n')

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
fn = len(j_actual) - tp
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

dist = np.linalg.norm(np.subtract(jarray_actual, jarray_predicted))

len_actual = np.linalg.norm(jarray_actual)
len_predicted = np.linalg.norm(jarray_predicted)
jarray_actual = jarray_actual / len_actual
jarray_predicted = jarray_predicted / len_predicted

correlation = np.correlate(jarray_actual, jarray_predicted)

#print 'Distance:\t' + str(dist)
print 'Correlation:\t' + str(correlation)


# plot all points
#print 'Plotting'
#matplotlib.pyplot.scatter(jarray_actual, jarray_predicted)
#matplotlib.pyplot.show()
print 'Done!'
