import string
import re
import numpy as np

'''
Read aligned reads from a file and compute the coverage vector for the chromosome
'''

# Genome chromosomes and lengths
# TODO: read this from a file
#chromosomes = ['2L', '2R', '3L', '3R', '4', 'M', 'X', '2LHet', '2RHet', '3LHet', '3RHet', 'XHet', 'YHet', 'U', 'Uextra']
#chr_lengths = [['2L', 23011544],
#               ['2R', 21146708],
#               ['3L', 24543557],
#               ['3R', 27905053],
#               ['4', 1351857],
#               ['M', 19517],
#               ['X', 22422827],
#               ['2LHet', 368872],
#               ['2RHet', 3288761],
#               ['3LHet', 2555491],
#               ['3RHet', 2517507],
#               ['XHet', 204112],
#               ['YHet', 347038],
#               ['U', 10049037],
#               ['Uextra', 29004656]]

chromosomes = ['2L']
chr_lengths = [['2L', 23011544]]


'''
Return the position in the 1d arry of the nt at the given position in a chromosome
'''
def getIndex(position, chromosome):
    index = 0
    for [a,b] in chr_lengths:
        if a == chromosome:
            return index + position
        else:
            index += b

'''
 Read the alignment data from a sam file and update the coverage arrays
'''
def readSAM(filename, cov):
    print 'Reading ' + filename

    hits = 0
    with open(filename, 'r') as tsv:
        rows = []
        for line in tsv:
            row = line.strip().split('\t')

            chr_name = row[2]
            if chr_name in chromosomes:
                start = int(row[3])

                pattern = row[5]
                sectionLens = []
                sectionOffsets = []
                currOffset = 0

                # Parse the cigar string
                # Split the string by non-numeric characters and
                #   add the offsets
                match = re.search("\D", pattern)
                while match:
                    index = match.start()
                    if pattern[index] == 'M':
                        sectionOffsets.append(currOffset)
                        sectionLens.append(int(pattern[:index]))
                    currOffset += int(''.join(pattern[:index]))
                    pattern = pattern[index+1:]
                    match = re.search("\D", pattern)

                # If the read maps to multiple places, weight each position to
                #   1 / # mappings
                wgt = 1
                for i in row[11:len(row)]:
                    if i[0:5] == 'NH:i:':
                        wgt = 1 / int(i[5:len(i)])
                
                for i in xrange(len(sectionLens)):
                    for j in xrange(sectionLens[i]):
                        index = getIndex(start+sectionOffsets[i]+j, chr_name)
                        cov[index] += wgt
                    hits += sectionLens[i];

    return cov,hits

'''
 Read the alignment data from a sam file and update the coverage arrays
'''
def readBED(filename, cov):
    print 'Reading ' + filename

    hits = 0
    with open(filename, 'r') as tsv:
        rows = []
        for line in tsv:
            row = line.strip().split('\t')
                
            chr_name = row[0]
            if chr_name == 'dmel_mitochondrion_genome':
                chr_name = 'M'
            if chr_name in chromosomes:
                start = int(row[1])
                    
                # Each read is represented as a list of substrand lengths and
                #   a corresponding list of their offsets from the read start
                sectionLens = [int(s) for s in string.split(row[10], ',')]
                sectionOffsets = [int(s) for s in string.split(row[11], ',')]

                for i in xrange(len(sectionLens)):
                    for j in xrange(sectionLens[i]):
                        index = getIndex(start+sectionOffsets[i]+j, chr_name)
                        cov[index] += 1
                    hits += sectionLens[i];

    return cov,hits;

