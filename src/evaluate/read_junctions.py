import string
import re

'''
Read junctions from a file in sam, bed, or gtf format
'''

# Genome chromosomes and lengths
# TODO: read this from a file
chromosomes = ['2L', '2R', '3L', '3R', '4', 'M', 'X', '2LHet', '2RHet', '3LHet', '3RHet', 'XHet', 'YHet', 'U', 'Uextra']
#chromosomes = ['chr' + el for el in chromosomes]
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
#chr_lengths = [['chr'+el[0], el[1]] for el in chr_lengths]

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
Read junctions from a SAM alignment file
Junctions are inferred from all spliced reads
'''
def readJunctionsSAM(filename, junctions, wgt):
    print 'Reading ' + filename

    with open(filename, 'r') as tsv:
        rows = []
        for line in tsv:
            row = line.strip().split('\t')

            chr_name = row[2]
            if chr_name in chromosomes:
                start = getIndex(int(row[3]), chr_name)

                pattern = row[5]
                junctionOffsets = []
                currOffset = 0

                # Parse the cigar string
                # Split the string by non-numeric characters and
                #   add the offsets
                origPattern = pattern
                match = re.search("\D", pattern)
                while match and match.start() < len(pattern)-1:
                    index = match.start()
                    currOffset += int(''.join(pattern[:index]))
                    pattern = pattern[index+1:]
                    junctionOffsets.append(currOffset)
                    match = re.search("\D", pattern)
                   
                for off in junctionOffsets:
                    i = off + start - 1
                    junctions[i] = 1

    return junctions

'''
Read junctions from a bed alignment file
Junctions are inferred from all spliced reads
'''
def readJunctionsBED(filename, junctions, wgt):
    print 'Reading ' + filename

    with open(filename, 'r') as tsv:
        rows = []
        for line in tsv:
            row = line.strip().split('\t')

            chr_name = row[0]
            if chr_name == 'dmel_mitochondrion_genome':
                chr_name = 'M'
            if chr_name in chromosomes:
                start = getIndex(int(row[1]), chr_name)
                
                # Each read is represented as a list of substrand lengths and
                #   a corresponding list of their offsets from the read start
                junctionLens = [int(s) for s in string.split(row[10], ',')]
                junctionOffsets = [int(s) for s in string.split(row[11], ',')]

                for i in xrange(len(junctionOffsets)):
                    off = junctionOffsets[i]
                    length = junctionLens[i]
                    if i > 0:
                        junctions[off+start] = 1
                    if i < (len(junctionOffsets)-1):
                        junctions[off+start+length] = 1

    return junctions

'''
Read all junctions from a gtf file
'''
def readJunctionsGTF(filename, junctions, wgt):
    print 'Reading ' + filename

    with open(filename, 'r') as tsv:
        rows = []
        for line in tsv:
            row = line.strip().split('\t')

            chr_name = row[0]
            if chr_name == 'dmel_mitochondrion_genome':
                chr_name = 'M'

            if chr_name in chromosomes:
                start = getIndex(int(row[3]) - 1, chr_name)
                end = getIndex(int(row[4]), chr_name)

                # Don't count endpoints of transcript
                # Each transcript is represented as a 'transcript' row
                #   and a list of 'exon' rows
                # Remove each of the transcript endpoints because they
                #   often can't be accurately matched
                if row[2] == 'transcript':
                    junctions[start] = 0
                    junctions[end] = 0

                if row[2] == 'exon':
                    if start in junctions:
                        junctions.pop(start, None)
                    else:
                        junctions[start] = 1
  
                    if end in junctions:
                        junctions.pop(end, None)
                    else:
                        junctions[end] = 1

                # To save all junctions including endpoints, comment out the
                #   section above and uncomment the lines below
#                junctions[start] = 1
#                junctions[end] = 1

    return junctions

