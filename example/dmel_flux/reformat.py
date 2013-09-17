import sys
import argparse
import os
import site
import fileinput

parser = argparse.ArgumentParser(description=\
                                     'Appends the file name to every read id in the fastq file.')
parser.add_argument(\
    '--input',type=str,default='',help='Path of input file to be modified in place')
parser.add_argument(\
    '--basename',type=str,default='',help='Path of input file to be modified in place')
args = parser.parse_args()

if __name__=="__main__":
    lineno = 0
    for line in fileinput.input(args.input, inplace=True):
        ln = line.rstrip()

        if lineno%4==0:
            print "@%s.%s" % (args.basename, ln[1:])
        else:
            print ln
        lineno+=1
