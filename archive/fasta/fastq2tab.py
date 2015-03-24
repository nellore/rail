'''
fastq2tab.py
(before align.py)

Converts all of the fastq reads into a tab delimited format

Input: FASTQ format
Output format:

 Format 1 (unpaired):                                                           
  1. Name                                                                       
  2. Nucleotide sequence                                                        
  3. Quality sequence
'''

import sys

index = 0
entry = ['']*4
for ln in sys.stdin:
    if index%4==0 and index!=0:
        line = entry[0].split(' ')[0]
        seqid,rdnm = line.split(".")
        seqid = seqid[1:]
        seq = entry[1]
        qual = entry[3]
        print "r_n%s;LB:%s\t%s\t%s"%(rdnm,seqid,seq,qual) 
    
    entry[index%4] = ln.rstrip()
    index+=1

line = entry[0].split(' ')[0]
seqid,rdnm = line.split(".")
seqid = seqid[1:]
seq = entry[1]
qual = entry[3]
print "r_n%s;LB:%s\t%s\t%s"%(rdnm,seqid,seq,qual) 


