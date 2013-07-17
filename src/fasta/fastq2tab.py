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
        seqid = entry[0].split(' ')[1]
        seq = entry[1]
        qual = entry[3]
        print "%s\t%s\t%s"%(seqid,seq,qual) 
    
    entry[index%4] = ln.rstrip()
    index+=1
seqid = entry[0].split(' ')[1]
seq = entry[1]
qual = entry[3]
print "%s\t%s\t%s"%(seqid,seq,qual) 


