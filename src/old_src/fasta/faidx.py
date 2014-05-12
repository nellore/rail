"""
faidx.py

Index a fasta file in the same manner as 'samtools faidx'
"""

import sys

outFa, outFaIdx = sys.argv[1:3]
assert outFa is not None
assert outFaIdx is not None

outFaH = open(outFa, 'w')
outFaIdxH = open(outFaIdx, 'w')

idx = []
nm, num, byteoff, myByteoff, nchar, ninc = None, None, 0, 0, 0, 0

while True:
    ln = sys.stdin.readline()
    lnlen = len(ln)
    outFaH.write(ln)
    if len(ln) == 0: break
    if ln[0] == '>':
        if nm is not None:
            outFaIdxH.write('\t'.join(map(str, [nm, num, myByteoff, nchar, ninc])))
            outFaIdxH.write('\n')
            num = 0
        myByteoff = byteoff + lnlen
        num = 0
        if ' ' in ln:
            nm = ln[1:ln.index(' ')].rstrip()
        else:
            nm = ln[1:].rstrip()
    else:
        ninc = max(ninc, len(ln))
        ln = ln.rstrip()
        nchar = max(nchar, len(ln))
        num += len(ln)
    byteoff += lnlen

if nm is not None:
    outFaIdxH.write('\t'.join(map(str, [nm, num, myByteoff, nchar, ninc])))
    outFaIdxH.write('\n')

outFaH.close()
outFaIdxH.close()
