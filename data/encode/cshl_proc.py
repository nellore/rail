"""
cshl_proc.py:

Take the MANIFEST.cshl file and preprocess it into a Myrna-style manifest file.
"""

import sys
import string

names = set()
paired = True
curName = None
rep = 1

for ln in sys.stdin:
    ln = ln.rstrip()
    if len(ln) == 0:
        continue
    if ln[0] == '#':
        toks = ln.split()
        if len(toks) == 1:
            continue # empty comment line
        curName = toks[-1]
        assert curName not in names, curName
        names.add(curName)
        curName = curName.replace('-', '_')
        rep = 1
        print >> sys.stderr, "Name: %s" % curName
    else:
        toks = string.split(ln, '\t')
        assert len(toks) == 2
        url = toks[0]
        urlToks = string.split(url, '/')
        srr = urlToks[-2]
        print >> sys.stderr, "  SRR: %s" % srr
        new_urls = [ 'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/%s/%s/%s_%d.fastq.gz' % (srr[:6], srr, srr, i) for i in [1, 2] ]
        print '\t'.join([new_urls[0], '0', new_urls[1], '0', '%s-1-%d' % (curName, rep)])
        rep += 1

print >> sys.stderr, "%d names" % len(names)

# # GSM758565: CshlLong_RnaSeq_HUVEC_nucleus_longPolyA
# ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR307/SRR307909/SRR307909.sra    SRR307909.sra
# ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR307/SRR307910/SRR307910.sra    SRR307910.sra
