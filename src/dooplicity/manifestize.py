#!/usr/bin/env python

"""
Author: Ben Langmead
  Date: 6/5/2017

Given a list of run accessions and corresponding run names on standard in, e.g.:

---
SRR32128090 SRP071321-PCells-16
SRR32128140 SRP071321-PCells-17
SRR32128150 SRP071321-PCells-18
---

...this script generates a manifest file with URLs for each, preferring ENA over
SRA where possible.  Prints manifest to standard out and a metadata table to
standard error.

TODO:
- Check for existence
- Get MD5s if available

For instructions on how to turn a run accession into a valid ENA URL, see:
http://www.ebi.ac.uk/ena/browse/read-download
"""

from __future__ import print_function
import subprocess
import os
import sys

query = []

for ln in sys.stdin:
    ln = ln.rstrip()
    toks = ln.split()
    srr = toks[0]
    query.append(srr)
    srr_pre = srr[:6]
    if len(srr) == 10:
        srr00 = '00' + srr[-1]
        url = 'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/%s/%s/%s/%s' % (srr_pre, srr00, srr, srr)
        print('\t'.join([url + '_1.fastq.gz', '0', url + '_2.fastq.gz', '0', toks[-1]]))
    else:
        url = 'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/%s/%s/%s' % (srr_pre, srr, srr)
        print('\t'.join([url + '_1.fastq.gz', '0', url + '_2.fastq.gz', '0', toks[-1]]))

sra_url = 'https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi'
cmd = "curl -s '%s?save=efetch&db=sra&rettype=runinfo&term=%s'" % (sra_url, '|'.join(query))
o, e = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True).communicate()
print(o, file=sys.stderr)
