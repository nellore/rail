#!/usr/bin/env python

"""
fastq2tab.py

Author: Ben Langmead
Date: 7/31/2013
Contact: langmea@cs.jhu.edu

Convert FASTQ to tab-delimited format.

TODO: Handle paired-end inputs.
"""

import sys
import gzip

for fn in sys.argv[1:]:
    fh = gzip.GzipFile(fn, 'r') if fn.endswith('.gz') else open(fn, 'r')
    print >> sys.stderr, "Processing '%s' ..." % fn
    while True:
        name1 = fh.readline().rstrip()
        if len(name1) == 0: break
        seq   = fh.readline().rstrip()
        name2 = fh.readline()
        qual  = fh.readline().rstrip()
        print '\t'.join([name1, seq, qual])
    fh.close()
