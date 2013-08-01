#!/usr/bin/env python

"""
fastq2tab.py

Author: Ben Langmead
Date: 7/31/2013
Contact: langmea@cs.jhu.edu

Convert FASTQ to tab-delimited format.

TODO: Handle paired-end inputs.
"""

import os
import gzip

def preprocess(fh, outfh, lab, filename=None):
    fh_is_file = True
    if not isinstance(fh, file):
        if filename is None:
            filename = os.path.basename(fh)
        fh_is_file = False
        fh = gzip.GzipFile(fh, 'r') if fh.endswith('.gz') else open(fh, 'r')
    assert filename is not None
    n = 0
    while True:
        name = fh.readline().rstrip()
        if len(name) == 0: break
        seq   = fh.readline().rstrip()
        fh.readline()
        qual  = fh.readline().rstrip()
        outfh.write('\t'.join([';'.join(["FN:" + filename.replace(';', '_'), "LB:" + lab.replace(';', '_'), str(n)]), seq, qual]))
        outfh.write('\n')
        n += 1
    if not fh_is_file:
        fh.close()

if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(description='Preprocess a FASTQ file into a file that can be used with Myrna & friends')
    
    parser.add_argument(\
        '--input', metavar='PATHS', type=str, nargs='+', required=True, help='Use this in FN: field in output read name')
    parser.add_argument(\
        '--output', metavar='PATHS', type=str, nargs='+', required=True, help='Write these preprocessed files')
    parser.add_argument(\
        '--filename', metavar='STR', type=str, help='Use this in FN: field in output read name')
    parser.add_argument(\
        '--label', metavar='STRS', type=str, nargs='+', required=True, help='Use this in LB: field in output read name')
    parser.add_argument(\
        '--test', action='store_const', const=True, default=False, help='Do unit tests')
    
    args = parser.parse_args()
    assert len(args.input) == len(args.label)
    assert len(args.input) == len(args.output)
    for fn, outfn, lab in zip(args.input, args.output, args.label):
        with open(outfn, 'w') as ofh:
            preprocess(fn, ofh, lab)
