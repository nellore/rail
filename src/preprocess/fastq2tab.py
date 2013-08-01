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
import site

base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
site.addsitedir(os.path.join(base_path, "util"))

import url

def preprocess(fh, outfh, lab, inputFormat="fastq", filename=None):
    fh_is_file = True
    if not isinstance(fh, file):
        if filename is None:
            filename = os.path.basename(fh)
        fh_is_file = False
        fh = gzip.GzipFile(fh, 'r') if fh.endswith('.gz') else open(fh, 'r')
    assert filename is not None
    n = 0
    pre = ';'.join([filename.replace(';', '_'), lab.replace(';', '_')])
    
    def handle(seq, qual=None):
        if qual is None: qual = 'I' * len(seq)
        outfh.write('\t'.join([';'.join([pre, str(n)]), seq, qual]))
        outfh.write('\n')
    
    if inputFormat == "fastq":
        while True:
            if len(fh.readline()) == 0:
                break
            seq = fh.readline().rstrip()
            fh.readline() # skip name line 2
            qual  = fh.readline().rstrip()
            handle(seq, qual, n)
            n += 1
    else:
        assert inputFormat == "fasta"
        lastLine = fh.readline()
        while len(lastLine) > 0:
            seqlines = []
            while True:
                ln = lastLine = fh.readline()
                if len(ln) == 0:
                    break
                if ln[0] == '>':
                    break
                seqlines.append(ln.rstrip())
            handle(''.join(seqlines), None, n)
            n += 1
    if not fh_is_file:
        fh.close()

if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(description='Preprocess a FASTQ file into a file that can be used with Myrna & friends')
    
    parser.add_argument(\
        '--input', metavar='PATHS', type=str, nargs='+', required=False, help='Use this in FN: field in output read name')
    parser.add_argument(\
        '--output', metavar='PATHS', type=str, nargs='+', required=True, help='Write these preprocessed files')
    parser.add_argument(\
        '--filename', metavar='STR', type=str, help='Use this in FN: field in output read name')
    parser.add_argument(\
        '--records-per-file', metavar='INT', type=int, default=500000, help='Allow a max of this many read records per output file')
    parser.add_argument(\
        '--label', metavar='STRS', type=str, nargs='+', required=True, help='Use this in LB: field in output read name')
    parser.add_argument(\
        '--gzip-output', action='store_const', const=True, default=False, help='Gzip compress output')
    parser.add_argument(\
        '--push', metavar='URL', type=str, required=False, help='Upload output files to this URL')
    parser.add_argument(\
        '--s3cfg', metavar='STR', type=str, required=False, help='s3cmd configuration file to use')
    parser.add_argument(\
        '--test', action='store_const', const=True, default=False, help='Do unit tests')
    
    args = parser.parse_args()
    
    inp, out, lab = args.input, args.output, args.label
    if args.input is None:
        # Take input from stdin
        pass
    
    assert len(args.input) == len(args.label)
    assert len(args.input) == len(args.output)
    myopen = gzip.open if args.gzip_output else open
    pushUrl = url.Url(args.push) if args.push is not None else None
    for fn, outfn, lab in zip(args.input, args.output, args.label):
        if args.gzip_output and not outfn.endswith('.gz'): outfn += ".gz"
        with myopen(outfn, 'w') as ofh: preprocess(fn, ofh, lab)
        if pushUrl is not None:
            if pushUrl.isS3():
                os.system("s3cmd %s put %s %s" % ("-c %s" % args.s3cfg if args.s3cfg is not None else "", outfn, pushUrl.toNonNativeUrl()))
            else:
                os.system("hadoop fs -put %s %s" % (outfn, pushUrl.toUrl()))
