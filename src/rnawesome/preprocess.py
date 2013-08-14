#!/usr/bin/env python

"""
preprocess.py

Author: Ben Langmead
Date: 7/31/2013
Contact: langmea@cs.jhu.edu
"""

desc = """
Convert FASTQ or FASTA to tab-delimited format.  Input can be remote to begin
with, in which case it's downloaded first.  Output can be remote, in which
case files are uploaded as they're completed.
"""

import os
import gzip
import site

base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
site.addsitedir(os.path.join(base_path, "util"))

import url
import filemover

class RecordHandler(object):
    """ Takes read records and handles the process of writing them to
        (rotating) files, compressing the files, uploading file, and deleting
        trash. """
    
    def __init__(self, outfn, pushDest=None, mover=None, maxnucs=120000000, gzip=False, keep=False):
        self.outfnPre, self.outfn = outfn, None
        self.gzip = gzip
        self.pushDest = pushDest
        self.mover = mover
        self.maxnucs = maxnucs
        self.n, self.nfile, self.nnucs = 0, 0, 0
        self.prefix = None
        self.ofh = None
        self.first = True
        self.keep = keep
    
    def setPrefix(self, pre):
        self.prefix = pre
    
    def push(self):
        """ Push the file we just finished writing """
        if self.pushDest is not None:
            print >> sys.stderr, "  Pushing '%s' to '%s' ..." % (self.outfn, self.pushDest.toUrl())
            self.mover.put(self.outfn, self.pushDest)
            if not self.keep:
                os.remove(self.outfn)
    
    def nextFile(self):
        self.finish()
        # Make new filename
        ofn = '.'.join([self.outfnPre, str(self.nfile)])
        self.nfile += 1
        if self.gzip: ofn += ".gz"
        self.outfn = ofn
        # Open (note: may be gzipped)
        self.ofh = gzip.open(ofn, 'w') if self.gzip else open(ofn, 'w')
    
    def add(self, seq1, seq2=None, qual1=None, qual2=None):
        """ Add a record, which might necessitate opening a new file, and
            possibly pushing the old one. """
        
        assert self.prefix is not None
        assert seq1 is not None
        
        if self.first or self.nnucs > self.maxnucs:
            self.n, self.nnucs = 0, 0
            self.nextFile()
            self.first = False
        
        if qual1 is None: qual1 = 'I' * len(seq1)
        if seq2 is not None and qual2 is None: qual2 = 'I' * len(seq2)
        pre, ofh = self.prefix, self.ofh
        if seq2 is None:
            ofh.write('\t'.join([';'.join([pre, str(self.n)]), seq1, qual1]))
        else:
            ofh.write('\t'.join([';'.join([pre, str(self.n)]), seq1, qual1, seq2, qual2]))
        ofh.write('\n')
        
        self.n += 1
        self.nnucs += len(seq1)
        if seq2 is not None:
            self.nnucs += len(seq2)
    
    def finish(self):
        # Close previous output filehandle
        if self.ofh is not None:
            self.ofh.close()
            self.push()

def hashize(s, n=10):
    ret = hex(hash(s))[-n:]
    while len(ret) < n:
        ret = '0' + ret
    return ret

def preprocess(handler, lab, fh1, fh2=None, inputFormat="fastq", filename=None, includeFilename=False):
    """ Preprocess an input file, which may be in either FASTA or FASTQ format,
        into a collection of 1 or more tab-delimited output files. """
    is_file = [True, True]
    fhs = [fh1, fh2]
    fullname = filename
    for i in xrange(2):
        if fhs[i] is not None and not isinstance(fhs[i], file):
            if filename is None:
                fullname = fhs[i]
                filename = os.path.basename(fhs[i])
            is_file[i] = False
            fhs[i] = gzip.GzipFile(fhs[i], 'r') if fhs[i].endswith('.gz') else open(fhs[i], 'r')
    fh1, fh2 = fhs
    
    fullnameHash = hashize(fullname)
    n = 0
    fields = []
    if includeFilename:
        # TOOD: should be careful not to pass on any undesirable chars
        fields.append("FN:" + filename.replace(';', '_'))
    fields.append("FH:" + fullnameHash)
    fields.append("LB:" + lab.replace(';', '_'))
    handler.setPrefix(';'.join(fields))
    
    if inputFormat == "fastq":
        while True:
            if len(fh1.readline()) == 0:
                break
            seq = fh1.readline().rstrip()
            fh1.readline() # skip name line 2
            qual  = fh1.readline().rstrip()
            if fh2 is not None:
                assert len(fh2.readline()) > 0
                seq2 = fh2.readline().rstrip()
                fh2.readline() # skip name line 2
                qual2 = fh2.readline().rstrip()
                handler.add(seq, seq2, qual, qual2)
            else:
                handler.add(seq, qual1=qual)
            n += 1
    else:
        assert inputFormat == "fasta"
        lastLine = fh1.readline()
        while len(lastLine) > 0:
            seqlines = []
            while True:
                ln = lastLine = fh1.readline()
                if len(ln) == 0: break
                if ln[0] == '>': break
                seqlines.append(ln.rstrip())
            if fh2 is not None:
                seqlines2 = [], []
                while True:
                    ln = lastLine = fh2.readline()
                    if len(ln) == 0: break
                    if ln[0] == '>': break
                    seqlines2.append(ln.rstrip())
                handler.add(''.join(seqlines), seq2=''.join(seqlines2))
            else:
                handler.add(''.join(seqlines))
            n += 1
    
    if not is_file[0]: fh1.close()
    if fh2 is not None and not is_file[1]: fh2.close()
    
    handler.finish()

def guessFasta(fn):
    for ext in [ '.fa', '.fasta', '.fna', '.fas' ]:
        for extext in [ '', '.gz' ]:
            if fn.endswith(ext + extext):
                return True
    return False

if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(description=desc)
    
    parser.add_argument(\
        '--input1', metavar='PATHS', type=str, nargs='+', required=False, help='Mate 1 files (or unpaired read files)')
    parser.add_argument(\
        '--input2', metavar='PATHS', type=str, nargs='+', required=False, help='Mate 2 files, parallel to those specified with --mate1')
    parser.add_argument(\
        '--output', metavar='PATHS', type=str, nargs='+', required=False, help='Write these preprocessed files')
    parser.add_argument(\
        '--filename', metavar='STR', type=str, help='Use this in FN: field in output read name')
    parser.add_argument(\
        '--nucs-per-file', metavar='INT', type=int, default=120000000, help='Allow a max of this many nucleotides per output file')
    parser.add_argument(\
        '--label', metavar='STRS', type=str, nargs='+', required=False, help='Use this in LB: field in output read name')
    parser.add_argument(\
        '--gzip-output', action='store_const', const=True, default=False, help='Gzip compress output')
    parser.add_argument(\
        '--push', metavar='URL', type=str, required=False, help='Upload output files to this URL')
    parser.add_argument(\
        '--ignore-first-token', action='store_const', const=True, default=False, help='Throw away first token of input; useful in Hadoop streaming context')
    parser.add_argument(\
        '--fasta', action='store_const', const=True, default=False, help='Force preprocessor to consider input to be FASTA')
    parser.add_argument(\
        '--keep', action='store_const', const=True, default=False, help='Keep input files that were downloaded; default is to delete them once preprocessed')
    parser.add_argument(\
        '--test', action='store_const', const=True, default=False, help='Do unit tests')
    
    filemover.addArgs(parser)
    
    args = parser.parse_args()
    
    inp1, inp2, out, lab = args.input1, args.input2, args.output, args.label
    if inp1 is None:
        # TODO: Should we handle these one-by-one instead of sucking them all in here?
        import string
        import sys
        inp1, inp2, out, lab = [], [], [], []
        for ln in sys.stdin:
            ln = ln.strip()
            if args.ignore_first_token:
                if '\t' in ln:
                    # Useful because Hadoop streaming sometimes (always?)
                    # seems to include a file offset as the first, leftmost
                    # token
                    ln = ln[ln.index('\t')+1:]
                else: continue
            if len(ln) == 0: continue
            if ln[0] == '#': continue
            toks = string.split(ln, '\t')
            if len(toks) != 3 and len(toks) != 5:
                raise RuntimeError("Malformed line; had # tokens != 3 or 5:\n" + ln)
            inp1.append(toks[0])
            lab.append(toks[-1])
            inp2.append(toks[2] if len(toks) == 5 else None)
            out.append(None)
    
    if inp2 is None: inp2 = [None] * len(inp1)
    assert len(inp1) == len(lab)
    assert len(inp1) == len(inp2)
    assert len(inp1) == len(out)
    if len(inp1) > 0:
        push = url.Url(args.push) if args.push is not None else None
        mover = filemover.FileMover(args=args)
        for fn1, fn2, outfn, lab in zip(inp1, inp2, out, lab):
            if fn2 is None:
                print >> sys.stderr, "Processing unpaired URL '%s' ..." % fn1
            else:
                print >> sys.stderr, "Processing paired URLs '%s', '%s' ..." % (fn1, fn2)
            # Determine whether input is fasta
            inputFormat = "fastq"
            if args.fasta or guessFasta(fn1): inputFormat = "fasta"
            # Does input file need to be pulled down?
            toDelete = []
            inpUrl1 = url.Url(fn1)
            if inpUrl1.isNotLocal():
                mover.get(inpUrl1)
                fn1 = os.path.basename(fn1)
                assert os.path.exists(fn1)
                if not args.keep: toDelete.append(fn1)
                if fn2 is not None:
                    inpUrl2 = url.Url(fn2)
                    assert inpUrl2.isNotLocal()
                    mover.get(inpUrl2)
                    fn2 = os.path.basename(fn2)
                    assert os.path.exists(fn2)
                    if not args.keep: toDelete.append(fn2)
            # Come up with an output filename
            if outfn is None: outfn = hashize(fn1)
            handler = RecordHandler(outfn, push, mover, args.nucs_per_file, args.gzip_output, args.keep)
            preprocess(handler, lab, fn1, fn2, inputFormat)
            for fn in toDelete: os.remove(fn)
