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
import sys
import gzip
import site

base_path = os.path.abspath(
                    os.path.dirname(os.path.dirname(os.path.dirname(
                        os.path.realpath(__file__)))
                    )
                )
utils_path = os.path.join(base_path, 'rna', 'utils')
site.addsitedir(utils_path)
site.addsitedir(base_path)

from dooplicity.ansibles import Url
import filemover


class RecordHandler(object):
    """ Takes read records and handles the process of writing them to
        (rotating) files, compressing the files, uploading file, and deleting
        trash. """

    def __init__(self, output_filename_prefix, push_destination=None, mover=None, maxnucs=120000000, gzip_output=False,
                 gzip_level=3, keep=False):
        self.output_filename_prefix, self.outfn = output_filename_prefix, None
        self.gzip_output = gzip_output
        self.gzip_level = gzip_level
        self.push_destination = push_destination
        self.mover = mover
        self.maxnucs = maxnucs
        self.n, self.nfile, self.nnucs = 0, 0, 0
        self.ofh = None
        self.first = True
        self.keep = keep
    
    def push(self):
        """ Push the file we just finished writing """
        if self.push_destination is not None:
            assert self.outfn is not None
            dest_url = self.push_destination.to_url()
            print >> sys.stderr, "  Pushing '%s' to '%s' ..." % (self.outfn, dest_url)
            if self.push_destination.isLocal():
                if not os.path.exists(dest_url):
                    print >> sys.stderr, "    Making local destination directory: %s" % dest_url
                    os.mkdir(dest_url)
            self.mover.put(self.outfn, self.push_destination)
            if not self.keep:
                os.remove(self.outfn)
    
    def next_file(self):
        self.finish()
        # Make new filename
        ofn = '.'.join([self.output_filename_prefix, str(self.nfile)])
        self.nfile += 1
        if self.gzip_output:
            ofn += ".gz"
        self.outfn = ofn
        # Open (note: may be gzipped)
        self.ofh = gzip.open(ofn, 'wb', self.gzip_level) if self.gzip_output else open(ofn, 'w')
    
    def add(self, prefix1, seq1, prefix2=None, seq2=None, qual1=None, qual2=None):
        """ Add a record, which might necessitate opening a new file, and
            possibly pushing the old one. """
        
        assert seq1 is not None
        
        if self.first or self.nnucs > self.maxnucs:
            self.n, self.nnucs = 0, 0
            self.next_file()
            self.first = False
        
        if qual1 is None:
            qual1 = 'I' * len(seq1)
        if seq2 is not None and qual2 is None:
            qual2 = 'I' * len(seq2)
        ofh = self.ofh
        if seq2 is None:
            ofh.write('\t'.join([';'.join([prefix1, str(self.n)]), seq1, qual1]))
        else:
            ofh.write('\t'.join([';'.join([prefix1, str(self.n)]), seq1, qual1, ';'.join([prefix2, str(self.n)]), seq2, qual2]))
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


def preprocess(handler, lab, fh1, fh2=None, input_format="fastq", filename=None, include_filename=False):
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
    
    n = 0
    # TODO: should be careful not to pass on any undesirable chars
    filename_read_name = "FN:" + filename.replace(';', '_') + "|"
    label = "LB:" + lab.replace(';', '_')
    
    if input_format == "fastq":
        while True:
            if len(fh1.readline()) == 0:
                break
            seq = fh1.readline().rstrip()
            original_name = fh1.readline()  # skip name line 2
            qual = fh1.readline().rstrip()
            if fh2 is not None:
                assert len(fh2.readline()) > 0
                seq2 = fh2.readline().rstrip()
                original_name_2 = fh2.readline()  # skip name line 2
                qual2 = fh2.readline().rstrip()
                read_name_hash = 'ID:' + hashize(original_name + original_name_2 + seq + qual + seq2 + qual2)
                handler.add((filename_read_name + ';' if include_filename else '')
                                + read_name_hash + '/1;' + label, seq,
                                prefix2=(filename_read_name + ';' if include_filename else '')
                                + read_name_hash + '/2;' + label, seq2=seq2, qual1=qual, qual2=qual2)
            else:
                read_name_hash = 'ID:' + hashize(original_name + seq + qual)
                handler.add(''.join([filename_read_name, read_name_hash]) + ";%s" % label, seq, qual1=qual)
            n += 1
    else:
        assert input_format == "fasta"
        last_line = fh1.readline()
        while len(last_line) > 0:
            seqlines = []
            while True:
                ln = last_line = fh1.readline()
                if len(ln) == 0 or ln[0] == '>':
                    break
                seqlines.append(ln.rstrip())
            seq = ''.join(seqlines)
            if fh2 is not None:
                seqlines2 = []
                while True:
                    ln = last_line = fh2.readline()
                    if len(ln) == 0 or ln[0] == '>':
                        break
                    seqlines2.append(ln.rstrip())
                seq2 = ''.join(seqlines2)
                read_name_hash = 'ID:' + hashize(fullname + ln + seq + seq2)
                handler.add((filename_read_name + ';' if include_filename else '')
                                + read_name_hash + '/1;' + label, seq,
                                prefix2=(filename_read_name + ';' if include_filename else '')
                                + read_name_hash + '/2;' + label, seq2=seq2, qual1=qual, qual2=qual2)
            else:
                read_name_hash = 'ID:' + hashize(fullname + ln + seq + qual)
                handler.add(''.join([filename_read_name, read_name_hash]) + ";%s" % label, seq, qual1=qual)
            n += 1
    
    if not is_file[0]:
        fh1.close()
    if fh2 is not None and not is_file[1]:
        fh2.close()
    
    handler.finish()


def guess_fasta(fn):
    for ext in ['.fa', '.fasta', '.fna', '.fas']:
        for extext in ['', '.gz']:
            if fn.endswith(ext + extext):
                return True
    return False


def go(args):
    inp1, inp2, out, lab = args.input1, args.input2, args.output, args.label
    if inp1 is None:
        inp1, inp2, out, lab = [], [], [], []
        for ln in sys.stdin:
            ln = ln.strip()
            if args.ignore_first_token:
                if '\t' in ln:
                    # Useful because Hadoop streaming sometimes (always?)
                    # seems to include a file offset as the first, leftmost
                    # token
                    ln = ln[ln.index('\t')+1:]
                else:
                    continue
            if len(ln) == 0 or ln[0] == '#':
                continue
            toks = ln.split()
            if len(toks) != 3 and len(toks) != 5:
                raise RuntimeError("Malformed line; had # tokens != 3 or 5:\n" + ln)
            inp1.append(toks[0])
            lab.append(toks[-1])
            inp2.append(toks[2] if len(toks) == 5 else None)
            out.append(None)

    if inp2 is None:
        inp2 = [None] * len(inp1)
    assert len(inp1) == len(lab)
    assert len(inp1) == len(inp2)
    assert len(inp1) == len(out)
    if len(inp1) > 0:
        push = args.push
        if push is not None:
            push = Url(push)
            if push.to_url()[-1] != '/':
                push = Url(push.to_url() + '/')
            assert push.to_url()[-1] == '/', push.to_url()
        mover = filemover.FileMover(args=args)
        for fn1, fn2, outfn, lab in zip(inp1, inp2, out, lab):
            if fn2 is None:
                print >> sys.stderr, "Processing unpaired URL '%s' ..." % fn1
            else:
                print >> sys.stderr, "Processing paired URLs '%s', '%s' ..." % (fn1, fn2)
            # Determine whether input is fasta
            input_format = "fastq"
            if args.fasta or guess_fasta(fn1):
                input_format = "fasta"
            # Does input file need to be pulled down?
            to_delete = []
            input_url_1 = Url(fn1)
            if input_url_1.is_local:
                mover.get(input_url_1)
                fn1 = os.path.basename(fn1)
                assert os.path.exists(fn1)
                if not args.keep:
                    to_delete.append(fn1)
                if fn2 is not None:
                    input_url_2 = Url(fn2)
                    assert not input_url_2.is_local
                    mover.get(input_url_2)
                    fn2 = os.path.basename(fn2)
                    assert os.path.exists(fn2)
                    if not args.keep:
                        to_delete.append(fn2)
            # Come up with an output filename
            if outfn is None:
                outfn = hashize(fn1)
            handler = RecordHandler(outfn, push, mover, args.nucs_per_file, args.gzip_output, args.gzip_level,
                                    args.keep)
            preprocess(handler, lab, fn1, fn2, input_format, include_filename=args.include_filename)
            for fn in to_delete:
                os.remove(fn)


if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(description=desc)

    parser.add_argument('--input1', metavar='PATHS', type=str, nargs='+', required=False,
                        help='Mate 1 files (or unpaired read files)')
    parser.add_argument('--input2', metavar='PATHS', type=str, nargs='+', required=False,
                        help='Mate 2 files, parallel to those specified with --mate1')
    parser.add_argument('--output', metavar='PATHS', type=str, nargs='+', required=False,
                        help='Write these preprocessed files')
    parser.add_argument('--filename', metavar='STR', type=str, help='Use this in FN: field in output read name')
    parser.add_argument('--nucs-per-file', metavar='INT', type=int, default=120000000,
                        help='Allow a max of this many nucleotides per output file')
    parser.add_argument('--label', metavar='STRS', type=str, nargs='+', required=False,
                        help='Use this in LB: field in output read name')
    parser.add_argument('--gzip-output', action='store_const', const=True, default=False, help='Gzip compress output')
    parser.add_argument('--gzip-level', metavar='INT', type=int, default=3, help='Level of gzip compression to use')
    parser.add_argument('--push', metavar='URL', type=str, required=False, help='Upload output files to this URL')
    parser.add_argument('--ignore-first-token', action='store_const', const=True, default=False,
                        help='Throw away first token of input; useful in Hadoop streaming context')
    parser.add_argument('--include-filename', action='store_const', const=True, default=False,
                        help='Put filename in FN: name field')
    parser.add_argument('--fasta', action='store_const', const=True, default=False,
                        help='Force preprocessor to consider input to be FASTA')
    parser.add_argument('--keep', action='store_const', const=True, default=False,
                        help='Keep input files that were downloaded; default is to delete them once preprocessed')
    parser.add_argument('--test', action='store_const', const=True, default=False, help='Do unit tests')
    parser.add_argument('--profile', action='store_const', const=True, default=False, help='Profile code')

    filemover.add_args(parser)
    master_args = parser.parse_args()
    if master_args.profile:
        import cProfile
        cProfile.run('go(master_args)')
    else:
        go(master_args)
