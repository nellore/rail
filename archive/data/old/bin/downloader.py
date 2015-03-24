"""
downloader.py

Download samples of reads from files on the Internet
"""

import os
import sys
import urllib
import urllib2
import string
import random

class ReservoirSampler(object):
    """ Simple reservoir sampler """
    def __init__(self, k):
        self.k = k
        self.r = []
        self.n = 0
    
    def add(self, obj):
        if self.n < self.k:
            self.r.append(obj)
        else:
            j = random.randint(0, self.n)
            if j < self.k:
                self.r[j] = obj
        self.n += 1
    
    def draw(self):
        return random.choice(self.r)
    
    def __len__(self):
        return self.n

def parseOptions(st):
    """ Parse option string, return options hash. """
    opts = {}
    for tok in string.split(st, ","):
        if tok == "gzip":
            opts["gzip"] = True
        elif tok.startswith("sample"):
            k, v = string.split(tok, "=")
            assert k == "sample"
            opts["sample"] = int(v)
        elif tok.startswith("format"):
            k, v = string.split(tok, "=")
            assert k == "format"
            opts["sample"] = v
    return opts

def processFile(url, loc):
    urllib.urlretrieve(url, loc)

def processReadFile(url, loc, opts, collector):
    fh = urllib2.urlopen(url)
    if "gzip" in opts:
        fh = gzip.GzipFile(fileobj=fh)
    if "format" not in opts or opts["format"] == "fastq":
        while True:
            head1 = fh.readline()
            if len(head1) == 0:
                break
            seq = fh.readline()
            head2 = fh.readline()
            qual = fh.readline()
            collector.add(''.join([head1, seq, head2, qual]))
    elif opts["format"] == "fasta":
        while True:
            head1 = fh.readline()
            if len(head1) == 0:
                break
            seq = fh.readline()
            collector.add(''.join([head1, seq, head2, qual]))
    fh.close()

class WriterCollector(object):
    """ Collector that, given a read, will immediately write it to file
        with given filename.  finalize method does nothing. """
    
    def __init__(self, ofn):
        self.fh = open(ofn, 'w')
    
    def __del__(self):
        self.fh.close()
    
    def add(self, o):
        self.fh.write(o)
    
    def finalize(self):
        pass

class SamplerCollector(object):
    """ Collector that, given a read, will stick it in its
        ReservoirSampler.  The finalize method then writes the reads in
        the reservoir out to the given filename. """
    
    def __init__(self, k, ofn):
        self.samp = ReservoirSampler(k)
        self.ofn = ofn
    
    def add(self, o):
        self.samp.add(o)
    
    def finalize(self):
        """ Write sampled reads in random order """
        with open(self.ofn, 'w') as fh:
            for rec in self.samp.r:
                fh.write(rec)

def processManifest(fn):
    """ Go through each URL in the manifest """
    with open(fn, 'r') as fh:
        for ln in fh:
            ln = ln.rstrip()
            if len(ln) == 0 or ln[0] == '#':
                continue
            toks = string.split(ln, '\t')
            opts = {}
            if len(toks) > 2:
                url, loc, optstr = toks
                opts = parseOptions(optstr)
            else:
                url, loc = toks
            d, b = os.path.split(loc)
            tmpWork = os.path.join(d, "." + b + ".working")
            tmpDone = os.path.join(d, "." + b + ".done")
            
            if not os.path.exists(tmpWork) and not os.path.exists(tmpDone):
                with open(tmpWork, 'w') as fh: fh.write('\n')
                
                if "reads" in opts:
                    if "%MATE%" in url:
                        assert '%MATE%' in loc
                        url1 = url.replace('%MATE%', '1')
                        url2 = url.replace('%MATE%', '2')
                        loc1 = loc.replace('%MATE%', '1')
                        loc2 = loc.replace('%MATE%', '2')
                        print >> sys.stderr, "Downloading paired-end reads files '%s'/'%s' ..." % (loc1, loc2)
                    else:
                        print >> sys.stderr, "Downloading unpaired reads file '%s' ..." % loc
                        if url.endswith(".gz"):
                            opts["gzip"] = True
                        elif url.endswith(".bz2"):
                            opts["bzip2"] = True
                        if "sample" in opts:
                            collector = SamplerCollector(opts["sample"], loc)
                        else:
                            collector = WriterCollector(loc)
                        processReadFile(url, loc, opts, collector)
                        collector.finalize()
                else:
                    print >> sys.stderr, "Downloading file '%s' ..." % loc
                    processFile(url, loc)
                
                with open(tmpDone, 'w') as fh: fh.write('\n')
                os.remove(tmpWork)
            elif os.path.exists(tmpWork):
                print >> sys.stderr, "'%s' already in process ..." % loc
            elif os.path.exists(tmpDone):
                print >> sys.stderr, "'%s' already done ..." % loc

def go(manifest):
    processManifest(manifest)

if __name__ == "__main__":
    import argparse
        
    parser = argparse.ArgumentParser(\
        description='Download and sample reads')
    
    parser.add_argument(\
        '--manifest', metavar='path', type=str, required=True,
        help='Manifest file to download.')
    parser.add_argument(\
        '--test', dest='test', action='store_const', const=True, default=False,
        help='Do unit tests')
    parser.add_argument(\
        '--profile', dest='profile', action='store_const', const=True,
        default=False, help='Profile the mode using python cProfile module')

    args = parser.parse_args()

    if not args.test:
        if args.profile:
            import cProfile
            cProfile.run('go(args.manifest)')
        else:
            go(args.manifest)
    
    if args.test:
        import unittest
        
        class TestCases(unittest.TestCase):
        
            def test_blah(self):
                pass
        
        unittest.main(argv=[sys.argv[0]])
        sys.exit()
