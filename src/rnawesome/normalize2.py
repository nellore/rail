'''
normalize2.py
(after walk_prenorm.py, before normalize_post.py)

Given position counts for a given sample, calculates some relevant per-sample
summary statistic, such as median, mean, total, or upper quartile.  Summary
statistics are passed along to normalize_post.py, which prints them all to a
file.

Tab-delimited input tuple columns:
 1. Sample label
 2. Chromosome name
 3. Genome Position
 4. Count (at some position)

Binning/sorting prior to this step:
 1. Binned by sample label
 2. Sorted by reference id
 3. Sorted by genome position

Other files:
 Sample files specified by sample name

Tab-delimited output tuple columns:
 1. Sample label
 2. Normalization factor

Todo: Push all bed and bigbed files to --out_dir on local mode

Note to self:  all of the reads are from the first sequence

'''

import sys
import argparse
import time
import subprocess
import os
import site
import string
timeSt = time.clock()

parser = argparse.ArgumentParser(description=\
    'Takes per-position counts for given samples and calculates a summary '
    'statistic such as median or 75% percentile.')

base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
site.addsitedir(os.path.join(base_path, "util"))

import url
import path
import filemover

parser.add_argument(\
    '--percentile', metavar='FRACTION', type=float, required=False, default=0.75,
    help='For a given sample, the per-position percentile to extract as the normalization factor')

parser.add_argument(\
    '--out_dir', type=str, required=True, default=".",
    help='The URL where all of the coverage vectors for each sample will be stored')

parser.add_argument(\
    '--hadoop_exe', type=str, required=False, default="",
    help='The location of the hadoop executable.')

parser.add_argument(\
    '--bigbed_exe', type=str, required=False, default="bedToBigBed",
    help='The location of the bigbed executable.')

parser.add_argument(\
    '--chrom_sizes', type=str, required=False,
    help='The location of chrom_sizes file required for bigbed conversion.')

parser.add_argument(\
    '--faidx', type=str, required=False,
    help='Path to a FASTA index that we can use instead of --chrom_sizes.')

parser.add_argument(\
    '--verbose', action='store_const', const=True, default=False, help='Prints out extra debugging statements')

filemover.addArgs(parser)

args = parser.parse_args()

if path.which(args.bigbed_exe) is None:
    raise RuntimeError("Could not find bedToBigBed exe; tried '%s'" % args.bigbed_exe)

ninp, nout = 0, 0          # # lines input/output so far
covhist = {}                   # coverage histogram

ofn, ofh = ".tmp.bed", None
last_refid = "\t"          # name of previous chromosome
last_samp = "\t"           # last sample last
last_cov = 0               # last coverage
last_pos = -1              # the last position 

def percentile(covhist, percentile=0.75):
    ''' Given a histogram (dictionary mapping keys to integers), return
        the element at the given percentile. '''
    cur = 0
    totnonz = sum(covhist.itervalues())
    for (k, v) in sorted(covhist.iteritems(), reverse=True):
        cur += v
        assert cur <= totnonz
        if cur > ((1.0 - percentile) * totnonz): return k
    raise RuntimeError("Should not reach this point")

def moveToHDFS(fin, fout):
    put_proc = subprocess.Popen([args.hadoop_exe, 'fs', '-put', fin, fout], stdout=sys.stderr)
    put_proc.wait()

mover = filemover.FileMover(args=args)
outUrl = url.Url(args.out_dir)

chromSizes = args.chrom_sizes
delChromSizes = False
if chromSizes is None:
    if args.faidx is None:
        raise RuntimeError("Must specify either --chrom_sizes or --faidx")
    import tempfile
    fh = tempfile.NamedTemporaryFile(mode='w', delete=False)
    with open(args.faidx, 'r') as ifh:
        for ln in ifh:
            toks = string.split(ln.rstrip(), '\t')
            print >>fh, " ".join([toks[0], toks[1]])
    fh.close()
    chromSizes = fh.name
    delChromSizes = True

# Parse the chrom_sizes file
reflens = {}
with open(chromSizes, 'r') as fh:
    for ln in fh:
        ln = ln.rstrip()
        if len(ln) == 0: continue
        ref, length = string.split(ln)
        length = int(length)
        reflens[ref] = length

assert os.path.exists(chromSizes)

def bedToBigBed(ifn, ofn, chromSizes):
    """ Run bedToBigBed on input file ifn, specifying chromSizes as file with
        reference lengths, and store output BigBed in ofn """
    assert os.path.exists(ifn)
    assert os.path.exists(chromSizes)
    assert not os.path.exists(ofn), "Already wrote '%s'" % ofn
    bigbed_cmd = [args.bigbed_exe, ifn, chromSizes, ofn]
    bigbed_proc = subprocess.Popen(bigbed_cmd, stdout=sys.stderr)
    ret = bigbed_proc.wait()
    if ret != 0:
        raise RuntimeError("bedToBigBed command '%s' returned exitlevel %d" % (' '.join(bigbed_cmd), ret))
    if args.verbose:
        assert os.path.exists(ofn)
        print >> sys.stderr, "bedToBigBed command '%s' succeeded" % ' '.join(bigbed_cmd)

def handleBigBed(bedfn, name, outUrl, chromSizes):
    assert os.path.exists(bedfn)
    bbfn = "%s.bb" % name if outUrl.isNotLocal() else "%s/%s.bb" % (outUrl.toUrl(), name)
    bedToBigBed(bedfn, bbfn, chromSizes)
    if outUrl.isNotLocal():
        assert os.path.exists(bbfn)
        mover.put(bbfn, outUrl.plus(bbfn))
    os.remove(bedfn)

ofh = open(ofn, 'w')
for ln in sys.stdin:
    
    toks = ln.rstrip().split('\t')
    assert len(toks) == 4
    samp, refid, pos, cov = toks[0], toks[1], int(toks[2]), int(toks[3])
    assert samp != last_samp or refid != last_refid or pos > last_pos
    assert refid in reflens, "Reference ID '%s' not amongst those in FASTA index: '%s'" % (refid, str(reflens.keys()))
    
    finishedSamp = samp != last_samp and last_samp != "\t"
    if finishedSamp:
        assert ofh is not None
        if last_cov != 0 and last_pos < reflens[last_refid]:
            st, en = last_pos, reflens[last_refid]
            ofh.write("%s\t%d\t%d\t%d\n" % (last_refid, st, en, last_cov))
            covhist[last_cov] = covhist.get(last_cov, 0) + (en - st)
        print "%s\t%d" % (last_samp, percentile(covhist, args.percentile)) 
        nout += 1
        covhist = {}
        ofh.close()
        handleBigBed(ofn, last_samp, outUrl, chromSizes)
        ofh = open(ofn, 'w')
    elif samp == last_samp:
        adjacentToPrev = refid == last_refid
        if adjacentToPrev:
            if last_cov != 0:
                covhist[last_cov] = covhist.get(last_cov, 0) + (pos - last_pos)
            if cov != last_cov and last_cov != 0:
                ofh.write("%s\t%d\t%d\t%d\n" % (refid, last_pos, pos, last_cov))
            elif last_cov != 0:
                assert cov == last_cov
                pos = last_pos # so next output interval extends back to prev pos
        elif last_cov != 0 and last_pos < reflens[last_refid]:
            st, en = last_pos, reflens[last_refid]
            ofh.write("%s\t%d\t%d\t%d\n" % (last_refid, st, en, last_cov))
            covhist[last_cov] = covhist.get(last_cov, 0) + (en - st)
    
    last_pos, last_refid, last_samp, last_cov = pos, refid, samp, cov
    ninp += 1

if last_samp != '\t':
    assert ofh is not None
    if last_cov != 0 and last_pos < reflens[last_refid]:
        st, en = last_pos, reflens[last_refid]
        ofh.write("%s\t%d\t%d\t%d\n" % (last_refid, st, en, last_cov))
        covhist[last_cov] = covhist.get(last_cov, 0) + (en - st)
    print "%s\t%d" % (last_samp, percentile(covhist, args.percentile)) 
    nout += 1
    ofh.close()
    handleBigBed(ofn, last_samp, outUrl, chromSizes)

if delChromSizes:
    os.remove(chromSizes)

# Done
timeEn = time.clock()
print >>sys.stderr, "DONE with normalize2.py; in/out = %d/%d; time=%0.3f secs" % (ninp, nout, timeEn-timeSt)
