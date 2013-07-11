'''
normalize.py
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
timeSt = time.clock()

parser = argparse.ArgumentParser(description=\
    'Takes per-position counts for given samples and calculates a summary '
    'statistic such as median or 75% percentile.')

parser.add_argument(\
    '--percentile', metavar='FRACTION', type=float, required=False, default=0.75,
    help='For a given sample, the per-position percentile to extract as the normalization factor')

parser.add_argument(\
    '--out_dir', type=str, required=False, default="",
    help='The directory where all of the coverage vectors for each sample will be stored')

parser.add_argument(\
    '--hadoop_exe', type=str, required=False, default="",
    help='The location of the hadoop executable.')

parser.add_argument(\
    '--bigbed_exe', type=str, required=False, default="",
    help='The location of the bigbed executable.')

parser.add_argument(\
    '--chrom_sizes', type=str, required=False, default="",
    help='The location of chrom_sizes file required for bigbed conversion.')

args = parser.parse_args()

ninp = 0                   # # lines input so far
nout = 0                   # # lines output so far
cov = dict()               # coverage histogram
totcov = 0                 # total coverage
totnonz = 0                # total positions with non-0 coverage

fname = "temp_file"        #temp file used to make a global file handle
samp_out = open(fname,'w')
last_chr = "\t"            # name of previous chromosome
last_samp = "\t"           # last sample last
frag_st  = -1              # start offset of current fragment
last_pos = -1              # the last position 
frag_dep = -1              # depth count of the current fragment

def percentile(cov):
    ''' Given a histogram (dictionary mapping keys to integers), return
        the upper quartile.  Same as 75th percentile.  Return a string.
        '''
    cur = 0
    for (k, v) in sorted(cov.iteritems(), reverse=True):
        cur += v
        assert cur <= totnonz
        if cur > ((1.0 - args.percentile) * totnonz):
            return str(k)
    raise RuntimeError("Should not reach this point")

#proc = subprocess.Popen(['/damsl/software/hadoop/hadoop-1.1.2/bin/hadoop', 'fs', '-put', '-', fname ], stdin=subprocess.PIPE)

def moveToHDFS(fin,fout):
    put_proc = subprocess.Popen([args.hadoop_exe, 'fs', '-put',fin, fout])
    put_proc.wait()


for ln in sys.stdin:
    ln = ln.rstrip()
    toks = ln.split('\t')
    assert len(toks) == 4
    samp, chr_name, pos, cv = toks[0],toks[1],int(toks[2]),int(toks[3])
    if samp != last_samp and last_samp != "\t":
        print "%s\t%s" % (last_samp, percentile(cov)) #just for testing now
        nout += 1
        cov = dict()
        totcov, totnonz = 0, 0
        
    if last_samp == "\t":
        samp_out.close()
        last_chr, frag_st, frag_dep, fname = chr_name, pos, cv, samp
        samp_out = open(fname,'w')        
    elif last_samp != samp:  #initialize file handles and convert to bigbed
        samp_out.close()
        bb_file = "%s.bb"%(last_samp) if args.hadoop_exe!="" else "%s/%s.bb"%(args.out_dir,last_samp)
        bigbed_proc = subprocess.Popen([args.bigbed_exe, fname, args.chrom_sizes, bb_file])
        bigbed_proc.wait()
        if args.hadoop_exe!="": #For hadoop mode - moves local file to HDFS
            out_fname = "%s/%s.bb"%(args.out_dir,last_samp)
            moveToHDFS(bb_file,out_fname)
        last_chr, frag_st, frag_dep, fname = chr_name, pos, cv, samp
        os.remove(fname)
        samp_out = open(fname,'w')
    elif last_chr!=chr_name:
        last_chr, frag_st, frag_dep, fname = chr_name, pos, cv, samp
    elif frag_dep!=cv and abs(last_pos-pos)==1: #record a new entry
        line = "%s\t%d\t%d\t%d\n"%(chr_name,frag_st,pos,frag_dep)
        samp_out.write(line)
        frag_st, frag_dep = pos, cv
    elif abs(last_pos-pos)>1:
        line = "%s\t%d\t%d\t%d\n"%(chr_name,frag_st,pos,0)
        samp_out.write(line)
        frag_st,frag_dep = pos, cv
    
    last_pos = pos
    last_samp = samp
    ninp += 1
    cov[cv] = cov.get(cv, 0) + 1
    totcov += cv
    totnonz += 1

if last_samp != "\t":
    print "%s\t%s" % (last_samp, percentile(cov))
    nout += 1

# Done
timeEn = time.clock()
print >>sys.stderr, "DONE with normalize.py; in/out = %d/%d; time=%0.3f secs" % (ninp, nout, timeEn-timeSt)

