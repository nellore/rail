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

'''

import sys
import argparse
import time
import subprocess
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

args = parser.parse_args()

ninp = 0                   # # lines input so far
nout = 0                   # # lines output so far
last_samp = "\t"           # id of previous partition
cov = dict()               # coverage histogram
totcov = 0                 # total coverage
totnonz = 0                # total positions with non-0 coverage

last_chr = "\t"            # name of previous chromosome
last_samp = "\t"           # last sample last
frag_st  = -1              # start offset of current fragment
frag_dep = -1              # depth count of the current fragment
#samp_out = open(fname,'w')    # file handle for samples

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

if args.hadoop_exe!="":
    fname = "%s/temp_file"%(args.out_dir)  #temp file used to make a global file handle
    print >>sys.stderr,fname
    proc = subprocess.Popen([args.hadoop_exe, 'fs', '-put', '-', fname ], stdin=subprocess.PIPE)
    proc.stdin.close()
    proc.wait()
else:
    fname = "%s/temp_file"%(args.out_dir)  #temp file used to make a global file handle
    print >>sys.stderr,fname
    samp_out = open(fname,'w')
#proc = open("temp.txt",'w')   
#proc.close()

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
    
    ##TODO: incorporate bigBED conversion here
    #This is for hadoop mode
    if args.hadoop_exe!="" and (last_samp =="\t" or last_samp != samp):  #initialize file handles
        proc.stdin.close()
        proc.wait()
        last_chr = chr_name
        frag_st = pos
        frag_dep = cv
        out_fname = "%s/%s"%(args.out_dir,samp)
        #samp_out = open(out_fname,'w')
        proc = subprocess.Popen([args.hadoop_exe, 'fs', '-put', '-', out_fname], stdin=subprocess.PIPE)
        #proc = open("temp.txt", "w")        
    elif args.hadoop_exe!="" and  frag_dep!=cv: #record a new entry
        line = "%s\t%d\t%d\t%d\n"%(chr_name,frag_st,pos,frag_dep)
        #samp_out.write(line)
        proc.stdin.write(line)
        #proc.write(line)
        frag_st = pos
        frag_dep = cv
    
    #This is for local mode
    if args.hadoop_exe=="" and (last_samp =="\t" or last_samp != samp):  #initialize file handles
        last_chr = chr_name
        frag_st = pos
        frag_dep = cv
        out_fname = "%s/%s"%(args.out_dir,samp)
        samp_out = open(out_fname,'w')
    elif args.hadoop_exe=="" and  frag_dep!=cv: #record a new entry
        line = "%s\t%d\t%d\t%d\n"%(chr_name,frag_st,pos,frag_dep)
        samp_out.write(line)
        frag_st = pos
        frag_dep = cv
    
    
    
        
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

