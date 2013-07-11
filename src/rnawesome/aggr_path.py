'''
aggr_path.py
(after hmm.py)

Collects the list of states from hmm and bins them into separate files
according to permutation and sample

Tab-delimited output tuple columns:
 1. Dataset id (0 for test, >0 for permutations)
 2. Reference ID
 3. Reference offset (0-based)
 4. Length of run of positions with this state
 5. HMM state

Other files:
 One file per permutation
 1. Reference ID
 2. Reference start position
 3. Reference end position
 4. HMM state

'''

import os
import sys
import argparse
import time
import subprocess
timeSt = time.clock()

parser = argparse.ArgumentParser(description=\
    'Takes per-position counts for given samples and calculates a summary '
    'statistic such as median or 75% percentile.')

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

fname, samp_out = None, None
last_perm = "\t"                 #last permutation id
ninp = 0

def moveToHDFS(fin,fout):
    put_proc = subprocess.Popen([args.hadoop_exe, 'fs', '-put',fin, fout])
    put_proc.wait()

for ln in sys.stdin:
    ninp += 1
    ln = ln.rstrip()
    toks = ln.split('\t')
    assert len(toks)==5
    data_id, ref_id, ref_off, ref_len, hmm_st = toks[0], toks[1], int(toks[2]), int(toks[3]), int(toks[4])
    #print toks
    if last_perm=="\t":
        if samp_out is not None:
            samp_out.close()
        fname = "perm%s.bed"%(data_id)
        print>>sys.stderr,"ID:",data_id
        samp_out = open(fname,'w')
    elif last_perm!=data_id:
        assert samp_out is not None
        samp_out.close()
        out_fname = "%s/perm%s.bb" % (args.out_dir, last_perm)
        bb_file = "perm%s.bb" % (last_perm) if args.hadoop_exe != "" else out_fname
        bigbed_proc = subprocess.Popen([args.bigbed_exe, fname, args.chrom_sizes, bb_file])
        bigbed_proc.wait()
        if args.hadoop_exe != "": #For hadoop mode - moves local file to HDFS
            moveToHDFS(bb_file,out_fname)
        os.remove(fname) # purge temporary bed file
        fname = "perm%s.bed"%(data_id)
        samp_out = open(fname,'w')
    line = "%s\t%d\t%d\t%d\n"%(ref_id,ref_off,ref_off+ref_len,hmm_st)
    samp_out.write(line)
    last_perm=data_id

#Can't forget about last file
samp_out.close()
out_fname =  "%s/perm%s.bb"%(args.out_dir,last_perm)
bb_file = "perm%s.bb" % (last_perm) if args.hadoop_exe != "" else out_fname
bigbed_proc = subprocess.Popen([args.bigbed_exe, fname, args.chrom_sizes, bb_file])
bigbed_proc.wait()
if args.hadoop_exe!="": #For hadoop mode - moves local file to HDFS
    moveToHDFS(bb_file,out_fname)
os.remove(fname) # purge temporary bed file

#Done
timeEn = time.clock()
print >>sys.stderr, "DONE with aggr_states.py; in = %d; time=%0.3f secs" % (ninp, timeEn-timeSt)
