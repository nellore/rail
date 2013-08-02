#!/usr/bin/env python

"""
align.py
(first step after preprocessing, before splice.py)

Alignment script for MapReduce pipelines.  Wraps Bowtie.  Has features for (1)
optionally extracting readlets (substrings) of configurable length, at a
configurable interval along the read, (2) optionally truncating reads or
omitting mates.  Each read or readlet is then written to the standard-in
filehandle of an open Bowtie process.  Output from the Bowtie process is parsed
and passed to the standard-out filehandle.  Alignments are in Bowtie format
(not SAM).

Tab-delimited input tuple columns; can be in any of 3 formats:
 Format 1 (unpaired):
  1. Name
  2. Nucleotide sequence
  3. Quality sequence
 Format 2 (paired, 5-column):
  1. Name
  2. Nucleotide sequence for mate 1
  3. Quality sequence for mate 1
  4. Nucleotide sequence for mate 2
  5. Quality sequence for mate 2
 Format 3 (paired, 6-column):
  1. Name for mate 1
  2. Nucleotide sequence for mate 1
  3. Quality sequence for mate 1
  4. Name for mate 2
  5. Nucleotide sequence for mate 2
  6. Quality sequence for mate 2

-Binning/sorting prior to this step:
 (none)

Exons:
Tab-delimited output tuple columns:
1. Partition ID for partition overlapped by interval
2. Interval start
3. Interval end (exclusive)
4. Reference ID
5. Sample label

Introns:
Tab-delimited output tuple columns:
1. Partition ID for partition overlapped by interval (includes strand information)
2. Interval start
3. Interval end (exclusive)
4. Reference ID
5. Sample label
6. Readlet Sequence on 5' site
6. Readlet Sequence on 3' site
"""

import sys
import os
import site
import argparse
import threading
import string
import time
import numpy

timeSt = time.clock()

base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
site.addsitedir(os.path.join(base_path, "bowtie"))
site.addsitedir(os.path.join(base_path, "read"))
site.addsitedir(os.path.join(base_path, "sample"))
site.addsitedir(os.path.join(base_path, "interval"))
site.addsitedir(os.path.join(base_path, "alignment"))
site.addsitedir(os.path.join(base_path, "fasta"))

import bowtie
import readlet
import sample
import interval
import partition
import eddist
import nw
import fasta

ninp = 0               # # lines input so far
nout = 0               # # lines output so far
pe = False
discardMate = None
lengths = dict()       # read lengths after truncation
rawLengths = dict()    # read legnths prior to truncation
qualCnts = dict()      # quality counts after adjustments
rawQualCnts = dict()   # quality counts before adjustments
qualAdd = None         # amt to add to qualities
truncateAmt = None     # amount to truncate reads
truncateTo = None      # amount to truncate reads

readletize = None      # if we're going to readletize,

xformReads = qualAdd is not None or truncateAmt is not None or truncateTo is not None

parser = argparse.ArgumentParser(description=\
    'Align reads using Bowtie, usually as the map step in a Hadoop program.')
parser.add_argument(\
    '--refseq', type=str, required=False,
    help='The fasta sequence of the reference genome. The fasta index of the '
         'reference genome is also required to be built via samtools')
parser.add_argument(\
    '--splice-overlap', type=int, default=10,
    help='The overlap length of spanning readlets when evaluating splice junctions')
parser.add_argument(\
    '--faidx', type=str, required=False, help='Fasta index file')

bowtie.addArgs(parser)
readlet.addArgs(parser)
partition.addArgs(parser)

parser.add_argument(\
    '--test', action='store_const', const=True, default=False, help='Run unit tests')
parser.add_argument(\
    '--profile', action='store_const', const=True, default=False, help='Profile the code')

# Collect the bowtie arguments first
argv = sys.argv
bowtieArgs = []
in_args = False
for i in xrange(1, len(sys.argv)):
    if in_args:
        bowtieArgs.append(sys.argv[i])
    if sys.argv[i] == '--':
        argv = sys.argv[:i]
        in_args = True

args = parser.parse_args(argv[1:])

def xformRead(seq, qual):
    # Possibly truncate and/or modify quality values
    # TODO: not implemented yet!
    newseq, newqual = "", ""
    if truncateAmt is not None:
        pass
    if truncateTo is not None:
        pass
    if qualAdd is not None:
        pass
    return newseq, newqual

_revcomp_trans = string.maketrans("ACGT", "TGCA")
def revcomp(s):
    return s[::-1].translate(_revcomp_trans)

bowtieOutDone = threading.Event()

"""
Applies Needleman Wunsch to correct splice junction gaps
"""
def correctSplice(read,ref_left,ref_right,fw):
    revread = revcomp(read)
    if not fw:
        ref_right = revcomp(ref_right)
        score1,leftDP  = nw.c_needlemanWunsch(ref_left, read, nw.matchCost)
        score2,rightDP = nw.c_needlemanWunsch(ref_right,revread, nw.matchCost)
        rightDP = numpy.fliplr(rightDP)
        rightDP = numpy.flipud(rightDP)
    else:
        ref_right = revcomp(ref_right)
        score1,leftDP  = nw.c_needlemanWunsch(ref_left, read, nw.matchCost)
        score2,rightDP = nw.c_needlemanWunsch(ref_right,revread, nw.matchCost)
        rightDP = numpy.fliplr(rightDP)
        rightDP = numpy.flipud(rightDP)

    total = leftDP+rightDP
    # print >> sys.stderr,"ref_left\t",ref_left
    # print >> sys.stderr,"ref_right\t",ref_right
    # print >> sys.stderr,"read    \t",read,"\n"

    # print >> sys.stderr,"Left Matrix\n",leftDP
    # print >> sys.stderr,"Right Matrix\n",rightDP
    # print >> sys.stderr,"Total Matrix\n",total
    index = numpy.argmax(total)
    max_  = numpy.max(total)
    n = len(read)+1
    r = index%n
    c = index/n
    return r,c,total[r,c]

def printIntrons(refid,rdseq,region_st,region_end,in_start,in_end,rdnm,fw):
    global nout
    offset = args.splice_overlap
    fw_char = "+" if fw else "-"
    if not fw:
        rdseq = revcomp(rdseq)

    left_st,left_end = region_st-offset,region_st
    right_st,right_end = region_end,region_end+offset

    left_flank = rdseq[left_st:left_end]
    left_overlap = rdseq[left_end:left_end+offset]
    right_overlap = rdseq[right_st-offset:right_st]
    right_flank = rdseq[right_st:right_end]


        # right_flank = rdseq[region_st-args.readletLen:region_st]
        # right_overlap = rdseq[region_st-args.readletLen-args.splice_overlap:region_st-args.readletLen]
        # left_overlap = rdseq[region_end+args.readletLen:region_end+args.readletLen+args.splice_overlap]
        # left_flank = rdseq[region_end:region_end+args.readletLen]

    # print >> sys.stderr,"Strand",fw_char
    # print >> sys.stderr,'left ',left_st,'\t',left_flank,'\t',left_end,'\t',left_overlap
    # print >> sys.stderr,'right',right_st,'\t',right_flank,'\t',right_end,'\t',right_overlap,'\n'

    for pt in iter(partition.partition(refid, in_start, in_end, binsz)):
        print "intron\t%s%s\t%012d\t%d\t%s\t%s\t%s\t%s\t%s\t%s" % (pt, fw_char, in_start, in_end, refid, sample.parseLab(rdnm),left_flank,left_overlap,right_flank,right_overlap)
        nout += 1


def handleIntron(k,in_start,in_end,rdseq,unmapped_st,unmapped_end,region_st,region_end,rdnm,fw,fnh):
    diff = unmapped_end-unmapped_st-1
    offset = args.splice_overlap  #offset of unmapped to region coordinate frames
    left_st, left_en = in_start - offset + 1, in_start - offset + diff + 1
    right_st, right_en = in_end + offset - diff, in_end + offset
    ref_left = fnh.fetch_sequence(k,left_st, left_en).upper()
    ref_right = fnh.fetch_sequence(k,right_st, right_en).upper()
    unmapped = rdseq[unmapped_st:unmapped_end]
    # print >> sys.stderr,"ref_left\t",ref_left
    # print >> sys.stderr,"ref_right\t",ref_right
    # print >> sys.stderr,"read    \t",unmapped,"\n"

    #print >> sys.stderr,"Forward Strand",fw
    #print >> sys.stderr,ref_left,len(ref_left),ref_right,len(ref_right),unmapped,len(unmapped)
    _, dj,_ = correctSplice(unmapped,ref_left,ref_right,fw)
    #print >> sys.stderr,"Before",region_st,region_end
    left_diff,right_diff = dj-offset, len(unmapped)-dj-offset
    #print >> sys.stderr,"Left Diff",left_diff,"Right Diff",right_diff
    region_st,region_end = region_st+left_diff,region_end-right_diff
    in_start,in_end = in_start+left_diff, in_end-right_diff
    #print >> sys.stderr,ref_left,ref_right,unmapped
    #left_diff,right_diff = (dj-region_st),(region_end-dj)
    # print >> sys.stderr,"Left Diff",left_diff,"Right Diff",right_diff
    #print >> sys.stderr,"After",region_st,region_end

    #region_st,region_end = region_st+left_diff,region_end-right_diff
    #in_start,in_end = in_start+left_diff,in_end-right_diff
    printIntrons(k,rdseq,region_st,region_end,in_start,in_end,rdnm,fw)

"""
Compares potential short intron with readlet
"""
def handleShortAlignment(k,in_start,in_end,rdseq,unmapped_st,unmapped_end,region_st,region_end,rdnm,fw,fnh):
    global nout
    refseq = fnh.fetch_sequence(k,in_start + 1, in_end + 1).upper() # Sequence from genome
    rdsubseq = rdseq[unmapped_st:unmapped_end]
    if not fw:
        rdsubseq = revcomp(rdsubseq)
    score,_ = nw.c_needlemanWunsch(refseq, rdsubseq, nw.inverseMatchCost)
    # TODO: redo this in terms of percent identity or some
    # other measure that adapts to length of the missing bit,
    # not just a raw score
    if score <= len(rdsubseq)/10:
        for pt in iter(partition.partition(k, in_start, in_end, binsz)):
            print "exon\t%s\t%012d\t%d\t%s\t%s" % (pt, in_start, in_end, k, sample.parseLab(rdnm))
            nout += 1
            #else:
            #handleIntron(k,in_start,in_end,rdseq,unmapped_st,unmapped_end,region_st,region_end,rdnm,fw,fnh)



def composeReadletAlignments(rdnm, rdals, rdseq):

    # TODO: We include strand info with introns, but not exons.  We might want
    # to include for both for the case where the RNA-seq protocol is stranded.

    global nout
    # Add this interval to the flattened interval collection for current read
    ivals = {}
    positions = dict()  #stores readlet number based keyed by position, strand and reference id
    for rdal in rdals:
        refid, fw, refoff0, seqlen, rlet_nm = rdal
        refoff0, seqlen, rlet_nm = int(refoff0), int(seqlen), int(rlet_nm)
        # Remember begin, end offsets for readlet w/r/t 5' end of the read
        #TODO: Consider stranded scenario.  These positions are always defined on 5' end
        if fw:
            positions[(refid, fw, refoff0)] = rlet_nm * args.readletIval
            positions[(refid, fw, refoff0 + seqlen)] = rlet_nm * args.readletIval + seqlen
        else:
            positions[(refid, fw, refoff0 + seqlen)] = rlet_nm * args.readletIval
            positions[(refid, fw, refoff0)] = rlet_nm * args.readletIval + seqlen


        if (refid, fw) not in ivals:
            ivals[(refid, fw)] = interval.FlatIntervals()
        ivals[(refid, fw)].add(interval.Interval(refoff0, refoff0 + seqlen))

    for kfw in ivals.iterkeys(): # for each chromosome covered by >= 1 readlet
        k, fw = kfw
        in_end, in_start = -1, -1
        for iv in sorted(iter(ivals[kfw])): # for each covered interval, left-to-right
            st, en = iv.start, iv.end
            assert en > st
            assert st >= 0 and en >= 0
            if in_end == -1 and in_start >= 0:
                in_end = st
            if in_start == -1:
                in_start = en
            if in_start >= 0 and in_end >= 0:
                region_st,region_end = min(positions[(k, fw, in_start)],positions[(k, fw, in_end)]), max(positions[(k, fw, in_start)],positions[(k, fw, in_end)])
                offset = args.splice_overlap
                unmapped_st,unmapped_end = region_st-offset,region_end+offset #Need to add on sequence to get more alignment info
                reflen,rdlet_len = in_end-in_start,region_end-region_st
                #print >> sys.stderr,"reflen",reflen,"unmappedlen",unmappedlen,"unmapped start",unmapped_st,"unmapped end",unmapped_end
                assert in_start<in_end
                # if in_start>in_end:
                #     print >> sys.stderr,"This should never happen!!!","ref_len",reflen,"<0"
                #     print >> sys.stderr,"In_start",in_start,"In_end",in_end
                if rdlet_len==0:
                    handleIntron(k,in_start,in_end,rdseq,unmapped_st,unmapped_end,region_st,region_end,rdnm,fw,fnh)
                elif abs(reflen-rdlet_len)/float(rdlet_len) < 0.05:
                    #Note: just a readlet missing due to error or variant
                    handleShortAlignment(k,in_start,in_end,rdseq,unmapped_st,unmapped_end,region_st,region_end,rdnm,fw,fnh)
                elif reflen > rdlet_len:
                    #print >> sys.stderr,"len>0","Region start",region_st,"Region end",region_end
                    handleIntron(k,in_start,in_end,rdseq,unmapped_st,unmapped_end,region_st,region_end,rdnm,fw,fnh)
                else:
                     print >> sys.stderr,"This should never happen!!!","ref_len",reflen,"<","rdlet_len",rdlet_len
                     print >> sys.stderr,"In_start",in_start,"In_end",in_end

                in_start, in_end = en, -1
            # Keep stringing rdid along because it contains the label string
            # Add a partition id that combines the ref id and some function of
            # the offsets
            for pt in iter(partition.partition(k, st, en, binsz)):
                print "exon\t%s\t%012d\t%d\t%s\t%s" % (pt, st, en, k, sample.parseLab(rdnm))
                nout += 1

def bowtieOutReadlets(st):
    ''' Process standard out (stdout) output from Bowtie.  Each line of output
        is another readlet alignment.  We *could* try to perform operations
        over all readlet alignments from the same read here.  Currently, we
        follow this step with an aggregation that bins by read id, then
        operate over bins in splice.py. '''
    global nout
    mem, cnt = {}, {}
    for line in st:
        if line[0] == '@':
            continue
        rdid, flags, refid, refoff1, _, _, _, _, _, seq, _, _ = string.split(line.rstrip(), '\t', 11)
        flags, refoff1 = int(flags), int(refoff1)
        seqlen = len(seq)
        toks = string.split(rdid, ';')
        rdnm = ';'.join(toks[:-3])
        rlet_nm = toks[2]
        cnt[rdnm] = cnt.get(rdnm, 0) + 1
        rdseq = toks[4]
        rdlet_n = int(toks[-2])
        if flags != 4:
            fw = (flags & 16) == 0
            if rdnm not in mem: mem[rdnm] = [ ]
            mem[rdnm].append((refid, fw, refoff1-1, seqlen,rlet_nm))
        if cnt[rdnm] == rdlet_n:
            # Last readlet
            if rdnm in mem:
                composeReadletAlignments(rdnm, mem[rdnm],rdseq)
                del mem[rdnm]
            del cnt[rdnm]
        nout += 1
    assert len(mem) == 0
    assert len(cnt) == 0
    bowtieOutDone.set()

def bowtieOut(st):
    ''' Process standard out (stdout) output from Bowtie.  Each line of output
        is another readlet alignment.  We *could* try to perform operations
        over all readlet alignments from the same read here.  Currently, we
        follow this step with an aggregation that bins by read id, then
        operate over bins in splice.py. '''
    global nout
    for line in st:
        sys.stdout.write(line)
        nout += 1
    bowtieOutDone.set()

def go():
    global ninp
    first = True
    for ln in sys.stdin:
        ln = ln.rstrip()
        toks = ln.split('\t')
        ninp += 1
        pair = False
        nm, seq, qual = None, None, None
        nm1, seq1, qual1 = None, None, None
        nm2, seq2, qual2 = None, None, None
        if len(toks) == 3:
            # Unpaired read
            nm, seq, qual = toks
            sample.hasLab(nm, mustHave=True) # check that label is present in name
        elif len(toks) == 5 or len(toks) == 6:
            # Paired-end read
            if len(toks) == 5:
                # 6-token version
                nm1, seq1, qual1, seq2, qual2 = toks
                nm2 = nm1
            else:
                # 5-token version
                nm1, seq1, qual1, nm2, seq2, qual2 = toks
            sample.hasLab(nm1, mustHave=True) # check that label is present in name
            if discardMate is not None:
                # We have been asked to discard one mate or the other
                if discardMate == 1:
                    nm, seq, qual = nm2, seq2, qual2 # discard mate 1
                else:
                    nm, seq, qual = nm1, seq1, qual1 # discard mate 2
            else:
                pair = True # paired-end read
        else:
            raise RuntimeError("Wrong number of tokens for line: " + ln)
        if pair:
            # Paired-end
            if xformReads:
                # Truncate and transform quality values
                seq1, qual1 = xformRead(seq1, qual1)
                seq2, qual2 = xformRead(seq2, qual2)
            if args.readletLen > 0:
                # Readletize
                rlets1 = readlet.readletize(args, nm1, seq1, qual1)
                for i in xrange(0, len(rlets1)):
                    nm_rlet, seq_rlet, qual_rlet = rlets1[i]
                    rdletStr = "%s;%d;%d;%s\t%s\t%s\n" % (nm_rlet, i, len(rlets1),seq1, seq_rlet, qual_rlet)
                    if first:
                        sys.stderr.write("First readlet: '%s'" % rdletStr)
                        first = False
                    proc.stdin.write(rdletStr)
                rlets2 = readlet.readletize(args, nm2, seq2, qual2)
                for i in xrange(0, len(rlets2)):
                    nm_rlet, seq_rlet, qual_rlet = rlets2[i]
                    rdletStr = "%s;%d;%d;%s\t%s\t%s\n" % (nm_rlet, i, len(rlets2),seq2, seq_rlet, qual_rlet)
                    if first:
                        sys.stderr.write("First readlet: '%s'" % rdletStr)
                        first = False
                    proc.stdin.write(rdletStr)
            else:
                rdStr = "%s\t%s\t%s\t%s\t%s\n" % (nm1, seq1, qual1, seq2, qual2)
                if first:
                    sys.stderr.write("First read: '%s'" % rdStr)
                    first = False
                proc.stdin.write(rdStr)
        else:
            # Unpaired
            if xformReads:
                # Truncate and transform quality values
                seq, qual = xformRead(seq, qual)
            if args.readletLen > 0:
                # Readletize
                rlets = readlet.readletize(args, nm, seq, qual)
                for i in xrange(0, len(rlets)):
                    nm_rlet, seq_rlet, qual_rlet = rlets[i]
                    rdletStr = "%s;%d;%d;%s\t%s\t%s\n" % (nm_rlet, i, len(rlets), seq, seq_rlet, qual_rlet)
                    if first:
                        sys.stderr.write("First readlet: '%s'" % rdletStr)
                        first = False
                    proc.stdin.write(rdletStr)
            else:
                rdStr = "%s\t%s\t%s\n" % (nm, seq, qual)
                if first:
                    sys.stderr.write("First read: '%s'" % rdStr)
                    first = False
                proc.stdin.write(rdStr)

    # Close and flush STDIN.
    proc.stdin.close()

    # Wait until the threads processing stdout and stderr are done.

    # Close stdout and stderr
    print >>sys.stderr, "Waiting for Bowtie stdout processing thread to finish"
    bowtieOutDone.wait()

    proc.stdout.close()

    # Done
    timeEn = time.clock()
    print >>sys.stderr, "DONE with align.py; in/out = %d/%d; time=%0.3f secs" % (ninp, nout, timeEn-timeSt)


def createTestFasta(fname,refid,refseq):
    fastaH = open(fname,'w')
    fastaIdx = open(fname+".fai",'w')
    fastaH.write(">%s\n%s\n"%(refid,refseq))
    fastaIdx.write("%s\t%d\t%d\t%d\t%d\n"%(refid,len(refseq),len(refid)+2,len(refseq),len(refseq)+1))
    fastaH.close()
    fastaIdx.close()

def test_fasta_create():
    rdseq = "ACGTACGT"
    refseq = "ACGTCCCCACGT"
    fname,refid = "test.fa","test"
    createTestFasta(fname,refid,refseq)
    fnh = fasta.fasta(fname)
    testseq = fnh.fetch_sequence(refid,1,len(refseq))
    assert testseq==refseq
    print >> sys.stderr,"Test Fasta Create Success!"
    os.remove(fname)
    os.remove(fname+".fai")


def test_short_alignment1():
    sys.stdout = open("test.out",'w')
    rdnm,fw = "0;LB:test",True
    #Mapped readlets: ACGT, ACGT
    rdseq,refseq,fname,refid = "GCACGTACGTCG","GCACGTCCCCCCCCCCCACGTCG","test.fa","test"
    createTestFasta(fname,refid,refseq)
    fnh = fasta.fasta(fname)
    region_st,region_end=6,6
    in_start,in_end=6,17
    #unmapped_st,unmapped_end = region_st-args.readletLen,region_end+args.readletLen
    offset = args.splice_overlap
    unmapped_st,unmapped_end = region_st-offset,region_end+offset
    printIntrons(refid,rdseq,region_st,region_end,in_start,in_end,rdnm,fw)
    handleIntron(refid,in_start,in_end,rdseq,unmapped_st,unmapped_end,region_st,region_end,rdnm,fw,fnh)
    sys.stdout.close()
    test_out = open("test.out",'r')
    line = test_out.readline().rstrip()
    testline = test_out.readline().rstrip()
    print >> sys.stderr, rdseq
    print >> sys.stderr, line,'\n',testline
    assert testline==line
    print >> sys.stderr,"Test Short Intron 1 Success!"
    os.remove(fname)
    os.remove(fname+".fai")
    os.remove("test.out")


def test_short_alignment2():
    sys.stdout = open("test.out",'w')
    rdnm,fw = "0;LB:test",False
    rdseq,refseq,fname,refid = "GCACGTACGTGC","GCACGTCCCCCCCCCCCACGTGC","test.fa","test"
    createTestFasta(fname,refid,refseq)
    fnh = fasta.fasta(fname)
    region_st,region_end=6,6
    in_start,in_end=6,17
    offset = args.splice_overlap
    unmapped_st,unmapped_end = region_st-offset,region_end+offset
    printIntrons(refid,rdseq,region_st,region_end,in_start,in_end,rdnm,fw)
    handleIntron(refid,in_start,in_end,rdseq,unmapped_st,unmapped_end,region_st,region_end,rdnm,fw,fnh)
    sys.stdout.close()
    test_out = open("test.out",'r')
    line = test_out.readline().rstrip()
    testline = test_out.readline().rstrip()
    print >> sys.stderr, rdseq
    print >> sys.stderr, line,'\n',testline
    assert testline==line
    print >> sys.stderr,"Test Short Intron 2 Success!"
    os.remove(fname)
    os.remove(fname+".fai")
    os.remove("test.out")


def test_short_alignment3():
    sys.stdout = open("test.out",'w')
    rdnm,fw = "0;LB:test",True
    #mapped reads: CGTA, TACG
    rdseq,refseq,fname,refid = "GCACGTACGTCG","GCACGTCCCCCCCCCCCACGTCG","test.fa","test"
    createTestFasta(fname,refid,refseq)
    fnh = fasta.fasta(fname)
    region_st,region_end=7,5
    in_start,in_end=7,16
    offset = args.splice_overlap
    unmapped_st,unmapped_end = region_st-offset,region_end+offset
    printIntrons(refid,rdseq,region_st-1,region_end+1,in_start-1,in_end+1,rdnm,fw)
    handleIntron(refid,in_start,in_end,rdseq,unmapped_st,unmapped_end,region_st,region_end,rdnm,fw,fnh)
    sys.stdout.close()
    test_out = open("test.out",'r')
    line = test_out.readline().rstrip()
    testline = test_out.readline().rstrip()
    print >> sys.stderr, rdseq
    print >> sys.stderr, line,'\n',testline
    assert testline==line
    print >> sys.stderr,"Test Short Intron 3 Success!"
    os.remove(fname)
    os.remove(fname+".fai")
    os.remove("test.out")


def test_correct_splice():
    left = "TTACGAAGGTTTGTA"
    right= "TAATTTAGATGGAGA"
    read = "TTACGAAGATGGAGA"
    # left = "AGTATCGAACCTGAAGCAAGTTACGAAGGTTTGTATAACAAAAATTATGTGAAAG"
    # right= "TAATATTTTCTTTTGAAATTTAATTTAGATGGAGAAATGGAAGCAGAGTGGCTAG"
    # read = "AGTATCGAACCTGAAGCAAGTTACGAAGATGGAGAAATGGAAGCAGAGTGGCTAG"
    fw = True
    r,c,score = correctSplice(read,left,right,fw)
    #print >> sys.stderr,read
    #print >> sys.stderr,left[:c],right[c:]
    assert left[:c]+right[c:] == read
    print >> sys.stderr,"Correct Splice Test Successful!!!"
def test():
    test_fasta_create()
    test_short_alignment1()
    test_short_alignment2()
    test_short_alignment3()
    test_correct_splice()

if args.test:
    binsz = 10000
    test()
elif args.profile:
    import cProfile
    cProfile.run('go()')
else:
    binsz = partition.binSize(args)
    if not os.path.exists(args.refseq):
        raise RuntimeError("No such --refseq file: '%s'" % args.refseq)
    if not os.path.exists(args.faidx):
        raise RuntimeError("No such --faidx file: '%s'" % args.faidx)

    fnh = fasta.fasta(args.refseq)

    proc = bowtie.proc(args, bowtieArgs=bowtieArgs, sam=True, outHandler=bowtieOutReadlets)
    go()
