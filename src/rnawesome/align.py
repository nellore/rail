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
7. Readlet Sequence on 3' site
8. Read name
"""

import sys
import os
import site
import argparse
import threading
import string
import numpy

base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
site.addsitedir(os.path.join(base_path, "bowtie"))
site.addsitedir(os.path.join(base_path, "read"))
site.addsitedir(os.path.join(base_path, "sample"))
site.addsitedir(os.path.join(base_path, "interval"))
site.addsitedir(os.path.join(base_path, "alignment"))
site.addsitedir(os.path.join(base_path, "fasta"))
site.addsitedir(os.path.join(base_path, "statsmath"))

import bowtie
import readlet
import sample
import interval
import partition
import needlemanWunsch
import fasta
import window

ninp = 0               # # lines input so far
nout = 0               # # lines output so far
pe = False             # What is this variable?
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
    '--splice-overlap', type=int, default=5,
    help='The overlap length of spanning readlets when evaluating splice junctions')
parser.add_argument(\
    '--faidx', type=str, required=False, help='Fasta index file')

bowtie.addArgs(parser)
readlet.addArgs(parser)
partition.addArgs(parser)

parser.add_argument(\
    '--serial', action='store_const', const=True, default=False, help="Run bowtie serially after rather than concurrently with the input-reading loop")
parser.add_argument(\
    '--keep-reads', action='store_const', const=True, default=False, help="Don't delete any temporary read file(s) created")
parser.add_argument(\
    '--write-reads', type=str, required=False, help='Write input reads to given tab-delimited file')
parser.add_argument(\
    '--test', action='store_const', const=True, default=False, help='Run unit tests')
parser.add_argument(\
    '--archive', metavar="PATH", type=str, help='Save input and command to a subdirectory (named using this process\'s PID) of PATH')
parser.add_argument(\
    '--profile', action='store_const', const=True, default=False, help='Profile the code')
parser.add_argument(\
    '--verbose', action='store_const', const=True, default=False, help='Prints out extra debugging statements')


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
Breaks ties in the Needleman-Wunsch algorithm by selecting the middle element in a set of maxes

e.g.  The 3rd row will be chosen
matrix = [ 1 2 3 4 5
           6 8 6 4 1
           1 2 8 4 5
           6 7 6 8 1
           1 2 3 4 5]
"""
def medianTieBreaker(dpmat,m_ind):
    m_elem = numpy.max(dpmat[:,m_ind])
    st = m_ind
    r,c = dpmat.shape
    while m_ind>=0 and m_ind<c and numpy.max(dpmat[:,m_ind])==m_elem:
        m_ind+=1
    end = m_ind
    return (st+end)/2

"""
Uses the sliding window algorithm to help place the flanking sequences close to the junction sites
"""
def windowTransform(read,site,cost):
    hist = [1]*( len(read) )
    scores = window.score(read, site, hist, cost)
    return scores

"""
Applies Needleman Wunsch to correct splice junction gaps
"""
def correctSplice(read,ref_left,ref_right,fw):
    revread = revcomp(read)
    """Needleman-Wunsch is a directional algorithm.  Since we are interested in scoring the 3' end of the right ref sequence,    we reverse complement the right ref sequence before applying the NW algorithm"""
    ref_right = revcomp(ref_right)
    score1,leftDP  = needlemanWunsch.needlemanWunsch(ref_left, read, needlemanWunsch.matchCost())
    score2,rightDP = needlemanWunsch.needlemanWunsch(ref_right,revread, needlemanWunsch.matchCost())

    #Once NW is applied, the right DP matrix must be transformed in the same coordinate frame as the left DP matrix
    rightDP = numpy.fliplr(rightDP)
    rightDP = numpy.flipud(rightDP)

    total = leftDP+rightDP
    index = numpy.argmax(total)
    max_  = numpy.max(total)

    n = len(read)+1
    r,c = index%n, index/n

    c = medianTieBreaker(total,c)
    r = numpy.argmax(total[:,c])
    return r,c,total[r,c],leftDP,rightDP,total

#Print all listed exons to stdout
def printExons(refid,in_start,in_end,rdnm):
    global nout
    for pt in iter(partition.partition(refid, in_start, in_end, binsz)):
        print "exon\t%s\t%012d\t%d\t%s\t%s" % (pt, in_start, in_end, refid, sample.parseLab(rdnm))
        nout += 1

"""
Returns a more human readable format of the string
"""
def readableFormat(s):
    return " ".join([s[i:i+10] for i in range(0,len(s),10)])
def formatList(s,l):
    return (" "*l).join( list( str(s) ) )

#Print all listed introns to stdout and the flanking sequences
def printIntrons(refid,rdseq,region_st,region_end,in_start,in_end,rdnm,fw,rdid,outhandle):
    global nout
    offset = args.splice_overlap
    fw_char = "+" if fw else "-"
    #Obtain coordinates for flanking coordinate frames
    left_st,left_end = region_st-offset,region_st
    right_st,right_end = region_end,region_end+offset

    flank1 = rdseq[left_st:left_end]
    overlap1 = rdseq[left_end:left_end+offset]
    overlap2 = rdseq[right_st-offset:right_st]
    flank2 = rdseq[right_st:right_end]
    if not fw:
        left_overlap = revcomp(flank1)
        left_flank = revcomp(overlap1)
        right_overlap = revcomp(flank2)
        right_flank = revcomp(overlap2)
        rdseq = revcomp(rdseq)
    else:
        left_flank = flank1
        left_overlap = overlap1
        right_flank = flank2
        right_overlap = overlap2
    """Since there is a possibility that one of the sequences may be out of boundaries (e.g. mapping error),
    the following checks to see if all of the sequences are valid"""
    if ( len(left_flank) == len(right_flank) and
         len(left_overlap) == len(right_overlap) and
         len(left_flank) == len(left_overlap)):
         for pt in iter(partition.partition(refid, in_start, in_end, binsz)):
            print >> outhandle,"intron\t%s%s\t%012d\t%d\t%s\t%s\t%s\t%s\t%s" % (pt, fw_char, in_start, in_end, refid, sample.parseLab(rdnm),left_flank,left_overlap,rdid)

            nout += 1

def handleIntron(k,in_start,in_end,rdseq,unmapped_st,unmapped_end,region_st,region_end,rdnm,fw,fnh,offset,rdid):
    diff = unmapped_end-unmapped_st-1
    #Obtain coordinates of flanking sequences surrounding splice site
    left_st,right_end = in_start-offset+1,in_end+offset
    left_end,right_st = left_st+diff,right_end-diff
    #Print directly to stdout if flanking sequences overlap too much
    if left_end<=left_st or right_end<=right_st:
        printIntrons(k,rdseq,region_st,region_end,in_start,in_end,rdnm,fw,rdid,sys.stdout)
    else:
        ref_left = fnh.fetch_sequence(k,left_st, left_end).upper()
        ref_right = fnh.fetch_sequence(k,right_st, right_end).upper()
        if not fw:
            readseq = revcomp(rdseq)
            unmapped = revcomp(rdseq[unmapped_st:unmapped_end])
        else:
            readseq = rdseq
            unmapped = rdseq[unmapped_st:unmapped_end]

        _, dj,score,leftDP,rightDP,total = correctSplice(unmapped,ref_left,ref_right,fw)
        left_diff,right_diff = dj, len(unmapped)-dj
        left_in_diff,right_in_diff = left_diff-offset,right_diff-offset
        if score>0:   #If crappy alignment, disregard corrections
            region_st,region_end = unmapped_st+left_diff,unmapped_end-right_diff
            in_start,in_end = in_start+left_in_diff,in_end-right_in_diff
        printIntrons(k,rdseq,region_st,region_end,in_start,in_end,rdnm,fw,rdid,sys.stdout)


"""
Compares potential short intron with readlet to check if it should be an exon instead
"""
def handleShortAlignment(k,in_start,in_end,rdseq,unmapped_st,unmapped_end,region_st,region_end,rdnm,fw,fnh):
    refseq = fnh.fetch_sequence(k,in_start + 1, in_end + 1).upper() # Sequence from genome
    rdsubseq = rdseq[unmapped_st:unmapped_end]
    if not fw:
        rdsubseq = revcomp(rdsubseq)
    score = needlemanWunsch.needlemanWunsch(refseq, rdsubseq, needlemanWunsch.matchCost())
    # TODO: redo this in terms of percent identity or some
    # other measure that adapts to length of the missing bit,
    # not just a raw score
    if score >= len(rdsubseq)*(9.0/10):
        printExons(k,in_start,in_end,rdnm)

def getIntervals(rdals,L):
    ivals = {}
    positions = dict()  #stores readlet number based keyed by position, strand and reference id
    for rdal in rdals:
        refid, fw, refoff0, seqlen, rlet_nm, rdid = rdal
        refoff0, seqlen, rlet_nm = int(refoff0), int(seqlen), int(rlet_nm)
        # Remember begin, end offsets for readlet w/r/t 5' end of the read
        if fw:
            positions[(refid, fw, refoff0)] = rlet_nm * args.readletIval
            positions[(refid, fw, refoff0 + seqlen)] = rlet_nm * args.readletIval + seqlen
        else:
            positions[(refid, fw, refoff0)] = rlet_nm * args.readletIval + seqlen
            positions[(refid, fw, refoff0 + seqlen)] = rlet_nm * args.readletIval
        if (refid, fw, rdid) not in ivals:
            ivals[(refid, fw, rdid)] = interval.FlatIntervals()
        ivals[(refid, fw, rdid)].add(interval.Interval(refoff0, refoff0 + seqlen))
    return ivals,positions

test_id = set(["r_n1688","r_n2805","r_n4526","r_n4833","r_n2020","r_n3512","r_n4446","r_n279","r_n1828","r_n5035","r_n1848","r_n2403","r_n3163","r_n3944","r_n4211"])
test_exons = open("test_exons.txt",'w')

def composeReadletAlignments(rdnm, rdals, rdseq):

    global nout
    L=len(rdseq)-1
    # Add this interval to the flattened interval collection for current read
    ivals,positions = getIntervals(rdals,L)
    for kfw in ivals.iterkeys(): # for each chromosome covered by >= 1 readlet
        k, fw, rdid = kfw
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
                if fw:
                    region_st,region_end = positions[(k, fw, in_start)],positions[(k, fw, in_end)]
                else:
                    region_end,region_st = positions[(k, fw, in_start)],positions[(k, fw, in_end)]

                offset = args.splice_overlap
                unmapped_st,unmapped_end = region_st-offset,region_end+offset
                reflen,rdlet_len = in_end-in_start, abs(region_end-region_st)

                assert in_start<in_end
                if abs(reflen-rdlet_len)/float(rdlet_len+1) < 0.05:
                    #Note: just a readlet missing due to sequencing error or variant
                    handleShortAlignment(k,in_start,in_end,rdseq,unmapped_st,unmapped_end,region_st,region_end,rdnm,fw,fnh)
                elif rdlet_len>reflen:
                    printExons(k,in_start,in_end,rdnm)
                else:
                    handleIntron(k,in_start,in_end,rdseq,unmapped_st,unmapped_end,region_st,region_end,rdnm,fw,fnh,offset,rd
                in_start, in_end = en, -1
            # Keep stringing rdid along because it contains the label string
            # Add a partition id that combines the ref id and some function of
            # the offsets
            for pt in iter(partition.partition(k, st, en, binsz)):
                print "exon\t%s\t%012d\t%d\t%s\t%s" % (pt, st, en, k, sample.parseLab(rdnm))
                nout += 1

def bowtieOutReadlets(st, reportMult=1.2):
    ''' Process standard out (stdout) output from Bowtie.  Each line of output
        is another readlet alignment.  We *could* try to perform operations
        over all readlet alignments from the same read here.  Currently, we
        follow this step with an aggregation that bins by read id, then
        operate over bins in splice.py. '''
    global nout
    mem, cnt = {}, {}
    report = 1
    line_nm = 0
    try:
        while True:
            line = st.readline()
            if len(line) == 0:
                break # no more output
            if line[0] == '@':
                continue # skip header
            nout += 1
            rdid, flags, refid, refoff1, _, _, _, _, _, seq, _, _ = string.split(line.rstrip(), '\t', 11)
            flags, refoff1 = int(flags), int(refoff1)
            if nout >= report:
                report *= reportMult
                if args.verbose:
                    print >> sys.stderr, "SAM output record %d: rdname='%s', flags=%d" % (nout, rdid, flags)
            seqlen = len(seq)
            toks = string.split(rdid, ';')
            rdnm = ';'.join(toks[:-3])
            rlet_nm = toks[2]
            cnt[rdnm] = cnt.get(rdnm, 0) + 1
            rd_name = toks[0]
            rdseq = toks[4]
            rdlet_n = int(toks[-2])
            if flags != 4:
                fw = (flags & 16) == 0
                if rdnm not in mem: mem[rdnm] = [ ]
                mem[rdnm].append((refid, fw, refoff1-1, seqlen,rlet_nm,rd_name))
            if cnt[rdnm] == rdlet_n:
                # Last readlet
                if rdnm in mem:
                    composeReadletAlignments(rdnm, mem[rdnm],rdseq)
                    del mem[rdnm]
                del cnt[rdnm]
            line_nm+=1
    except IOError as e:
        print >> sys.stderr, "I/O error while reading output from Bowtie ({0}): {1}".format(e.errno, e.strerror)
        sys.exit(20)
    except ValueError as e:
        print >> sys.stderr, "Value error while reading output from Bowtie: " + str(e)
        sys.exit(30)
    except TypeError as e:
        print >> sys.stderr, "Type error while reading output from Bowtie: " + str(e)
        raise
        sys.exit(35)
    except:
        print >> sys.stderr, "Unexpected error while reading output from Bowtie:%s" % (sys.exc_info()[0])
        raise
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

def writeReads(fhs, reportMult=1.2):
    """ Parse input reads, optionally transform them and/or turn them into
        readlets. """
    global ninp
    report = 1.0
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
                # 5-token version
                nm1, seq1, qual1, seq2, qual2 = toks
                nm2 = nm1
            else:
                # 6-token version
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
                    if ninp >= report and i == 0:
                        report *= reportMult
                        if args.verbose:
                            print >> sys.stderr, "First readlet from read %d: '%s'" % (ninp, rdletStr.rstrip())
                    for fh in fhs: fh.write(rdletStr)
                rlets2 = readlet.readletize(args, nm2, seq2, qual2)
                for i in xrange(0, len(rlets2)):
                    nm_rlet, seq_rlet, qual_rlet = rlets2[i]
                    rdletStr = "%s;%d;%d;%s\t%s\t%s\n" % (nm_rlet, i, len(rlets2),seq2, seq_rlet, qual_rlet)
                    if ninp >= report and i == 0:
                        report *= reportMult
                        if args.verbose:
                            print >> sys.stderr, "First readlet from read %d: '%s'" % (ninp, rdletStr.rstrip())
                    for fh in fhs: fh.write(rdletStr)
            else:
                rdStr = "%s\t%s\t%s\t%s\t%s\n" % (nm1, seq1, qual1, seq2, qual2)
                if ninp >= report:
                    report *= reportMult
                    if args.verbose:
                        print >> sys.stderr, "Read %d: '%s'" % (ninp, rdStr.rstrip())
                for fh in fhs: fh.write(rdStr)
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
                    if ninp >= report and i == 0:
                        report *= reportMult
                        if args.verbose:
                            print >> sys.stderr, "First readlet from read %d: '%s'" % (ninp, rdletStr.rstrip())
                    for fh in fhs: fh.write(rdletStr)
            else:
                rdStr = "%s\t%s\t%s\n" % (nm, seq, qual)
                if ninp >= report:
                    report *= reportMult
                    if args.verbose:
                        print >> sys.stderr, "Read %d: '%s'" % (ninp, rdStr.rstrip())
                for fh in fhs: fh.write(rdStr)

def go():

    import time
    timeSt = time.clock()

    archiveFh, archiveDir = None, None
    if args.archive is not None:
        archiveDir = os.path.join(args.archive, str(os.getpid()))
        if args.verbose:
            print >> sys.stderr, "Putting --archive reads and command in '%s'" % archiveDir
        if not os.path.exists(archiveDir):
            os.makedirs(archiveDir)
        archiveFh = open(os.path.join(archiveDir, "reads.tab5"), 'w')

    if args.serial:
        # Reads are written to a file, then Bowtie reads them from the file
        import tempfile
        if args.write_reads is None:
            tmpdir = tempfile.mkdtemp()
            readFn = os.path.join(tmpdir, 'reads.tab5')
        else:
            readFn = args.write_reads
        with open(readFn, 'w') as fh:
            fhs = [fh]
            if args.archive is not None: fhs.append(archiveFh)
            writeReads(fhs)
        assert os.path.exists(readFn)
        proc, mycmd = bowtie.proc(args, readFn=readFn, bowtieArgs=bowtieArgs, sam=True, outHandler=bowtieOutReadlets, stdinPipe=False)
    else:
        # Reads are written to Bowtie process's stdin directly
        proc, mycmd = bowtie.proc(args, readFn=None, bowtieArgs=bowtieArgs, sam=True, outHandler=bowtieOutReadlets, stdinPipe=True)
        fhs = [proc.stdin]
        if args.archive is not None: fhs.append(archiveFh)
        writeReads(fhs)
        proc.stdin.close()

    if args.archive is not None:
        archiveFh.close()
        with open(os.path.join(archiveDir, "bowtie_cmd.sh"), 'w') as ofh:
            ofh.write(mycmd + '\n')
        with open(os.path.join(archiveDir, "align_py_cmd.sh"), 'w') as ofh:
            ofh.write(' '.join(sys.argv) + '\n')
    if args.verbose:
        print >>sys.stderr, "Waiting for Bowtie to finish"
    bowtieOutDone.wait()
    proc.stdout.close()
    if args.verbose:
        print >>sys.stderr, "Bowtie finished"

    # Remove any temporary reads files created
    if args.serial and args.write_reads is None and not args.keep_reads:
        print >>sys.stderr, "Cleaning up temporary files"
        import shutil
        shutil.rmtree(tmpdir)

    timeEn = time.clock()
    print >>sys.stderr, "DONE with align.py; in/out = %d/%d; time=%0.3f secs" % (ninp, nout, timeEn-timeSt)

#Only used for testing
def createTestFasta(fname,refid,refseq):
    fastaH = open(fname,'w')
    fastaIdx = open(fname+".fai",'w')
    fastaH.write(">%s\n%s\n"%(refid,refseq))
    fastaIdx.write("%s\t%d\t%d\t%d\t%d\n"%(refid,len(refseq),len(refid)+2,len(refseq),len(refseq)+1))
    fastaH.close()
    fastaIdx.close()


if not args.test and not args.profile:
    binsz = partition.binSize(args)
    if not os.path.exists(args.refseq):
        raise RuntimeError("No such --refseq file: '%s'" % args.refseq)
    if not os.path.exists(args.faidx):
        raise RuntimeError("No such --faidx file: '%s'" % args.faidx)
    fnh = fasta.fasta(args.refseq)
    go()
elif args.profile:
    binsz = partition.binSize(args)
    if not os.path.exists(args.refseq):
        raise RuntimeError("No such --refseq file: '%s'" % args.refseq)
    if not os.path.exists(args.faidx):
        raise RuntimeError("No such --faidx file: '%s'" % args.faidx)
    import cProfile
    cProfile.run('go()')
else:
    del sys.argv[1:]
    import unittest
    binsz = 10000
    #test()

    class TestAlignFunctions1(unittest.TestCase):
        def setUp(self):
            #A visual representation of the reference sequence and the read
            """ACGAAGGACT GCTTGACATC GGCCACGATA AC                                                                                                                 AACCT TTTTTGCGCC AATCTTAAGA GCCTTCT"""
            """ACGAAGGACT GCTTGACATC GGCCACGATA ACCTGAGTCG ATAGGACGAA ACAAGTATAT ATTCGAAAAT TAATTAATTC CGAAATTTCA ATTTCATCCG ACATGTATCT ACATATGCCA CACTTCTGGT TGGACAACCT TTTTTGCGCC A"""
            self.rdseq  = "ACGAAGGACTGCTTGACATCGGCCACGATAACAACCTTTTTTGCGCCAATCTTAAGAGCCTTCT"
            self.refseq = "ACGAAGGACTGCTTGACATCGGCCACGATAACCTGAGTCGATAGGACGAAACAAGTATATATTCGAAAATTAATTAATTCCGAAATTTCAATTTCATCCGACATGTATCTACATATGCCACACTTCTGGTTGGACAACCTTTTTTGCGCCA"

            self.fasta = "test.fa"
            self.faidx = "test.fa.fai"
            createTestFasta(self.fasta,"test",self.refseq)

        def tearDown(self):
            os.remove(self.fasta)
            os.remove(self.faidx)

        def test_short_alignment1(self):

            sys.stdout = open("test.out",'w')
            rdnm,fw,refid = "0;LB:test",True,"test"
            #mapped reads:
            #ref st,end = (0,33),(135,166)
            #region_st,region_end=32,33
            #in_start,in_end=31,134
            fnh = fasta.fasta(self.fasta)
            region_st,region_end=32,33
            in_start,in_end=31,132
            offset = 5
            rdid = "testid"
            unmapped_st,unmapped_end = region_st-offset,region_end+offset
            printIntrons(refid,self.rdseq,region_st,region_end,31,132,rdnm,fw,rdid,sys.stdout)
            handleIntron(refid,in_start,in_end,self.rdseq,unmapped_st,unmapped_end,region_st,region_end,rdnm,fw,fnh,offset,"testid")
            sys.stdout.close()
            test_out = open("test.out",'r')
            line = test_out.readline().rstrip()
            testline = test_out.readline().rstrip()
            print >> sys.stderr, readableFormat(self.rdseq)
            print >> sys.stderr, line,'\n',testline
            self.assertEquals(line,testline)




        def test_correct_splice1(self):
            left = "TTACGAAGGTTTGTA"
            right= "TAATTTAGATGGAGA"
            read = "TTACGAAGATGGAGA"
            # left = "AGTATCGAACCTGAAGCAAGTTACGAAGGTTTGTATAACAAAAATTATGTGAAAG"
            # right= "TAATATTTTCTTTTGAAATTTAATTTAGATGGAGAAATGGAAGCAGAGTGGCTAG"
            # read = "AGTATCGAACCTGAAGCAAGTTACGAAGATGGAGAAATGGAAGCAGAGTGGCTAG"
            fw = True
            r,c,score,_,_,_ = correctSplice(read,left,right,fw)
            #print >> sys.stderr,read
            #print >> sys.stderr,left[:c],right[c:]
            assert left[:c]+right[c:] == read
            print >> sys.stderr,"Correct Splice Test Successful!!!"

        def test_correct_splice2(self):

            read = "ACGATAACCTTTTTT"
            left = "ACGATAACCTGAGTC"
            right= "TGGACAACCTTTTTT"
            fw = True
            r,c,score,_,_,_ = correctSplice(read,left,right,fw)
            #print >> sys.stderr,read
            #print >> sys.stderr,left[:c],right[c:]
            assert left[:c]+right[c:] == read
            print >> sys.stderr,"Correct Splice Test Successful!!!"

        def test_windowTransform(self):
            left = "ACGATAACCTGAGTCG"
            right= "GGTTGGACAACCTTTT"
            read     = "ACGATAACAACCTTTT"
            ref_left,ref_right = left,right
            revread = revcomp(read)
            """Needleman-Wunsch is a directional algorithm.  Since we are interested in scoring the 3' end of the right ref sequence,    we reverse complement the right ref sequence before applying the NW algorithm"""
            ref_right = revcomp(ref_right)
            score1,leftDP  = needlemanWunsch.needlemanWunsch(ref_left, read, needlemanWunsch.matchCost())
            score2,rightDP = needlemanWunsch.needlemanWunsch(ref_right,revread, needlemanWunsch.matchCost())

            #Once NW is applied, the right DP matrix must be transformed in the same coordinate frame as the left DP matrix
            rightDP = numpy.fliplr(rightDP)
            rightDP = numpy.flipud(rightDP)

            #Apply window transform to find proper sites for flanking sequences
            left_site,right_site = ("CT","AC")
            """
            Need to offset windows
            ref  --------------GT.........AG-----------
                              ^<          >>^
            """
            cost = -2
            win_left = numpy.matrix([0]+windowTransform(left,left_site,cost)+[0])
            win_right = numpy.matrix([0]+[0]+windowTransform(right,right_site,cost))

            print >> sys.stderr,"Before window transform"
            print >> sys.stderr,"read ",read
            print >> sys.stderr,"left ",left
            print >> sys.stderr,"right",right
            print >> sys.stderr,"leftDP\n",leftDP
            print >> sys.stderr,"rightDP\n",rightDP

            leftDP = leftDP+win_left
            rightDP = rightDP+win_right
            total = leftDP+rightDP
            win_right[0,16] = -10
            win_left[0,16] = -10
            print >> sys.stderr,"After window transform"
            print >> sys.stderr,"read\n   ",formatList(read,3)
            print >> sys.stderr,"left\n   ",formatList(left,3)
            print >> sys.stderr,win_left
            print >> sys.stderr,"leftDP\n",leftDP
            print >> sys.stderr,"right\n   ",formatList(right,3)
            print >> sys.stderr,win_right
            print >> sys.stderr,"rightDP\n",rightDP
            print >> sys.stderr,"read\n   ",formatList(read,3)
            print >> sys.stderr,"total\n",total


    unittest.main()

