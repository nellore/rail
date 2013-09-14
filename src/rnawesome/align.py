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
    '--splice-overlap', type=int, default=10,
    help='The overlap length of spanning readlets when evaluating splice junctions')
parser.add_argument(\
    '--faidx', type=str, required=False, help='Fasta index file')
parser.add_argument(\
    '--max-intron-length', type=int, required=False,default=1000000,
    help='Filters out all potential introns longer than this length')
parser.add_argument(\
    '--exon-differentials', action='store_const', const=True, default=False,
    help='Print exon differentials (+1s and -1s)')
parser.add_argument(\
    '--exon-intervals', action='store_const', const=True, default=False,
    help='Print exon intervals')

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
    '--intron-partition-overlap', type=int, required=False, default=20,
    help='Amount that partitions overlap their left and right neighbors by.')
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

# Print all listed exons to stdout
def printExons(refid,in_start,in_end,rdid):
    global nout
    lab = sample.parseLab(rdid)
    if args.exon_differentials:
        for pt, pt_st, pt_en in iter(partition.partition(refid, in_start, in_end, binsz)):
            # Print increment at interval start
            assert in_start < pt_en
            print "exon_diff\t%s\t%s\t%012d\t1" % (pt, lab, max(pt_st, in_start))
            nout += 1
            # Possibly print decrement at interval end
            assert in_end > pt_st
            if in_end < pt_en:
                print "exon_diff\t%s\t%s\t%012d\t-1" % (pt, lab, in_end)
                nout += 1
    if args.exon_intervals:
        for pt, _, _ in iter(partition.partition(refid, in_start, in_end, binsz)):
            print "exon_ival\t%s\t%012d\t%d\t%s" % (pt, in_start, in_end, lab)
            nout += 1

"""
Returns a more human readable format of the string
"""
def readableFormat(s):
    return " ".join([s[i:i+10] for i in range(0,len(s),10)])
def formatList(s,l):
    return (" "*l).join( list( str(s) ) )


# Print all listed introns to stdout and the flanking sequences
def printIntrons(refid,rdseq,regionSt,regionEnd,intronSt,intronEnd,rdid,fw,outhandle):
    if abs(intronEnd-intronSt) > args.max_intron_length:
        if args.verbose or args.test: print >> sys.stderr,"Huge candidate intron filtered at %s:%d-%d"%(refid,intronSt,intronEnd)
        return

    global nout
    offset = args.splice_overlap
    fw_char = "+" if fw else "-"
    #Obtain coordinates for flanking coordinate frames
    left_st,left_end = regionSt-offset,regionSt
    right_st,right_end = regionEnd,regionEnd+offset
    left_flank  = rdseq[left_st:left_end]   if fw else revcomp(rdseq[left_st:left_end])
    right_flank = rdseq[right_st:right_end] if fw else revcomp(rdseq[right_st:right_end])
    assert len(left_flank)==len(right_flank),"Bad flanks %s %s found in %s:%d-%d"%(left_flank,right_flank,refid,intronSt,intronEnd)
    lab = sample.parseLab(rdid)
    for pt, _, _ in iter(partition.partitionStartOverlaps(refid, intronSt, intronEnd, binsz, fudge=args.intron_partition_overlap)):
        print >> outhandle, "intron\t%s%s\t%012d\t%d\t%s\t%s\t%s\t%s" % (pt, fw_char, intronSt, intronEnd, lab,left_flank,right_flank,rdid)
        nout += 1

def handleUnmappedReadlets(refid,intronSt,intronEnd,rdseq,region_st,region_end,rdid,fw,fnh,offset):
    """Remaps unmapped portions of the original read between mapped readlets"""
    """
    Flanks                      ====    ====
    Read          |=================----====================|
    Ref           |=============------------------==========|
    Mapped Flanks           ====----          ----====
    """
    assert region_st<region_end
    uLen = region_end-region_st #unmapped portion
    print >> sys.stderr,"uLen",uLen
    leftSt,  leftEnd  = intronSt, intronSt+uLen
    rightSt, rightEnd = intronEnd-uLen, intronEnd

    ref_left = fnh.fetch_sequence(refid,leftSt+1, leftEnd).upper()
    ref_right = fnh.fetch_sequence(refid,rightSt+1, rightEnd).upper()

    unmapped = rdseq[region_st:region_end] if fw else revcomp(rdseq[unmapped_st:unmapped_end])
    _, diffpos, score, leftDP, rightDP, total = correctSplice(unmapped,ref_left,ref_right,fw)
    left_diff, right_diff    = diffpos, len(unmapped)-diffpos
    if score<uLen and (args.verbose or args.test): print >> sys.stderr,"Bad Needleman-Wunsch realignment"
    region_st,  region_end = region_st+left_diff,  region_end-right_diff
    intronSt,   intronEnd     = intronSt+left_diff,  intronEnd-right_diff

    printIntrons(refid,rdseq,region_st,region_end,intronSt,intronEnd,rdid,fw,sys.stdout)

def handleOverlappingFlanks(refid,intronSt,intronEnd,rdseq,region_st,region_end,rdid,fw,fnh,offset):
    """Remaps unmapped portions of the original read between mapped readlets"""
    """
    Left Flank                    ===|=
    Right Flank                     =|===
    Read          |==================|===================|
    Ref           |=============------------------==========|
    Mapped Flanks            ====                ====
    """
    assert region_st>region_end
    regLen = region_st-region_end
    region_st,region_end = region_end,region_st
    intronSt, intronEnd = intronSt-regLen, intronEnd+regLen
    handleUnmappedReadlets(refid,intronSt,intronEnd,rdseq,region_st,region_end,rdid,fw,fnh,offset)

def handleIntron(refid,intronSt,intronEnd,rdseq,region_st,region_end,rdid,fw,fnh,offset):
    """
    intronSt: reference offset for beginning of intron
    intronEnd: reference offset for end of intron
    region_st: offset from 5' end of read of LHS of splice
    region_end: offset from 5' end of read of RHS of splice
    """

    if region_st==region_end:
        """
        Scenario 1: The perfect scenario - the flanking sequences already have a good estimate
        Read   |============================================|
                                      /\
        Genome |======================--=====================|
        Flanks                   ^===^  ^===^
        """
        printIntrons(refid,rdseq,region_st,region_end,intronSt,intronEnd,rdid,fw,sys.stdout)
        return
    elif region_st<region_end:
        """
        Scenario 2: Unmapped region: Need Needleman-Wunsch to remap unmapped portions
                                       ^   ^ (Unmapped region)
        Read   |=======================-----======================|
                                      /     \
        Genome |======================-------=====================|
        Flanks                   ^===^       ^===^
        """
        handleUnmappedReadlets(refid,intronSt,intronEnd,rdseq,region_st,region_end,rdid,fw,fnh,offset)
    else:
        """
        Scenario 3: Overlapping flanking sequences - flanking sequences will overlap in the original read
                                      ^^
        Read   |=============================================|
                                      /      \
        Genome |======================-------=====================|
        Flanks                      ^===^  ^===^
        """
        handleOverlappingFlanks(refid,intronSt,intronEnd,rdseq,region_st,region_end,rdid,fw,fnh,offset)

    # ulen = unmapped_end - unmapped_st - 1 #length of the unmapped region
    # # Obtain reference genome coordinates of flanking sequences surrounding splice site
    # #TODO: Something funky is going on with this calculation
    # leftSt,  rightEnd = intronSt-offset+1,   intronEnd+offset
    # leftEnd, rightSt  = leftSt+ulen,        rightEnd-ulen

    #Scenario 3:
    # if leftEnd<=leftSt or rightEnd<=rightSt:
    #     if args.verbose or args.test: print >> sys.stderr,"Overlapping Flanking sequence scenario"
    #     printIntrons(refid,rdseq,region_st,region_end,intronSt,intronEnd,rdid,fw,sys.stdout)

    # else:
    # ref_left = fnh.fetch_sequence(refid,leftSt, leftEnd).upper()
    # ref_right = fnh.fetch_sequence(refid,rightSt, rightEnd).upper()

    # unmapped = rdseq[unmapped_st:unmapped_end] if fw else revcomp(rdseq[unmapped_st:unmapped_end])

    # _, diffpos, score, leftDP, rightDP, total = correctSplice(unmapped,ref_left,ref_right,fw)
    # left_diff,    right_diff    = diffpos,          len(unmapped)-diffpos
    # left_in_diff, right_in_diff = left_diff-offset, right_diff-offset
    # if score>0:   #If crappy alignment, disregard corrections
    #     if args.verbose or args.test: print >> sys.stderr,"Bad Needleman-Wunsch realignment"
    #     region_st,  region_end = unmapped_st+left_diff,  unmapped_end-right_diff
    #     intronSt,   intronEnd     = intronSt+left_in_diff,  intronEnd-right_in_diff

    #printIntrons(refid,rdseq,region_st,region_end,intronSt,intronEnd,rdid,fw,sys.stdout)


"""
Compares potential short intron with readlet to check if it should be an exon instead
"""
def handleShortAlignment(k,in_start,in_end,rdseq,unmapped_st,unmapped_end,region_st,region_end,rdid,fw,fnh):
    refseq = fnh.fetch_sequence(k,in_start + 1, in_end + 1).upper() # Sequence from genome
    rdsubseq = rdseq[unmapped_st:unmapped_end]
    if not fw:
        rdsubseq = revcomp(rdsubseq)
    score = needlemanWunsch.needlemanWunsch(refseq, rdsubseq, needlemanWunsch.matchCost())
    # TODO: redo this in terms of percent identity or some
    # other measure that adapts to length of the missing bit,
    # not just a raw score
    if score >= len(rdsubseq)*(9.0/10):
        printExons(k,in_start,in_end,rdid)

def getIntervals(rdals):
    ivals, positions = {}, {}
    for rdal in rdals:
        refid, fw, refoff0, seqlen, rlet_nm, _ = rdal
        refoff0, seqlen, rlet_nm = int(refoff0), int(seqlen), int(rlet_nm)
        # Remember begin, end offsets for readlet w/r/t 5' end of the read
        l, r = rlet_nm * args.readletIval, rlet_nm * args.readletIval + seqlen
        if not fw: l, r = r, l
        positions[(refid, fw, refoff0)] = l
        positions[(refid, fw, refoff0 + seqlen)] = r
        if (refid, fw) not in ivals:
            ivals[(refid, fw)] = interval.FlatIntervals()
        ivals[(refid, fw)].add(interval.Interval(refoff0, refoff0 + seqlen))
    return ivals, positions


def composeReadletAlignments(rdid, rdals, rdseq):
    global nout
    # Add this interval to the flattened interval collection for current read
    ivals, positions = getIntervals(rdals)
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
                region_st, region_end = positions[(k, fw, in_start)], positions[(k, fw, in_end)]
                if not fw: region_st, region_end = region_end, region_st
                offset = args.splice_overlap
                unmapped_st,unmapped_end = region_st-offset,region_end+offset
                reflen,rdlet_len = in_end-in_start, abs(region_end-region_st)

                assert in_start < in_end
                if abs(reflen-rdlet_len)/float(rdlet_len+1) < 0.05:
                    #Note: just a readlet missing due to sequencing error or variant
                    handleShortAlignment(k,in_start,in_end,rdseq,unmapped_st,unmapped_end,region_st,region_end,rdid,fw,fnh)
                elif rdlet_len>reflen:
                    printExons(k,in_start,in_end,rdid)
                else:
                    handleIntron(k,in_start,in_end,rdseq,region_st,region_end,rdid,fw,fnh,offset)
                in_start, in_end = en, -1
            # Keep stringing rdid along because it contains the label string
            # Add a partition id that combines the ref id and some function of
            # the offsets
            printExons(k, st, en, rdid)

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
            # Make sure mate id is part of the read name for now, so we don't
            # accidentally investigate the gap between the mates as though it's
            # a spliced alignment
            rdid = ';'.join(toks[:-3])
            cnt[rdid] = cnt.get(rdid, 0) + 1
            rd_name, rlet_nm, rdseq = toks[0], toks[3], toks[5]
            rdlet_n = int(toks[-2])
            if flags != 4:
                fw = (flags & 16) == 0
                if rdid not in mem: mem[rdid] = [ ]
                mem[rdid].append((refid, fw, refoff1-1, seqlen, rlet_nm, rd_name))
            if cnt[rdid] == rdlet_n:
                # Last readlet
                if rdid in mem:
                    # Remove mate ID so that rest of pipeline sees same rdid
                    # for both mates
                    composeReadletAlignments(rdid[:-2], mem[rdid],rdseq)
                    del mem[rdid]
                del cnt[rdid]
            line_nm += 1
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
    sys.stdout.flush()
    bowtieOutDone.set()

def writeReads(fhs, reportMult=1.2):
    """ Parse input reads, optionally transform them and/or turn them into
        readlets. """
    global ninp
    report = 1.0
    for ln in sys.stdin:
        toks = ln.rstrip().split('\t')
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
                rlets1 = readlet.readletize(args, nm1 + ';1', seq1, qual1)
                for i in xrange(0, len(rlets1)):
                    nm_rlet, seq_rlet, qual_rlet = rlets1[i]
                    rdletStr = "%s;%d;%d;%s\t%s\t%s\n" % (nm_rlet, i, len(rlets1), seq1, seq_rlet, qual_rlet)
                    if ninp >= report and i == 0:
                        report *= reportMult
                        if args.verbose: print >> sys.stderr, "First readlet from read %d: '%s'" % (ninp, rdletStr.rstrip())
                    for fh in fhs: fh.write(rdletStr)
                rlets2 = readlet.readletize(args, nm2 + ';2', seq2, qual2)
                for i in xrange(0, len(rlets2)):
                    nm_rlet, seq_rlet, qual_rlet = rlets2[i]
                    rdletStr = "%s;%d;%d;%s\t%s\t%s\n" % (nm_rlet, i, len(rlets2), seq2, seq_rlet, qual_rlet)
                    if ninp >= report and i == 0:
                        report *= reportMult
                        if args.verbose: print >> sys.stderr, "First readlet from read %d: '%s'" % (ninp, rdletStr.rstrip())
                    for fh in fhs: fh.write(rdletStr)
            else:
                rdStr = "%s\t%s\t%s\t%s\t%s\n" % (nm1, seq1, qual1, seq2, qual2)
                if ninp >= report:
                    report *= reportMult
                    if args.verbose: print >> sys.stderr, "Read %d: '%s'" % (ninp, rdStr.rstrip())
                for fh in fhs: fh.write(rdStr)
        else:
            # Unpaired
            if xformReads:
                # Truncate and transform quality values
                seq, qual = xformRead(seq, qual)
            if args.readletLen > 0:
                # Readletize
                rlets = readlet.readletize(args, nm + ';0', seq, qual)
                for i in xrange(0, len(rlets)):
                    nm_rlet, seq_rlet, qual_rlet = rlets[i]
                    rdletStr = "%s;%d;%d;%s\t%s\t%s\n" % (nm_rlet, i, len(rlets), seq, seq_rlet, qual_rlet)
                    if ninp >= report and i == 0:
                        report *= reportMult
                        if args.verbose: print >> sys.stderr, "First readlet from read %d: '%s'" % (ninp, rdletStr.rstrip())
                    for fh in fhs: fh.write(rdletStr)
            else:
                rdStr = "%s\t%s\t%s\n" % (nm, seq, qual)
                if ninp >= report:
                    report *= reportMult
                    if args.verbose: print >> sys.stderr, "Read %d: '%s'" % (ninp, rdStr.rstrip())
                for fh in fhs: fh.write(rdStr)

def go():

    import time
    timeSt = time.time()

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
        proc, mycmd, threads = bowtie.proc(args, readFn=readFn,
                                           bowtieArgs=bowtieArgs, sam=True,
                                           outHandler=bowtieOutReadlets,
                                           stdinPipe=False)
    else:
        # Reads are written to Bowtie process's stdin directly
        proc, mycmd, threads = bowtie.proc(args, readFn=None,
                                           bowtieArgs=bowtieArgs, sam=True,
                                           outHandler=bowtieOutReadlets,
                                           stdinPipe=True)
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
    if args.verbose: print >>sys.stderr, "Waiting for Bowtie to finish"
    bowtieOutDone.wait()
    for thread in threads:
        if args.verbose: print >> sys.stderr, "  Joining a thread..."
        thread.join()
        if args.verbose: print >> sys.stderr, "    ...joined!"
    sys.stdout.flush()
    if args.verbose:
        print >>sys.stderr, "Bowtie finished"

    # Remove any temporary reads files created
    if args.serial and args.write_reads is None and not args.keep_reads:
        print >>sys.stderr, "Cleaning up temporary files"
        import shutil
        shutil.rmtree(tmpdir)

    print >> sys.stderr, "DONE with align.py; in/out = %d/%d; time=%0.3f secs" % (ninp, nout, time.time()-timeSt)

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
        ###Big Note:  We are going to assume base-0 indexing for everything
        def setUp(self):
            #A visual representation of the reference sequence and the read
            """Read"""
            """ACGAAGGACT GCTTGACATC GGCCACGATA ACAACCTTTT TTGCGCCAAT CTTAAGAGCC TTCT"""
            #             ^10        ^20        ^30        ^40        ^50        ^60
            """Genome"""
            """ACGAAGGACT GCTTGACATC GGCCACGATA ACCTGAGTCG ATAGGACGAA ACAAGTATAT ATTCGAAAAT TAATTAATTC CGAAATTTCA ATTTCATCCG ACATGTATCT ACATATGCCA CACTTCTGGT TGGACAACCT TTTTTGCGCC A"""
            """ACGAAGGACT GCTTGACATC GGCCACGATA AC                                                                                                                 AACCT TTTTTGCGCC AATCTTAAGA GCCTTCT"""
            #             ^10        ^20        ^30        ^40        ^50        ^60        ^70        ^80        ^90        ^100       ^110       ^120       ^130       ^140       ^150
            self.rdseq  = "ACGAAGGACTGCTTGACATCGGCCACGATAACAACCTTTTTTGCGCCAATCTTAAGAGCCTTCT"
            self.refseq = "ACGAAGGACTGCTTGACATCGGCCACGATAACCTGAGTCGATAGGACGAAACAAGTATATATTCGAAAATTAATTAATTCCGAAATTTCAATTTCATCCGACATGTATCTACATATGCCACACTTCTGGTTGGACAACCTTTTTTGCGCCA"
            self.testDump = "test.out"
            self.fasta = "test.fa"
            self.faidx = "test.fa.fai"
            createTestFasta(self.fasta,"test",self.refseq)
            open(self.testDump,'w') #Just to initialize file
        def tearDown(self):
            os.remove(self.fasta)
            os.remove(self.faidx)
            os.remove(self.testDump)


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

            # print >> sys.stderr,"Before window transform"
            # print >> sys.stderr,"read ",read
            # print >> sys.stderr,"left ",left
            # print >> sys.stderr,"right",right
            # print >> sys.stderr,"leftDP\n",leftDP
            # print >> sys.stderr,"rightDP\n",rightDP

            leftDP = leftDP+win_left
            rightDP = rightDP+win_right
            total = leftDP+rightDP
            win_right[0,16] = -10
            win_left[0,16] = -10
            # print >> sys.stderr,"After window transform"
            # print >> sys.stderr,"read\n   ",formatList(read,3)
            # print >> sys.stderr,"left\n   ",formatList(left,3)
            # print >> sys.stderr,win_left
            # print >> sys.stderr,"leftDP\n",leftDP
            # print >> sys.stderr,"right\n   ",formatList(right,3)
            # print >> sys.stderr,win_right
            # print >> sys.stderr,"rightDP\n",rightDP
            # print >> sys.stderr,"read\n   ",formatList(read,3)
            # print >> sys.stderr,"total\n",total
        """
        Scenario 1:
        Read   |============================================|
                                      /\
        Genome |======================--=====================|
        Flanks                   ^===^  ^===^
        """
        def testScenario1(self):
            sys.stdout = open(self.testDump,'w')
            rdid,fw,refid = "0;LB:test",True,"test"
            #leftSt,leftEnd = 21,37  #left coords
            #rightSt,rightEnd = 143,154  #left coords

            iSt,iEnd = 32,135 #intron coords
            rSt,rEnd = 32,32  #region coords
            offset = 10
            fnh = fasta.fasta(self.fasta)
            handleIntron(refid,iSt,iEnd,self.rdseq,
                         rSt,rEnd,rdid,fw,fnh,offset)
            sys.stdout.close()
            test_out = open(self.testDump,'r')
            testLine = test_out.readline().rstrip()
            toks = testLine.split("\t")
            st,end,lab,leftFlank,rightFlank,rdid = int(toks[2]), int(toks[3]), toks[4], toks[5], toks[6], toks[7]
            self.assertEquals(st,32)
            self.assertEquals(end,135)
            self.assertEquals(leftFlank,"CCACGATAAC")
            self.assertEquals(rightFlank,"AACCTTTTTT")

        """
        Scenario 2: Unmapped region
                                       ^   ^ (Unmapped region)
        Read   |=======================-----======================|
                                      /     \
        Genome |======================-------=====================|
        Flanks                   ^===^       ^===^
        """
        def testScenario2(self):
            sys.stdout = open(self.testDump,'w')
            rdid,fw,refid = "0;LB:test",True,"test"
            #leftSt,leftEnd = 21,37  #left coords
            #rightSt,rightEnd = 143,154  #left coords

            iSt,iEnd = 28,139 #intron coords
            rSt,rEnd = 28,36  #region coords
            offset = 10
            fnh = fasta.fasta(self.fasta)
            handleIntron(refid,iSt,iEnd,self.rdseq,
                         rSt,rEnd,rdid,fw,fnh,offset)
            sys.stdout.close()
            test_out = open(self.testDump,'r')
            testLine = test_out.readline().rstrip()
            toks = testLine.split("\t")
            st,end,lab,leftFlank,rightFlank,rdid = int(toks[2]), int(toks[3]), toks[4], toks[5], toks[6], toks[7]
            self.assertTrue( abs(st-32) < 4)
            self.assertTrue( abs(end-135) < 4)
            scoreLeft,_  = needlemanWunsch.needlemanWunsch(leftFlank,"CCACGATAAC" , needlemanWunsch.matchCost())
            scoreRight,_ = needlemanWunsch.needlemanWunsch(rightFlank,"AACCTTTTTT" , needlemanWunsch.matchCost())
            print >> sys.stderr,"Scores",scoreLeft,scoreRight
            self.assertTrue(scoreLeft > 4)
            self.assertTrue(scoreRight > 4)
        """
            Scenario 3: Overlapping flanking sequences - flanking sequences will overlap in the original read
                                      ^^
        Read   |=============================================|
                                      /      \
        Genome |======================-------=====================|
        Flanks                      ^===^  ^===^
        """
        def testScenario3(self):
            sys.stdout = open(self.testDump,'w')
            rdid,fw,refid = "0;LB:test",True,"test"
            #leftSt,leftEnd = 21,37  #left coords
            #rightSt,rightEnd = 143,154  #left coords

            iSt,iEnd = 36,131 #intron coords
            rSt,rEnd = 36,28  #region coords
            offset = 10
            fnh = fasta.fasta(self.fasta)
            handleIntron(refid,iSt,iEnd,self.rdseq,
                         rSt,rEnd,rdid,fw,fnh,offset)
            sys.stdout.close()
            test_out = open(self.testDump,'r')
            testLine = test_out.readline().rstrip()
            toks = testLine.split("\t")
            st,end,lab,leftFlank,rightFlank,rdid = int(toks[2]), int(toks[3]), toks[4], toks[5], toks[6], toks[7]

            self.assertTrue( abs(st-32) < 4)
            self.assertTrue( abs(end-135) < 4)
            scoreLeft,_  = needlemanWunsch.needlemanWunsch(leftFlank,"CCACGATAAC" , needlemanWunsch.matchCost())
            scoreRight,_ = needlemanWunsch.needlemanWunsch(rightFlank,"AACCTTTTTT" , needlemanWunsch.matchCost())
            print >> sys.stderr,"Scores",scoreLeft,scoreRight
            self.assertTrue(scoreLeft > 4)
            self.assertTrue(scoreRight > 4)


    unittest.main()

