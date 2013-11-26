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

An input stream can be a mix of the above formats.

Binning/sorting prior to this step:
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
from Queue import Queue
import string
import numpy
import shutil

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
import tempfile

ninp = 0               # # lines input so far
nout = 0               # # lines output so far
pe = False             # What is this variable? (Probably stands for "print exons" -AN)
discardMate = -1       # don't discard mates by default
lengths = dict()       # read lengths after truncation
rawLengths = dict()    # read lengths prior to truncation
qualCnts = dict()      # quality counts after adjustments
rawQualCnts = dict()   # quality counts before adjustments
qualAdd = None         # amt to add to qualities
truncateAmt = None     # amount to truncate reads
truncateTo = None      # amount to truncate reads
readletize = None      # if we're going to readletize

# Review with Ben: do we want to transform reads?
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
    '--min-intron-length', type=int, required=False,default=10,
    help='Filters out all potential introns shorter than this length')
parser.add_argument(\
    '--exon-differentials', action='store_const', const=True, default=False,
    help='Print exon differentials (+1s and -1s)')
parser.add_argument(\
    '--exon-intervals', action='store_const', const=True, default=False,
    help='Print exon intervals')

bowtie.addArgs(parser)
readlet.addArgs(parser)
partition.addArgs(parser)

#parser.add_argument(\
#    '--serial', action='store_const', const=True, default=False, help="Run bowtie serially after rather than concurrently with the input-reading loop")
parser.add_argument(\
    '--only-readletize', action='store_const', const=True, default=False, help="Don't run first pass of Bowtie, where full reads are aligned")
# Since args.serial is always true now, keep_reads is equivalent to archive; review with Ben
#parser.add_argument(\
#    '--keep-reads', action='store_const', const=True, default=False, help="Don't delete any temporary read file(s) created")
# args.write_reads would need to accommodate two different filenames (for each pass of Bowtie) now, so kill it; review with Ben
#parser.add_argument(\
#    '--write-reads', type=str, required=False, help='Write input reads to given tab-delimited file')
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

# Global args variable below takes input
# from command line, but args is still
# an argument of the go() function
args = parser.parse_args(argv[1:])

'''class RenamedTemporaryFile(object):
    """For creating a temporary file object (using tempfile) that will be renamed
    on exit. From
    http://stackoverflow.com/questions/12003805/threadsafe-and-fault-tolerant-file-writes .
    This class facilitates threadsafe writes. Review with Ben; even when only single thread (OutputThread) 
    is writing to file, vanilla file I/O causes hiccups (missed characters at beginnings of lines)
    once every ~10k lines. Confirmed that output buffers were always flushed."""

    def __init__(self, final_path, **kwargs):
        tmpfile_dir = kwargs.pop('dir', None)
        if tmpfile_dir is None:
            tmpfile_dir = os.path.dirname(final_path)

        self.tmpfile = tempfile.NamedTemporaryFile(dir=tmpfile_dir, **kwargs)
        self.final_path = final_path

    def __getattr__(self, attr):
        return getattr(self.tmpfile, attr)

    def __enter__(self):
        self.tmpfile.__enter__()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if exc_type is None:
            self.tmpfile.delete = False
            result = self.tmpfile.__exit__(exc_type, exc_val, exc_tb)
            os.rename(self.tmpfile.name, self.final_path)
        else:
            result = self.tmpfile.__exit__(exc_type, exc_val, exc_tb)

        return result'''

def xformRead(seq, qual):
    # Possibly truncate and/or modify quality values
    # TODO: not implemented yet!
    # Review with Ben: if we'll never transform reads, we should kill this
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

"""
Applies Needleman-Wunsch to correct splice junction gaps
"""
def correctSplice(read,ref_left,ref_right,fw):

    revread = revcomp(read)

    # Needleman-Wunsch is a directional algorithm.  Since we are interested in scoring the 3' end of the right ref sequence, we reverse complement the right ref sequence before applying the NW algorithm
    ref_right = revcomp(ref_right)
    _, leftDP  = needlemanWunsch.needlemanWunsch(ref_left, read, needlemanWunsch.matchCost())
    _, rightDP = needlemanWunsch.needlemanWunsch(ref_right, revread, needlemanWunsch.matchCost())

    # Once NW is applied, right DP matrix must be transformed in the same
    # coordinate frame as the left DP matrix
    rightDP = numpy.fliplr(rightDP)
    rightDP = numpy.flipud(rightDP)

    total = leftDP+rightDP
    index = numpy.argmax(total)

    n = len(read) + 1
    _, c = index % n, index / n

    c = needlemanWunsch.medianTieBreaker(total,c)
    r = numpy.argmax(total[:,c])

    return r,c,total[r,c],leftDP,rightDP,total

# Print all listed exons to stdout
def printExons(refid,in_start,in_end,rdid,exon_differentials,exon_intervals):
    global nout
    lab = sample.parseLab(rdid)
    if exon_differentials:
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
    if exon_intervals:
        for pt, _, _ in iter(partition.partition(refid, in_start, in_end, binsz)):
            print "exon_ival\t%s\t%012d\t%d\t%s" % (pt, in_start, in_end, lab)
            nout += 1

# Print all listed introns to stdout and the flanking sequences
def printIntrons(args, refid,rdseq,regionSt,regionEnd,intronSt,intronEnd,rdid,fw,outhandle):
    if abs(intronEnd-intronSt) > args.max_intron_length:
        if args.verbose or args.test: print >> sys.stderr,"Huge candidate intron filtered at %s:%d-%d"%(refid,intronSt,intronEnd)
        return

    global nout
    fw_char = "+" if fw else "-"
    lab = sample.parseLab(rdid)
    for pt, _, _ in iter(partition.partitionStartOverlaps(refid, intronSt, intronEnd, binsz, fudge=args.intron_partition_overlap)):
        print >> outhandle, "intron\t%s%s\t%012d\t%d\t%s\t%s" % (pt, fw_char, intronSt, intronEnd, lab,rdid)
        nout += 1

def handleUnmappedReadlets(args, refid,intronSt,intronEnd,rdseq,regionSt,regionEnd,rdid,fw,fnh,offset):
    """Remaps unmapped portions of the original read between mapped readlets"""
    """
    Flanks                      ====    ====
    Read          |=================----====================|
    Ref           |=============------------------==========|
    Mapped Flanks           ====----          ----====
    """
    assert regionSt<regionEnd
    #uLen = regionEnd-regionSt #unmapped portion
    #readlet length may be shorter than expected, so define uLen as follows
    uLen = len(rdseq[regionSt:regionEnd])
    leftSt,  leftEnd  = intronSt, intronSt+uLen
    rightSt, rightEnd = intronEnd-uLen, intronEnd

    ref_left = fnh.fetch_sequence(refid,leftSt+1, leftEnd).upper()
    ref_right = fnh.fetch_sequence(refid,rightSt+1, rightEnd).upper()

    unmapped = rdseq[regionSt:regionEnd] if fw else revcomp(rdseq[regionSt:regionEnd])
    _, diffpos, score, _, _, _ = correctSplice(unmapped,ref_left,ref_right,fw)
    left_diff, right_diff  = diffpos, len(unmapped)-diffpos
    if score<uLen and (args.verbose or args.test): print >> sys.stderr,"Bad Needleman-Wunsch realignment with %s = %s + %s found intron %s:%d-%d and read region %s-%s\n"%(unmapped,ref_left,ref_right,refid,intronSt,intronEnd,regionSt,regionEnd)
    regionSt,  regionEnd = regionSt+left_diff-1,  regionEnd-right_diff
    intronSt,   intronEnd     = intronSt+left_diff,  intronEnd-right_diff

    printIntrons(args, refid,rdseq,regionSt,regionEnd,intronSt,intronEnd,rdid,fw,sys.stdout)

def handleOverlappingFlanks(args, refid,intronSt,intronEnd,rdseq,regionSt,regionEnd,rdid,fw,fnh,offset):
    """Remaps unmapped portions of the original read between mapped readlets"""
    """
    Left Flank                    ===|=
    Right Flank                     =|===
    Read          |==================|===================|
    Ref           |=============------------------==========|
    Mapped Flanks            ====                ====
    """
    assert regionSt>regionEnd
    regLen = regionSt-regionEnd
    regionSt,regionEnd = regionEnd,regionSt
    intronSt, intronEnd = intronSt-regLen, intronEnd+regLen #read just intron boundaries
    handleUnmappedReadlets(args, refid,intronSt,intronEnd,rdseq,regionSt,regionEnd,rdid,fw,fnh,offset)

def handleIntron(args, refid,intronSt,intronEnd,rdseq,regionSt,regionEnd,rdid,fw,fnh,offset):
    """
    intronSt: reference offset for beginning of intron
    intronEnd: reference offset for end of intron
    regionSt: offset from 5' end of read of LHS of splice
    regionEnd: offset from 5' end of read of RHS of splice
    """
    if regionSt==regionEnd:
        """
        Scenario 1: The perfect scenario - the flanking sequences already have a good estimate
        Read   |============================================|
                                      /\
        Genome |======================--=====================|
        Flanks                   ^===^  ^===^
        """
        printIntrons(args, refid,rdseq,regionSt,regionEnd,intronSt,intronEnd,rdid,fw,sys.stdout)
        return
    elif regionSt<regionEnd:
        """
        Scenario 2: Unmapped region - Need Needleman-Wunsch to remap unmapped portions
                                       ^   ^ (Unmapped region)
        Read   |=======================-----======================|
                                      /     \
        Genome |======================-------=====================|
        Flanks                   ^===^       ^===^
        """
        handleUnmappedReadlets(args, refid,intronSt,intronEnd,rdseq,regionSt,regionEnd,rdid,fw,fnh,offset)
    else:
        """
        Scenario 3: Overlapping flanking sequences - flanking sequences will overlap in the original read
                                      ^^
        Read   |=============================================|
                                      /      \
        Genome |======================-------=====================|
        Flanks                      ^===^  ^===^
        """
        handleOverlappingFlanks(args, refid,intronSt,intronEnd,rdseq,regionSt,regionEnd,rdid,fw,fnh,offset)

def handleShortAlignment(args, k,intronSt,intronEnd,rdseq,regionSt,regionEnd,rdid,fw,fnh,exon_differentials,exon_intervals):
    """
    When the intervals of unaligned read & reference characters are similar in
    length, we say there's no intron but we also check to see if there's
    enough similarity in there to say that the middle is also part of the
    exon.
    """
    if regionSt>=regionEnd:
        #Reject because this exon interval has already been evaluated by the surrounding readlets
        return
    #Note: fetch_sequence is base 1 indexed ainclusive
    refseq = fnh.fetch_sequence(k, intronSt + 1, intronEnd).upper() # Sequence from genome
    rdsubseq = rdseq[regionSt:regionEnd]
    if not fw:
        rdsubseq = revcomp(rdsubseq)
    score,scoreMat = needlemanWunsch.needlemanWunsch(refseq, rdsubseq, needlemanWunsch.matchCost())

    # TODO: redo this in terms of percent identity or some
    # other measure that adapts to length of the missing bit,
    # not just a raw score (This is done below, no? -AN)
    if score >= len(rdsubseq)*(5.0/10):
        printExons(k,intronSt,intronEnd,rdid,exon_differentials,exon_intervals)

def getIntervals(args, rdals):
    """ Given a collection of readlet-alignment tuples, turn them into a set
        of interval.Interval objects, one per chromosome overlapped at least
        once. """
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

def composeReadletAlignments(args, rdid, rdals, rdseq, fnh):
    rdseq = rdseq.upper()
    ivals, positions = getIntervals(args, rdals)
    for kfw in ivals.iterkeys(): # for each chromosome covered by >= 1 readlet
        k, fw = kfw
        lastEn = None
        for iv in sorted(iter(ivals[kfw])): # for each covered interval, left-to-right
            st, en = iv.start, iv.end
            assert en > st and st >= 0 and en >= 0
            if lastEn is not None:
                intronSt, intronEnd = lastEn, st
                regionSt, regionEnd = positions[(k, fw, intronSt)], positions[(k, fw, intronEnd)]
                if not fw: regionSt, regionEnd = regionEnd, regionSt
                # Now regionSt, regionEn are w/r/t alignment's "left" end
                offset = args.splice_overlap
                reflen, readlen = intronEnd - intronSt, regionEnd - regionSt
                potentialIntronLen = reflen - readlen
                if potentialIntronLen < args.min_intron_length:
                    handleShortAlignment(args, k, intronSt, intronEnd, rdseq, regionSt, regionEnd, rdid, fw, fnh,args.exon_differentials,args.exon_intervals)
                else:
                    handleIntron(args, k, intronSt, intronEnd, rdseq, regionSt, regionEnd, rdid, fw, fnh, offset)
                intronSt, intronEnd = en, -1
            # Keep stringing rdid along because it contains the label string
            # Add a partition id that combines the ref id and some function of
            # the offsets
            printExons(k, st, en, rdid, args.exon_differentials, args.exon_intervals)
            lastEn = en

class OutputThread(threading.Thread):

    """ A worker thread that examines SAM output from Bowtie and emits
        appropriate tuples for exons and introns.  Each line of output
        is another SAM readlet alignment. """
    def __init__(self, args, st, done, outq, tmpdir, fnh, firstPass=True, reportMult=1.2):
        super(OutputThread, self).__init__()
        self.st = st
        self.firstPass = firstPass
        self.reportMult = 1.2 # shouldn't this be self.reportMult = reportMult? -AN
        self.nout = 0
        self.done = done
        self.outq = outq
        self.fnh = fnh
        self.tmpdir = tmpdir
        self.args = args
        self.stoprequest = threading.Event()

    def run(self):
        """ Main driver method for the output thread """
        mem, cnt = {}, {}
        report = 1
        line_nm = 0
        st = self.st
        exc = None
        nsam = 0
        exitlevel = 0
        unmappedFn = os.path.join(self.tmpdir, 'unmappedreads.tsv')
        try:
            with open(unmappedFn, 'w') as unmappedfh:
                # This puts unmapped reads in the temporary directory
                # If it's the second pass, the with above is redundant
                # but does not affect performance
                while not self.stoprequest.isSet():
                    line = st.readline()
                    if len(line) == 0:
                        break # no more output
                    if line[0] == '@':
                        continue # skip header
                    nsam += 1
                    rdid, flags, refid, refoff1, _, _, _, _, _, seq, qual, _ = string.split(line.rstrip(), '\t', 11)
                    flags, refoff1 = int(flags), int(refoff1)
                    if nsam >= report:
                        report *= self.reportMult
                        if self.args.verbose:
                            print >>sys.stderr, "SAM output record %d: rdname='%s', flags=%d" % (nsam, rdid, flags)
                    seqlen = len(seq)
                    toks = rdid.split(';')
                    # Make sure mate id is part of the read name for now, so we don't
                    # accidentally investigate the gap between the mates as though it's
                    # a spliced alignment
                    rdid = ';'.join(toks[:-3])
                    cnt[rdid] = cnt.get(rdid, 0) + 1
                    rd_name, rlet_nm, rdseq = toks[0], toks[4], toks[6]
                    rdlet_n = int(toks[-2])
                    if flags != 4:
                        fw = (flags & 16) == 0
                        if rdid not in mem: mem[rdid] = [ ]
                        mem[rdid].append((refid, fw, refoff1-1, seqlen, rlet_nm, rd_name))
                    elif self.firstPass:
                        if len(("%s\t%s\t%s" % (rdid[:-2], seq, qual)).split('\t')) == 3 and len(seq) == len(qual):
                            print >>unmappedfh, "%s\t%s\t%s" % (rdid[:-2], seq, qual)
                    if cnt[rdid] == rdlet_n:
                        # Last readlet
                        if rdid in mem:
                            # Remove mate ID so that rest of pipeline sees same rdid
                            # for both mates
                            composeReadletAlignments(self.args, rdid[:-2], mem[rdid], rdseq, self.fnh)
                            del mem[rdid]
                        del cnt[rdid]
                    line_nm += 1
        except Exception as e:
            print >> sys.stderr, "Exception while reading output from Bowtie"
            exitlevel = 1
            exc = e
        sys.stdout.flush()
        self.outq.put(exitlevel)
        self.done.set()
        while len(st.readline()) > 0:
            pass
        if exc is not None:
            import traceback
            traceback.print_exc()
        else:
            assert len(mem) == 0
            assert len(cnt) == 0

class OutputThread2(threading.Thread):

    """ A worker thread that examines SAM output from Bowtie and emits
        appropriate tuples for exons and introns.  Each line of output
        is another SAM readlet alignment. """
    def __init__(self, args, st, done, outq, unmappedst, fnh, reportMult=1.2):
        super(OutputThread2, self).__init__()
        self.st = st
        self.reportMult = 1.2 # shouldn't this be self.reportMult = reportMult? -AN
        self.nout = 0
        self.done = done
        self.outq = outq
        self.fnh = fnh
        self.args = args
        self.unmappedst = unmappedst
        self.stoprequest = threading.Event()

    def run(self):
        """ Main driver method for the output thread """
        mem, cnt = {}, {}
        report = 1
        line_nm = 0
        exc = None
        nsam = 0
        exitlevel = 0
        try:
            # This puts unmapped reads in the temporary directory
            # If it's the second pass, the with above is redundant
            # but does not affect performance
            while not self.stoprequest.isSet():
                line = self.st.readline()
                if len(line) == 0:
                    break # no more output
                if line[0] == '@':
                    continue # skip header
                nsam += 1
                rdid, flags, refid, refoff1, _, _, _, _, _, seq, qual, _ = string.split(line.rstrip(), '\t', 11)
                flags, refoff1 = int(flags), int(refoff1)
                if nsam >= report:
                    report *= self.reportMult
                    if self.args.verbose:
                        print >>sys.stderr, "SAM output record %d: rdname='%s', flags=%d" % (nsam, rdid, flags)
                seqlen = len(seq)
                toks = rdid.split(';')
                # Make sure mate id is part of the read name for now, so we don't
                # accidentally investigate the gap between the mates as though it's
                # a spliced alignment
                rdid = ';'.join(toks[:-3])
                cnt[rdid] = cnt.get(rdid, 0) + 1
                rd_name, rlet_nm, rdseq = toks[0], toks[4], toks[6]
                rdlet_n = int(toks[-2])
                if flags != 4:
                    fw = (flags & 16) == 0
                    if rdid not in mem: mem[rdid] = [ ]
                    mem[rdid].append((refid, fw, refoff1-1, seqlen, rlet_nm, rd_name))
                else:
                    if len(("%s\t%s\t%s" % (rdid[:-2], seq, qual)).split('\t')) == 3 and len(seq) == len(qual):
                        print >>self.unmappedst, "%s\t%s\t%s" % (rdid[:-2], seq, qual)
                if cnt[rdid] == rdlet_n:
                    # Last readlet
                    if rdid in mem:
                        # Remove mate ID so that rest of pipeline sees same rdid
                        # for both mates
                        composeReadletAlignments(self.args, rdid[:-2], mem[rdid], rdseq, self.fnh)
                        del mem[rdid]
                    del cnt[rdid]
                line_nm += 1
        except Exception as e:
            print >> sys.stderr, "Exception while reading output from Bowtie"
            exitlevel = 1
            exc = e
        sys.stdout.flush()
        self.outq.put(exitlevel)
        self.done.set()
        while len(self.st.readline()) > 0:
            pass
        if exc is not None:
            import traceback
            traceback.print_exc()
        else:
            assert len(mem) == 0
            assert len(cnt) == 0

def writeReads(args, fh, bowtieOutDone, inst=None, tmpdir=None, firstPass=True, reportMult=1.2):
    """ Parse input reads, optionally transform them and/or turn them into
        readlets. """
    global ninp
    report = 1.0
    assert firstPass == (tmpdir is None) # tmpdir is provided only on second pass
    if inst is None:
        if firstPass or args.only_readletize:
            inst = sys.stdin
        else:
            unmappedfh = open(os.path.join(tmpdir, 'unmappedreads.tsv'), 'r')
            inst = unmappedfh
    for ln in inst:
        if bowtieOutDone.is_set():
            return
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
                if (discardMate & 1) != 0:
                    nm, seq, qual = nm2, seq2, qual2 # discard mate 1
                if (discardMate & 2) != 0:
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
            if args.readletLen > 0 and not firstPass:
                # Readletize
                rlets1 = readlet.readletize(args, nm1 + ';1', seq1, qual1)
                for i in xrange(len(rlets1)):
                    nm_rlet, seq_rlet, qual_rlet = rlets1[i]
                    rdletStr = "%s;%d;%d;%s\t%s\t%s" % (nm_rlet, i, len(rlets1), seq1, seq_rlet, qual_rlet)
                    if ninp >= report and i == 0:
                        report *= reportMult
                        if args.verbose: print >> sys.stderr, "First readlet from read %d: '%s'" % (ninp, rdletStr.rstrip())
                    print >>fh, rdletStr
                rlets2 = readlet.readletize(args, nm2 + ';2', seq2, qual2)
                for i in xrange(len(rlets2)):
                    nm_rlet, seq_rlet, qual_rlet = rlets2[i]
                    rdletStr = "%s;%d;%d;%s\t%s\t%s" % (nm_rlet, i, len(rlets2), seq2, seq_rlet, qual_rlet)
                    if ninp >= report and i == 0:
                        report *= reportMult
                        if args.verbose: print >> sys.stderr, "First readlet from read %d: '%s'" % (ninp, rdletStr.rstrip())
                    print >>fh, rdletStr
            else:
                # Align entire reads first, but use same convention for rdStr as for readlets
                # Review with Ben; this fixes bug in OutputThread
                rdStr = "%s;%d;%d;%s\t%s\t%s" % (nm + ';0', 0, 1, seq, seq, qual)
                if ninp >= report:
                    report *= reportMult
                    if args.verbose: print >> sys.stderr, "First read %d: '%s'" % (ninp, rdStr.rstrip())
                print >>fh, rdStr
            # Code below was not compatible with that in OutputThread
            '''else:
                rdStr = "%s\t%s\t%s\t%s\t%s\n" % (nm1, seq1, qual1, seq2, qual2)
                if ninp >= report:
                    report *= reportMult
                    if args.verbose: print >> sys.stderr, "Read %d: '%s'" % (ninp, rdStr.rstrip())
                for fh in fhs: fh.write(rdStr)'''
        else:
            # Unpaired
            if xformReads:
                # Truncate and transform quality values
                seq, qual = xformRead(seq, qual)
            if args.readletLen > 0 and not firstPass:
                # Readletize
                rlets = readlet.readletize(args, nm + ';0', seq, qual)
                for i in xrange(len(rlets)):
                    nm_rlet, seq_rlet, qual_rlet = rlets[i]
                    rdletStr = "%s;%d;%d;%s\t%s\t%s" % (nm_rlet, i, len(rlets), seq, seq_rlet, qual_rlet)
                    if ninp >= report and i == 0:
                        report *= reportMult
                        if args.verbose: print >> sys.stderr, "First readlet from read %d: '%s'" % (ninp, rdletStr.rstrip())
                    print >>fh, rdletStr
            else:
                rdStr = "%s;%d;%d;%s\t%s\t%s" % (nm + ';0', 0, 1, seq, seq, qual)
                if ninp >= report:
                    report *= reportMult
                    if args.verbose: print >> sys.stderr, "First read %d: '%s'" % (ninp, rdStr.rstrip())
                print >>fh, rdStr

    fh.flush()

import itertools
def erat2():
    # generate primes; uses sieve of eratosthenes
    # this is used in performBowtieIteration's code to obtain readlets,
    # ensuring that readlet "breaks" don't overlap
    # from Python Cookbook
    D = {}
    yield 2
    for q in itertools.islice(itertools.count(3), 0, None, 2):
        p = D.pop(q, None)
        if p is None:
            D[q*q] = q
            yield q
        else:
            x = p + q
            while x in D or not (x&1):
                x += p
            D[x] = p

def performBowtieIteration(args, bowtieArgs, tmpdir, fnh, passNumber=0, readletDepth=3, readLabelRoot="readlets", reportMult=1.2, inst=sys.stdin):

    global ninp
    itninp = 0
    bowtieOutDone = threading.Event()
    threads = []
    outq = Queue()
    report = 1.0
    assert readletDepth >= 1

    # Prepare readlets
    readletFn = os.path.join(tmpdir, "%s_pass_%d.tsv" % (readLabelRoot, passNumber))
    with open(readletFn, 'w') as readletfh:

        if not passNumber:
            # if it's the first pass of Bowtie, don't split reads
            for ln in inst:
                toks = ln.rstrip().split('\t')
                ninp += 1
                itninp += 1
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
                        if (discardMate & 1) != 0:
                            nm, seq, qual = nm2, seq2, qual2 # discard mate 1
                        if (discardMate & 2) != 0:
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
                    else:
                        # Align entire reads first, but use same convention for rdStr as for readlets
                        # Review with Ben; this fixes bug in OutputThread
                        rdStr = "%s;%d;%d;%s\t%s\t%s" % (nm + ';0', 0, 1, seq, seq, qual)
                        if itninp >= report:
                            report *= reportMult
                            if args.verbose: print >> sys.stderr, "First read %d: '%s'" % (itninp, rdStr)
                        print >>readletfh, rdStr
                else:
                    # Unpaired
                    if xformReads:
                        # Truncate and transform quality values
                        seq, qual = xformRead(seq, qual)
                    else:
                        rdStr = "%s;%d;%d;%s\t%s\t%s" % (nm + ';0', 0, 1, seq, seq, qual)
                        if itninp >= report:
                            report *= reportMult
                            if args.verbose: print >> sys.stderr, "First read %d: '%s'" % (itninp, rdStr.rstrip())
                        print >>readletfh, rdStr
        else:
            primes = []
            nextPrimeGen = erat2()
            nextPrime = nextPrimeGen.next()
            while nextPrime <= readletDepth:
                primes.append(nextPrime)
                nextPrime = nextPrimeGen.next()
            for ln in inst:
                ninp += 1
                itninp += 1
                toks = ln.rstrip().split('\t')
                # every performBowtieIteration() call after first should be fed data in three-token format
                assert len(toks) == 3
                nm, seq, qual = toks
                seqlen = len(seq)
                assert seqlen == len(qual)
                readlets = []
                for currentDepth in reversed(primes):
                    readletLength = seqlen / currentDepth
                    if readletLength < args.readletLen:
                        continue # don't output alignments if the readletLength is below some minimum
                    readlets = readlets + [(nm + ';0', seq[i:i+readletLength], qual[i:i+readletLength]) for i in range(0, seqlen, readletLength)]
                for i in xrange(len(readlets)):
                    nm_rlet, seq_rlet, qual_rlet = readlets[i]
                    rdletStr = "%s;%d;%d;%s\t%s\t%s" % (nm_rlet, i, len(readlets), seq, seq_rlet, qual_rlet)
                    if itninp >= report and i == 0:
                        report *= reportMult
                        if args.verbose: print >> sys.stderr, "First readlet from read %d: '%s'" % (ninp, rdletStr.rstrip())
                    print >>readletfh, rdletStr
            if not itninp: return False # stop bowtying if there are no more readlets to bowtie

    proc, mycmd, threads = bowtie.proc(args, readFn=readletFn,
                                           bowtieArgs=bowtieArgs, sam=True,
                                           stdoutPipe=True, stdinPipe=False)

    unmappedFn = os.path.join(tmpdir, "unmapped_%s_pass_%d.tsv" % (readLabelRoot, passNumber))
    with open(unmappedFn, 'w') as unmappedst:
        bowtieOutDone.clear()
        outThread = OutputThread2(args, proc.stdout, bowtieOutDone, outq, unmappedst, fnh)
        threads.append(outThread)
        outThread.start()
        if args.verbose: print >>sys.stderr, "Waiting for pass #%d of Bowtie to finish" % passNumber
        bowtieOutDone.wait()
        exitlevel = outq.get()
        if exitlevel != 0:
            raise RuntimeError('Bad exitlevel from output thread: %d' % exitlevel)
        for thread in threads:
            if args.verbose: print >> sys.stderr, "  Joining a thread..."
            thread.join()
            if args.verbose: print >> sys.stderr, "    ...joined!"
        sys.stdout.flush()
        if args.verbose: print >>sys.stderr, "Pass #%d of Bowtie is finished" % passNumber

    return unmappedFn

def go2(args, bowtieArgs, tmpdir, fnh, inst):

    import time
    timeSt = time.time()

    tmpdir = tempfile.mkdtemp()
    nextfn = performBowtieIteration(args, bowtieArgs, tmpdir, fnh, inst=inst)
    passNumber = 1
    while nextfn:
        with open(nextfn) as nextst:
            nextfn = performBowtieIteration(args, bowtieArgs, tmpdir, fnh, passNumber=passNumber, inst=nextst)
            passNumber += 1

    print >> sys.stderr, "DONE with align.py; in/out = %d/%d; time=%0.3f secs" % (ninp, nout, time.time()-timeSt)

def go(args, bowtieArgs, tmpdir, fnh, inst, rdlab="reads.tsv", rdletlab="readlets.tsv"):

    import time
    timeSt = time.time()

    bowtieOutDone = threading.Event()

    # For quick test:
    # args.only_readletize = True
    
    # So that all appropriate output directories are at least created, even if
    # they perhaps end up empty
    # if args.exon_differentials: print 'exon_diff\tDUMMY'
    # if args.exon_intervals: print 'exon_ival\tDUMMY'
    # print 'intron\tDUMMY'

    # Instead of using the approach below,
    # rename temporary directory to archive
    '''archiveFh, archiveDir = None, None
    if args.archive is not None:
        archiveDir = os.path.join(args.archive, str(os.getpid()))
        if args.verbose:
            print >> sys.stderr, "Putting --archive reads and command in '%s'" % archiveDir
        if not os.path.exists(archiveDir):
            os.makedirs(archiveDir)
        if not args.only_readletize:
            archiveFh = open(os.path.join(archiveDir, rdlab), 'w')'''

    outq = Queue()

    if not args.only_readletize:
        threads = []

        # Since reads are written to a file before being read
        # in Bowtie's second pass anyway, ALWAYS run serially
        # if args.serial:
        # Reads are written to a file, then Bowtie reads them from the file
        # import tempfile
        # if args.write_reads is None:
        readFn = os.path.join(tmpdir, rdlab)
        #else:
        #    readFn = args.write_reads
        with open(readFn, 'w') as fh:
            # if args.archive is not None: fhs.append(archiveFh)
            writeReads(args, fh, bowtieOutDone, inst=inst, firstPass=True)
        assert os.path.exists(readFn)
        proc, mycmd, threads = bowtie.proc(args, readFn=readFn,
                                           bowtieArgs=bowtieArgs, sam=True,
                                           stdoutPipe=True, stdinPipe=False)
        outThread = OutputThread(args, proc.stdout, bowtieOutDone, outq, tmpdir, fnh)
        threads.append(outThread)
        outThread.start()
        '''else:
            # Reads are written to Bowtie process's stdin directly
            proc, mycmd, threads = bowtie.proc(args, readFn=None,
                                           bowtieArgs=bowtieArgs, sam=True,
                                           stdoutPipe=True, stdinPipe=True)
            outThread = OutputThread(proc.stdout, bowtieOutDone, outq)
            threads.append(outThread)
            outThread.start()
            fhs = [proc.stdin]
            if args.archive is not None: fhs.append(archiveFh)
            writeReads(fhs, firstPass=True)
            proc.stdin.close()'''

        if args.archive is not None:
            with open(os.path.join(tmpdir, "bowtie_cmd%s.sh" % ("" if args.readletLen <= 0 else "_pass1")), 'w') as ofh:
                ofh.write(mycmd + '\n')
            with open(os.path.join(tmpdir, "align_py_cmd%s.sh" % ("" if args.readletLen <= 0 else "_pass1")), 'w') as ofh:
                ofh.write(' '.join(sys.argv) + '\n')
        if args.verbose: 
            print >>sys.stderr, "Waiting for Bowtie%s to finish" % ("" if args.readletLen <= 0 else "'s first pass")
        bowtieOutDone.wait()

        exitlevel = outq.get()
        if exitlevel != 0:
            raise RuntimeError('Bad exitlevel from output thread: %d' % exitlevel)
        for thread in threads:
            if args.verbose: print >> sys.stderr, "  Joining a thread..."
            thread.join()
            if args.verbose: print >> sys.stderr, "    ...joined!"
        sys.stdout.flush()
        if args.verbose:
            print >>sys.stderr, "Bowtie%s finished" % ("" if args.readletLen <= 0 else "'s first pass")

        # Remove temporary read files from first pass of Bowtie
        '''if args.write_reads is None and not args.keep_reads:
            print >>sys.stderr, "Cleaning up temporary files"
            shutil.rmtree(tmpdir)'''

        bowtieOutDone.clear()

    if args.readletLen > 0:
        # Perform Bowtie's second pass only
        # if readletizing is requested
        #if args.archive is not None:
        #    archiveFh = open(os.path.join(archiveDir, rdletlab), 'w')

        # outq = Queue()
        threads = []

        # if args.serial:
        # Reads are written to a file, then Bowtie reads them from the file
        # if args.write_reads is None:
        readFn = os.path.join(tmpdir, rdletlab)
        #else:
        #    readFn = args.write_reads
        with open(readFn, 'w') as fh:
            # fhs = [fh]
            # if args.archive is not None: fhs.append(archiveFh)
            writeReads(args, fh, bowtieOutDone, tmpdir=tmpdir, firstPass=False)
        assert os.path.exists(readFn)
        proc, mycmd, threads = bowtie.proc(args, readFn=readFn,
                                           bowtieArgs=bowtieArgs, sam=True,
                                           stdoutPipe=True, stdinPipe=False)
        outThread = OutputThread(args, proc.stdout, bowtieOutDone, outq, tmpdir, fnh, firstPass=False)
        threads.append(outThread)
        outThread.start()
        '''else:
            # Reads are written to Bowtie process's stdin directly
            proc, mycmd, threads = bowtie.proc(args, readFn=None,
                                           bowtieArgs=bowtieArgs, sam=True,
                                           stdoutPipe=True, stdinPipe=True)
            outThread = OutputThread(proc.stdout, bowtieOutDone, outq)
            threads.append(outThread)
            outThread.start()
            fhs = [proc.stdin]
            if args.archive is not None: fhs.append(archiveFh)
            writeReads(fhs, firstPass=False)
            proc.stdin.close()'''

        if args.archive is not None:
            archiveFh.close()
            with open(os.path.join(archiveDir, "bowtie_cmd%s.sh" % ("" if args.only_readletize else "_pass2")), 'w') as ofh:
                ofh.write(mycmd + '\n')
            with open(os.path.join(archiveDir, "align_py_cmd%s.sh" % ("" if args.only_readletize else "_pass2")), 'w') as ofh:
                ofh.write(' '.join(sys.argv) + '\n')
        if args.verbose: print >>sys.stderr, "Waiting for Bowtie%s to finish" % ("" if args.only_readletize else "'s second pass")
        bowtieOutDone.wait()

        exitlevel = outq.get()
        if exitlevel != 0:
            raise RuntimeError('Bad exitlevel from output thread: %d' % exitlevel)
        for thread in threads:
            if args.verbose: print >> sys.stderr, "  Joining a thread..."
            thread.join()
            if args.verbose: print >> sys.stderr, "    ...joined!"
        sys.stdout.flush()
        if args.verbose:
            print >>sys.stderr, "Bowtie%s finished" % ("" if args.only_readletize else "'s second pass")

        print >> sys.stderr, "DONE with align.py; in/out = %d/%d; time=%0.3f secs" % (ninp, nout, time.time()-timeSt)
        # If 
        # Remove temporary read files from first and second passes of Bowtie

if not args.test:
    binsz = partition.binSize(args)
    fnh = fasta.fasta(args.refseq)
    tmpdir = tempfile.mkdtemp()
    if not os.path.exists(args.refseq):
        raise RuntimeError("No such --refseq file: '%s'" % args.refseq)
    if not os.path.exists(args.faidx):
        raise RuntimeError("No such --faidx file: '%s'" % args.faidx)
    if args.profile:
        import cProfile
        cProfile.run('go(args, bowtieArgs, tmpdir, fnh, sys.stdin)')
    else:
        go(args, bowtieArgs, tmpdir, fnh, sys.stdin) # Replace this with go2() to test alternative readletizing scheme
        if args.archive is not None:
            archiveDir = os.path.join(args.archive, str(os.getpid()))
            os.rename(tmpdir, archiveDir)
        #elif args.write_reads is None and not args.keep_reads:
        else:
            print >>sys.stderr, "Cleaning up temporary files"
            import shutil
            shutil.rmtree(tmpdir)
else:
    del sys.argv[1:]
    import unittest
    binsz = 10000

    # Only used for testing
    def createTestFasta(fname,refid,refseq):
        fastaH = open(fname,'w')
        fastaIdx = open(fname+".fai",'w')
        fastaH.write(">%s\n%s\n"%(refid,refseq))
        fastaIdx.write("%s\t%d\t%d\t%d\t%d\n"%(refid,len(refseq),len(refid)+2,len(refseq),len(refseq)+1))
        fastaH.close()
        fastaIdx.close()

    class TestAlignFunctions1(unittest.TestCase):
        ###Big Note:  We are going to assume base-0 indexing for everything
        def setUp(self):
            #A visual representation of the reference sequence and the read
            #Test1
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
            testargs = ["--refseq", self.fasta, "--exon-differentials"]
            self.args = parser.parse_args(testargs)
            createTestFasta(self.fasta,"test",self.refseq)
            open(self.testDump,'w') #Just to initialize file

        def tearDown(self):
            os.remove(self.fasta)
            os.remove(self.faidx)
            os.remove(self.testDump)

        def test_short_alignment1(self):
            sys.stdout = open(self.testDump,'w')
            rdid,fw,refid = "0;LB:test",True,"test"
            iSt,iEnd = 10,19  #intron coords
            rSt,rEnd = 10,19  #region coords
            offset = 10
            fnh = fasta.fasta(self.fasta)
            handleShortAlignment(self.args,refid,iSt,iEnd,self.rdseq,
                                 rSt,rEnd,rdid,fw,fnh,0,True)
            sys.stdout.close()
            test_out = open(self.testDump,'r')
            testLine = test_out.readline().rstrip()
            toks = testLine.split("\t")
            st,end = int(toks[2]), int(toks[3])
            self.assertEquals(st,10)
            self.assertEquals(end,19)

        def test_correct_splice1(self):
            left = "TTACGAAGGTTTGTA"
            right= "TAATTTAGATGGAGA"
            read = "TTACGAAGATGGAGA"
            # left = "AGTATCGAACCTGAAGCAAGTTACGAAGGTTTGTATAACAAAAATTATGTGAAAG"
            # right= "TAATATTTTCTTTTGAAATTTAATTTAGATGGAGAAATGGAAGCAGAGTGGCTAG"
            # read = "AGTATCGAACCTGAAGCAAGTTACGAAGATGGAGAAATGGAAGCAGAGTGGCTAG"
            fw = True
            _, c, _, _, _, _ = correctSplice(read,left,right,fw)
            assert left[:c]+right[c:] == read

        def test_correct_splice2(self):
            read = "ACGATAACCTTTTTT"
            left = "ACGATAACCTGAGTC"
            right= "TGGACAACCTTTTTT"
            fw = True
            _,c,_,_,_,_ = correctSplice(read,left,right,fw)
            assert left[:c]+right[c:] == read

        def test1Scenario1(self):
            sys.stdout = open(self.testDump,'w')
            rdid,fw,refid = "0;LB:test",True,"test"
            #leftSt,leftEnd = 21,37  #left coords
            #rightSt,rightEnd = 143,154  #left coords
            iSt,iEnd = 32,135 #intron coords
            rSt,rEnd = 32,32  #region coords
            offset = 10
            fnh = fasta.fasta(self.fasta)
            handleIntron(self.args, refid,iSt,iEnd,self.rdseq,
                         rSt,rEnd,rdid,fw,fnh,offset)
            sys.stdout.close()
            test_out = open(self.testDump,'r')
            testLine = test_out.readline().rstrip()
            toks = testLine.split("\t")
            st,end,rdid = int(toks[2]), int(toks[3]), toks[4]
            self.assertEquals(st,32)
            self.assertEquals(end,135)

        """
        Scenario 2: Unmapped region
                                       ^   ^ (Unmapped region)
        Read   |=======================-----======================|
                                      /     \
        Genome |======================-------=====================|
        Flanks                   ^===^       ^===^
        """
        def test1Scenario2(self):
            sys.stdout = open(self.testDump,'w')
            rdid,fw,refid = "0;LB:test",True,"test"
            #leftSt,leftEnd = 21,37  #left coords
            #rightSt,rightEnd = 143,154  #left coords

            iSt,iEnd = 28,139 #intron coords
            rSt,rEnd = 28,36  #region coords
            offset = 10
            fnh = fasta.fasta(self.fasta)
            handleIntron(self.args, refid,iSt,iEnd,self.rdseq,
                         rSt,rEnd,rdid,fw,fnh,offset)
            sys.stdout.close()
            test_out = open(self.testDump,'r')
            testLine = test_out.readline().rstrip()
            toks = testLine.split("\t")
            st,end,rdid = int(toks[2]), int(toks[3]), toks[4]
            self.assertTrue( abs(st-32)<4 )
            self.assertTrue( abs(end-135)<4 )

        """
            Scenario 3: Overlapping flanking sequences - flanking sequences will overlap in the original read
                                      ^^
        Read   |=============================================|
                                      /      \
        Genome |======================-------=====================|
        Flanks                      ^===^  ^===^
        """
        def test1Scenario3(self):
            sys.stdout = open(self.testDump,'w')
            rdid,fw,refid = "0;LB:test",True,"test"

            iSt,iEnd = 36,131 #intron coords
            rSt,rEnd = 36,28  #region coords
            offset = 10
            fnh = fasta.fasta(self.fasta)
            handleIntron(self.args, refid,iSt,iEnd,self.rdseq,
                         rSt,rEnd,rdid,fw,fnh,offset)
            sys.stdout.close()
            test_out = open(self.testDump,'r')
            testLine = test_out.readline().rstrip()
            toks = testLine.split("\t")
            st,end,rdid = int(toks[2]), int(toks[3]), toks[4]
            self.assertTrue( abs(st-32)<4 )
            self.assertTrue( abs(end-135)<4 )

    class TestAlignFunctions2(unittest.TestCase):
        ###Big Note:  We are going to assume base-0 indexing for everything
        def setUp(self):
            #A visual representation of the reference sequence and the read

            #Test2
            """Read"""
            """ACGAAGGACT GCTTGACATC GGCCACGATA ACAACCTTTT TTGCGCCAAT CTTAAGAGCC TTCT"""
            #             ^10        ^20        ^30        ^40        ^50        ^60
            """Genome"""
            """ACGAAGGACT GCTTGACATC GGCCAAAAAA AACTGAGTCG ATAGGACGAA ACAAGTATAT ATTCGAAAAT TAATTAATTC CGAAATTTCA ATTTCATCCG ACATGTATCT ACATATGCCA CACTTCTGGT TGGACTTTTT TTTTTGCGCC A"""
            """ACGAAGGACT GCTTGACATC GGCCAAAAAA AA                                                                                                                 TTTTT TTTTTGCGCC AATCTTAAGA GCCTTCT"""
            #             ^10        ^20        ^30        ^40        ^50        ^60        ^70        ^80        ^90        ^100       ^110       ^120       ^130       ^140       ^150

            self.rdseq  = "ACGAAGGACTGCTTGACATCGGCCAAAAAAAATTTTTTTTTTGCGCCAATCTTAAGAGCCTTCT"
            self.refseq = "ACGAAGGACTGCTTGACATCGGCCAAAAAAAACTGAGTCGATAGGACGAAACAAGTATATATTCGAAAATTAATTAATTCCGAAATTTCAATTTCATCCGACATGTATCTACATATGCCACACTTCTGGTTGGACTTTTTTTTTTGCGCCA"
            self.testDump = "test.out"
            self.fasta = "test.fa"
            self.faidx = "test.fa.fai"
            testargs = ["--refseq", self.fasta, "--exon-differentials"]
            self.args = parser.parse_args(testargs)
            createTestFasta(self.fasta,"test",self.refseq)
            open(self.testDump,'w') #Just to initialize file
        def tearDown(self):
            os.remove(self.fasta)
            os.remove(self.faidx)
            os.remove(self.testDump)


        """
        Scenario 2: Unmapped region
                                       ^   ^ (Unmapped region)
        Read   |=======================-----======================|
                                      /     \
        Genome |======================-------=====================|
        Flanks                   ^===^       ^===^
        """
        def test2Scenario2(self):
            sys.stdout = open(self.testDump,'w')
            rdid,fw,refid = "0;LB:test",True,"test"
            #leftSt,leftEnd = 21,37  #left coords
            #rightSt,rightEnd = 143,154  #left coords

            iSt,iEnd = 28,139 #intron coords
            rSt,rEnd = 28,36  #region coords
            offset = 10
            fnh = fasta.fasta(self.fasta)
            handleIntron(self.args, refid,iSt,iEnd,self.rdseq,
                         rSt,rEnd,rdid,fw,fnh,offset)
            sys.stdout.close()
            test_out = open(self.testDump,'r')
            testLine = test_out.readline().rstrip()
            toks = testLine.split("\t")
            st,end,rdid = int(toks[2]), int(toks[3]), toks[4]

            self.assertEquals( st,32 )
            self.assertEquals( end,135 )

        """
            Scenario 3: Overlapping flanking sequences - flanking sequences will overlap in the original read
                                      ^^
        Read   |=============================================|
                                      /      \
        Genome |======================-------=====================|
        Flanks                      ^===^  ^===^
        """
        def test2Scenario3(self):
            sys.stdout = open(self.testDump,'w')
            rdid,fw,refid = "0;LB:test",True,"test"
            #leftSt,leftEnd = 21,37  #left coords
            #rightSt,rightEnd = 143,154  #left coords

            iSt,iEnd = 36,131 #intron coords
            rSt,rEnd = 36,28  #region coords
            offset = 10
            fnh = fasta.fasta(self.fasta)
            handleIntron(self.args, refid,iSt,iEnd,self.rdseq,
                         rSt,rEnd,rdid,fw,fnh,offset)
            sys.stdout.close()
            test_out = open(self.testDump,'r')
            testLine = test_out.readline().rstrip()
            toks = testLine.split("\t")
            st,end,rdid = int(toks[2]), int(toks[3]), toks[4]

            self.assertEquals( st,32 )
            self.assertEquals( end,135 )
            self.assertEquals( st,32 )
            self.assertEquals( end,135 )

    class TestAlignFunctions3(unittest.TestCase):
        ###Big Note:  We are going to assume base-0 indexing for everything
        def setUp(self):
            #A visual representation of the reference sequence and the read

            #Test2
            """Read"""
            """ACGAAGGACT GCTTGACATC GGCCACGATA ACAACCTTTT TTGCGCCAAT CTTAAGAGCC TTCT"""
            #             ^10        ^20        ^30        ^40        ^50        ^60
            """Genome"""
            """ACGAAGGACT GCTTGACATC GGCCAAAAAA AACTGAGTCG ATAGGACGAA ACAAGTATAT ATTCGAAAAT TAATTAATTC CGAAATTTCA ATTTCATCCG ACATGTATCT ACATATGCCA CACTTCTGGT TGGACTTTTT TTTTTGCGCC A"""
            """ACGAAGGACT GCTTGACATC GGCCAAAAAA AA                                                                                                                 TTTTT TTTTTGCGCC AATCTTAAGA GCCTTCT"""
            #             ^10        ^20        ^30        ^40        ^50        ^60        ^70        ^80        ^90        ^100       ^110       ^120       ^130       ^140       ^150

            self.rdseq  = "ACGATGGACTGCTTGACTCGGCCAAAAAAAATTTTTTTTTTGCGCCAATCTTAAGAGCCTTCT"
            self.refseq = "ACGAAGGACTGCTTGACATCGGCCAAAAAAAACTGAGTCGATAGGACGAAACAAGTATATATTCGAAAATTAATTAATTCCGAAATTTCAATTTCATCCGACATGTATCTACATATGCCACACTTCTGGTTGGACTTTTTTTTTTGCGCCA"
            self.testDump = "test.out"
            self.fasta = "test.fa"
            self.faidx = "test.fa.fai"
            testargs = ["--refseq", self.fasta, "--exon-differentials"]
            self.args = parser.parse_args(testargs)
            createTestFasta(self.fasta,"test",self.refseq)
            open(self.testDump,'w') #Just to initialize file
        def tearDown(self):
            os.remove(self.fasta)
            os.remove(self.faidx)
            os.remove(self.testDump)


        def test_short_alignment1(self):
            """Tests mismatches: ~ score=7"""
            sys.stdout = open(self.testDump,'w')
            rdid,fw,refid = "0;LB:test",True,"test"
            iSt,iEnd = 0,9  #intron coords
            rSt,rEnd = 0,9  #region coords
            offset = 10
            fnh = fasta.fasta(self.fasta)
            handleShortAlignment(self.args,refid,iSt,iEnd,self.rdseq,
                                 rSt,rEnd,rdid,fw,fnh,0,True)
            sys.stdout.close()
            test_out = open(self.testDump,'r')
            testLine = test_out.readline().rstrip()
            toks = testLine.split("\t")
            st,end = int(toks[2]), int(toks[3])
            self.assertEquals(st,0)
            self.assertEquals(end,9)

        def test_short_alignment2(self):
            """Tests mismatches: ~ score=5"""
            sys.stdout = open(self.testDump,'w')
            rdid,fw,refid = "0;LB:test",True,"test"
            iSt,iEnd = 10,19  #intron coords
            rSt,rEnd = 10,19  #region coords
            offset = 10
            fnh = fasta.fasta(self.fasta)
            handleShortAlignment(self.args,refid,iSt,iEnd,self.rdseq,
                                 rSt,rEnd,rdid,fw,fnh,0,True)
            sys.stdout.close()
            test_out = open(self.testDump,'r')
            testLine = test_out.readline().rstrip()
            toks = testLine.split("\t")
            st,end = int(toks[2]), int(toks[3])
            self.assertEquals(st,10)
            self.assertEquals(end,19)

    class TestExactMatch(unittest.TestCase):
        ###Big Note:  We are going to assume base-0 indexing for everything
        def setUp(self):

            # TestExactMatch -- Does an exact match get aligned on Bowtie's first pass?
            # -- Other tests on exact matches can be written as test methods like test_first_pass_bowtie
            """Read"""
            """ACGAAGGACT GCTTGACATC GGCCAAAAAA AACTGAGTCG ATAGGACGAA ACAAGTATAT ATT"""
            #             ^10        ^20        ^30        ^40        ^50        ^60
            """Genome"""
            """ACGAAGGACT GCTTGACATC GGCCAAAAAA AACTGAGTCG ATAGGACGAA ACAAGTATAT ATTCGAAAAT TAATTAATTC CGAAATTTCA ATTTCATCCG ACATGTATCT ACATATGCCA CACTTCTGGT TGGACTTTTT TTTTTGCGCC A"""
            """ACGAAGGACT GCTTGACATC GGCCAAAAAA AA                                                                                                                 TTTTT TTTTTGCGCC AATCTTAAGA GCCTTCT"""
            #             ^10        ^20        ^30        ^40        ^50        ^60        ^70        ^80        ^90        ^100       ^110       ^120       ^130       ^140       ^150

            rdseq  = "ACGAAGGACTGCTTGACATCGGCCAAAAAAAACTGAGTCGATAGGACGAAACAAGTATATATT"
            refseq = "ACGAAGGACTGCTTGACATCGGCCAAAAAAAACTGAGTCGATAGGACGAAACAAGTATATATTCGAAAATTAATTAATTCCGAAATTTCAATTTCATCCGACATGTATCTACATATGCCACACTTCTGGTTGGACTTTTTTTTTTGCGCCA"
            self.tmpdir = tempfile.mkdtemp()
            self.testDump = os.path.join(self.tmpdir, "test.out")
            fastaf = os.path.join(self.tmpdir, "test.fa")
            faidx = os.path.join(self.tmpdir, "test.fa.fai")
            tabinp = os.path.join(self.tmpdir, "test.tsv")
            self.rdlab = "reads.tsv"
            self.rdletlab = "readlets.tsv"
            idxroot = os.path.join(self.tmpdir, "testgenome")
            bowtieExe = "bowtie"
            bowtieBuildExe = "bowtie-build"
            readletLen = 15 # Note readlet length is less than 25 bp!
            readletIval = 5
            partitionLen = 30000
            testargs = ["--refseq", fastaf, "--bowtieIdx", idxroot, "--bowtieExe", bowtieExe, "--readletLen", \
                str(readletLen), "--readletIval", str(readletIval), "--partition-len", str(partitionLen), "--exon-differentials"]
            self.args = parser.parse_args(testargs)
            createTestFasta(fastaf,"test",refseq)
            self.fnh = fasta.fasta(fastaf)
            # write 3-token tab-delimited file
            with open(tabinp, 'w') as oh:
                oh.write("FH:44d2170ded;LB:testlabel;0\t%s\t%s\n" % (rdseq, "~"*len(rdseq))) # make base quality scores as high as possible
            self.inst = open(tabinp, 'r')
            # create Bowtie index files from reference fasta
            import subprocess
            with open(os.devnull, 'w') as nullst:
                proc = subprocess.Popen([bowtieBuildExe, fastaf, idxroot], stdout=nullst)
                proc.wait()
            open(self.testDump,'w') # Just to initialize file

        def tearDown(self):
            # Remove all temporary files
            import shutil
            shutil.rmtree(self.tmpdir)

        def testFirstPassBowtie(self):
            # An unusually large unit test
            # Essentially mirrors go(), but archiving is killed
            # May serve as prototype for better-organized go()
            # Lines with args.archive and args.verbose from go() are not included

            sys.stdout = open(self.testDump, 'w')

            go(self.args, [], self.tmpdir, self.fnh, inst=self.inst, rdlab=self.rdlab, rdletlab=self.rdletlab)

            sys.stdout.close()

            # reads.tsv should have data (this file catches first-pass Bowtie output),
            # while readlets.tsv should be empty

            # reads5.tab is just os.path.join(self.tmpdir, self.rdlab)

            self.assertNotEqual(os.stat(os.path.join(self.tmpdir, self.rdlab)).st_size, 0)
            self.assertEqual(os.stat(os.path.join(self.tmpdir, self.rdletlab)).st_size, 0)

    class Test30bpIntron(unittest.TestCase):
        ###Big Note:  We are going to assume base-0 indexing for everything
        def setUp(self):

            # TestExactMatch -- Does an exact match get aligned on Bowtie's first pass?
            # -- Other tests on exact matches can be written as test methods like test_first_pass_bowtie
            """Read"""
            """ACGAAGGACT GCTTGACATC                                  ACAAGTATAT ATTCGAAAAT TAATTAATTC"""
            #             ^10        ^20        ^30        ^40        ^50        ^60
            """Genome"""
            """ACGAAGGACT GCTTGACATC GGCCAAAAAA AACTGAGTCG ATAGGACGAA ACAAGTATAT ATTCGAAAAT TAATTAATTC CGAAATTTCA ATTTCATCCG ACATGTATCT ACATATGCCA CACTTCTGGT TGGACTTTTT TTTTTGCGCC A"""
            """ACGAAGGACT GCTTGACATC GGCCAAAAAA AA                                                                                                                 TTTTT TTTTTGCGCC AATCTTAAGA GCCTTCT"""
            #             ^10        ^20        ^30        ^40        ^50        ^60        ^70        ^80        ^90        ^100       ^110       ^120       ^130       ^140       ^150

            rdseq  = "ACGAAGGACTGCTTGACATCACAAGTATATATTCGAAAATTAATTAATTC"
            refseq = "ACGAAGGACTGCTTGACATCGGCCAAAAAAAACTGAGTCGATAGGACGAAACAAGTATATATTCGAAAATTAATTAATTCCGAAATTTCAATTTCATCCGACATGTATCTACATATGCCACACTTCTGGTTGGACTTTTTTTTTTGCGCCA"
            self.tmpdir = tempfile.mkdtemp()
            self.testDump = os.path.join(self.tmpdir, "test.out")
            fastaf = os.path.join(self.tmpdir, "test.fa")
            faidx = os.path.join(self.tmpdir, "test.fa.fai")
            tabinp = os.path.join(self.tmpdir, "test.tsv")
            self.rdlab = "reads.tsv"
            self.rdletlab = "readlets.tsv"
            idxroot = os.path.join(self.tmpdir, "testgenome")
            bowtieExe = "bowtie"
            bowtieBuildExe = "bowtie-build"
            readletLen = 15 # Note readlet length is less than 25 bp here!
            readletIval = 5
            partitionLen = 30000
            testargs = ["--refseq", fastaf, "--bowtieIdx", idxroot, "--bowtieExe", bowtieExe, "--readletLen", \
                str(readletLen), "--readletIval", str(readletIval), "--partition-len", str(partitionLen), "--exon-differentials"]
            self.args = parser.parse_args(testargs)
            createTestFasta(fastaf,"test",refseq)
            self.fnh = fasta.fasta(fastaf)
            # write 3-token tab-delimited file
            with open(tabinp, 'w') as oh:
                oh.write("FH:44d2170ded;LB:testlabel;0\t%s\t%s\n" % (rdseq, "~"*len(rdseq))) # make base quality scores as high as possible
            self.inst = open(tabinp, 'r')
            # create Bowtie index files from reference fasta
            import subprocess
            with open(os.devnull, 'w') as nullst:
                proc = subprocess.Popen([bowtieBuildExe, fastaf, idxroot], stdout=nullst)
                proc.wait()
            open(self.testDump,'w') # Just to initialize file

        def tearDown(self):
            # Remove all temporary files
            import shutil
            shutil.rmtree(self.tmpdir)

        def testSecondPassBowtie(self):
            # An unusually large unit test
            # Essentially mirrors go(), but archiving is killed
            # May serve as prototype for better-organized go()
            # Lines with args.archive and args.verbose from go() are not included
            sys.stdout = open(self.testDump, 'w')

            go(self.args, [], self.tmpdir, self.fnh, inst=self.inst, rdlab=self.rdlab, rdletlab=self.rdletlab)

            sys.stdout.close()

            # reads.tsv should have data (this file catches first-pass Bowtie output),
            # and readlets.tsv should ALSO have data (since there are readlets here too)
            self.assertNotEqual(os.stat(os.path.join(self.tmpdir, self.rdlab)).st_size, 0)
            self.assertNotEqual(os.stat(os.path.join(self.tmpdir, self.rdletlab)).st_size, 0)

    class TestMultimappedRead(unittest.TestCase):
        ###Big Note:  We are going to assume base-0 indexing for everything
        def setUp(self):

            # TestExactMatch -- Does an exact match get aligned on Bowtie's first pass?
            # -- Other tests on exact matches can be written as test methods like test_first_pass_bowtie
            """Read"""
            """ACGAAGGACT GCTTGACATC GGCCAAAAAA AACTGAGTCG ATAGGACGAA ACAAGTATAT """
            #             ^10        ^20        ^30        ^40        ^50        
            """Genome"""
            """ACGAAGGACT GCTTGACATC GGCCAAAAAA AACTGAGTCG ATAGGACGAA ACAAGTATAT ACGAAGGACT GCTTGACATC GGCCAAAAAA AACTGAGTCG ATAGGACGAA ACAAGTATAT ACGAAGGACT GCTTGACATC GGCCAAAAAA AACTGAGTCG ATAGGACGAA ACAAGTATAT"""
            #             ^10        ^20        ^30        ^40        ^50        ^60        ^70        ^80        ^90        ^100       ^110       ^120       ^130       ^140       ^150       ^160       ^170

            rdseq  = "ACGAAGGACTGCTTGACATCGGCCAAAAAAAACTGAGTCGATAGGACGAAACAAGTATATATT"
            refseq = "ACGAAGGACTGCTTGACATCGGCCAAAAAAAACTGAGTCGATAGGACGAAACAAGTATATATTACGAAGGACTGCTTGACATCGGCCAAAAAAAACTGAGTCGATAGGACGAAACAAGTATATATTACGAAGGACTGCTTGACATCGGCCAAAAAAAACTGAGTCGATAGGACGAAACAAGTATATATT"
            self.tmpdir = tempfile.mkdtemp()
            self.testDump = os.path.join(self.tmpdir, "test.out")
            fastaf = os.path.join(self.tmpdir, "test.fa")
            faidx = os.path.join(self.tmpdir, "test.fa.fai")
            tabinp = os.path.join(self.tmpdir, "test.tsv")
            self.rdlab = "reads.tsv"
            self.rdletlab = "readlets.tsv"
            idxroot = os.path.join(self.tmpdir, "testgenome")
            bowtieExe = "bowtie"
            bowtieBuildExe = "bowtie-build"
            readletLen = 15 # Note readlet length is less than 25 bp here!
            readletIval = 5
            partitionLen = 30000
            testargs = ["--refseq", fastaf, "--bowtieIdx", idxroot, "--bowtieExe", bowtieExe, "--readletLen", \
                str(readletLen), "--readletIval", str(readletIval), "--partition-len", str(partitionLen), "--exon-differentials"]
            self.args = parser.parse_args(testargs)
            createTestFasta(fastaf,"test",refseq)
            self.fnh = fasta.fasta(fastaf)
            # write 3-token tab-delimited file
            with open(tabinp, 'w') as oh:
                oh.write("FH:44d2170ded;LB:testlabel;0\t%s\t%s\n" % (rdseq, "~"*len(rdseq))) # make base quality scores as high as possible
            self.inst = open(tabinp, 'r')
            # create Bowtie index files from reference fasta
            import subprocess
            with open(os.devnull, 'w') as nullst:
                proc = subprocess.Popen([bowtieBuildExe, fastaf, idxroot], stdout=nullst)
                proc.wait()
            open(self.testDump,'w') # Just to initialize file

        def tearDown(self):
            # Remove all temporary files
            import shutil
            shutil.rmtree(self.tmpdir)

        def testMultimappedBowtie(self):
            # An unusually large unit test
            # Essentially mirrors go(), but archiving is killed
            # May serve as prototype for better-organized go()
            # Lines with args.archive and args.verbose from go() are not included
            sys.stdout = open(self.testDump, 'w')

            go(self.args, [], self.tmpdir, self.fnh, inst=self.inst, rdlab=self.rdlab, rdletlab=self.rdletlab)

            sys.stdout.close()

            # reads.tsv should have caught alignments, but readlets.tsv should not
            self.assertNotEqual(os.stat(os.path.join(self.tmpdir, self.rdlab)).st_size, 0)
            self.assertEqual(os.stat(os.path.join(self.tmpdir, self.rdletlab)).st_size, 0)

    unittest.main()

