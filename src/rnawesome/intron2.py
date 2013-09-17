"""
Tab-delimited input tuple columns:
1. Partition ID for partition overlapped by interval (also includes strand information)
2. Interval start
3. Interval end (exclusive)
4. Sample label
5. Readlet Sequence before 5' site
6. Readlet Sequence after 5' site
7. Readlet Sequence before 3' site
8. Readlet Sequence after 3' site

Tab-delimited splice site output tuple columns:
1. Reference ID
2. 5' start
3. 3' start
4. Sample label
5. Read frequency (number of times sample read overlapped junction)

Tab-delimited cooccurence output tuple columns:
1. rdid
2. refID
3. left_site
4. right_site

Strandedness
============

If the input reads are *not* stranded, then we have to allow for all relevant
splice junction (donor/acceptor) motifs, including e.g. both GT-AG and CT-AC.
If the input reads are stranded, we might wish to ignore motifs that are
inconsistent with the strand of the flanking aligned sequences.  I added the
--stranded option to do this.

"""

import os
import sys
import argparse
import site
from collections import defaultdict
import time
timeSt = time.time()

base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for subdir in [ 'fasta', 'splice', 'read', 'alignment', 'statsmath', 'util', 'interval' ]:
    site.addsitedir(os.path.join(base_path, subdir))

import partition
import chrsizes
import fasta
import readlet
import needlemanWunsch
import histogram
import window
from strip import Strip, cluster
from pick import SiteSelector
import motifs

parser = argparse.ArgumentParser(description='Reports splice junction information')
parser.add_argument(\
    '--test', action='store_const', const=True, default=False, help='Run unit tests')
parser.add_argument(\
    '--verbose', action='store_const', const=True, default=False, help='Print more warnings and informative messages')
parser.add_argument(\
    '--refseq', type=str, required=False,
    help='The fasta sequence of the reference genome. The fasta index of the reference genome is also required')
parser.add_argument(\
    '--cluster-radius', type=int, required=False, default=10,
    help='The radius of the clustering algorithm and the radius of the sliding window algorithm.')
parser.add_argument(\
    '--intron-partition-overlap', type=int, required=False, default=20,
    help='Amount that partitions overlap their left and right neighbors by.')
parser.add_argument(\
    '--per-site', action='store_const', const=True, default=False,
    help='Output one record for every splice site, giving information about number of times a read from each label spans the site')
parser.add_argument(\
    '--per-span', action='store_const', const=True, default=False,
    help='Output one record for every instance where a read spans a splice site')
parser.add_argument(\
    '--stranded', action='store_const', const=True, default=False,
    help='Assume input reads come from the sense strand')

readlet.addArgs(parser)
partition.addArgs(parser)

args = parser.parse_args()

ninp, nout = 0, 0 # # lines input/output so far

def go(ifh, ofh, verbose=False, refseq=None):

    global ninp
    
    starts, ends, labs, rdids = [], [], [], []
    last_pt, last_refid = '\t', '\t'
    last_strand = None
    last_pt_st, last_pt_en = None, None
    fnh = fasta.fasta(refseq)
    binsz = partition.binSize(args)
    fudge = args.intron_partition_overlap
    mat = {}
    sel, selwat, selcri = \
        SiteSelector(fnh, motifs.sa), SiteSelector(fnh, motifs.saf), \
        SiteSelector(fnh, motifs.sar)
    
    def binIntervals(clustering):
        bins = {}
        for st, en, lab, rdid in zip(starts, ends, labs, rdids):
            i, j = en - st, st
            if (i, j) not in clustering.cmap:
                continue
            idx = clustering.cmap[(i, j)]
            if idx in bins:
                bins[idx].append((st, en, lab, rdid))
            else:
                bins[idx] = [(st, en, lab, rdid)]
        return bins
    
    def handlePartition():
        assert last_strand is not None
        assert last_pt != '\t'
        
        if verbose:
            print >> sys.stderr, "For partition %s:[%d, %d)" % (last_pt, last_pt_st, last_pt_en)
        if len(mat) == 0: return
        
        # Step 1. Make a strip
        strip = Strip.fromHistogram(mat)
        nelts = len(strip)
        
        # Step 2. Cluster splice sites within strip
        clustering = cluster(strip, N=args.cluster_radius)
        nclustsPre = len(clustering)
        if verbose:
            print >> sys.stderr, "  %d possible splice junction clustered down to %d" % (nelts, nclustsPre)
        
        clustering.limitTo(last_pt_st, last_pt_en)
        nclustsPost = len(clustering)
        if verbose:
            print >> sys.stderr, "  %d after overlap removal" % (nclustsPost)
        
        # Step 3. Build a bins object
        bins = binIntervals(clustering)
        
        # Step 4: Apply sliding windows to find splice junction locations
        for introns in bins.itervalues():
            sts, ens, lab, rdid = zip(*introns)
            if args.stranded:
                assert last_strand == '-' or last_strand == '+'
                if last_strand == '+':
                    ssites = selwat.handleCluster(last_refid, sts, ens)
                else:
                    ssites = selcri.handleCluster(last_refid, sts, ens)
            else:
                ssites = sel.handleCluster(last_refid, sts, ens)
            if len(ssites) == 0:
                if verbose:
                    print >> sys.stderr, "Warning: A cluster with %d intronic chunks had no splice site" % len(introns)
                continue
            _, _, motif, st, en = ssites[0] # just take highest-scoring
            _, motifl, motifr = motif
            if args.per_span:
                for intron in introns:
                    _, _, l, rdid = intron
                    ofh.write("span\t%s\t%s\t%012d\t%d\t%s\t%s\t%s\n" % (rdid, last_refid, st, en+2, motifl, motifr, l))
            if args.per_site:
                d = defaultdict(int)
                for l in lab: d[l] += 1
                for l, c in d.iteritems():
                    ofh.write("site\t%s\t%012d\t%d\t%s\t%s\t%s\t%d\n" % (last_refid, st, en+2, motifl, motifr, l, c))
    
    for ln in ifh:
        # Parse next read
        toks = ln.rstrip().split('\t')
        assert len(toks) >= 5
        pt, st, en, lab, rdid = \
            toks[0], int(toks[1]), int(toks[2]), toks[3], toks[4]
        assert en > st
        refid, pt_st, pt_en = partition.parse(pt[:-1], binsz)
        strand = pt[-1]
        assert st >= pt_st - fudge and st < pt_en + fudge, \
            "Intron start %d not in partition [%d, %d), partition id=%s" % (st, pt_st, pt_en, pt)
        
        if last_pt != pt and last_pt != '\t':
            handlePartition()
            starts, ends, labs, rdids = [], [], [], []
            mat = {}
        
        i, j = en - st, st
        mat[(i, j)] = mat.get((i, j), 0) + 1
        
        starts.append(st)
        ends.append(en)
        labs.append(lab)
        rdids.append(rdid)
        last_pt, last_strand, last_refid = pt, strand, refid
        last_pt_st, last_pt_en = pt_st, pt_en
        ninp += 1
    
    if last_pt != '\t': handlePartition()
    
    timeEn = time.time()
    print >>sys.stderr, "DONE with intron2.py; in/out = %d/%d; time=%0.3f secs" % (ninp, nout, timeEn-timeSt)

if not args.test:
    go(sys.stdin, sys.stdout, args.verbose, args.refseq)
else:
    del sys.argv[1:]
    import unittest
    from cStringIO import StringIO
    
    class TestIntronFunctions1(unittest.TestCase):
        
        def test1(self):
            fafn, _ = fasta.writeIndexedFasta("ref1", "G" * 100 + "GT" + "T" * 96 + "AG")
            inp = StringIO('\t'.join(['ref1;0', '100', '200', 'A', 'ACGT', 'TGCA', 'read1']))
            out = StringIO()
            go(inp, out, False, fafn)
    
    unittest.main()
