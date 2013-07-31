"""
tornado_config.py

Parse and organize parameters controlling the behavior of the Tornado pipeline.
"""

import os
import site

base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
site.addsitedir(os.path.join(base_path, "util"))

import path

def addArgs(parser):
    parser.add_argument(\
        '--quality', metavar='STR', type=str, default="phred33", help='Quality string format.')
    parser.add_argument(\
        '--readlet-length', metavar='INT', type=int, default="25", help='Substring length to extract from read.')
    parser.add_argument(\
        '--readlet-interval', metavar='INT', type=int, default="5", help='Distance between substrings to extract from read.')
    parser.add_argument(\
        '--partition-length', metavar='INT', type=int, default=10000, help='Size of genome partitions to use.')
    parser.add_argument(\
        '--cluster-radius', metavar='INT', type=int, default="2", help='For clustering candidate introns into junctions.')
    parser.add_argument(\
        '--bowtie-exe', metavar='STR', type=str, help='Bowtie executable to use.  Must exist at this path on all the cluster nodes.')
    parser.add_argument(\
        '--bowtie-args', metavar='STR', type=str, help='Arguments to pass to Bowtie.')
    parser.add_argument(\
        '--downsample-reads', metavar='FRACTION', type=float, default=1.0, help='Fraction of reads to randomly downsample to.')
    parser.add_argument(\
        '--truncate-reads', metavar='INT', type=int, default=None, help='Truncate reads to INT base pairs.')
    parser.add_argument(\
        '--log-add', metavar='FRACTION', type=float, default=32.0, help='Add this to counts before taking log2 to obtain log count.')
    parser.add_argument(\
        '--normalize-percentile', metavar='FRACTION', type=float, default=0.75, help='Percentile non-zero coverage to use for normalization (0.75).')
    parser.add_argument(\
        '--permutations', metavar='INT', type=int, default=5, help='# permutations to try.')
    parser.add_argument(\
        '--discard-mate1', action='store_const', const=True, help='Discard mate 1s for paired-end reads.')
    parser.add_argument(\
        '--discard-mate2', action='store_const', const=True, help='Discard mate 2s for paired-end reads.')
    parser.add_argument(\
        '--pool-tech-replicates', action='store_const', const=True, help='Pool together technical replicates instead of treating them separately (not recommended).')
    parser.add_argument(\
        '--pool-bio-replicates', action='store_const', const=True, help='Pool together biological and technical replicates instead of treating them separately (not recommended).')
    parser.add_argument(\
        '--hmm-overlap', metavar='INT', type=int, default=100, help='Number of positions of overlap between adjacent emission strings.')

class TornadoConfig(object):
    
    def __init__(self, args):
        q = self.quality = args.quality
        if q != "phred64" and q != "phred33" and q != "solexa64":
            raise RuntimeError("Unknown argument for --quality: '%s'" % q)
        self.bowtieExe = args.bowtie_exe
        l = self.readletLen = args.readlet_length
        if l < 4:
            raise RuntimeError("Argument for --readlet-length must be >= 4; was %d" % l)
        i = self.readletIval = args.readlet_interval
        if i < 1:
            raise RuntimeError("Argument for --readlet-interval must be >= 1; was %d" % i)
        r = self.clusterRadius = args.cluster_radius
        if r < 0:
            raise RuntimeError("Argument for --cluster-radius must be >= 0; was %d" % r)
        p = self.partitionLen = args.partition_length
        if p < 100:
            raise RuntimeError("Argument for --partition-length must be >= 100; was %d" % p)
        self._bowtieArgs = args.bowtie_args or "-v 1 -m"
        d = self.downsampleReads = args.downsample_reads
        if d <= 0.0 or d >= 1.00001:
            raise RuntimeError("Argument for --downsample-reads must be in (0, 1]; was %f" % d)
        t = self.truncateReads = args.truncate_reads
        if t is not None and t < 4:
            raise RuntimeError("Argument for --truncate-reads must be >= 4; was %d" % t)
        l = self.logAdd = args.log_add
        if l <= 0.0:
            raise RuntimeError("Argument for --log-add must be in > 0.0; was %f" % l)
        p = self.normPercentile = args.normalize_percentile
        if p <= 0.0 or p >= 1.0:
            raise RuntimeError("Argument for --normalize-percentile must be in (0, 1); was %f" % p)
        p = self.numPermutations = args.permutations or 5
        if p < 0:
            raise RuntimeError("Argument for --permutations must be >= 0; was %d" % p)
        self.discardMate1 = args.discard_mate1
        self.discardMate2 = args.discard_mate2
        self.poolTech = args.pool_tech_replicates
        self.poolBio = args.pool_bio_replicates
        o = self.hmmOlap = args.hmm_overlap
        if o < 0:
            raise RuntimeError("Argument for --hmm-overlap must be >= 0; was %d" % o)
    
    def bowtieArgs(self):
        return ' '.join([self._bowtieArgs, "--mm", "-t"])
