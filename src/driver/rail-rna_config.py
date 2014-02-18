"""
rail-rna_config.py

Parse and organize parameters controlling the behavior of the Rail-RNA pipeline.
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
        '--readlet-interval', metavar='INT', type=int, default="4", help='Distance between substrings to extract from read.')
    parser.add_argument(\
        '--capping-fraction', metavar='FRACTION', type=float, default="0.85", help='Successive capping readlets on a given end of a read are tapered in size exponentially with this fractional base.')
    parser.add_argument(\
        '--partition-length', metavar='INT', type=int, default=30000, help='Size of genome partitions to use.')
    parser.add_argument(\
        '--cluster-radius', metavar='INT', type=int, default="50", help='For clustering candidate introns into junctions.')
    parser.add_argument(\
        '--intron-partition-overlap', metavar='INT', type=int, default="20", help='# of nucleotides of overlap between intron-finding partitions.')
    parser.add_argument(\
        '--bowtie-exe', metavar='STR', type=str, help='Bowtie executable to use. Must exist at this path on all the cluster nodes.')
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
        '--stranded', action='store_const', const=True, help='RNA-seq data is stranded?')
    parser.add_argument(\
        '--hmm-overlap', metavar='INT', type=int, default=100, help='Number of positions of overlap between adjacent emission strings.')
    parser.add_argument(\
        '--output-bam-by-chromosome', action='store_const', const=True, default=True,
        help='Split SAM/BAM output files up by RNAME (typically chromosome) if '
             'True; otherwise output single SAM/BAM file consolidating all '
             'spliced alignments')
    parser.add_argument(\
        '--bam-basename', type=str, required=False, 
        default='spliced_alignments',
        help='The basename (including path) of all SAM/BAM output. Basename is '
             'followed by ".bam" or ".sam" if --by-chromosome=False; otherwise, '
             'basename is followed by ".RNAME.bam" or ".RNAME.sam" for each '
             'RNAME. Ignored if --out is not specified '
             '(that is, if --out is stdout)')
    parser.add_argument(\
        '--output-bed-by-chromosome', action='store_const', const=True, default=False,
        help='Split BED output files up by chrom (typically chromosome) if '
             'True; otherwise output single BED file consolidating all '
             'spliced alignments')
    parser.add_argument(\
        '--bed-basename', type=str, required=False, 
        default='junctions',
        help='The basename (including path) of all BED output. Basename is '
             'followed by ".bed" if --by-chromosome=False; otherwise, basename is '
             'followed by ".[chrom].bed" for each [chrom]. Ignored if --out is '
             'not specified (that is, if --out is stdout)')
    parser.add_argument(\
        '--samtools-exe', metavar='EXE', type=str, required=False,
        help='Path to executable for samtools. Must exist at this path on all the cluster nodes.')
    parser.add_argument(\
        '--output-sam', action='store_const', const=True, default=False, 
        help='Output SAM files if True; otherwise output BAM files')
    parser.add_argument('--do-not-search_for_caps',
        action='store_const',
        const=True,
        default=False,
        help='Ordinarily, reference is search for the segment of a read (a '
             'cap) that precedes the first EC and the cap that follows the '
             'last EC. Such caps are subsequently added as ECs themselves. '
             'Use this command-line parameter to turn the feature off')
    parser.add_argument('--min-cap-query-size', type=int, required=False,
        default=8,
        help='The reference is not searched for a segment of a read that '
             'precedes the first EC or follows the last EC smaller than this '
             'size')
    parser.add_argument('--cap-search-window-size', type=int, required=False,
        default=1000,
        help='The size (in bp) of the reference subsequence in which to '
             'search for a cap --- i.e., a segment of a read that follows '
             'the last EC or precedes the first EC.')

class Rail_RNAConfig(object):
    
    def __init__(self, args):
        q = self.quality = args.quality
        if q != "phred64" and q != "phred33" and q != "solexa64":
            raise RuntimeError("Unknown argument for --quality: '%s'" % q)
        self.bowtieExe = args.bowtie_exe
        l = self.readletLen = args.readlet_length
        if l < 4:
            raise RuntimeError("Argument for --readlet-length must be >= 4; was %d" % l)
        c = self.capping_fraction = args.capping_fraction
        if not 0 <= c < 1:
            raise RuntimeError("Argument for --capping-fraction must be on [0, 1); was %d" % c)
        i = self.readletIval = args.readlet_interval
        if i < 1:
            raise RuntimeError("Argument for --readlet-interval must be >= 1; was %d" % i)
        r = self.clusterRadius = args.cluster_radius
        if r < 0:
            raise RuntimeError("Argument for --cluster-radius must be >= 0; was %d" % r)
        i = self.intronPartitionOlap = args.intron_partition_overlap
        if i < 0:
            raise RuntimeError("Argument for --intron-partition-overlap must be >= 0; was %d" % i)
        p = self.partitionLen = args.partition_length
        if p < 100:
            raise RuntimeError("Argument for --partition-length must be >= 100; was %d" % p)
        self._bowtieArgs = args.bowtie_args or "-v 1 -a -m 40"
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
        self.stranded = args.stranded
        self.samtoolsExe = args.samtools_exe
        self.output_bed_by_chromosome = args.output_bed_by_chromosome
        self.bed_basename = args.bed_basename
        self.output_bam_by_chromosome = args.output_bam_by_chromosome
        self.bam_basename = args.bam_basename
        self.output_sam = args.output_sam
        self.poolTech = args.pool_tech_replicates
        self.poolBio = args.pool_bio_replicates
        self.do_not_search_for_caps = args.do_not_search_for_caps
        if args.min_cap_query_size < 0:
            raise RuntimeError("Argument for --min-cap-query-size must be in > 0; was %d" % args.min_cap_query_size)
        self.min_cap_query_size = args.min_cap_query_size
        if args.cap_search_window_size < 0:
            raise RuntimeError("Argument for --cap-search-window-size must be in > 0; was %d" % args.cap_search_window_size)
        self.cap_search_window_size = args.cap_search_window_size
        o = self.hmmOlap = args.hmm_overlap
        if o < 0:
            raise RuntimeError("Argument for --hmm-overlap must be >= 0; was %d" % o)
    
    def bowtieArgs(self):
        return ' '.join([self._bowtieArgs, "-t", "--sam-nohead", "--startverbose"])
