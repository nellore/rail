"""
get_junctions.py
Abhi Nellore / July 15, 2014

Extracts junctions from transcripts in GTF file (read from stdin). Considers
only 'exon' features. Tailored for Gencode v12 annotation available at
ftp://ftp.sanger.ac.uk/pub/gencode/release_12/gencode.v12.annotation.gtf.gz.

Writes output to stdout exclusively for STAR. Each line of output represents
an intron and takes the form
<rname><TAB><1-based start position><TAB><1-based end-position>.

Also writes stderr output with intron length distribution and stats on splice
site motifs if Bowtie index is specified.

Might take a bit of memory if the GTF is massive.
"""
import sys
from collections import defaultdict

if __name__ == '__main__':
    import argparse
    # Print file's docstring if -h is invoked
    parser = argparse.ArgumentParser(description=__doc__, 
                formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--bowtie2-idx', type=str, required=False,
            default=None,
            help=('Path to Bowtie 2 genome index used to determine '
                  'whether junctions are canonical. Canonical/'
                  'noncanonical status and intron length '
                  'distribution are output only if this '
                  'is specified.')
        )
    args = parser.parse_args()
    if args.bowtie2_idx is not None:
        from count_introns import BowtieIndexReference
        reference_index = BowtieIndexReference(args.bowtie2_idx)
        intron_lengths = defaultdict(int)
        motif_counts = defaultdict(int)
    exons = defaultdict(set)
    for line in sys.stdin:
        if line[0] == '#': continue
        tokens = line.strip().split('\t')
        if tokens[2].lower() != 'exon': continue
        '''key: transcript_id
           value: (rname, exon start (1-based), exon end (1-based))

        transcript_id in token 12 is decked with " on the left and "; on the
        right; kill them in key below.
        '''
        attribute = tokens[-1].split(';')
        id_index = [i for i, name in enumerate(attribute)
                    if 'transcript_id' in name]
        assert len(id_index) == 1, ('More than one transcript ID specified; '
                                    'offending line: %s') % line 
        id_index = id_index[0]
        attribute[id_index] = attribute[id_index].strip()
        quote_index = attribute[id_index].index('"')
        exons[attribute[id_index][quote_index+1:-1]].add(
                (tokens[0], int(tokens[3]), int(tokens[4]))
            )

    for transcript_id in exons:
        exons_from_transcript = sorted(list(exons[transcript_id]))
        # Recall that GTF is end-inclusive, and so is STAR's junctions.txt
        for i in xrange(1, len(exons_from_transcript)):
            if exons_from_transcript[i][0] == exons_from_transcript[i-1][0]:
                # Kill any introns 4 bases or smaller
                if (exons_from_transcript[i][1]
                    - exons_from_transcript[i-1][2]) < 5:
                    continue
                intron = (exons_from_transcript[i][0],
                            exons_from_transcript[i-1][2] + 1,
                            exons_from_transcript[i][1] - 1)
                print '\t'.join(intron)
                if args.bowtie2_idx is not None:
                    length = intron[2] - intron[1] + 1
                    motif = (reference_index.get_stretch(intron[0],
                                                            intron[1] - 1,
                                                            2),
                             reference_index.get_stretch(intron[0],
                                                            intron[1]
                                                            + length - 2,
                                                            2))
                    intron_lengths[length] += 1
                    motif_counts[motif_count] += 1
    if args.bowtie2_idx is not None:
        for length, frequency in sorted(intron_lengths.items()):
            print >>sys.stderr, '%d\t%d' % (length, frequency)
        for motif, frequency in motif_counts:
            print >>sys.stderr, '%s\t%s' % (motif, frequency)