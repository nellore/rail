"""
get_junctions.py
Abhi Nellore / July 15, 2014

Extracts junctions from transcripts in GTF file (read from stdin). Considers
only 'exon' features. Tailored for Gencode v12 annotation available at
ftp://ftp.sanger.ac.uk/pub/gencode/release_12/gencode.v12.annotation.gtf.gz.

Writes output to stdout exclusively for STAR. Each line of output represents
an intron and takes the form
<rname><TAB><1-based start position><TAB><1-based end-position>.

Might take a bit of memory if the GTF is massive.
"""
import sys
from collections import defaultdict

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
            if exons_from_transcript[i][1] - exons_from_transcript[i-1][2] < 5:
                continue
            print '\t'.join([
                        exons_from_transcript[i][0],
                        str(exons_from_transcript[i-1][2] + 1),
                        str(exons_from_transcript[i][1] - 1)
                    ])