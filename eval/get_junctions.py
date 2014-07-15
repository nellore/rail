"""
get_junctions.py
Abhi Nellore / July 15, 2014

Extracts junctions from transcripts in GTF file (read from stdin). Considers
only 'exon' features. Tailored for Gencode v12 annotation available at
ftp://ftp.sanger.ac.uk/pub/gencode/release_12/gencode.v12.annotation.gtf.gz.

Writes output to stdout exclusively for STAR. Each line of output represents
an intron and takes the form
<rname><TAB><1-based start position><TAB><1-based end-position>.
"""
import sys
from collections import defaultdict

exons = defaultdict(list)
for line in sys.stdin:
	tokens = line.strip().split('\t')