"""
gtf.py
"""

import string
import argparse
import re

parser = argparse.ArgumentParser(description='Parse a GTF file.')

parser.add_argument(\
    '--fasta', metavar='path', type=str, nargs='+',
    help='FASTA file(s) containing reference genome sequences')
parser.add_argument(\
    '--gtf', metavar='path', type=str, nargs='+', required=True,
    help='GTF file(s) containing gene annotations')
parser.add_argument(\
    '--introns-out', metavar='PATH', type=str, required=False,
    help='File to write intron information to')
parser.add_argument(\
    '--exons-out', metavar='PATH', type=str, required=False,
    help='File to write exon information to')

args = parser.parse_args()

gene_id_re = re.compile("gene_id \"([^\"]+)\"")
xscript_id_re = re.compile("transcript_id \"([^\"]+)\"")

class Annot(object):
    def __init__(self, refid, st0, en0, orient, feature, score, frame, attrs):
        self.refid = refid
        self.st0 = st0 # 0-based
        self.en0 = en0 # 0-based, exclusive
        self.orient = orient
        self.feature = feature
        self.score = score
        self.frame = frame
        self.attrs = attrs

class Transcript(object):
    """ Right now I ignore all associations with TSSs and promoters """
    def __init__(self, exons, start=None, stop=None, tssId=None):
        self.exons = exons # list of exon Annots
        for exon in exons:
            assert exon.feature == "exon"
        self.start = start # start coden Annot
        self.stop = stop
        self.tssId = tssId

for gtf in args.gtf:
    with open(gtf, 'r') as gtfFh:
        last_gene_id = None
        xscripts = {}
        for line in gtfFh:
            gene_id, xscript_id = None, None
            refid, source, feature, start1, end1, score, orient, frame, attr = string.split(line, '\t')
            if feature == "CDS":
                continue
            if feature == "start_codon":
                continue
            if feature == "stop_codon":
                continue
            assert feature == "exon"
            start1, end1 = int(start1), int(end1)
            assert end1 >= start1
            annot = Annot(refid, start1-1, end1)
            ma = gene_id_re.search(attr)
            if ma is None:
                raise RuntimeError("No gene_id!:\n" + line)
            gene_id = ma.group(1)
            ma = xscript_id_re.search(attr)
            if ma is None:
                raise RuntimeError("No transcript_id!:\n" + line)
            xscript_id = ma.group(1)
            if gene_id != last_gene_id:
                last_gene_id = gene_id
