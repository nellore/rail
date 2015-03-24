"""
motifs.py

Catalog of the most common splice motifs.

Turunen JJ, Niemela EH, Verma B, Frilander MJ. The significant other: splicing
by the minor spliceosome. Wiley Interdiscip Rev RNA. 2013 Jan-Feb;4(1):61-76.
doi: 10.1002/wrna.1141.

http://onlinelibrary.wiley.com/doi/10.1002/wrna.1141/full

Thanaraj TA, Clark F. Human GC-AG alternative intron isoforms with weak donor 
sites show enhanced consensus at acceptor exon positions. Nucleic Acids Res. 2001
Jun 15;29(12):2581-93


"""

s1f = (0, 'GT', 'AG')
s1r = (0, 'CT', 'AC')
s1 = (s1f, s1r)

s2f = (1, 'GC', 'AG')
s2r = (1, 'CT', 'GC')
s2 = (s2f, s2r)

s3f = (2, 'AT', 'AC')
s3r = (2, 'GT', 'AT')
s3 = (s3f, s3r)

saf = [s1f, s2f, s3f]
sar = [s1r, s2r, s3r]
sa = sorted(saf + sar)
