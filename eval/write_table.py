#!/user/bin/env python
"""
write_table.py

Writes table of raw coverages of introns from output of grab_introns.py.
"""
import os
import sys

a = {}
samples = set()
index_to_sample_name = {}

for counter, i in enumerate(os.listdir(sys.argv[1])):
	sample_name = i.rpartition('.')[0].partition('.')[2]
	samples.add(sample_name)
	index_to_sample_name[counter] = sample_name
	print >>sys.stderr, i
	with open(os.path.join(sys.argv[1], i)) as f:
		for line in f:
			line = line.strip().split('\t')
			intron = tuple(line[:-1])
			if intron not in a:
				a[intron] = [0]*700
			coverage = int(line[-1])
			a[intron][counter] = coverage

for i in xrange(len(samples)):
	sys.stdout.write('\t' + index_to_sample_name[i])
sys.stdout.write('\n')

for intron in a:
	sys.stdout.write(';'.join([intron[0] + intron[-1], intron[1], intron[2]]))
	for i in xrange(len(samples)):
		sys.stdout.write('\t' + str(a[intron][i]))
	sys.stdout.write('\n')