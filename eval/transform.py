"""
transform.py

This script just changes the URLs from a GEUVADIS manifest file to local resources.
"""

import sys

for line in sys.stdin:
	tokens = line.strip().split('\t')
	if not (len(tokens) == 5):
		sys.stdout.write(line)
		continue
	for i in [0, 2]:
		tokens[i] = '/scratch0/langmead-fs1/data/big_public_datasets/geuvadis/' + tokens[i].rpartition('/')[2]
	print '\t'.join(tokens)

sys.stdout.flush()