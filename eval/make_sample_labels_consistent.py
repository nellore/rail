"""
make_sample_labels_consistent.py

Makes sample labels from 56-/28- manifest files agree
with sample labels from 112-sample manifest file.
"""

import argparse
import sys

# Print file's docstring if -h is invoked
parser = argparse.ArgumentParser(description=__doc__, 
            formatter_class=argparse.RawDescriptionHelpFormatter)
# Add command-line arguments
parser.add_argument(
        '--samples', type=str, required=False,
        default='GEUVADIS_112.manifest',
        help='Manifest file with 112 samples'
    )
parser.add_argument(
        '--subset', type=str, required=False,
        default='GEUVADIS_56.manifest',
        help='Manifest file with subset of 112 samples'
    )

args = parser.parse_args()

mappings = {}
with open(args.samples) as source_stream:
    for line in source_stream:
        if (not line.strip()) or line[0] == '#': continue
        tokens = line.strip().split('\t')
        mappings[tokens[0]] = tokens[-1]

with open(args.subset) as subset_stream:
    for line in subset_stream:
        if (not line.strip()) or line[0] == '#':
            sys.stdout.write(line)
            continue
        tokens = line.strip().split('\t')
        print '\t'.join(tokens[:-1] + [mappings[tokens[0]]])

sys.stdout.flush()