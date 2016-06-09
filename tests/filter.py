#!/usr/bin/env python
"""
filter.py

When obtaining a coverage vector from BAM, bedtools only either:
  1) treats Ns and Ds in CIGAR strings as uncovered (with -split)
  2) treats Ns and Ds in CIGAR strings as covered

By default, however, Rail treats Ns as uncovered but Ds as covered. So to test
whether Rail's default output is consistent with the output of bedtools
operating on BAM, one must first transform the BAM so deletions are interpreted
as matches. This is what the --deletions-to-matches option does.
--uniques filters out multimapping reads.

SAM is read from stdin; processed SAM is written to stdout.
"""

import sys
import re

def is_int(value):
    try:
      int(value)
    except (ValueError, TypeError):
      return False
    return True

if __name__ == '__main__':
    import argparse
    # Print file's docstring if -h is invoked
    parser = argparse.ArgumentParser(description=__doc__, 
                formatter_class=argparse.RawDescriptionHelpFormatter)
    # Add command-line arguments
    parser.add_argument('--deletions-to-matches', action='store_const',
            const=True, default=False,
            help='transforms SAM so deletions are interpreted as matches'
        )
    parser.add_argument('--uniques', action='store_const', const=True,
            default=False,
            help='filters out multimapping reads'
        )
    for line in sys.stdin:
        if line[0] == '@':
            sys.stdout.write(line)
            continue
        tokens = line.strip().split('\t')
        if args.uniques and 'NH:i:1' not in tokens: continue
        if args.deletions_to_matches:
            to_add = 0
            split_cigar = re.split(r'([MINDS])', tokens[5])[:-1]
            for i, k in enumerate(split_cigar):
                if k == 'D':
                    to_add += int(split_cigar[i-1])
            tokens[5] = tokens[5].replace('D', 'M')
            tokens[9] = tokens[9] + 'A'*to_add
            tokens[10] = tokens[10] + '#'*to_add
        print '\t'.join(tokens)
