"""
validate.py

Reads in the bed file containing the estimated splice sites and the pickle file containing the transcripts and provides the following statistics
1.  Number of exactly correct splice junctions
2.  Number of splice junctions within 1 radius (specified by user)
3.  Number of splice junctions completely off
4.  Plot of the distribution of error
"""
import re
import os
import site
import argparse
import sys
import math
import pickle
import bisect
import copy
from collections import Counter
from collections import defaultdict
base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
site.addsitedir(os.path.join(base_path, "annotation"))
site.addsitedir(os.path.join(base_path, "struct"))
site.addsitedir(os.path.join(base_path, "fasta"))
site.addsitedir(os.path.join(base_path, "statsmath"))
site.addsitedir(os.path.join(base_path, "util"))

import gtf
import search
import fasta
import window
import display
import counter
parser = argparse.ArgumentParser(description=\
                                     'Splice junction validator')
parser.add_argument(\
    '--bed-file', metavar='path', type=str, required=False, default=""
    help='Path of the estimated splice sites bed file')
parser.add_argument(\
    '--radius', type=int, required=False,default=10,
    help='The radius of tolerance for identifying splice site neighborhoods')
parser.add_argument(\
    '--window-radius', type=int, required=False,default=50,
    help='The radius of display window')
parser.add_argument(\
    '--refseq', type=str, required=False, default=""
    help='The reference sequence')
parser.add_argument(\
    '--lib-file', type=str, required=False, default=""
    help='The library file containing all of the correct positions of the fragments')
parser.add_argument(\
    '--profile', action='store_const', const=True, default=False,
    help='Profile simulation generation')
parser.add_argument(\
    '--test', action='store_const', const=True, default=False,
    help='Run unittests')


display.addArgs(parser)

args = parser.parse_args()


def go():

    #Step 1: Isolate all read in .lib file that span splice junctions.  These are the annotated splice junctions
    
    #Step 2: Compare detected splice junctions to annotated splice junctions

    #Step 3: Output two files: false positive regions and false negative regions

    
if __name__=="__main__":
    if args.profile:
        import cProfile
        cProfile.run('go()')
    else:
        go_flux()
