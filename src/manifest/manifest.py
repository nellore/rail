'''
manifest.py
'''

import os

def addArgs(parser):
    parser.add_argument(\
        '--manifest', metavar='FILE', type=str, required=False,
         help='Manifest file')

def labels(args):
    labs = []
    if not os.path.isfile(args.manifest):
        raise RuntimeError("No such --manifest file '%s'" % args.manifest)
    fh = open(args.manifest)
    for line in fh:
        line = line.rstrip()
        toks = line.split('\t')
        if len(toks) == 3:
            _, _, lab = toks
            labs.append(lab)
        elif len(toks) == 5:
            _, _, _, _, lab = toks
            labs.append(lab)
        else:
            raise RuntimeError("Wrong number of tokens in line: " + line)
    return labs
