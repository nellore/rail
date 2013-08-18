'''
manifest.py
'''

import os

def addArgs(parser):
    parser.add_argument(\
        '--manifest', metavar='PATH', type=str, required=False, help='Manifest file')

def labels(args):
    """ Parse the manifest file and return a list of all the labels in the
        file """
    labs = []
    if not os.path.isfile(args.manifest):
        raise RuntimeError("No such --manifest file '%s'" % args.manifest)
    with open(args.manifest, 'r') as fh:
        for line in fh:
            line = line.strip()
            if len(line) == 0 or line[0] == '#':
                continue
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
