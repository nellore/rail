"""
manifest.py
Part of Rail-RNA

Contains some helpful tools for parsing the file manifest, which lists the
locations of files containing raw reads and their corresponding sample names.

There is some redundancy here to accommodate legacy code.
"""

import os

def add_args(parser):
    parser.add_argument(\
        '--manifest', metavar='PATH', type=str, required=False,
        help='Manifest file')

def labels(args):
    """ Parse the manifest file and return a list of all the labels in the
        file.

        args: command-line arguments

        Return value: list of all labels in file.
    """
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
                raise RuntimeError('Wrong number of tokens in line: ' + line)
    return labs

class LabelsAndIndices:
    """ Parses the manifest file to create manifest dictionary. """
    def __init__(self, manifest_file):
        self.label_to_index = {}
        self.index_to_label = {}
        with open(manifest_file) as manifest_stream:
            for i, line in enumerate(manifest_stream):
                tokens = line.rstrip().split('\t')
                assert len(tokens) <= 5, ('Manifest file has invalid line: %s' 
                                          % line)
                index_string = str(i)
                self.label_to_index[tokens[-1]] = index_string
                self.index_to_label[index_string] = tokens[-1]