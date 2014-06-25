"""
manifest.py
Part of Rail-RNA

Contains some helpful tools for parsing the file manifest, which lists the
locations of files containing raw reads and their corresponding sample names.
"""

import os
from collections import deque

def add_args(parser):
    parser.add_argument(\
        '--manifest', metavar='PATH', type=str, required=False,
        help='Manifest file')

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

_charset = '0123456789abcdefghijklmnopqrstuvwxyz'

def string_from_int(integer):
    """ Converts NONNEGATIVE integer to base-36 representation.

        To invert this operation, use int(str, 36).
        Based on http://stackoverflow.com/questions/2063425/
        python-elegant-inverse-function-of-intstring-base.

        Used to encode sample indexes.

        Return value: string encoding base-36 representation of integer.
    """
    integer, remainder = divmod(integer, 36)
    to_return = deque()
    while integer:
        to_return.appendleft(_charset[remainder])
        integer, remainder = divmod(integer, 36)
    to_return.appendleft(_charset[remainder])
    return ''.join(to_return)
