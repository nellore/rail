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
            i = 0
            for line in manifest_stream:
                line = line.strip()
                if line[0] == '#' or not line: continue
                tokens = line.split('\t')
                assert len(tokens) <= 5, ('Manifest file has invalid line: %s' 
                                          % line)
                index_string = str(i)
                self.label_to_index[tokens[-1]] = index_string
                self.index_to_label[index_string] = tokens[-1]
                i += 1

class LabelsAndIndicesWithClusters:
    """ Parses the manifest file to create manifest dictionary. """
    def __init__(self, manifest_file, by='group'):
        """ Parses sample file to map indexes to samples and sample clusters.

            A sample cluster is specified in 'by' as either 'group' or
            'biorep'.

            manifest_file: Myrna-style manifest file with the hash twist
                mentioned above.
            by: 'group' or 'biorep'
        """
        self.label_to_index = {}
        self.index_to_label = {}
        self.index_to_cluster = {}
        self.cluster_to_indices = defaultdict(int)
        assert by in ['group', 'biorep']
        with open(manifest_file) as manifest_stream:
            i = 0
            for line in manifest_stream:
                line = line.strip()
                if line[0] == '#' or not line: continue
                tokens = line.split('\t')
                assert len(tokens) <= 5, ('Manifest file has invalid line: %s' 
                                          % line)
                index_string = str(i)
                self.label_to_index[tokens[-1]] = index_string
                self.index_to_label[index_string] = tokens[-1]
                divisions = tokens[-1].split('-')
                assert len(divisions) == 3, (('Manifest file has invalid '
                                              'sample label: %s') % tokens[-1])
                if by == 'group':
                    self.index_to_cluster[i] = divisions[0]
                else:
                    # by == 'biorep'
                    self.index_to_cluster[i] = '-'.join(divisions[:2])
                i += 1
        for index in self.index_to_cluster:
            self.cluster_to_indices[self.index_to_cluster[index]] += 2**index

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
