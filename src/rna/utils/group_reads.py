"""
group_reads.py
Part of Rail-RNA

Contains a class for grouping read sequences to be aligned to the same
transcriptome Bowtie 2 index.
"""
import hashlib

def add_args(parser):
    parser.add_argument(
        '--index-count', type=int, required=False,
        default=50,
        help=('Number of transcriptome Bowtie 2 indexes to which reads '
              'are to be assigned.')
    )

class IndexGroup(object):
    """ Includes method for assigning reads to indexes. """

    def __init__(self, index_count):
        self.index_count = index_count

    def index_group(self, seq):
        return ('%012d'
                 % (int(hashlib.md5(seq).hexdigest(), 16) % self.index_count))