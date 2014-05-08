"""
partition.py

Contains a generator for partitioning genome into bins.
"""

def addArgs(parser):
    parser.add_argument(\
        '--partition-length', metavar='LEN', type=int, required=False,
        help='Length of a single genome partition')

def partition(rname, pos, end_pos, bin_size):
    """ Assigns the interval rname:[pos, end_pos) to one or more partitions.

        rname: RNAME on which interval lies
        pos: start position of interval (inclusive) AND 1-BASED
        end_pos: end position if interval (exclusive) AND 1-BASED
        bin_size: number of bases spanned by partition

        Yield value: Tuple (rname + ';' + partition number starting at 0,
                                start position of partition (1-BASED),
                                end position of partition (1-BASED))
    """
    first_bin = (pos - 1) / bin_size
    last_bin = (end_pos - 2) / bin_size
    for bin_number in xrange(first_bin, last_bin + 1):
        bin_pos = bin_number * bin_size + 1
        bin_end_pos = bin_pos + bin_size
        yield ';'.join([rname, str(bin_number)]), bin_pos, bin_end_pos

if __name__ == '__main__':
    import unittest

    class TestPartition(unittest.TestCase):
        """ Tests partition(); needs not fixture. """

        def test_interval_within_partition(self):
            """ Fails if first partition not returned. """
            self.assertEquals(list(partition('chr1', 27, 5001, 5000)),
                            [('chr1;0', 1, 5001)])

        def test_interval_across_three_partitions(self):
            """ Fails if three partitions are not returned. """
            self.assertEquals(list(partition('chr1', 20, 10002, 5000)),
                            [('chr1;0', 1, 5001), ('chr1;1', 5001, 10001),
                             ('chr1;2', 10001, 15001)])
        
    unittest.main()