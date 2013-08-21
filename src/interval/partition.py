'''
partition.py
'''

import math

def addArgs(parser):
    parser.add_argument(\
        '--ntasks', metavar='NUMTASKS', type=int, required=False,
        help='Number of reduce tasks')
    parser.add_argument(\
        '--genomeLen', metavar='LEN', type=int, required=False,
        help='Total length of the genome; required so that we can accurately calculate bin sizes')
    parser.add_argument(\
        '--partition-length', metavar='LEN', type=int, required=False,
        help='Length of a single genome partition')

def binSize(args):
    if args.partition_length is not None:
        return args.partition_length
    elif args.genomeLen is not None:
        return int(math.ceil(1.0 * args.genomeLen / args.ntasks))
    else:
        return 10000

def partition(refid, st, en, binSize):
    ''' Assign the interval refid:[st, en) to one or more partitions
        based on partition bin size and the interval's start and end
        positions. '''
    binid_st = int(st / binSize)
    binid_en = int((en-1) / binSize)
    return [ (";".join([refid, str(i)]), i * binSize, (i+1) * binSize) for i in xrange(binid_st, binid_en+1) ]

def partitionOverlaps(refid, st, en, binSize, fudge=0):
    ''' Assign the interval refid:[st, en) to one or more partitions based on
        which partitions are overlapped by the interval.  The interval is
        (optionally) extended on both ends by some fudge factor before
        overlapping. '''
    binid_st = int((st - fudge    ) / binSize)
    binid_en = int((en + fudge - 1) / binSize)
    return [ (";".join([refid, str(i)]), i * binSize, (i+1) * binSize) for i in xrange(binid_st, binid_en+1) ]

def partitionStartOverlaps(refid, st, en, binSize, fudge=0):
    ''' Assign the interval refid:[st, en) to one or more partitions based on
        which partitions are overlapped by the interval.  The interval is
        (optionally) extended on both ends by some fudge factor before
        overlapping. '''
    binid_st = int((st - fudge) / binSize)
    binid_en = int((st + fudge) / binSize)
    return [ (";".join([refid, str(i)]), i * binSize, (i+1) * binSize) for i in xrange(binid_st, binid_en+1) ]

def parse(st, binSz):
    ''' Parse a partition id. '''
    toks = st.split(";")
    if len(toks) != 2:
        raise RuntimeError("Expected two tokens separated by ;, got %d: '%s'" % (len(toks), st))
    refid, i = toks[0], int(toks[1])
    b = i * binSz
    return refid, b, b + binSz

if __name__ == '__main__':
    import unittest

    class TestFlatIntervals(unittest.TestCase):

        def test1(self):
            pt = partition("blah", 27, 37, 10)
            self.assertEqual(2, len(pt))
            refid, st, en = parse(pt[0][0], 10)
            self.assertEqual("blah", refid)
            self.assertEqual(20, st)
            self.assertEqual(30, en)
            refid, st, en = parse(pt[1][0], 10)
            self.assertEqual("blah", refid)
            self.assertEqual(30, st)
            self.assertEqual(40, en)
        
        def test2(self):
            ps = partitionOverlaps('p', 99, 100, 100, fudge=0)
            self.assertEquals(1, len(ps))
        
        def test3(self):
            ps = partitionOverlaps('p', 99, 100, 100, fudge=1)
            self.assertEquals(2, len(ps))
        
        def test4(self):
            ps = partitionStartOverlaps('p', 90, 95, 100, fudge=8)
            self.assertEquals(1, len(ps))
        
        def test5(self):
            ps = partitionStartOverlaps('p', 90, 95, 100, fudge=11)
            self.assertEquals(2, len(ps))
        
        def test6(self):
            ps = partitionOverlaps('p', 90, 95, 100, fudge=8)
            self.assertEquals(2, len(ps))
        
    unittest.main()
