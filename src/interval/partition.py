'''
partition.py
'''

import math

def addArgs(parser):
    parser.add_argument(\
        '--ntasks', metavar='NUMTASKS', type=int, required=True,
        help='Number of reduce tasks')
    parser.add_argument(\
        '--genomeLen', metavar='LEN', type=int, required=True,
        help='Total length of the genome; required so that we can accurately calculate bin sizes')

def binSize(args):
    return int(math.ceil(1.0 * args.genomeLen / args.ntasks))

def partition(refid, st, en, binSize):
    ''' Assign the interval refid:[st, en) to one or more partitions
        based on partition bin size and the interval's start and end
        positions. '''
    binid_st = int(st / binSize)
    binid_en = int((en-1) / binSize)
    return [ ";".join([refid, str(i)]) for i in xrange(binid_st, binid_en+1) ]

def parse(st, binSz):
    ''' Parse a partition id. '''
    toks = st.split(";")
    if len(toks) != 2:
        raise RuntimeError("Expected two tokens separated by underscore, got %d: '%s'" % (len(toks), st))
    refid, i = toks[0], int(toks[1])
    b = i * binSz
    return refid, b, b + binSz

if __name__ == '__main__':
    import unittest

    class TestFlatIntervals(unittest.TestCase):

        def test1(self):
            pt = partition("blah", 27, 37, 10)
            self.assertEqual(2, len(pt))
            refid, st, en = parse(pt[0], 10)
            self.assertEqual("blah", refid)
            self.assertEqual(20, st)
            self.assertEqual(30, en)
            refid, st, en = parse(pt[1], 10)
            self.assertEqual("blah", refid)
            self.assertEqual(30, st)
            self.assertEqual(40, en)

    unittest.main()
