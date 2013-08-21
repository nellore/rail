"""
strip.py

One strategy for clustering splice sites.  Think of potential splice sites as
elements in a 2D grid, where x axis is LHS of the splice site and y axis is
RHS.  Because of minimum and maximum intron lengths, only a stretch of
diagonals in the grid's upper triangle contains possible splice sites.  We
subdivide this stretch of diagonals into vertical strips (parallelograms) and
search for clusters within each strip.

I can't think of any clustering scheme s.t. we can look for clusters in
different strips in parallel and get the correct answer without iterating
and/or using unreasonably large partitions.  This implementation gives up on
finding a clustering that is globally correct, with the hope that we can
correct globally inconsistent results later.
"""

from splice_clustering import SpliceClustering

class Strip(object):
    """ This is a vertical strip from the stretch of diagonals that can hold
        legal splice sites.  The strip is represented as a matrix; this is a
        little tricky because the strip is really a parallelogram that we've
        "collapsed" into a rectangle.  If the lower-left corner of the
        parallelogram is at (I, J), element (i, j) of the rectangle is element
        (i+j, j) of the parallelogram.
        
          o
         oo
        ooo => ooo
        oo     ooo
        o      ooo
        
    """
    
    def __init__(self, elts, nrow, ncol):
        """ elts is a list of (count, i, j) tuples """
        self.nrow, self.ncol = nrow, ncol
        self.elts = elts
        self.elts.sort(reverse=True)
        self.eltset = set()
        for _, i, j in self.elts:
            self.eltset.add((i, j))
    
    @classmethod
    def fromDenseMatrix(cls, mat):
        """ Create a strip from a dense 2D matrix with elts = counts """
        nrow, ncol = len(mat), len(mat[0])
        elts = []
        for i in xrange(0, nrow):
            for j in xrange(0, ncol):
                if mat[i][j] > 0:
                    elts.append((mat[i][j], i, j))
        return cls(elts, nrow, ncol)
    
    @classmethod
    def fromHistogram(cls, hist, nrow, ncol):
        """ Create a strip from a histogram of (st, en): count """
        elts = [ (cnt, c[0], c[1]) for c, cnt in hist.iteritems() ]
        return cls(elts, nrow, ncol)

def cluster(strip, N=5):
    assert strip.nrow > 0
    assert strip.ncol > 0
    clusters = []
    covered = set()
    for _, i, j in strip.elts:
        if (i, j) not in covered:
            covered.add((i, j))
            clusters.append([(i, j)])
            for j2 in xrange(max(j-N, 0), min(j+N+1, strip.ncol)):
                if (i, j2) in strip.eltset and (i, j2) not in covered:
                    clusters[-1].append((i, j2))
                    covered.add((i, j2))
            assert len(covered) <= len(strip.elts)
            if len(covered) == len(strip.elts):
                return SpliceClustering.fromClusterList(clusters)
    raise RuntimeError("Should not get here!")

def clusterExtra(strip, npad, N=5):
    """ Handle a strip where 'npad' columns on the left and right sides are
        """
    clusters = cluster(strip, N=N)
    newclusters = []
    assert strip.ncol > 2*npad
    newncol = strip.ncol - 2*npad
    for clust in clusters.clist:
        newclusters.append([])
        leftmost = strip.ncol
        for i, j in clust:
            leftmost = min(leftmost, j-npad)
            newclusters[-1].append((i, j-npad))
        if len(newclusters[-1]) == 0 or leftmost < 0 or leftmost >= newncol:
            newclusters.pop()
    return SpliceClustering.fromClusterList(newclusters)

if __name__ == '__main__':
    import unittest

    class TestCircularCountBuffer(unittest.TestCase):
        
        def test1(self):
            mat = [[0, 0, 0, 0, 1, 0],
                   [0, 0, 0, 0, 0, 0],
                   [0, 0, 0, 2, 0, 0],
                   [0, 0, 0, 0, 0, 0],
                   [0, 0, 3, 0, 3, 0],
                   [0, 0, 0, 0, 0, 0],
                   [0, 0, 0, 0, 0, 0],
                   [4, 4, 0, 0, 0, 0]]
            strip = Strip.fromDenseMatrix(mat)
            self.assertEqual(6, len(strip.elts))
            self.assertEqual(4, strip.elts[0][0])
            self.assertEqual(4, strip.elts[1][0])
        
        def test2(self):
            mat = [[0, 0, 0, 0, 1, 0],
                   [0, 0, 0, 0, 0, 0],
                   [0, 0, 0, 2, 0, 0],
                   [0, 0, 0, 0, 0, 0],
                   [0, 0, 3, 0, 3, 0],
                   [0, 0, 0, 0, 0, 0],
                   [0, 0, 0, 0, 0, 0],
                   [4, 4, 0, 0, 0, 0]]
            clusters = cluster(Strip.fromDenseMatrix(mat), 1)
            self.assertEqual(5, len(clusters))
        
        def test3(self):
            mat = [[0, 0, 0, 0, 1, 0],
                   [0, 0, 0, 0, 0, 0],
                   [0, 0, 0, 2, 0, 0],
                   [0, 0, 0, 0, 0, 0],
                   [0, 0, 3, 0, 3, 0],
                   [0, 0, 0, 0, 0, 0],
                   [0, 0, 0, 0, 0, 0],
                   [4, 4, 0, 0, 0, 0]]
            clusters = cluster(Strip.fromDenseMatrix(mat), 2)
            self.assertEqual(4, len(clusters))
        
        def test4(self):
            mat = [[0, 0, 0, 0, 1, 0],
                   [0, 0, 0, 2, 0, 0],
                   [0, 0, 3, 0, 3, 0],
                   [5, 4, 4, 4, 0, 0]]
            clusters = cluster(Strip.fromDenseMatrix(mat), 2)
            clist = clusters.clist
            self.assertEqual(5, len(clist))
            self.assertTrue([(3, 0), (3, 1), (3, 2)] in clist)
            self.assertTrue([(3, 3)] in clist)
            self.assertTrue([(2, 4), (2, 2)] in clist)
            self.assertTrue([(1, 3)] in clist)
            self.assertTrue([(0, 4)] in clist)
        
        def test5(self):
            mat = [ (1, 0, 4), (2, 1, 3), (3, 2, 2), (3, 2, 4),
                    (5, 3, 0), (4, 3, 1), (4, 3, 2), (4, 3, 3) ]
            clusters = cluster(Strip(mat, 4, 6), 2)
            clist = clusters.clist
            self.assertEqual(5, len(clist))
            self.assertTrue([(3, 0), (3, 1), (3, 2)] in clist)
            self.assertTrue([(3, 3)] in clist)
            self.assertTrue([(2, 4), (2, 2)] in clist)
            self.assertTrue([(1, 3)] in clist)
            self.assertTrue([(0, 4)] in clist)
        
        def test6(self):
            mat = [ (1, 0, 4), (2, 1, 3), (3, 2, 2), (3, 2, 4),
                    (5, 3, 0), (4, 3, 1), (4, 3, 2), (4, 3, 3) ]
            clusters = clusterExtra(Strip(mat, 4, 6), 1, N=2)
            clist = clusters.clist
            self.assertEqual(4, len(clist))
            self.assertTrue([(3, 2)] in clist)
            self.assertTrue([(2, 3), (2, 1)] in clist)
            self.assertTrue([(1, 2)] in clist)
            self.assertTrue([(0, 3)] in clist)
        
        def test7(self):
            mat = [ (1, 0, 4), (2, 1, 3), (3, 2, 2), (3, 2, 4),
                    (5, 3, 0), (4, 3, 1), (4, 3, 2), (4, 3, 3) ]
            clusters = clusterExtra(Strip(mat, 4, 6), 2, N=2)
            clist = clusters.clist
            self.assertEqual(3, len(clist))
            self.assertTrue([(2, 2), (2, 0)] in clist)
            self.assertTrue([(1, 1)] in clist)
            self.assertTrue([(3, 1)] in clist)
        
        def test8(self):
            hist = { (0, 4) : 1, (1, 3) : 2, (2, 2) : 3, (2, 4) : 3,
                     (3, 0) : 5, (3, 1) : 4, (3, 2) : 4, (3, 3) : 4 }
            clusters = clusterExtra(Strip.fromHistogram(hist, 4, 6), 2, N=2)
            clist = clusters.clist
            self.assertEqual(3, len(clist))
            self.assertTrue([(2, 2), (2, 0)] in clist)
            self.assertTrue([(1, 1)] in clist)
            self.assertTrue([(3, 1)] in clist)
    
    unittest.main()
