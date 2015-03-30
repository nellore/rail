"""
splice_clustering.py

Encapsulates a clustering of splice sites in (st, en) space.
"""

class SpliceClustering(object):
    """ A clustering of splice sites in (st, en) space """
    
    def __init__(self, clist, cmap):
        self.clist, self.cmap = clist, cmap
    
    def __len__(self):
        return len(self.clist)
    
    @staticmethod
    def listToMap(clist):
        mp = {}
        for i in xrange(len(clist)):
            for elt in clist[i]:
                assert elt not in mp
                mp[elt] = i
        return mp
    
    @staticmethod
    def mapToList(cmap):
        assert len(cmap) > 0
        nelt = max(cmap.itervalues()) + 1
        ls = [ [] for _ in xrange(nelt) ]
        for elt, clusti in cmap.iteritems():
            ls[clusti].append(elt)
        return ls
    
    @classmethod
    def fromClusterList(cls, clist):
        cmap = SpliceClustering.listToMap(clist)
        return cls(clist, cmap)
    
    @classmethod
    def fromClusterMap(cls, cmap):
        clist = SpliceClustering.mapToList(cmap)
        return cls(clist, cmap)
    
    def limitTo(self, sti, stf):
        """ Modify the clustering so that it only contains clusters whose
            leftmost element is within the given range of start positions """
        newclist = []
        for clust in self.clist:
            leftmost = min([ j for _, j in clust ])
            if leftmost >= sti and leftmost < stf:
                newclist.append(clust)
        self.cmap = SpliceClustering.listToMap(newclist)
        self.clist = newclist

if __name__ == '__main__':
    import unittest
    from strip import Strip, cluster

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
            clustering = cluster(strip, N=2)
            self.assertEqual(4, len(clustering))
            clustering.limitTo(1, 5)
            self.assertEqual(3, len(clustering))
            self.assertIn([(4, 4), (4, 2)], clustering.clist)
    
    unittest.main()
