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
