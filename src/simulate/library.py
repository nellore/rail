"""
Parser/container for flux simulator .lib files
"""

base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
site.addsitedir(os.path.join(base_path, "struct"))
site.addsitedir(os.path.join(base_path, "annotation"))

import search
import defaultdict
import gtf

class fragment(object):
    def __init__(self,st1,en1,fragid,copies):
        self.st1=st1
        self.en1=en1
        self.fragid = fragid
        self.copies = copies

class library(object):
    def __init__(self,fname,xscripts):
        self.fragments = defaultdict( list )
        self.fname = fname
        self.junctions = []
        self.findAnnotatedJunctions(xscripts)
        
    def findAnnotatedJunctions(self,xscripts):
        """
        Finds all splice sites that are spanned by a library fragment
        """
        with open(fname,'r') as fh:
            for ln in fh:
                toks = ln.split("\t")
                assert len(toks)==4
                st1,en1,fragid,copies = int(toks[0]), int(toks[1]), toks[2], int(toks[3])
                refid = fragid.split("\t")[0]
                xid = fragid.split("@")[1]
                sites = xscripts[xid].getOverlapSites(st1,en1)
                self.junctions += sites


    def find(self,refid,guess_site):
        result = search.find_tuple(self.junctions[refid],guess_site)
        return fragment(result[0],result[1],result[2],
