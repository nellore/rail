"""
url.py

For parsing, handling and manipulating URLs which might be on local
filesystem, in HDFS, or in either of the S3 filesystems (native or block).
"""

class Url(object):
    def __init__(self, s):
        if ':' in s:
            cidx = s.index(':')
            preColon = s[:cidx].lower()
            if preColon[:2] == "s3n":
                self.type = "s3n"
            elif preColon[:2] == "s3":
                self.type = "s3"
            elif preColon[:4] == "hdfs":
                self.type = "hdfs"
            elif preColon[:4] == "http":
                self.type = "http"
            elif preColon[:4] == "ftp":
                self.type = "ftp"
            elif preColon[:5] == "local":
                self.type = "local"
            else:
                raise RuntimeError("Unrecognized URL; not S3, HDFS or local: '%s'" % s)
            self.rest = s[cidx+1:]
        else:
            self.type = "local"
            self.rest = s
    
    def isS3(self):
        return self.type[:2] == 's3'
    
    def isWgettable(self):
        return self.type in ['ftp', 'http']
    
    def isNotLocal(self):
        return self.type != 'local'
    
    def toUrl(self):
        if self.type == 'local':
            return self.rest
        else:
            return self.type + ':' + self.rest
    
    def toUpperUrl(self):
        """ Useful for hiding protocol names from Elastic MapReduce so it
            doesn't mistake a URL passed as a mapper argument as an input
            URL. """
        if self.type == 'local':
            return self.rest
        else:
            return self.type.upper() + ':' + self.rest
    
    def toNonNativeUrl(self):
        if self.type[:2] == 's3':
            return 's3:' + self.rest
        else:
            return self.type + ':' + self.rest

