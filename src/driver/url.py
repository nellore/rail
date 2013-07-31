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
            elif preColon[:5] == "local":
                self.type = "local"
            else:
                raise RuntimeError("Unrecognized URL; not S3, HDFS or local: '%s'" % s)
            self.rest = s[cidx+1:]
        else:
            self.type = "local"
            self.rest = s
    
    def isS3(self):
        return self.type[:2] == "s3"
    
    def toUrl(self):
        if self.type == "local":
            return self.rest
        else:
            return self.type + ':' + self.rest
    
    def toUpperUrl(self):
        if self.type == "local":
            return self.rest
        else:
            return self.type.upper() + ':' + self.rest
    
    def toNonNativeUrl(self):
        if self.type[:2] == "s3":
            return 's3:' + self.rest
        else:
            return self.type + ':' + self.rest

