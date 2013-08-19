"""
filemover.py

Utilities for moving files around amongst local, hdfs, s3, http and ftp
filesystems.
"""

import os
import sys
import subprocess

def addArgs(parser):
    """ Set up arguments related to moving files around """
    parser.add_argument(\
        '--s3cfg', metavar='STR', type=str, required=False,
        help='s3cmd configuration file to use (only relevant if some output is being pushed to S3)')
    parser.add_argument(\
        '--acl-public', action='store_const', const=True, default=False,
        help='Make files uploaded to S3 publicly-readable (only relevant if some output is being pushed to S3)')

class FileMover(object):
    """ Responsible for details on how to move files to and from URLs. """
    
    def __init__(self, args=None, s3cred=None, s3public=False):
        if args is not None:
            self.s3cred, self.s3public = args.s3cfg, args.acl_public
        else:
            self.s3cred, self.s3public = s3cred, s3public
    
    def put(self, fn, url):
        """ Upload a local file to a url """
        assert os.path.exists(fn)
        if url.isS3():
            cmdl = ['s3cmd']
            if self.s3cred is not None:
                cmdl.append('-c')
                cmdl.append(self.s3cred)
            cmdl.append('sync')
            if self.s3public:
                cmdl.append("--acl-public")
            cmdl.append(fn)
            cmdl.append(url.toNonNativeUrl())
        elif url.isWgettable():
            raise RuntimeError("I don't know how to upload to http/ftp URLs")
        else:
            cmdl = ['hadoop', 'fs', '-put']
            cmdl.append(fn)
            cmdl.append('/'.join([url.toUrl(), os.path.basename(fn)]))
        cmd = ' '.join(cmdl)
        print >> sys.stderr, "  Push command: '%s'" % cmd
        extl = subprocess.Popen(cmdl, stdout=sys.stderr).wait()
        print >> sys.stderr, "    Exitlevel: %d" % extl
        if extl > 0:
            raise RuntimeError("Non-zero exitlevel %d from push command '%s'" % (extl, cmd))
    
    def get(self, url, dest="."):
        """ Get a file to local directory """
        if url.isS3():
            cmdl = ["s3cmd"]
            if self.s3cred is not None:
                cmdl.append("-c")
                cmdl.append(self.s3cred)
            cmdl.append("get")
            cmdl.append(url.toNonNativeUrl())
            cmdl.append(dest)
            cmd = ' '.join(cmdl)
            extl = subprocess.Popen(cmdl, stdout=sys.stderr).wait()
            if extl > 0:
                raise RuntimeError("Non-zero exitlevel %d from s3cmd get command '%s'" % (extl, cmd))
        elif url.isWgettable():
            oldp = os.getcwd()
            os.chdir(dest)
            cmdl = ['wget', '-t', '4', '-T', '20', '-w', '25']
            cmdl.append(url.toUrl())
            cmd = ' '.join(cmdl)
            extl = subprocess.Popen(cmdl, stdout=sys.stderr).wait()
            os.chdir(oldp)
            if extl > 0:
                raise RuntimeError("Non-zero exitlevel %d from wget command '%s'" % (extl, cmd))
        else:
            cmdl = ["hadoop", "fs", "-get"]
            cmdl.append(url.toUrl())
            cmdl.append(dest)
            cmd = ' '.join(cmdl)
            extl = subprocess.Popen(cmdl, stdout=sys.stderr).wait()
            if extl > 0:
                raise RuntimeError("Non-zero exitlevel %d from hadoop fs -get command '%s'" % (extl, cmd))
