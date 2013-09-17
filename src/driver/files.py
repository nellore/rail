"""
files.py

Responsible for replacing file placeholders (like %MANIFEST%) in mapper and
reducer commands.
"""

import os

def addFileArgs(parser):
    """ Arguments that refer to files that have to be on the local filesystem,
        or make their way there. """
    parser.add_argument(\
        '--manifest', metavar='PATH', type=str, required=False, help='URL for manifest file')

class FileConfigEmr(object):
    """ When in emr mode, this class's config method is the way to replace
        file-related placeholders. """
    def __init__(self, emrLocalDir):
        self.emrLocalDir = os.path.abspath(emrLocalDir)
    
    def config(self, s):
        return s.replace('%MANIFEST%', os.path.join(self.emrLocalDir, 'MANIFEST'))

class FileConfigHadoop(object):
    """ When in hadoop mode, this class's config method is the way to replace
        application-related placeholders. """
    def __init__(self, manifest, checkLocal):
        assert manifest.isLocal()
        self.manifest = os.path.abspath(manifest.toUrl())
        self.checkLocal = checkLocal
    
    def config(self, s):
        if self.checkLocal and not os.path.isfile(self.manifest):
            raise RuntimeError('No such --manifest file: "%s"' % self.manifest)
        return s.replace('%MANIFEST%', self.manifest)

class FileConfigLocal(object):
    """ When in local mode, this class's config method is the way to replace
        file-related placeholders. """
    def __init__(self, manifest, checkLocal):
        assert manifest.isLocal()
        self.manifest = os.path.abspath(manifest.toUrl())
        self.checkLocal = checkLocal
    
    def config(self, s):
        if self.checkLocal and not os.path.isfile(self.manifest):
            raise RuntimeError('No such --manifest file: "%s"' % self.manifest)
        return s.replace('%MANIFEST%', self.manifest)
