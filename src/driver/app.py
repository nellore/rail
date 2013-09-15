"""
app.py

Responsible for replacing application-related placeholders (like %BASE%) in
mapper and reducer commands.
"""

import os

base_path = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

def checkApp(appName, base):
    # Check existence of application directory
    if not os.path.exists(base) or not os.path.isdir(base):
        raise RuntimeError('No such base directory for app "%s" file: "%s"' % (appName, base))

class AppConfigEmr(object):
    """ When in emr mode, this class's config method is the way to replace
        application-related placeholders. """
    def __init__(self, emrLocalDir):
        self.emrLocalDir = os.path.abspath(emrLocalDir)
    
    def config(self, s):
        return s.replace('%BASE%', self.emrLocalDir)

class AppConfigHadoop(object):
    """ When in hadoop mode, this class's config method is the way to replace
        application-related placeholders. """
    def __init__(self, appName, checkLocal):
        self.appName = appName
        self.checkLocal = checkLocal
    
    def config(self, s):
        if self.checkLocal:
            checkApp(self.appName, base_path)
        return s.replace('%BASE%', base_path)

class AppConfigLocal(object):
    """ When in local mode, this class's config method is the way to replace
        reference-related placeholders. """
    def __init__(self, appName, checkLocal):
        self.appName = appName
        self.checkLocal = checkLocal
    
    def config(self, s):
        if self.checkLocal:
            checkApp(self.appName, base_path)
        return s.replace('%BASE%', base_path)
