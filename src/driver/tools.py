"""
tools.py

Routines for ensuring and/or checking that appropriate helpers tools are
installed.

Specifically, we're responsible for replacing the tool placeholders in mapper
and reducer commands.

The bootstrapTool function generates bootstrap-action configuration strings,
but should perhaps be moved to an EMR-specific module.
"""

import os
import site

base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
site.addsitedir(os.path.join(base_path, "util"))

import path

toolDescs = { 'bowtie'        : ['BOWTIE'],
              'bowtie-build'  : ['BOWTIE_BUILD'],
              'bowtie2'       : ['BOWTIE2'],
              'bowtie2-build' : ['BOWTIE2_BUILD'],
              'sra-toolkit'   : ['SRA_TOOLKIT'],
              'R'             : ['R'],
              'bedToBigBed'   : ['BEDTOBIGBED', 'KENTTOOLS'],
              'samtools'      : ['SAMTOOLS'],
              'pypy'          : ['PYPY']}

# The following tools are installed in EMR mode with apt-get, and therefore
# can be found in the PATH on the EMR machines.  All other tools get placed in
# the EMR local dir by the bootstrap actions.
emrPathTools = set(['bowtie', 'bowtie-build', 'bowtie2', 'bowtie2-build', 'samtools', 'sra-toolkit', 'R', 'pypy'])

def checkLocalTool(appName, exe, envs, check=True):
    """ Check for the existence of an executable in the local environment.
        First check if the user specified a location in an environment
        variable. If so, try to use that. Otherwise, look in the PATH. """
    assert appName is not None
    assert exe is not None
    for env in envs:
        var = '_'.join([appName.upper(), env.upper()])
        if var in os.environ:
            envExe = os.environ[var]
            if check and not path.is_exe(envExe):
                raise RuntimeError('Value in "%s" environment variable does not point to executable: %s' % (var, envExe))
            return envExe
    # Favor PyPy over Python
    if exe == 'pypy':
        pathwhich = path.which(exe)
        if pathwhich is None:
            pathwhich = path.which('python')
    else:
        pathwhich = path.which(exe)
    if pathwhich is None:
        raise RuntimeError('Could not find path for %s' % exe + '.')
    return pathwhich

def configureTools(appName, s, check=True):
    """ Given name of application and a string that might have %TOOLNAME%
        placeholders in it, try to replace all the placeholders with paths to
        binaries. """
    for k, vl in toolDescs.iteritems():
        tok = '%' + k.upper() + '%'
        if tok in s:
            exe = checkLocalTool(appName, k, vl, check=check)
            s = s.replace(tok, exe)
    return s

class ToolConfigEmr(object):
    """ When in emr mode, this class's config method is the way to replace
        tool placeholders. """
    def __init__(self, emrLocalDir):
        self.emrLocalDir = emrLocalDir
    
    def config(self, s):
        for k in toolDescs.iterkeys():
            tok = '%' + k.upper() + '%'
            if tok in s:
                if k in emrPathTools:
                    # apt-get installed it in the PATH
                    s = s.replace(tok, k)
                else:
                    # a bootstrap action installed it in /mnt
                    s = s.replace(tok, '/'.join([self.emrLocalDir, 'bin', k]))
        return s

class ToolConfigHadoop(object):
    """ When in hadoop mode, this class's config method is the way to replace
        tool placeholders. """
    def __init__(self, appName, checkLocal):
        self.appName = appName
        self.checkLocal = checkLocal
    
    def config(self, s):
        return configureTools(self.appName, s, check=self.checkLocal)

class ToolConfigLocal(object):
    """ When in local mode, this class's config method is the way to replace
        tool placeholders. """
    def __init__(self, appName, checkLocal):
        self.appName = appName
        self.checkLocal = checkLocal
    
    def config(self, s):
        return configureTools(self.appName, s, check=self.checkLocal)

# Move this to an EMR-specific module?
def bootstrapTool(name, src=None, dest=None):
    """ Create a bootstrap action to install the given tool.  There must be a
        bootstrap action called install-<tool-name>.sh in the
        s3://rail-emr/bootstrap directory and the script must not need any
        arguments. """
    ret = ['--bootstrap-action s3://rail-emr/bootstrap/install-%s.sh' % name,
           '--bootstrap-name "%s"' % name ]
    if src is not None and dest is not None:
        ret.append('--args "%s,%s"' % (src.toNonNativeUrl(), dest))
    elif src is not None:
        ret.append('--args "%s"' % src.toNonNativeUrl())
    elif dest is not None:
        ret.append('--args "%s"' % dest)
    return ' '.join(ret)
