"""
tools.py

Ensure certain executables are installed, either via bootstrap action or by
downloading and installing.
"""

def bootstrapTool(name, dest=None):
    """ Create a bootstrap action to install the given tool.  There must be a
        bootstrap action called install-<tool-name>.sh in the
        s3://tornado-emr/bootstrap directory and the script must not need any
        arguments. """
    ret = ['--bootstrap-action s3://tornado-emr/bootstrap/install-%s.sh' % name,
           '--bootstrap-name "%s"' % name ]
    if dest is not None:
        ret.append('--args "%s"' % dest)
    return ' '.join(ret)
