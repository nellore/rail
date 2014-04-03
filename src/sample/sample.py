'''
sample.py

Helper functions for checking for and parsing sample names. 
'''

import re

parseLabReStr = ".*;LB:([^;]*)"
parseIDReStr = ".*;ID:([^;]*)"
parseLabRe = re.compile(parseLabReStr)
parseIDRe = re.compile(parseIDReStr)

def parseLab(nm):
    """Parse sample label from read name.  Label is at the end of the read name like: <name>;LB:<label>"""
    m = parseLabRe.match(nm)
    if m.group(1) is None:
        raise RuntimeError("Read name did not match regular expression '%s': %s" % (parseLabReStr, nm))
    return m.group(1)

def parseID(nm):
    """Parse sample label from read name.  Label is at the end of the read name like: <name>;LB:<label>"""
    m = parseIDRe.match(nm)
    if m.group(1) is None:
        raise RuntimeError("Read name did not match regular expression '%s': %s" % (parseIDReStr, nm))
    return m.group(1)

def hasLab(nm, mustHave=False):
    """Return true iff read name has label at the end."""
    ret = parseLabRe.match(nm)
    if not ret and mustHave:
        raise RuntimeError("Read name did not match regular expression '%s': %s" % (parseLabReStr, nm))
