'''
sample.py
'''

import re

parseLabReStr = ".*;LB:(.*)$"
parseLabRe = re.compile(parseLabReStr)

def parseLab(nm):
    ''' Parse sample label from read name.  Label is at the end of the
        read name like: <name>;LB:<label> '''
    m = parseLabRe.match(nm)
    if m.group(1) is None:
        raise RuntimeError("Read name did not match regular expression '%s': %s" % (parseLabReStr, nm))
    return m.group(1)

def hasLab(nm, mustHave=False):
    ''' Return true iff read name has label at the end. '''
    ret = parseLabRe.match(nm)
    if not ret and mustHave:
        raise RuntimeError("Read name did not match regular expression '%s': %s" % (parseLabReStr, nm))
