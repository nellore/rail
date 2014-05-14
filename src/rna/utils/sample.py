"""
sample.py
Part of Rail-RNA

Contains helper functions for checking for and parsing sample names. 
"""

import re

parse_label_pattern = ".*;LB:([^;]*)"
parse_id_pattern = ".*;ID:([^;]*)"
parse_label_compiled = re.compile(parse_label_pattern)
parse_id_compiled = re.compile(parse_id_pattern)

def parse_label(qname):
    """ Parse sample label from read's QNAME. 

        The label is at the end of the read name: <name>;LB:<label>.

        qname: QNAME in SAM format

        Return value: sample label
    """
    name_match = parse_label_compiled.match(qname)
    if name_match.group(1) is None:
        raise RuntimeError(('Read name did not match regular expression '
                            '"%s": %s') % (parse_label_pattern, qname))
    return name_match.group(1)

def parse_id(qname):
    """ Parse ID from read name.

        The ID is embedded in the read name: <name>;ID:<label>.

        qname: QNAME in SAM format

        Return value: ID
    """
    name_match = parse_id_compiled.match(qname)
    if name_match.group(1) is None:
        raise RuntimeError(('Read name did not match regular expression '
                            '"%s": %s') % (parse_id_pattern, qname))
    return name_match.group(1)

def has_label(qname):
    """ Return True iff read name has label at the end.

        qname: QNAME in SAM format

        Return value: True iff read name has a label.
    """
    if not parse_label_compiled.match(qname):
        return False
    return True
