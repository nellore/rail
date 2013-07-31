"""
hadoop.py

Routines that help to abstract some things about Hadoop that vary from version
to version.
"""

def confArg(v):
    if v[0] == 0 and v[1] < 19:
        return "-jobconf"
    return "-D"

def keyFields(v, n):
    if v[0] == 0 and v[1] < 19:
        return "num.key.fields.for.partition=%d" % n
    else:
        return "mapred.text.key.partitioner.options=-k1,%d" % n
