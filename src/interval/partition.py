'''
partition.py
'''

import math

def addArgs(parser):
    parser.add_argument(\
        '--ntasks', metavar='NUMTASKS', type=int, required=True,
        help='Number of reduce tasks')
    parser.add_argument(\
        '--genomeLen', metavar='LEN', type=int, required=True,
        help='Total length of the genome; required so that we can accurately calculate bin sizes')

def binSize(args):
    return math.ceil(1.0 * args.genomeLen / args.ntasks)

def partition(refid, st, en, binSize):
    ''' Assign to one or more partitions based on partition bin size
        and start and end positions '''
    binid_st = int(st / binSize)
    binid_en = int((en-1) / binSize)
    return [ "_".join([refid, str(i)]) for i in xrange(binid_st, binid_en+1) ]
