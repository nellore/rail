"""
readlet.py
Description: contains the function readletize, which splits a read into readlets. 
    Adds capping readlets and accounts for uncovered edges.
    See doc/spliced_alignment.pdf for more information.
"""

def addArgs(parser):
    parser.add_argument(\
        '--readletLen', metavar='LEN', type=int, required=False, default=0,
         help='If readlets are desired, length of readlets')
    parser.add_argument(\
        '--readletIval', metavar='IVAL', type=int, required=False,
        help='If readlets are desired, interval between readlet starts')

# readletize code minus capping and accounting for uncovered edges (as of 11/5/2013)
'''
def readletize(args, nm, seq, qual):
    st, en = 0, args.readletLen
    seqlen = len(seq)
    rlets = []
    while en <= seqlen:
        rlet_seq = seq[st:en]
        rlet_qual = qual[st:en]
        assert len(rlet_seq) == args.readletLen
        assert len(rlet_qual) == args.readletLen
        rlets.append((nm, rlet_seq, rlet_qual))
        st += args.readletIval
        en += args.readletIval
    return rlets
'''

# readletize code with capping and accounting for uncovered edges
def readletize(args, nm, seq, qual):
    st, en = args.readletIval, args.readletLen + args.readletIval
    seqlen = len(seq)
    rlets = []
    while en <= seqlen:
        rlet_seq = seq[st:en]
        rlet_qual = qual[st:en]
        assert len(rlet_seq) == args.readletLen
        assert len(rlet_qual) == args.readletLen
        rlets.append((nm, rlet_seq, rlet_qual))
        st += args.readletIval
        en += args.readletIval
    # add caps
    en = args.readletLen
    while en >= 8:
        rlets.append((nm, seq[:en], qual[:en]))
        rlets.append((nm, seq[-en:], qual[-en:]))
        en -= args.readletIval
    return rlets
