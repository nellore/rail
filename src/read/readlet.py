"""
readlet.py
Description: contains the function readletize, which splits a read into readlets. 
    Adds capping readlets and accounts for uncovered edges.
    See doc/spliced_alignment.pdf for more information.
"""

def addArgs(parser):
    parser.add_argument(\
        '--readletLen', metavar='LEN', type=int, required=False, default=0,
         help='If readlets are desired, length of all readlets but capping readlets')
    parser.add_argument(\
        '--readletIval', metavar='IVAL', type=int, required=False,
        help='If readlets are desired, interval between readlet starts')
    parser.add_argument(\
        '--cappingFraction', metavar='CFRAC', type=float, required=False, default=.75,
        help='If readlets are desired, the length of each successive capping readlet is multiplied by this fraction')
    parser.add_argument(\
        '--minReadletLen', metavar='MINLEN', type=int, required=False, default=5,
         help='If readlets are desired, minimum length of capping readlets')

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
        rlets.append((nm, rlet_seq, rlet_qual, st))
        st += args.readletIval
        en += args.readletIval
    # add caps
    en = args.readletLen
    while en >= args.minReadletLen:
        rlets.append((nm, seq[:en], qual[:en], 0))
        rlets.append((nm, seq[-en:], qual[-en:], seqlen-en-1))
        en *= args.cappingFraction
        en = int(round(en))
    return rlets
