'''
readlet.py
'''

def addArgs(parser):
    parser.add_argument(\
        '--readletLen', metavar='LEN', type=int, required=False, default=0,
         help='If readlets are desired, length of readlets')
    parser.add_argument(\
        '--readletIval', metavar='IVAL', type=int, required=False,
        help='If readlets are desired, interval between readlet starts')

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
