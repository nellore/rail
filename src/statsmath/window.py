"""
Scores a set of windows based off of splice site
"""
def score(seq, site, hist,cost):
    wsize = len(site) # window size
    nwins = len(seq)-wsize+1
    wins = [0]*nwins

    for i in range(0,nwins):
        for j in range(0,len(site)):
            s = 1 if site[j]==seq[i+j] else cost
            wins[i]+=s*hist[i+j]
    return wins
