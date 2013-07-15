


class Sequence(object):
    def __init__(self,st,en,seq):
        self.start = st
        self.end = en
        self.seq = seq
        
    def __str__(self):
        return "%d\t%s\t%d"%(self.start,self.seq,self.end)
    
    def reverse_complement(self,genome_len):
        self.start = genome_len-self.start
        self.end = genome_len-self.end
        i,j = 0,len(self.seq)-1
    
        def complement(base):
            if base=="A":
                return "T"
            elif base=="T":
                return "A"
            elif base=="C":
                return "G"
            elif base=="G":
                return "C"
            else: #Either N, R, or Y
                return base
        newseq = ""
        for i in range(0,len(self.seq)):
            j = len(self.seq)-i-1
            tmp = complement(self.seq[j])
            newseq+=tmp
        self.seq = newseq
            
