

lines = open("../hmm.out",'r').readlines()
out = open("../state_seq.out",'w')
for l in lines:
    line = l.rstrip()
    tabs = line.split("\t")
    state = int(tabs[3])
    length = int(tabs[2])
    schar = "0dED"[state]
    out.write(schar*length)
    
