"""
Utilities to handle samtools faidx indexed fasta
 
Copyright (c) 2013, Allen Yu.  All rights reserved.
 
Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:
 
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
 
    * Redistributions in binary form must reproduce the above
      copyright notice, this list of conditions and the following
      disclaimer in the documentation and/or other materials provided
      with the distribution.
 
    * The names of its contributors may not be used to endorse or 
      promote products derived from this software without specific 
      prior written permission.
 
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""
 
 
class fasta:
    def __init__(self, fasta_file):        
        self.faidx = {}
         
        self.fasta_file=fasta_file
        try:
            self.fasta_handle=open(fasta_file)
        except IOError:
            print "Reference sequence doesn't exist"
 
        try:
            self.faidx_handle=open(fasta_file+".fai")
        except IOError:
            print "samtools faidx file doesn't exist for reference"   
        self.load_faidx()
         
    #Function to cache fasta index in dictionary
    #faidx format contains the following columns:
    ##.the name of the sequence
    ##.the length of the sequence
    ##.the offset of the first base in the file
    ##.the number of bases in each fasta line
    ##.the number of bytes in each fasta line    
    def load_faidx(self):
        for line in self.faidx_handle:
            line=line.strip()
            cols=line.split('\t')
            chrom = cols[0]
            slen,offset,blen,bytelen=[int(i) for i in cols[1:]]
            self.faidx[chrom]=(slen,offset,blen,bytelen)
     
    #Function to fetch sequence from an indexed fasta
    #*chrom--Chromosome name (str)
    #*start--Start position (1-based) (int)
    #*end--End position (1-based) (int)
    #*keepN--Keep ambiguous bases in sequence (boolean)        
    def fetch_sequence(self, chrom, start, end):
        #Fetch a sequence from start to end in 1-based coordinates
        seq=""
         
        if not self.faidx.has_key(chrom):
            raise ValueError('Chromosome %s not found in reference' % chrom)
        slen,offset,blen,bytelen=self.faidx[chrom]
        start = start-1 #To 0-base
        #Sanity check of start and end position
        if start<0:
            raise IndexError('Sequence window out of bound--Chr: %s\tStart:%d\tEnd:%s' % (chrom,start+1,end))
        elif start>=end:
            raise ValueError('Start position %d is larger than end position %d' % (start+1,end))
        elif end>slen and start-(end-slen)>=0: #end is out of bound, adjust the window towards start 
            end=slen
            start=start-(end-slen)
        elif end>slen:
            raise IndexError('Sequence window out of bound--Chr: %s\tStart:%d\tEnd:%s' % (chrom,start+1,end))
         
        self.fasta_handle.seek(offset+start/blen*bytelen+start%blen)
         
        while len(seq)<end-start:
            line=self.fasta_handle.readline()
            line=line.rstrip() #Remove newline symbols
            seq=seq+line
             
        #chomp off extra bases
        return seq[:end-start]
     
    def __exit__(self, type, value, traceback):
        self.fasta_handle.close()
        self.faidx_handle.close()
