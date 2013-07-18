"""
Generates the run_hadoop.sh header 
"""
import sys
import argparse
import os
import site
base_path = os.path.dirname(os.path.dirname((os.path.abspath(__file__))))
site.addsitedir(os.path.join(base_path, "src/fasta"))

import chrsizes

parser = argparse.ArgumentParser(description=\
                                     'Generates a run script.')
parser.add_argument(\
    '--ntasks',type=int,default=20,help='Number of tasks to be divided up')
parser.add_argument(\
    '--hmm_overlap',type=int,default=30,help='Overlap length of Viterbi sequences')
parser.add_argument(\
    '--permutations',type=int,default=5,help='Number permutations required to conduct permutation test')
parser.add_argument(\
    '--readlet_len',type=int,default=30,help='Length of readlets'
)
parser.add_argument(\
    '--readlet_ival',type=int,default=5,help='Sampling interval length of readlets'
)
parser.add_argument(\
    '--igenome',type=str,required=True,help='Path of Illumina iGenome package'
)
parser.add_argument(\
    '--rnaseq',type=str,required=True,help='Path of input RNA sequence data'
)
parser.add_argument(\
    '--manifest',type=str,required=True,help='Path of manifest file'
)
parser.add_argument(\
    '--intermediate',type=str,default="$PWD/intermediate",help='Path of the intermediate folder'
)

args = parser.parse_args()

def print_parameters(ntasks,hmm_overlap,permutations,readlet_len,readlet_ival):
    print "NTASKS=%d"%(ntasks)
    print "HMM_OVERLAP=%d"%(hmm_overlap)
    print "PERMUTATIONS=%d"%(permutations)
    print "READLET_LEN=%d"%(readlet_len)
    print "READLET_IVAL=%d"%(readlet_ival)
    
def print_paths(igenome,rnaseq,manifest,intermediate):
    cwd = os.getcwd()
    bowtie_idx = "%s/Sequence/BowtieIndex/genome"%igenome
    fasta_idx = "%s/Sequence/WholeGenomeFasta/genome.fa.fai"%igenome
    genome_length = chrsizes.totalLength(fasta_idx)
    print "TORNADO=%s/.."%cwd
    print "IGENOME=%s"%igenome
    print "BOWTIE_IDX=%s"%bowtie_idx
    print "GENOME=%s/Sequence/WholeGenomeFasta/genome.fa"%igenome
    print "FASTA_IDX=%s"%fasta_idx
    print "RNASEQ=%s"%rnaseq
    print "INDEX1=%s.1.ebwt"%(bowtie_idx)
    print "INDEX2=%s.2.ebwt"%(bowtie_idx)
    print "INDEX3=%s.3.ebwt"%(bowtie_idx)
    print "INDEX4=%s.4.ebwt"%(bowtie_idx)
    print "INDEX5=%s.rev.1.ebwt"%(bowtie_idx)
    print "INDEX6=%s.rev.2.ebwt"%(bowtie_idx)
    print "GENOME_LEN=%d"%(genome_length)
    print "MANIFEST_FN=%s"%manifest
    print "INTERMEDIATE_DIR=%s"%(intermediate)

    
if __name__=="__main__":
    print_parameters(args.ntasks,
                     args.hmm_overlap,
                     args.permutations,
                     args.readlet_len,
                     args.readlet_ival)
    print_paths(args.igenome,args.rnaseq,args.manifest,args.intermediate)
    
