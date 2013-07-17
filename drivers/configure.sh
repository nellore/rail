##Modify the parameters below
NTASKS=20
HMM_OVERLAP=30
PERMUTATIONS=5
READLET_LEN=30
READLET_IVAL=5
IGENOME="/home/jamie/Documents/summer2013/workspace/tornado/example/drosophila/Drosophila_melanogaster/UCSC/dm3"
RNASEQ="/damsl/projects/myrna2/langmead/tornado/data/trapnell/*.fastq"
MANIFEST="fly.manifest"
INTERMEDIATE_DIR=$PWD
MODE="local"

#Don't modify anything below here!!!
if [ "hadoop" = "$MODE" ]; then
    python gen_script.py --ntasks=$NTASKS --hmm_overlap=$HMM_OVERLAP --permutations=$PERMUTATIONS --readlet_len=$READLET_LEN --readlet_ival=$READLET_IVAL --igenome=$IGENOME --rnaseq=$RNASEQ --manifest=$MANIFEST --intermediate=$INTERMEDIATE_DIR> run.sh
    cat hadoop_base.sh >> run.sh
elif [ "local" = "$MODE" ]; then
    python gen_script.py --ntasks=$NTASKS --hmm_overlap=$HMM_OVERLAP --permutations=$PERMUTATIONS --readlet_len=$READLET_LEN --readlet_ival=$READLET_IVAL --igenome=$IGENOME --rnaseq=$RNASEQ --manifest=$MANIFEST --intermediate=$INTERMEDIATE_DIR> run.sh
    cat local_base.sh >> run.sh
else
    echo "That option doesn't exist!"
fi
