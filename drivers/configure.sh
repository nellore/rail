##Modify the parameters below
NTASKS=20
HMM_OVERLAP=30
PERMUTATIONS=5
READLET_LEN=30
READLET_IVAL=5
RADIUS=10

IGENOME="/damsl/projects/myrna2/software/tornado/example/drosophila/Drosophila_melanogaster/UCSC/dm3"
RNASEQ="/damsl/projects/myrna2/software/tornado/example/sim_splice/*.tab"
MANIFEST="/damsl/projects/myrna2/software/tornado/example/sim_splice/fly.manifest"
INTERMEDIATE_DIR="/damsl/projects/myrna2/software/tornado/example/sim_splice/intermediate"
HDFS_DIR="/user/jmorton/sim_splice"

# IGENOME="/damsl/projects/myrna2/langmead/igenomes/Homo_sapiens/UCSC/hg19"
# RNASEQ="/damsl/projects/myrna2/software/tornado/example/sim_human/*.fastq"
# MANIFEST="/damsl/projects/myrna2/software/tornado/example/sim_human/human.manifest"
# INTERMEDIATE_DIR="/damsl/projects/myrna2/software/tornado/example/sim_human/intermediate"
# HDFS_DIR="/user/jmorton/sim_human"

MODE="hadoop"      #Note: make sure that HADOOP_HOME are already set in .bashrc

#System variables that hopefully shouldn't be tinkered with
PYTHONPATH=$PYTHONPATH:/home/jmorton/.local
PYTHONUSERBASE=/home/jmorton/.local
PYTHONUSERSITE=/home/jmorton/.local/lib/python2.7/site-packages
PYTHONCOMPILED=$PWD


#Don't modify anything below here!!!
if [ "hadoop" = "$MODE" ]; then
    echo 'HADOOP_FILES='$HDFS_DIR'' > run.sh
    echo 'PYTHONPATH='$PYTHONPATH'' >> run.sh
    echo 'PYTHONUSERBASE='$PYTHONUSERBASE'' >> run.sh
    echo 'PYTHONUSERSITE='$PYTHONUSERSITE'' >> run.sh
    echo 'PYTHONCOMPILED='$PYTHONCOMPILED'' >> run.sh
    echo 'RADIUS='$RADIUS'' >> run.sh
    python gen_script.py --ntasks=$NTASKS --hmm_overlap=$HMM_OVERLAP --permutations=$PERMUTATIONS --readlet_len=$READLET_LEN --readlet_ival=$READLET_IVAL --igenome=$IGENOME --rnaseq=$RNASEQ --manifest=$MANIFEST --intermediate=$INTERMEDIATE_DIR >> run.sh
    cat hadoop_base.sh >> run.sh
    chmod +x run.sh
elif [ "local" = "$MODE" ]; then
    echo 'RADIUS='$RADIUS'' >> run.sh
    python gen_script.py --ntasks=$NTASKS --hmm_overlap=$HMM_OVERLAP --permutations=$PERMUTATIONS --readlet_len=$READLET_LEN --readlet_ival=$READLET_IVAL --igenome=$IGENOME --rnaseq=$RNASEQ --manifest=$MANIFEST --intermediate=$INTERMEDIATE_DIR> run.sh
    cat local_base.sh >> run.sh
    chmod +x run.sh
else
    echo "That option doesn't exist!"
fi
