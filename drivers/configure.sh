##Modify the parameters below
NTASKS=20
HMM_OVERLAP=30
PERMUTATIONS=5
READLET_LEN=30
READLET_IVAL=5
IGENOME="/scratch0/igm2/langmead/igenomes/Homo_sapiens/UCSC/hg19"
RNASEQ="/home/jmorton/data/human/*.fastq"
MANIFEST="/home/jmorton/workspace/tornado/drivers/../example/sim_splice/fly.manifest"
INTERMEDIATE_DIR="/home/jmorton/workspace/tornado/drivers/../example/sim_human/intermediate"
PYTHONPATH=$PYTHONPATH:/home/jmorton/.local/bin:/home/jmorton/.local/lib
PYTHONUSERBASE=/home/jmorton/.local
PYTHONUSERSITE=/home/jmorton/.local/lib/python2.7/site-packages
MODE="hadoop"      #Note: make sure that HADOOP_HOME are already set in .bashrc
HDFS_DIR="/user/jmorton/sim_human"

#Don't modify anything below here!!!
if [ "hadoop" = "$MODE" ]; then
    echo 'HADOOP_FILES='$HDFS_DIR'' > run.sh
    echo 'PYTHONPATH='$PYTHONPATH'' >> run.sh
    echo 'PYTHONUSERBASE='$PYTHONUSERBASE'' >> run.sh
    echo 'PYTHONUSERSITE='$PYTHONUSERSITE'' >> run.sh
    python gen_script.py --ntasks=$NTASKS --hmm_overlap=$HMM_OVERLAP --permutations=$PERMUTATIONS --readlet_len=$READLET_LEN --readlet_ival=$READLET_IVAL --igenome=$IGENOME --rnaseq=$RNASEQ --manifest=$MANIFEST --intermediate=$INTERMEDIATE_DIR >> run.sh
    cat hadoop_base.sh >> run.sh
    chmod +x run.sh
elif [ "local" = "$MODE" ]; then
    python gen_script.py --ntasks=$NTASKS --hmm_overlap=$HMM_OVERLAP --permutations=$PERMUTATIONS --readlet_len=$READLET_LEN --readlet_ival=$READLET_IVAL --igenome=$IGENOME --rnaseq=$RNASEQ --manifest=$MANIFEST --intermediate=$INTERMEDIATE_DIR> run.sh
    cat local_base.sh >> run.sh
    chmod +x run.sh
else
    echo "That option doesn't exist!"
fi
