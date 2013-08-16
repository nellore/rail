#Modify the parameters below
NTASKS=10000
HMM_OVERLAP=30
PERMUTATIONS=5
READLET_LEN=35
READLET_IVAL=5
RADIUS=10
SPLICE_OVERLAP=10
MAIN="$PWD/../example/sim_splice"
IGENOME="$PWD/../example/drosophila/Drosophila_melanogaster/UCSC/dm3"
RNASEQ="$MAIN/*.tab"
MANIFEST="$MAIN/fly.manifest"
INTERMEDIATE_DIR="$MAIN/intermediate"
HDFS_DIR="/user/jmorton/sim_splice"

MODE="local"      #Note: make sure that HADOOP_HOME are already set in .bashrc

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
    echo 'SPLICE_OVERLAP'=$SPLICE_OVERLAP'' >> run.sh
    echo 'SITES_FILE='$SITES_FILE'' >> run.sh
   
    python gen_script.py --ntasks=$NTASKS --hmm_overlap=$HMM_OVERLAP --permutations=$PERMUTATIONS --readlet_len=$READLET_LEN --readlet_ival=$READLET_IVAL --igenome=$IGENOME --rnaseq=$RNASEQ --manifest=$MANIFEST --intermediate=$INTERMEDIATE_DIR >> run.sh
    cat hadoop_base.sh >> run.sh
    chmod +x run.sh
elif [ "local" = "$MODE" ]; then
    python gen_script.py --ntasks=$NTASKS --hmm_overlap=$HMM_OVERLAP --permutations=$PERMUTATIONS --readlet_len=$READLET_LEN --readlet_ival=$READLET_IVAL --igenome=$IGENOME --rnaseq=$RNASEQ --manifest=$MANIFEST --intermediate=$INTERMEDIATE_DIR> run.sh
    echo 'SITES_FILE='$SITES_FILE'' >> run.sh
    echo 'SPLICE_OVERLAP'=$SPLICE_OVERLAP'' >> run.sh
    echo 'RADIUS='$RADIUS'' >> run.sh
    cat local_base.sh >> run.sh
    chmod +x run.sh
else
    echo "That option doesn't exist!"
fi
