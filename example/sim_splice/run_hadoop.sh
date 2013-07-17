# Parameters
GENOME_LEN=124078213
NTASKS=20
HMM_OVERLAP=100
PERMUTATIONS=5
READLET_LEN=25
READLET_IVAL=4

BOWTIE_EXE=$BOWTIE_HOME/bowtie
BOWTIE_IDX=/damsl/projects/myrna2/software/tornado/example/drosophila/Drosophila_melanogaster/UCSC/dm3/Sequence/BowtieIndex/genome
FASTA_IDX=/damsl/projects/myrna2/software/tornado/example/drosophila/Drosophila_melanogaster/UCSC/dm3/Sequence/WholeGenomeFasta/genome.fa.fai
GENOME=/damsl/projects/myrna2/software/tornado/example/drosophila/Drosophila_melanogaster/UCSC/dm3/Sequence/WholeGenomeFasta/genome.fa
INPUT=*.tab6

#Bowtie indexs
INDEX1=$BOWTIE_IDX.1.ebwt
INDEX2=$BOWTIE_IDX.2.ebwt  
INDEX3=$BOWTIE_IDX.3.ebwt  
INDEX4=$BOWTIE_IDX.4.ebwt  
INDEX5=$BOWTIE_IDX.rev.1.ebwt  
INDEX6=$BOWTIE_IDX.rev.2.ebwt

#Hadoop file locations
HADOOP_FILES=/user/hduser/sim_splice
STREAMING=$HADOOP_HOME/contrib/streaming/hadoop-streaming*.jar
HADOOP_EXE=$HADOOP_HOME/bin/hadoop

#Note: Must copy over the entire tornado directory
TORNADO=/damsl/projects/myrna2/software/tornado
EXAMPLE=$TORNADO/example
SCR_DIR=$TORNADO/src
UTIL=$SCR_DIR/util
RNAWESOME=$SCR_DIR/rnawesome
UCSC_TOOLS=$TORNADO/tools/ucsc_tools
#UCSC tools location
BIGBED_EXE=`which bedToBigBed`
SIM_SIMPLE=$EXAMPLE/sim_splice
MANIFEST_FN=$SIM_SIMPLE/fly.manifest

if [ $? -ne 0 ] ; then
        echo "bedToBigBed must be in PATH"
	exit 1
fi

CHROM_SIZES=$PWD/chrom.sizes
cat $FASTA_IDX | cut -f -2 > $CHROM_SIZES

INTERMEDIATE_DIR=$SIM_SIMPLE/intermediate

# Step 1: Readletize input reads and use Bowtie to align the readlets 
ALIGN_AGGR="cat"
ALIGN="python $RNAWESOME/align.py"
ALIGN_ARGS=''$ALIGN' --bowtieArgs '\''-v 2 -m 1 -p 6'\'' --bowtieExe '$BOWTIE_EXE' --bowtieIdx='$BOWTIE_IDX' --readletLen '$READLET_LEN' --readletIval '$READLET_IVAL' --manifest '$MANIFEST_FN' --refseq='$GENOME' --v2 --ntasks='$NTASKS' --genomeLen='$GENOME_LEN''
ALIGN_OUT=$HADOOP_FILES/align_output

# Step 2: Collapse identical intervals from same sample
MERGE_AGGR1="sort -n -k2,2"
MERGE_AGGR2="sort -s -k1,1"
MERGE="python $RNAWESOME/merge.py"
MERGE_OUT=$HADOOP_FILES/merge_output

# Step 3: Walk over genome windows and emit per-sample, per-position
#         coverage tuples
WALK_PRENORM_AGGR1="sort -n -k2,2"
WALK_PRENORM_AGGR2="sort -s -k1,1"
WALK_PRENORM="python $RNAWESOME/walk_prenorm.py"
WALK_PRENORM_OUT=$HADOOP_FILES/walk_prenorm_output
WALK_ARGS=''$WALK_PRENORM' --ntasks='$NTASKS' --genomeLen='$GENOME_LEN''

# Step 4: For all samples, take all coverage tuples for the sample and
#         from them calculate a normalization factor
SAMPLE_OUT=$HADOOP_FILES/sample_output
NORMALIZE_AGGR="sort -k1,1"
NORMALIZE="python $RNAWESOME/normalize.py"
NORMALIZE_OUT=$HADOOP_FILES/normalize_output
NORMALIZE_ARGS=''$NORMALIZE' --percentile 0.75 --out_dir='$SAMPLE_OUT' --bigbed_exe='$BIGBED_EXE' --chrom_sizes='$CHROM_SIZES' --hadoop_exe='$HADOOP_EXE''
BIGBED_EXE="$UCSC_TOOLS/bedToBigBed"

# Step 5a: Collect all the norm factors together and write to file
NORMALIZE_POST_AGGR="cat"
NORMALIZE_POST="python $RNAWESOME/normalize_post.py"
NORMALIZE_POST_ARGS=''$NORMALIZE_POST' --manifest '$MANIFEST_FN''
NORMALIZE_POST_OUT=$HADOOP_FILES/normalize_post_output

# Step 5b: Analyze introns
INTRON="python $RNAWESOME/intron.py"
INTRON_ARGS=''$INTRON' --refseq='$GENOME' --readletIval '$READLET_IVAL''

# Step 5c: Formats all splice sites into bed files
SITE2BED="python $UTIL/site2bed.py"
INTRON_OUT=$HADOOP_FILES/intron_output
SITEBED_OUT=$HADOOP_FILES/sitebed_output

# Step 6: Walk over genome windows again (taking output from Step 2)
#         but this time, calculate per-position coverage vectors and
#         fit a linear model to each
WALK_FIT="python $RNAWESOME/walk_fit.py"
WALK_FIT_ARGS=''$WALK_FIT' --ntasks='$NTASKS' --genomeLen='$GENOME_LEN' --seed=777 --normals '${INTERMEDIATE_DIR}/norm.tsv' --permutations='$PERMUTATIONS' --permutations-out='${INTERMEDIATE_DIR}/norm.tsv''
#WALK_FIT_ARGS=''$WALK_FIT' --ntasks='$NTASKS' --genomeLen='$GENOME_LEN' --seed=777 --normals '$NORMALIZE_POST_OUT/part-00000''
WALK_FIT_OUT=$HADOOP_FILES/walkfit_output 

# Step 7: Given all the t-statistics, moderate them and emit moderated
#         t-stats
EBAYES_AGGR="cat"
EBAYES="python $RNAWESOME/ebayes.py"
EBAYES_ARGS=''$EBAYES' --ntasks='$NTASKS' --genomeLen='$GENOME_LEN' --hmm-overlap='$HMM_OVERLAP''
EBAYES_OUT=$HADOOP_FILES/ebayes_output

# Step 8: Given all the moderated t-statistics, calculate the HMM
#         parameters to use in the next step
HMM_PARAMS_AGGR="cat"
HMM_PARAMS="python $RNAWESOME/hmm_params.py"
HMM_PARAMS_ARGS=''$HMM_PARAMS' --null'
HMM_PARAMS_OUT=$HADOOP_FILES/hmm_params_output

# Step 9: Given sorted bins of moderated t-statistics, and HMM
#         parameters, run the HMM
HMM_AGGR1="sort -n -k2,2"
HMM_AGGR2="sort -s -k1,1"
HMM="python $RNAWESOME/hmm.py"
HMM_ARGS=''$HMM' --ntasks='$NTASKS' --genomeLen='$GENOME_LEN' --params '${INTERMEDIATE_DIR}/hmm_params.tsv' --hmm-overlap='$HMM_OVERLAP''

HMM_OUT=$HADOOP_FILES/hmm_output

# Step 10: Output HMM path into permutation files
AGGR_PATH="python $RNAWESOME/aggr_path.py"
PERM_OUT=$HADOOP_FILES/permutation_output
NULL_OUT=$HADOOP_FILES/null 
AGGR_PATH_ARGS=''$AGGR_PATH' --out_dir='$PERM_OUT' --bigbed_exe='$BIGBED_EXE' --chrom_sizes='$CHROM_SIZES' --hadoop_exe='$HADOOP_EXE''

# Temporary files so we can form a DAG
WALK_IN_TMP=${TMPDIR}walk_in.tsv
HMM_IN_TMP=${INTERMEDIATE_DIR}/hmm_in.tsv

echo "Temporary file or walk_fit.py input is '$WALK_IN_TMP'"
echo "Temporary file for hmm.py input is '$HMM_IN_TMP'"

#Remove and Create files
rm -r intermediate
mkdir -p intermediate
hadoop dfs -rmr $SAMPLE_OUT
hadoop dfs -mkdir $SAMPLE_OUT
hadoop dfs -rmr $PERM_OUT
hadoop dfs -mkdir $PERM_OUT

#copy files over to hdfs
hadoop dfs -mkdir $HADOOP_FILES
hadoop dfs -copyFromLocal $INPUT $HADOOP_FILES

#Step 1 ALIGN
hadoop dfs -rmr $ALIGN_OUT
time hadoop jar $STREAMING \
    -D mapred.text.key.partitioner.options=-k1,1 \
    -D stream.num.map.output.key.fields=1 \
    -libjars multiplefiles.jar \
    -partitioner org.apache.hadoop.mapred.lib.KeyFieldBasedPartitioner \
    -file "$SCR_DIR/rnawesome/align.py" \
    -file "$SCR_DIR/bowtie/bowtie.py" \
    -file "$SCR_DIR/read/readlet.py" \
    -file "$SCR_DIR/read/truncate.py" \
    -file "$SCR_DIR/interval/interval.py" \
    -file "$SCR_DIR/interval/partition.py" \
    -file "$SCR_DIR/manifest/manifest.py" \
    -file "$SCR_DIR/sample/sample.py" \
    -file "$INDEX1"\
    -file "$INDEX2"\
    -file "$INDEX3"\
    -file "$INDEX4"\
    -file "$INDEX5"\
    -file "$INDEX6"\
    -file "$GENOME"\
    -file "$BOWTIE_EXE" \
    -outputformat org.myorg.MultipleOutputFormat \
    -mapper "$ALIGN_ARGS" \
    -reducer cat \
    -input $HADOOP_FILES/$INPUT -output $ALIGN_OUT


##Check $? after completion.  If not 0, then print an error message and quit
if [ $? -ne 0 ] ; then
    echo "ALIGN step failed, now exiting"
    exit 1
fi


# #Sort MERGE output
hadoop dfs -rmr $MERGE_OUT
time hadoop jar $STREAMING \
    -D mapred.text.key.partitioner.options=-k1,1 \
    -D stream.num.map.output.key.fields=2 \
    -partitioner org.apache.hadoop.mapred.lib.KeyFieldBasedPartitioner \
    -mapper 'cut -f 2-' \
    -reducer "$MERGE" \
    -input $ALIGN_OUT/exon/part* -output $MERGE_OUT

if [ $? -ne 0 ] ; then
    echo "MERGE step failed, now exiting"
    exit 1
fi

#Step 3 WALK_PRENORM
hadoop dfs -rmr $WALK_PRENORM_OUT
time hadoop jar $STREAMING \
    -D mapred.reduce.tasks=0 \
    -file "$SCR_DIR/rnawesome/walk_prenorm.py" \
    -file "$SCR_DIR/interval/partition.py" \
    -file "$SCR_DIR/manifest/manifest.py" \
    -file "$SCR_DIR/struct/circular.py" \
    -mapper "$WALK_ARGS" \
    -input $MERGE_OUT/*part* -output $WALK_PRENORM_OUT

if [ $? -ne 0 ] ; then
    echo "WALK_PRENORM step failed, now exiting"
    exit 1
fi


#Step 4 NORMALIZE
hadoop dfs -rmr $NORMALIZE_OUT
time hadoop jar $STREAMING \
    -D mapred.text.key.partitioner.options=-k1,1 \
    -D stream.num.map.output.key.fields=3 \
    -partitioner org.apache.hadoop.mapred.lib.KeyFieldBasedPartitioner \
    -file "$SCR_DIR/rnawesome/normalize.py" \
    -file "$SCR_DIR/interval/partition.py" \
    -file "$SCR_DIR/manifest/manifest.py" \
    -file "$SCR_DIR/struct/circular.py" \
    -file "$CHROM_SIZES" \
    -file "$BIGBED_EXE" \
    -mapper cat \
    -reducer "$NORMALIZE_ARGS" \
    -input $WALK_PRENORM_OUT/*part* -output $NORMALIZE_OUT

if [ $? -ne 0 ] ; then
    echo "NORMALIZE step failed, now exiting"
    exit 1
fi

#Step 5 NORMALIZE_POST
hadoop dfs -rmr $NORMALIZE_POST_OUT
time hadoop jar $STREAMING \
    -file "$SCR_DIR/rnawesome/normalize_post.py" \
    -file "$SCR_DIR/manifest/manifest.py" \
    -mapper cat \
    -reducer "$NORMALIZE_POST_ARGS" \
    -input $NORMALIZE_OUT/*part* -output $NORMALIZE_POST_OUT

if [ $? -ne 0 ] ; then
    echo "NORMALIZE_POST step failed, now exiting"
    exit 1
fi

#Step 6 INTRON
hadoop dfs -rmr $INTRON_OUT
time hadoop jar $STREAMING \
    -D mapred.text.key.partitioner.options=-k1,1 \
    -D stream.num.map.output.key.fields=3 \
    -partitioner org.apache.hadoop.mapred.lib.KeyFieldBasedPartitioner \
    -file "$SCR_DIR/rnawesome/intron.py" \
    -mapper 'cut -f 2-' \
    -reducer "$INTRON_ARGS" \
    -input $ALIGN_OUT/intron/*part* -output $INTRON_OUT

if [ $? -ne 0 ] ; then
    echo "INTRON step failed, now exiting"
    exit 1
fi


#Step 5c
hadoop dfs -rmr $SITEBED_OUT
time hadoop jar $STREAMING \
    -D mapred.text.key.partitioner.options=-k1,1 \
    -D stream.num.map.output.key.fields=3 \
    -partitioner org.apache.hadoop.mapred.lib.KeyFieldBasedPartitioner \
    -file "$SCR_DIR/check/site2bed.py" \
    -mapper cat \
    -reducer "$SITE2BED" \
    -input $INTRON_OUT -output $SITEBED_OUT

if [ $? -ne 0 ] ; then
    echo "SITEBED step failed, now exiting"
    exit 1
fi


# #Copy files to local machine
# rm ${INTERMEDIATE_DIR}/norm.tsv
# hadoop dfs -copyToLocal $NORMALIZE_POST_OUT/part* ${INTERMEDIATE_DIR}/norm.tsv
# rm ${INTERMEDIATE_DIR}/walk_in_input.tsv
# hadoop dfs -copyToLocal $MERGE_OUT/part* ${INTERMEDIATE_DIR}/walk_in_input.tsv
# rm ${INTERMEDIATE_DIR}/splice_sites.bed
# hadoop dfs -copyToLocal $SITEBED_OUT/part* ${INTERMEDIATE_DIR}/splice_sites.bed


# #Step 6 WALK_FIT 
# hadoop dfs -rmr $WALK_FIT_OUT
# hadoop jar $STREAMING \
#     -D mapred.text.key.partitioner.options=-k1,1 \
#     -D stream.num.map.output.key.fields=2 \
#     -partitioner org.apache.hadoop.mapred.lib.KeyFieldBasedPartitioner \
#     -file "$SCR_DIR/rnawesome/walk_fit.py" \
#     -file "$SCR_DIR/rnawesome/ebayes.py" \
#     -file "$SCR_DIR/interval/partition.py" \
#     -file "$SCR_DIR/struct/circular.py" \
#     -file "${INTERMEDIATE_DIR}/norm.tsv" \
#     -mapper cat \
#     -reducer "$WALK_FIT_ARGS" \
#     -input $MERGE_OUT/*part* -output $WALK_FIT_OUT

# #Step 7 EBAYES
# hadoop dfs -rmr $EBAYES_OUT
# hadoop jar $STREAMING \
#     -D mapred.text.key.partitioner.options=-k1,1 \
#     -D stream.num.map.output.key.fields=2 \
#     -partitioner org.apache.hadoop.mapred.lib.KeyFieldBasedPartitioner \
#     -file "$SCR_DIR/rnawesome/ebayes.py" \
#     -file "$SCR_DIR/interval/partition.py" \
#     -file "$SCR_DIR/struct/circular.py" \
#     -mapper cat \
#     -reducer "$EBAYES_ARGS" \
#     -input $WALK_FIT_OUT/*part* -output $EBAYES_OUT

# # # #Step 8 HMM_PARAMS
# hadoop dfs -rmr $HMM_PARAMS_OUT
# hadoop jar $STREAMING \
#     -file "$SCR_DIR/rnawesome/hmm_params.py" \
#     -mapper cat \
#     -reducer "$HMM_PARAMS_ARGS" \
#     -input $EBAYES_OUT/*part* -output $HMM_PARAMS_OUT

# # #Copy HMM params to local machine
# hadoop dfs -copyToLocal $HMM_PARAMS_OUT/part* ${INTERMEDIATE_DIR}/hmm_params.tsv

# #Step 9 HMM
# hadoop dfs -rmr $HMM_OUT
# hadoop jar $STREAMING \
#     -D mapred.text.key.partitioner.options=-k1,1 \
#     -D stream.num.map.output.key.fields=2 \
#     -file "$SCR_DIR/rnawesome/hmm.py" \
#     -file "$SCR_DIR/interval/partition.py" \
#     -file "$INTERMEDIATE_DIR/hmm_params.tsv" \
#     -partitioner org.apache.hadoop.mapred.lib.KeyFieldBasedPartitioner \
#     -mapper cat \
#     -reducer "$HMM_ARGS" \
#     -input $EBAYES_OUT/*part* -output $HMM_OUT

# #Step 10 AGGR_PATH
# hadoop dfs -rmr $NULL_OUT
# hadoop jar $STREAMING \
#     -D mapred.text.key.partitioner.options=-k1,1 \
#     -D stream.num.map.output.key.fields=3 \
#     -partitioner org.apache.hadoop.mapred.lib.KeyFieldBasedPartitioner \
#     -file "$SCR_DIR/rnawesome/aggr_path.py" \
#     -file "$CHROM_SIZES" \
#     -file "$BIGBED_EXE" \
#     -mapper cat \
#     -reducer "$AGGR_PATH_ARGS" \
#     -input $HMM_OUT/*part* -output $NULL_OUT


# # # #Check $? after completion.  If not 0, then print an error message and quit
# # # # if [$? -ne 0]
# # # # then
# # # #     echo "Splice step failed, now exiting"
# # # #     exit 0
# # # # fi
