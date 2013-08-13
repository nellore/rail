BIGBED_EXE=`which bedToBigBed`
if [ $? -ne 0 ] ; then
        echo "bedToBigBed must be in PATH"
	exit 1
fi

BOWTIE_EXE=`which bowtie`
if [ $? -ne 0 ] ; then
        echo "bowtie must be in PATH"
	exit 1
fi

HADOOP_EXE=`which hadoop`
if [ $? -ne 0 ] ; then
        echo "hadoop must be in PATH"
	exit 1
fi

PYTHON=`which python`

STREAMING=$HADOOP_HOME/contrib/streaming/hadoop-streaming*.jar
EXAMPLE=$TORNADO/example
SCR_DIR=$TORNADO/src
UTIL=$SCR_DIR/util
RNAWESOME=$SCR_DIR/rnawesome

CHROM_SIZES=$PWD/chrom.sizes
cat $FASTA_IDX | cut -f -2 > $CHROM_SIZES

# Step 0: Format Fastq files into tab delimited format
ALIGN_IN=$HADOOP_FILES/align_input
FASTQ2TAB="$PYTHON $SCR_DIR/fasta/fastq2tab.py"

# Step 1: Readletize input reads and use Bowtie to align the readlets 
ALIGN="$PYTHON $RNAWESOME/align.py"
ALIGN_ARGS=''$ALIGN' --bowtieExe '$BOWTIE_EXE' --bowtieIdx='$BOWTIE_IDX' --readletLen '$READLET_LEN' --readletIval '$READLET_IVAL' --refseq='$GENOME' --ntasks='$NTASKS' --genomeLen='$GENOME_LEN' --faidx='$FASTA_IDX' --splice-overlap='$SPLICE_OVERLAP' -- -v 2 -m 1 --mm'
ALIGN_OUT=$HADOOP_FILES/align_output

# Step 2: Collapse identical intervals from same sample
MERGE="$PYTHON $RNAWESOME/merge.py"
MERGE_OUT=$HADOOP_FILES/merge_output

# Step 3: Walk over genome windows and emit per-sample, per-position
#         coverage tuples
WALK_PRENORM="$PYTHON $RNAWESOME/walk_prenorm.py"
WALK_PRENORM_OUT=$HADOOP_FILES/walk_prenorm_output
WALK_ARGS=''$WALK_PRENORM' --ntasks='$NTASKS' --genomeLen='$GENOME_LEN''

# Step 4: For all samples, take all coverage tuples for the sample and
#         from them calculate a normalization factor
SAMPLE_OUT=$HADOOP_FILES/sample_output
NORMALIZE="$PYTHON $RNAWESOME/normalize.py"
NORMALIZE_OUT=$HADOOP_FILES/normalize_output
NORMALIZE_ARGS=''$NORMALIZE' --percentile 0.75 --out_dir='$SAMPLE_OUT' --bigbed_exe='$BIGBED_EXE' --chrom_sizes='$CHROM_SIZES' --hadoop_exe='$HADOOP_EXE''

# Step 5a: Collect all the norm factors together and write to file
NORMALIZE_POST="$PYTHON $RNAWESOME/normalize_post.py"
NORMALIZE_POST_ARGS=''$NORMALIZE_POST' --manifest '$MANIFEST_FN''
NORMALIZE_POST_OUT=$HADOOP_FILES/normalize_post_output

# Step 5b: Analyze introns
INTRON="$PYTHON $RNAWESOME/intron.py"
#INTRON_ARGS=''$INTRON' --refseq='$GENOME' --readletIval '$READLET_IVAL' --readletLen '$READLET_LEN' --radius='$RADIUS' --sites-file='$SITES_FILE''
INTRON_ARGS=''$INTRON' --refseq='$GENOME' --readletIval '$READLET_IVAL' --readletLen '$READLET_LEN' --radius='$RADIUS''

# Step 5c: Formats all splice sites into bed files
SITE2BED="$PYTHON $UTIL/site2bed.py"
INTRON_OUT=$HADOOP_FILES/intron_output
SITEBED_OUT=$HADOOP_FILES/sitebed_output

# Step 6: Walk over genome windows again (taking output from Step 2)
#         but this time, calculate per-position coverage vectors and
#         fit a linear model to each
WALK_FIT="$PYTHON $RNAWESOME/walk_fit.py"
WALK_FIT_ARGS=''$WALK_FIT' --ntasks='$NTASKS' --genomeLen='$GENOME_LEN' --seed=777 --normals '${INTERMEDIATE_DIR}/norm.tsv' --permutations='$PERMUTATIONS' --permutations-out='${INTERMEDIATE_DIR}/norm.tsv''
WALK_FIT_OUT=$HADOOP_FILES/walkfit_output 

# Step 7: Given all the t-statistics, moderate them and emit moderated
#         t-stats
EBAYES_AGGR="cat"
EBAYES="$PYTHON $RNAWESOME/ebayes.py"
EBAYES_ARGS=''$EBAYES' --ntasks='$NTASKS' --genomeLen='$GENOME_LEN' --hmm-overlap='$HMM_OVERLAP''
EBAYES_OUT=$HADOOP_FILES/ebayes_output

# Step 8: Given all the moderated t-statistics, calculate the HMM
#         parameters to use in the next step
HMM_PARAMS_AGGR="cat"
HMM_PARAMS="$PYTHON $RNAWESOME/hmm_params.py"
HMM_PARAMS_ARGS=''$HMM_PARAMS' --null'
HMM_PARAMS_OUT=$HADOOP_FILES/hmm_params_output

# Step 9: Given sorted bins of moderated t-statistics, and HMM
#         parameters, run the HMM
HMM_AGGR1="sort -n -k2,2"
HMM_AGGR2="sort -s -k1,1"
HMM="$PYTHON $RNAWESOME/hmm.py"
HMM_ARGS=''$HMM' --ntasks='$NTASKS' --genomeLen='$GENOME_LEN' --params '${INTERMEDIATE_DIR}/hmm_params.tsv' --hmm-overlap='$HMM_OVERLAP''

HMM_OUT=$HADOOP_FILES/hmm_output

# Step 10: Output HMM path into permutation files
AGGR_PATH="$PYTHON $RNAWESOME/aggr_path.py"
PERM_OUT=$HADOOP_FILES/permutation_output
NULL_OUT=$HADOOP_FILES/null 
AGGR_PATH_ARGS=''$AGGR_PATH' --out_dir='$PERM_OUT' --bigbed_exe='$BIGBED_EXE' --chrom_sizes='$CHROM_SIZES' --hadoop_exe='$HADOOP_EXE''

# Temporary files so we can form a DAG
WALK_IN_TMP=${TMPDIR}walk_in.tsv
HMM_IN_TMP=${INTERMEDIATE_DIR}/hmm_in.tsv

echo "Temporary file or walk_fit.py input is '$WALK_IN_TMP'"
echo "Temporary file for hmm.py input is '$HMM_IN_TMP'"

#Remove and Create files
rm -r ${INTERMEDIATE_DIR}
mkdir -p ${INTERMEDIATE_DIR}
hadoop dfs -rmr $SAMPLE_OUT
hadoop dfs -mkdir $SAMPLE_OUT
hadoop dfs -rmr $PERM_OUT
hadoop dfs -mkdir $PERM_OUT

# hadoop dfs -rmr $ALIGN_IN
# hadoop dfs -mkdir $ALIGN_IN
# for FASTQ in $RNASEQ
# do
#   FILE=`basename $FASTQ`  
#   cat <<EOF > $FILE.sh 
#   cat $FASTQ | $FASTQ2TAB > $FILE.tab
# EOF
#   sh $FILE.sh &
# done
# wait

# for FASTQsh in *fastq.sh
# do
#   rm $FASTQsh
# done

# for FASTQ in *.tab
# do
#   FILE=`basename $FASTQ`  
#   cat $FILE | hadoop dfs -put - "$ALIGN_IN/$FILE"
#   rm $FILE
# done

# #Step 1 ALIGN
# hadoop dfs -rmr $ALIGN_OUT
# time hadoop jar $STREAMING \
#     -D mapred.text.key.partitioner.options=-k1,1 \
#     -D stream.num.map.output.key.fields=1 \
#     -D mapred.reduce.tasks=0 \
#     -libjars multiplefiles.jar \
#     -cmdenv PYTHONPATH=$PYTHONPATH \
#     -cmdenv PYTHONUSERBASE=$PYTHONUSERBASE \
#     -cmdenv PYTHONUSERSITE=$PYTHONUSERSITE \
#     -partitioner org.apache.hadoop.mapred.lib.KeyFieldBasedPartitioner \
#     -file "$SCR_DIR/rnawesome/align.py" \
#     -file "$SCR_DIR/bowtie/bowtie.py" \
#     -file "$SCR_DIR/read/readlet.py" \
#     -file "$SCR_DIR/read/truncate.py" \
#     -file "$SCR_DIR/interval/interval.py" \
#     -file "$SCR_DIR/interval/partition.py" \
#     -file "$SCR_DIR/manifest/manifest.py" \
#     -file "$SCR_DIR/sample/sample.py" \
#     -file "$SCR_DIR/util/path.py" \
#     -file "$GENOME" \
#     -file "$FASTA_IDX" \
#     -file "$INDEX1"\
#     -file "$INDEX2"\
#     -file "$INDEX3"\
#     -file "$INDEX4"\
#     -file "$INDEX5"\
#     -file "$INDEX6"\
#     -file "$GENOME"\
#     -file "$BOWTIE_EXE" \
#     -outputformat edu.jhu.cs.MultipleOutputFormat \
#     -mapper "$ALIGN_ARGS" \
#     -input $ALIGN_IN/*.tab -output $ALIGN_OUT


# ##Check $? after completion.  If not 0, then print an error message and quit
# if [ $? -ne 0 ] ; then
#     echo "ALIGN step failed, now exiting"
#     exit 1
# fi


# # #Sort MERGE output
# hadoop dfs -rmr $MERGE_OUT
# time hadoop jar $STREAMING \
#     -D mapred.reduce.tasks=32 \
#     -D mapred.text.key.partitioner.options=-k1,1 \
#     -D stream.num.map.output.key.fields=2 \
#     -cmdenv PYTHONPATH=$PYTHONPATH \
#     -cmdenv PYTHONUSERBASE=$PYTHONUSERBASE \
#     -cmdenv PYTHONUSERSITE=$PYTHONUSERSITE \
#     -partitioner org.apache.hadoop.mapred.lib.KeyFieldBasedPartitioner \
#     -mapper 'cat' \
#     -reducer "$MERGE" \
#     -input $ALIGN_OUT/exon/part* -output $MERGE_OUT

# if [ $? -ne 0 ] ; then
#     echo "MERGE step failed, now exiting"
#     exit 1
# fi

# #Step 3 WALK_PRENORM
# hadoop dfs -rmr $WALK_PRENORM_OUT
# time hadoop jar $STREAMING \
#     -D mapred.reduce.tasks=32 \
#     -D mapred.text.key.partitioner.options=-k1,1 \
#     -D stream.num.map.output.key.fields=3 \
#     -cmdenv PYTHONPATH=$PYTHONPATH \
#     -cmdenv PYTHONUSERBASE=$PYTHONUSERBASE \
#     -cmdenv PYTHONUSERSITE=$PYTHONUSERSITE \
#     -file "$SCR_DIR/rnawesome/walk_prenorm.py" \
#     -file "$SCR_DIR/interval/partition.py" \
#     -file "$SCR_DIR/manifest/manifest.py" \
#     -file "$SCR_DIR/struct/circular.py" \
#     -mapper cat \
#     -reducer "$WALK_ARGS" \
#     -input $MERGE_OUT/*part* -output $WALK_PRENORM_OUT

# if [ $? -ne 0 ] ; then
#     echo "WALK_PRENORM step failed, now exiting"
#     exit 1
# fi


# #Step 4 NORMALIZE
# hadoop dfs -rmr $NORMALIZE_OUT
# time hadoop jar $STREAMING \
#     -D mapred.reduce.tasks=32 \
#     -D mapred.text.key.partitioner.options=-k1,1 \
#     -D stream.num.map.output.key.fields=2 \
#     -cmdenv PYTHONPATH=$PYTHONPATH \
#     -cmdenv PYTHONUSERBASE=$PYTHONUSERBASE \
#     -cmdenv PYTHONUSERSITE=$PYTHONUSERSITE \
#     -partitioner org.apache.hadoop.mapred.lib.KeyFieldBasedPartitioner \
#     -file "$SCR_DIR/rnawesome/normalize.py" \
#     -file "$SCR_DIR/interval/partition.py" \
#     -file "$SCR_DIR/manifest/manifest.py" \
#     -file "$SCR_DIR/struct/circular.py" \
#     -file "$CHROM_SIZES" \
#     -file "$BIGBED_EXE" \
#     -mapper cat \
#     -reducer "$NORMALIZE_ARGS" \
#     -input $WALK_PRENORM_OUT/*part* -output $NORMALIZE_OUT

# if [ $? -ne 0 ] ; then
#     echo "NORMALIZE step failed, now exiting"
#     exit 1
# fi

# #Step 5 NORMALIZE_POST
# hadoop dfs -rmr $NORMALIZE_POST_OUT
# time hadoop jar $STREAMING \
#     -D mapred.reduce.tasks=32 \
#     -cmdenv PYTHONPATH=$PYTHONPATH \
#     -cmdenv PYTHONUSERBASE=$PYTHONUSERBASE \
#     -cmdenv PYTHONUSERSITE=$PYTHONUSERSITE \
#     -file "$SCR_DIR/rnawesome/normalize_post.py" \
#     -file "$SCR_DIR/manifest/manifest.py" \
#     -mapper cat \
#     -reducer "$NORMALIZE_POST_ARGS" \
#     -input $NORMALIZE_OUT/*part* -output $NORMALIZE_POST_OUT

# if [ $? -ne 0 ] ; then
#     echo "NORMALIZE_POST step failed, now exiting"
#     exit 1
# fi

#Step 6 INTRON
hadoop dfs -rmr $INTRON_OUT
time hadoop jar $STREAMING \
    -D mapred.reduce.tasks=32 \
    -D mapred.text.key.partitioner.options=-k1,1 \
    -D stream.num.map.output.key.fields=3 \
    -libjars multiplefiles.jar \
    -cmdenv PYTHONPATH=$PYTHONPATH \
    -cmdenv PYTHONUSERBASE=$PYTHONUSERBASE \
    -cmdenv PYTHONUSERSITE=$PYTHONUSERSITE \
    -partitioner org.apache.hadoop.mapred.lib.KeyFieldBasedPartitioner \
    -outputformat edu.jhu.cs.MultipleOutputFormat \
    -file "$SCR_DIR/rnawesome/intron.py" \
    -mapper 'cat' \
    -reducer "$INTRON_ARGS" \
    -input $ALIGN_OUT/intron/*part* -output $INTRON_OUT

if [ $? -ne 0 ] ; then
    echo "INTRON step failed, now exiting"
    exit 1
fi


#Step 5c
hadoop dfs -rmr $SITEBED_OUT
time hadoop jar $STREAMING \
    -D mapred.reduce.tasks=32 \
    -D mapred.text.key.partitioner.options=-k1,1 \
    -D stream.num.map.output.key.fields=3 \
    -cmdenv PYTHONPATH=$PYTHONPATH \
    -cmdenv PYTHONUSERBASE=$PYTHONUSERBASE \
    -cmdenv PYTHONUSERSITE=$PYTHONUSERSITE \
    -cmdenv PYTHONCOMPILED=$PYTHONCOMPILED \
    -partitioner org.apache.hadoop.mapred.lib.KeyFieldBasedPartitioner \
    -file "$SCR_DIR/util/site2bed.py" \
    -mapper cat \
    -reducer "$SITE2BED" \
    -input $INTRON_OUT/site/*part* -output $SITEBED_OUT

if [ $? -ne 0 ] ; then
    echo "SITEBED step failed, now exiting"
    exit 1
fi


#Copy files to local machine
rm ${INTERMEDIATE_DIR}/normalize
mkdir ${INTERMEDIATE_DIR}/normalize
hadoop dfs -copyToLocal $NORMALIZE_POST_OUT/part* ${INTERMEDIATE_DIR}/normalize
cat ${INTERMEDIATE_DIR}/normalize/* > ${INTERMEDIATE_DIR}/norm.tsv
rm ${INTERMEDIATE_DIR}/walk_in_input
mkdir ${INTERMEDIATE_DIR}/walk_in_input
hadoop dfs -copyToLocal $MERGE_OUT/part* ${INTERMEDIATE_DIR}/walk_in_input
rm ${INTERMEDIATE_DIR}/splice_sites
mkdir ${INTERMEDIATE_DIR}/splice_sites
hadoop dfs -copyToLocal $SITEBED_OUT/part* ${INTERMEDIATE_DIR}/splice_sites


#Step 6 WALK_FIT 
hadoop dfs -rmr $WALK_FIT_OUT
hadoop jar $STREAMING \
    -D mapred.reduce.tasks=32 \
    -D mapred.text.key.partitioner.options=-k1,1 \
    -D stream.num.map.output.key.fields=2 \
    -cmdenv PYTHONPATH=$PYTHONPATH \
    -cmdenv PYTHONUSERBASE=$PYTHONUSERBASE \
    -cmdenv PYTHONUSERSITE=$PYTHONUSERSITE \
    -partitioner org.apache.hadoop.mapred.lib.KeyFieldBasedPartitioner \
    -file "$SCR_DIR/rnawesome/walk_fit.py" \
    -file "$SCR_DIR/rnawesome/ebayes.py" \
    -file "$SCR_DIR/interval/partition.py" \
    -file "$SCR_DIR/struct/circular.py" \
    -file "${INTERMEDIATE_DIR}/norm.tsv" \
    -mapper cat \
    -reducer "$WALK_FIT_ARGS" \
    -input $MERGE_OUT/*part* -output $WALK_FIT_OUT

#Step 7 EBAYES
hadoop dfs -rmr $EBAYES_OUT
hadoop jar $STREAMING \
    -D mapred.reduce.tasks=32 \
    -D mapred.text.key.partitioner.options=-k1,1 \
    -D stream.num.map.output.key.fields=2 \
    -cmdenv PYTHONPATH=$PYTHONPATH \
    -cmdenv PYTHONUSERBASE=$PYTHONUSERBASE \
    -cmdenv PYTHONUSERSITE=$PYTHONUSERSITE \
    -partitioner org.apache.hadoop.mapred.lib.KeyFieldBasedPartitioner \
    -file "$SCR_DIR/rnawesome/ebayes.py" \
    -file "$SCR_DIR/interval/partition.py" \
    -file "$SCR_DIR/struct/circular.py" \
    -mapper cat \
    -reducer "$EBAYES_ARGS" \
    -input $WALK_FIT_OUT/*part* -output $EBAYES_OUT

# # #Step 8 HMM_PARAMS
hadoop dfs -rmr $HMM_PARAMS_OUT
hadoop jar $STREAMING \
    -D mapred.reduce.tasks=32 \
    -cmdenv PYTHONPATH=$PYTHONPATH \
    -cmdenv PYTHONUSERBASE=$PYTHONUSERBASE \
    -cmdenv PYTHONUSERSITE=$PYTHONUSERSITE \
    -file "$SCR_DIR/rnawesome/hmm_params.py" \
    -mapper cat \
    -reducer "$HMM_PARAMS_ARGS" \
    -input $EBAYES_OUT/*part* -output $HMM_PARAMS_OUT

# #Copy HMM params to local machine
hadoop dfs -copyToLocal $HMM_PARAMS_OUT/part* ${INTERMEDIATE_DIR}/hmm_params.tsv

#Step 9 HMM
hadoop dfs -rmr $HMM_OUT
hadoop jar $STREAMING \
    -D mapred.reduce.tasks=32 \
    -D mapred.text.key.partitioner.options=-k1,1 \
    -D stream.num.map.output.key.fields=2 \
    -cmdenv PYTHONPATH=$PYTHONPATH \
    -cmdenv PYTHONUSERBASE=$PYTHONUSERBASE \
    -cmdenv PYTHONUSERSITE=$PYTHONUSERSITE \
    -file "$SCR_DIR/rnawesome/hmm.py" \
    -file "$SCR_DIR/interval/partition.py" \
    -file "$INTERMEDIATE_DIR/hmm_params.tsv" \
    -partitioner org.apache.hadoop.mapred.lib.KeyFieldBasedPartitioner \
    -mapper cat \
    -reducer "$HMM_ARGS" \
    -input $EBAYES_OUT/*part* -output $HMM_OUT

#Step 10 AGGR_PATH
hadoop dfs -rmr $NULL_OUT
hadoop jar $STREAMING \
    -D mapred.reduce.tasks=32 \
    -D mapred.text.key.partitioner.options=-k1,1 \
    -D stream.num.map.output.key.fields=3 \
    -cmdenv PYTHONPATH=$PYTHONPATH \
    -cmdenv PYTHONUSERBASE=$PYTHONUSERBASE \
    -cmdenv PYTHONUSERSITE=$PYTHONUSERSITE \
    -partitioner org.apache.hadoop.mapred.lib.KeyFieldBasedPartitioner \
    -file "$SCR_DIR/rnawesome/aggr_path.py" \
    -file "$CHROM_SIZES" \
    -file "$BIGBED_EXE" \
    -mapper cat \
    -reducer "$AGGR_PATH_ARGS" \
    -input $HMM_OUT/*part* -output $NULL_OUT


