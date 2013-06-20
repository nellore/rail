# Parameters
GENOME_LEN=1000
NTASKS=10
HMM_OVERLAP=30

#Hadoop file locations
HADOOP_FILES=/user/hduser/sim_simple
STREAMING=$HADOOP_HOME/contrib/streaming/hadoop-streaming*.jar

#Note: Must copy over the entire tornado directory
TORNADO=/home/hduser/workspace/tornado
EXAMPLE=$TORNADO/example


#Bowtie locations
BOWTIE=/home/hduser/software/bowtie-0.12.8/bowtie
INDEXS=$EXAMPLE/fasta
BOWTIE_IDX=$INDEXS/lambda_virus

SCR_DIR=/home/hduser/workspace/tornado/src


SIM_SIMPLE=$EXAMPLE/sim_simple
MANIFEST_FN=$SIM_SIMPLE/eg1.manifest

mkdir -p intermediate
INTERMEDIATE_DIR=intermediate/

# Step 1a: Readletize input reads and use Bowtie to align the readlets 
ALIGN_AGGR="cat"
ALIGN="python $SCR_DIR/rnawesome/align.py"
ALIGN_ARGS=''$ALIGN' --bowtieArgs '\''-v 2 -m 1 -p 6'\'' --bowtieExe '$BOWTIE' --bowtieIdx='$BOWTIE_IDX' --readletLen 20 --readletIval 2 --manifest '$MANIFEST_FN''
ALIGN_OUT=/user/hduser/sim_simple/align_output

# Step 1b: Bring together all the readlet alignment intervals for a given read and 
SPLICE_AGGR="cat"
SPLICE="python $SCR_DIR/rnawesome/splice.py"
SPLICE_ARGS=''$SPLICE' --ntasks='$NTASKS' --genomeLen='$GENOME_LEN' --manifest '$MANIFEST_FN''
SPLICE_OUT=/user/hduser/sim_simple/splice_output

# Steps 1a and 1b happen together in the same map step in practice

# Step 2: Collapse identical intervals from same sample
MERGE_AGGR1="sort -n -k2,2"
MERGE_AGGR2="sort -s -k1,1"
MERGE="python $SCR_DIR/rnawesome/merge.py"
MERGE_OUT=/user/hduser/sim_simple/merge_output

# Step 3: Walk over genome windows and emit per-sample, per-position
#         coverage tuples
WALK_PRENORM_AGGR1="sort -n -k2,2"
WALK_PRENORM_AGGR2="sort -s -k1,1"
WALK_PRENORM="python $SCR_DIR/rnawesome/walk_prenorm.py"
WALK_OUT=/user/hduser/sim_simple/walk_output
WALK_ARGS=''$WALK_PRENORM' --ntasks='$NTASKS' --genomeLen='$GENOME_LEN' --manifest '$MANIFEST_FN''

# Step 4: For all samples, take all coverage tuples for the sample and
#         from them calculate a normalization factor
NORMALIZE_AGGR="sort -k1,1"
NORMALIZE="python $SCR_DIR/normalize.py"
NORMALIZE_OUT=/user/hduser/sim_simple/normalize_output

# Step 5: Collect all the norm factors together and write to file
NORMALIZE_POST_AGGR="cat"
NORMALIZE_POST="python $SCR_DIR/normalize_post.py"
NORMALIZE_POST_OUT=/user/hduser/sim_simple/normalize_post_output

# Step 6: Walk over genome windows again (taking output from Step 2)
#         but this time, calculate per-position coverage vectors and
#         fit a linear model to each
WALK_FIT="python $SCR_DIR/walk_fit.py"
WALK_FIT=/user/hduser/sim_simple/walkfit_output

# Step 7: Given all the t-statistics, moderate them and emit moderated
#         t-stats
EBAYES_AGGR="cat"
EBAYES="python $SCR_DIR/ebayes.py"
EBAYES_OUT=/user/hduser/sim_simple/ebayes_output

# Step 8: Given all the moderated t-statistics, calculate the HMM
#         parameters to use in the next step
HMM_PARAMS_AGGR="cat"
HMM_PARAMS="python $SCR_DIR/hmm_params.py"
#HMM_PARAMS_OUT=/user/hduser/sim_simple/hmm_params_output

# Step 9: Given sorted bins of moderated t-statistics, and HMM
#         parameters, run the HMM
HMM_AGGR1="sort -n -k2,2"
HMM_AGGR2="sort -s -k1,1"
HMM="python $SCR_DIR/hmm.py"
#HMM_OUT=/user/hduser/sim_simple/hmm_output

# Temporary files so we can form a DAG
WALK_IN_TMP=${TMPDIR}walk_in.tsv
HMM_IN_TMP=${INTERMEDIATE_DIR}hmm_in.tsv

echo "Temporary file for walk_fit.py input is '$WALK_IN_TMP'"
echo "Temporary file for hmm.py input is '$HMM_IN_TMP'"

# echo $ALIGN_ARGS
# echo $TORNADO
# #FILES=`find $TORNADO -name \*`
# FILES=`ls`
# echo $FILES
# FILES=$(echo $FILES|sed '/s a -file')
# echo $FILES

#copy files over to hdfs
#hadoop dfs -copyFromLocal *.tab5 $HADOOP_FILES
#hadoop dfs -mkdir $HADOOP_FILES/output

#Step 1 ALIGN
# hadoop dfs -rmr $ALIGN_OUT
# hadoop jar $STREAMING \
# -D mapred.reduce.tasks=0 \
# -file "$SCR_DIR/rnawesome/align.py" \
# -file "$SCR_DIR/bowtie/bowtie.py" \
# -file "$SCR_DIR/read/readlet.py" \
# -file "$SCR_DIR/read/truncate.py" \
# -file "$SCR_DIR/interval/interval.py" \
# -file "$SCR_DIR/interval/partition.py" \
# -file "$SCR_DIR/manifest/manifest.py" \
# -file "$SCR_DIR/sample/sample.py" \
# -file "$INDEXS/lambda_virus.rev.1.ebwt" \
# -file "$INDEXS/lambda_virus.rev.2.ebwt" \
# -file "$INDEXS/lambda_virus.1.ebwt" \
# -file "$INDEXS/lambda_virus.2.ebwt" \
# -file "$INDEXS/lambda_virus.3.ebwt" \
# -file "$INDEXS/lambda_virus.4.ebwt" \
# -file "$INDEXS/lambda_virus.fa" \
# -file "$BOWTIE" \
# -mapper "$ALIGN_ARGS" \
# -input $HADOOP_FILES/*.tab5 -output $ALIGN_OUT

#Check $? after completion.  If not 0, then print an error message and quit
# if [$? -ne 0]
# then
#     echo "Alignment step failed, now exiting"
#     exit 0
# fi

#Step 2 SPLICE and MERGE
hadoop dfs -rmr $SPLICE_OUT
hadoop jar $STREAMING \
    -file "$SCR_DIR/rnawesome/splice.py" \
    -file "$SCR_DIR/interval/interval.py" \
    -file "$SCR_DIR/interval/partition.py" \
    -file "$SCR_DIR/manifest/manifest.py" \
    -file "$SCR_DIR/sample/sample.py" \
    -file "$SCR_DIR/rnawesome/merge.py" \
    -jobconf num.key.fields.for.partition=1 \
    -jobconf stream.num.map.output.key.fields=2 \
    -partitioner org.apache.hadoop.mapred.lib.KeyFieldBasedPartitioner \
    -mapper "$SPLICE_ARGS" \
    -reducer "$MERGE" \
    -input $ALIGN_OUT/*part* -output $SPLICE_OUT


#Step 3 WALK_PRENORM (sort)
hadoop dfs -rmr $WALK_OUT
hadoop jar $STREAMING \
    -jobconf num.key.fields.for.partition=1 \
    -jobconf stream.num.map.output.key.fields=2 \
    -partitioner org.apache.hadoop.mapred.lib.KeyFieldBasedPartitioner \
    -mapper org.apache.hadoop.mapred.lib.IdentityMapper \
    -reducer org.apache.hadoop.mapred.lib.IdentityReducer \
    -input $SPLICE_OUT/*part* -output $WALK_OUT


#Step 4 NORMALIZE
hadoop dfs -rmr $SPLICE_OUT
hadoop jar $STREAMING \
    -file "$SCR_DIR/rnawesome/walk_prenorm.py" \
    -file "$SCR_DIR/interval/partition.py" \
    -file "$SCR_DIR/manifest/manifest.py" \
    -file "$SCR_DIR/struct/circular.py" \
    -jobconf num.key.fields.for.partition=1 \
    -jobconf stream.num.map.output.key.fields=2 \
    -partitioner org.apache.hadoop.mapred.lib.KeyFieldBasedPartitioner \
    -mapper org.apache.hadoop.mapred.lib.IdentityMapper \
    -reducer org.apache.hadoop.mapred.lib.IdentityReducer \
    -mapper "$WALK_ARGS" \
    -reducer "$MERGE" \
    -input $SPLICE_OUT/*part* -output $WALK_PRENORM




    # -partitioner org.apache.hadoop.mapred.lib.KeyFieldBasedPartitioner \
    # -D stream.num.map.output.key.fields=2 \
    # -D mapred.output.key.comparator.class=org.apache.hadoop.mapred.lib.KeyFieldBasedComparator \
    # -D mapred.text.key.partitioner.options=-k1,1 \
    # -D mapred.text.key.comparator.options=-k2,2n \


#Check $? after completion.  If not 0, then print an error message and quit
# if [$? -ne 0]
# then
#     echo "Splice step failed, now exiting"
#     exit 0
# fi


# hadoop dfs -rmr $TMPOUT
# hadoop jar $STREAMING \
# -D mapred.reduce.tasks=0 \
# -mapper $ALIGN \
# -file "find $TORNADO -name \*"
# -input $HADOOP_FILES/*.tab5 -output $TMPOUT


# bin/hadoop jar contrib/streaming/hadoop-*streaming*.jar \
# -file /home/hduser/mapper.py    -mapper /home/hduser/mapper.py \
# -file /home/hduser/reducer.py   -reducer /home/hduser/reducer.py \
# -input /user/hduser/gutenberg/* -output /user/hduser/gutenberg-output


# cat *.tab5 \
# 	| $ALIGN_AGGR | $ALIGN \
# 		--bowtieArgs '-v 2 -m 1 -p 6' \
# 		--bowtieExe /home/hduser/software/bowtie-0.12.8/bowtie \
# 		--bowtieIdx=../fasta/lambda_virus \
# 		--readletLen 20 \
# 		--readletIval 2 \
# 		--manifest $MANIFEST_FN \
# 	| $SPLICE_AGGR | $SPLICE \
# 		--ntasks=$NTASKS \
# 	--genomeLen=$GENOME_LEN \
# 		--manifest $MANIFEST_FN \
# 	| $MERGE_AGGR1 | $MERGE_AGGR2 | $MERGE \
# 	| $WALK_PRENORM_AGGR1 | $WALK_PRENORM_AGGR2 | tee $WALK_IN_TMP | $WALK_PRENORM \
# 		--manifest $MANIFEST_FN \
# 		--ntasks=$NTASKS \
# 		--genomeLen=$GENOME_LEN \
# 	| $NORMALIZE_AGGR | $NORMALIZE \
# 		--percentile 0.75 \
# 	| $NORMALIZE_POST_AGGR | $NORMALIZE_POST \
# 		--manifest $MANIFEST_FN > ${INTERMEDIATE_DIR}norm.tsv

# cp $WALK_IN_TMP ${INTERMEDIATE_DIR}walk_in_input.tsv

# cat $WALK_IN_TMP \
# 	| $WALK_FIT \
# 		--ntasks=$NTASKS \
# 		--genomeLen=$GENOME_LEN \
# 		--seed=777 \
# 		--normals ${INTERMEDIATE_DIR}norm.tsv \
# 	| $EBAYES_AGGR | $EBAYES \
# 		--ntasks=$NTASKS \
# 		--genomeLen=$GENOME_LEN \
# 		--hmm-overlap=$HMM_OVERLAP \
# 	| tee $HMM_IN_TMP | $HMM_PARAMS_AGGR | $HMM_PARAMS \
# 		--null \
# 		--out ${INTERMEDIATE_DIR}hmm_params.tsv 

# cat $HMM_IN_TMP \
# 	| $HMM_AGGR1 | $HMM_AGGR2 | $HMM \
# 		--ntasks=$NTASKS \
# 		--genomeLen=$GENOME_LEN \
# 		--params ${INTERMEDIATE_DIR}hmm_params.tsv \
# 		--hmm-overlap=$HMM_OVERLAP

# echo DONE

# echo "Normalization file:"
# cat ${INTERMEDIATE_DIR}norm.tsv

# echo "HMM parameter file:"
# cat ${INTERMEDIATE_DIR}hmm_params.tsv
