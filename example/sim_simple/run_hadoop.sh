# Parameters
GENOME_LEN=48502
NTASKS=10
HMM_OVERLAP=30

#Hadoop file locations
HADOOP_FILES=/user/hduser/sim_simple
STREAMING=$HADOOP_HOME/contrib/streaming/hadoop-streaming*.jar

#Note: Must copy over the entire tornado directory
TORNADO=/damsl/projects/myrna2/software/tornado
EXAMPLE=$TORNADO/example
SCR_DIR=$TORNADO/src
RNAWESOME=$SCR_DIR/rnawesome

#Bowtie locations
BOWTIE=$TORNADO/tools/bowtie-0.12.8/bowtie
INDEXS=$EXAMPLE/fasta
BOWTIE_IDX=$INDEXS/lambda_virus


SIM_SIMPLE=$EXAMPLE/sim_simple
MANIFEST_FN=$SIM_SIMPLE/eg1.manifest

INTERMEDIATE_DIR=$SIM_SIMPLE/intermediate
#rm -r intermediate
#mkdir -p intermediate

# Step 1a: Readletize input reads and use Bowtie to align the readlets 
ALIGN_AGGR="cat"
ALIGN="python $RNAWESOME/align.py"
ALIGN_ARGS=''$ALIGN' --bowtieArgs '\''-v 2 -m 1 -p 6'\'' --bowtieExe '$BOWTIE' --bowtieIdx='$BOWTIE_IDX' --readletLen 20 --readletIval 2 --manifest '$MANIFEST_FN''
ALIGN_OUT=/user/hduser/sim_simple/align_output

# Step 1b: Bring together all the readlet alignment intervals for a given read and 
SPLICE_AGGR="cat"
SPLICE="python $RNAWESOME/splice.py"
SPLICE_ARGS=''$SPLICE' --ntasks='$NTASKS' --genomeLen='$GENOME_LEN' --manifest '$MANIFEST_FN''
SPLICE_OUT=/user/hduser/sim_simple/splice_output

# Steps 1a and 1b happen together in the same map step in practice

# Step 2: Collapse identical intervals from same sample
MERGE_AGGR1="sort -n -k2,2"
MERGE_AGGR2="sort -s -k1,1"
MERGE="python $RNAWESOME/merge.py"
MERGE_OUT=/user/hduser/sim_simple/merge_output

# Step 3: Walk over genome windows and emit per-sample, per-position
#         coverage tuples
WALK_PRENORM_AGGR1="sort -n -k2,2"
WALK_PRENORM_AGGR2="sort -s -k1,1"
WALK_PRENORM="python $RNAWESOME/walk_prenorm.py"
WALK_PRENORM_OUT=/user/hduser/sim_simple/walk_prenorm_output
WALK_ARGS=''$WALK_PRENORM' --ntasks='$NTASKS' --genomeLen='$GENOME_LEN' --manifest '$MANIFEST_FN''

# Step 4: For all samples, take all coverage tuples for the sample and
#         from them calculate a normalization factor
NORMALIZE_AGGR="sort -k1,1"
NORMALIZE="python $RNAWESOME/normalize.py"
NORMALIZE_OUT=/user/hduser/sim_simple/normalize_output

# Step 5: Collect all the norm factors together and write to file
NORMALIZE_POST_AGGR="cat"
NORMALIZE_POST="python $RNAWESOME/normalize_post.py"
NORMALIZE_POST_ARGS=''$NORMALIZE_POST' --manifest '$MANIFEST_FN''
NORMALIZE_POST_OUT=/user/hduser/sim_simple/normalize_post_output

# Step 6: Walk over genome windows again (taking output from Step 2)
#         but this time, calculate per-position coverage vectors and
#         fit a linear model to each
WALK_FIT="python $RNAWESOME/walk_fit.py"
WALK_FIT_ARGS=''$WALK_FIT' --ntasks='$NTASKS' --genomeLen='$GENOME_LEN' --seed=777 --normals '${INTERMEDIATE_DIR}/norm.tsv''
#WALK_FIT_ARGS=''$WALK_FIT' --ntasks='$NTASKS' --genomeLen='$GENOME_LEN' --seed=777 --normals '$NORMALIZE_POST_OUT/part-00000''
WALK_FIT_OUT=/user/hduser/sim_simple/walkfit_output 
# Step 7: Given all the t-statistics, moderate them and emit moderated
#         t-stats
EBAYES_AGGR="cat"
EBAYES="python $RNAWESOME/ebayes.py"
EBAYES_ARGS=''$EBAYES' --ntasks='$NTASKS' --genomeLen='$GENOME_LEN' --hmm-overlap='$HMM_OVERLAP''
EBAYES_OUT=/user/hduser/sim_simple/ebayes_output

# Step 8: Given all the moderated t-statistics, calculate the HMM
#         parameters to use in the next step
HMM_PARAMS_AGGR="cat"
HMM_PARAMS="python $RNAWESOME/hmm_params.py"
HMM_PARAMS_ARGS=''$HMM_PARAMS' --null'
HMM_PARAMS_OUT=/user/hduser/sim_simple/hmm_params_output

# Step 9: Given sorted bins of moderated t-statistics, and HMM
#         parameters, run the HMM
HMM_AGGR1="sort -n -k2,2"
HMM_AGGR2="sort -s -k1,1"
HMM="python $RNAWESOME/hmm.py"
HMM_ARGS=''$HMM' --ntasks='$NTASKS' --genomeLen='$GENOME_LEN' --params '${INTERMEDIATE_DIR}/hmm_params.tsv' --hmm-overlap='$HMM_OVERLAP''

HMM_OUT=/user/hduser/sim_simple/hmm_output

# Temporary files so we can form a DAG
WALK_IN_TMP=${TMPDIR}walk_in.tsv
HMM_IN_TMP=${INTERMEDIATE_DIR}/hmm_in.tsv

echo "Temporary file or walk_fit.py input is '$WALK_IN_TMP'"
echo "Temporary file for hmm.py input is '$HMM_IN_TMP'"

# echo $ALIGN_ARGS
# echo $TORNADO
# #FILES=`find $TORNADO -name \*`
# FILES=`ls`
# echo $FILES
# FILES=$(echo $FILES|sed '/s a -file')
# echo $FILES

#copy files over to hdfs
# hadoop dfs -mkdir $HADOOP_FILES
# hadoop dfs -copyFromLocal *.tab5 $HADOOP_FILES

#Step 1 ALIGN
hadoop dfs -rmr $ALIGN_OUT
hadoop jar $STREAMING \
    -D mapred.reduce.tasks=0 \
    -file "$SCR_DIR/rnawesome/align.py" \
    -file "$SCR_DIR/bowtie/bowtie.py" \
    -file "$SCR_DIR/read/readlet.py" \
    -file "$SCR_DIR/read/truncate.py" \
    -file "$SCR_DIR/interval/interval.py" \
    -file "$SCR_DIR/interval/partition.py" \
    -file "$SCR_DIR/manifest/manifest.py" \
    -file "$SCR_DIR/sample/sample.py" \
    -file "$INDEXS/lambda_virus.rev.1.ebwt" \
    -file "$INDEXS/lambda_virus.rev.2.ebwt" \
    -file "$INDEXS/lambda_virus.1.ebwt" \
    -file "$INDEXS/lambda_virus.2.ebwt" \
    -file "$INDEXS/lambda_virus.3.ebwt" \
    -file "$INDEXS/lambda_virus.4.ebwt" \
    -file "$INDEXS/lambda_virus.fa" \
    -file "$BOWTIE" \
    -mapper "$ALIGN_ARGS" \
    -input $HADOOP_FILES/*.tab5 -output $ALIGN_OUT

# #Check $? after completion.  If not 0, then print an error message and quit
# # if [$? -ne 0]
# # then
# #     echo "Alignment step failed, now exiting"
# #     exit 0
# # fi

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


#Sort MERGE output
hadoop dfs -rmr $MERGE_OUT
hadoop jar $STREAMING \
    -D mapred.text.key.partitioner.options=-k1,1 \
    -D stream.num.map.output.key.fields=2 \
    -partitioner org.apache.hadoop.mapred.lib.KeyFieldBasedPartitioner \
    -mapper cat \
    -reducer org.apache.hadoop.mapred.lib.IdentityReducer \
    -input $SPLICE_OUT/*part* -output $MERGE_OUT

# #Step 3 and 4  and NORMALIZE
# hadoop dfs -rmr $NORMALIZE_OUT
# hadoop jar $STREAMING \
#     -D mapred.text.key.partitioner.options=-k1,1 \
#     -D stream.num.map.output.key.fields=2 \
#     -file "$SCR_DIR/rnawesome/normalize.py" \
#     -file "$SCR_DIR/rnawesome/walk_prenorm.py" \
#     -file "$SCR_DIR/interval/partition.py" \
#     -file "$SCR_DIR/manifest/manifest.py" \
#     -file "$SCR_DIR/struct/circular.py" \
#     -partitioner org.apache.hadoop.mapred.lib.KeyFieldBasedPartitioner \
#     -reducer "$WALK_ARGS" \
#     -mapper cat \
#     -input $MERGE_OUT/*part* -output $NORMALIZE_OUT

#Step 3 WALK_PRENORM
hadoop dfs -rmr $WALK_PRENORM_OUT
hadoop jar $STREAMING \
    -D mapred.reduce.tasks=0 \
    -file "$SCR_DIR/rnawesome/normalize.py" \
    -file "$SCR_DIR/rnawesome/walk_prenorm.py" \
    -file "$SCR_DIR/interval/partition.py" \
    -file "$SCR_DIR/manifest/manifest.py" \
    -file "$SCR_DIR/struct/circular.py" \
    -mapper "$WALK_ARGS" \
    -input $MERGE_OUT/*part* -output $WALK_PRENORM_OUT


#Step 4 NORMALIZE
hadoop dfs -rmr $NORMALIZE_OUT
hadoop jar $STREAMING \
    -file "$SCR_DIR/rnawesome/normalize.py" \
    -file "$SCR_DIR/interval/partition.py" \
    -file "$SCR_DIR/manifest/manifest.py" \
    -file "$SCR_DIR/struct/circular.py" \
    -mapper cat \
    -reducer "$NORMALIZE" \
    -input $WALK_PRENORM_OUT/*part* -output $NORMALIZE_OUT


#Step 5 NORMALIZE_POST
hadoop dfs -rmr $NORMALIZE_POST_OUT
hadoop jar $STREAMING \
    -file "$SCR_DIR/rnawesome/normalize_post.py" \
    -file "$SCR_DIR/manifest/manifest.py" \
    -mapper cat \
    -reducer "$NORMALIZE_POST_ARGS" \
    -input $NORMALIZE_OUT/*part* -output $NORMALIZE_POST_OUT

#Copy files to local machine
hadoop dfs -copyToLocal $NORMALIZE_POST_OUT/part* ${INTERMEDIATE_DIR}/norm.tsv
hadoop dfs -copyToLocal $MERGE_OUT/part* ${INTERMEDIATE_DIR}/walk_in_input.tsv


#Step 6 WALK_FIT 
hadoop dfs -rmr $WALK_FIT_OUT
hadoop jar $STREAMING \
    -D mapred.text.key.partitioner.options=-k1,1 \
    -D stream.num.map.output.key.fields=2 \
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
    -D mapred.text.key.partitioner.options=-k1,1 \
    -D stream.num.map.output.key.fields=2 \
    -partitioner org.apache.hadoop.mapred.lib.KeyFieldBasedPartitioner \
    -file "$SCR_DIR/rnawesome/ebayes.py" \
    -file "$SCR_DIR/interval/partition.py" \
    -file "$SCR_DIR/struct/circular.py" \
    -mapper cat \
    -reducer "$EBAYES_ARGS" \
    -input $WALK_FIT_OUT/*part* -output $EBAYES_OUT


# # #Check $? after completion.  If not 0, then print an error message and quit
# # # if [$? -ne 0]
# # # then
# # #     echo "Splice step failed, now exiting"
# # #     exit 0
# # # fi

# # #Step 8 HMM_PARAMS
hadoop dfs -rmr $HMM_PARAMS_OUT
hadoop jar $STREAMING \
    -file "$SCR_DIR/rnawesome/hmm_params.py" \
    -mapper cat \
    -reducer "$HMM_PARAMS_ARGS" \
    -input $EBAYES_OUT/*part* -output $HMM_PARAMS_OUT

# #Copy HMM params to local machine
hadoop dfs -copyToLocal $HMM_PARAMS_OUT/part* ${INTERMEDIATE_DIR}/hmm_params.tsv

#Step 9 HMM
hadoop dfs -rmr $HMM_OUT
hadoop jar $STREAMING \
    -D mapred.text.key.partitioner.options=-k1,1 \
    -D stream.num.map.output.key.fields=2 \
    -file "$SCR_DIR/rnawesome/hmm.py" \
    -file "$SCR_DIR/interval/partition.py" \
    -file "$INTERMEDIATE_DIR/hmm_params.tsv" \
    -partitioner org.apache.hadoop.mapred.lib.KeyFieldBasedPartitioner \
    -mapper cat \
    -reducer "$HMM_ARGS" \
    -input $EBAYES_OUT/*part* -output $HMM_OUT


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
