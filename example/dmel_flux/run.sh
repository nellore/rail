NTASKS=10000
HMM_OVERLAP=30
PERMUTATIONS=5
READLET_LEN=35
READLET_IVAL=5
IGENOME=$IGENOMES_HOME/Drosophila_melanogaster/UCSC/dm3
BOWTIE_IDX=$IGENOMES_HOME/Drosophila_melanogaster/UCSC/dm3/Sequence/BowtieIndex/genome
GENOME=$IGENOMES_HOME/Drosophila_melanogaster/UCSC/dm3/Sequence/WholeGenomeFasta/genome.fa
FASTA_IDX=$IGENOMES_HOME/Drosophila_melanogaster/UCSC/dm3/Sequence/WholeGenomeFasta/genome.fa.fai
INDEX1=$IGENOMES_HOME/Drosophila_melanogaster/UCSC/dm3/Sequence/BowtieIndex/genome.1.ebwt
INDEX2=$IGENOMES_HOME/Drosophila_melanogaster/UCSC/dm3/Sequence/BowtieIndex/genome.2.ebwt
INDEX3=$IGENOMES_HOME/Drosophila_melanogaster/UCSC/dm3/Sequence/BowtieIndex/genome.3.ebwt
INDEX4=$IGENOMES_HOME/Drosophila_melanogaster/UCSC/dm3/Sequence/BowtieIndex/genome.4.ebwt
INDEX5=$IGENOMES_HOME/Drosophila_melanogaster/UCSC/dm3/Sequence/BowtieIndex/genome.rev.1.ebwt
INDEX6=$IGENOMES_HOME/Drosophila_melanogaster/UCSC/dm3/Sequence/BowtieIndex/genome.rev.2.ebwt
INPUT_DIR=preprocessed_reads
GENOME_LEN=129682844
MANIFEST_FN=dmel_flux.manifest
INTERMEDIATE_DIR=intermediate
SITES_FILE=
SPLICE_OVERLAP=10
RADIUS=10
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

SCR_DIR=$RAIL_HOME/src/rail-rna

# Step 1: Readletize input reads and use Bowtie to align the readlets 
ALIGN_AGGR="cat"
ALIGN="python $SCR_DIR/align.py"

# Step 2: Collapse identical intervals from same sample
# In Hadoop, we want to partition by first field then sort by second
# -partitioner org.apache.hadoop.mapred.lib.KeyFieldBasedPartitioner
# -D stream.num.map.output.key.fields=2
# -D mapred.text.key.partitioner.options=-k1,1
MERGE_AGGR1="grep '^exon'"
MERGE_AGGR2="cut -f 2-"
MERGE_AGGR3="sort -n -k2,2"
MERGE_AGGR4="sort -s -k1,1"
MERGE="python $SCR_DIR/merge.py"

NORMALIZE_PRE_AGGR="sort -n -k1,3"
NORMALIZE_PRE="python $SCR_DIR/normalize_pre.py"

INTRON_AGGR1="grep '^intron'"
INTRON_AGGR2="cut -f 2-"
INTRON_AGGR3="sort -n -k2,2"
INTRON_AGGR4="sort -s -k1,1"
INTRON="python $SCR_DIR/intron2.py"
UTIL=$RAIL_HOME/src/util
SITE2BED="python $UTIL/site2bed2.py"
EXONS2BED="python $UTIL/exons2bed.py"
INTRONS2BED="python $UTIL/introns2bed.py"
JUNCS2BED="python $UTIL/junc2site.py"

SITE_AGGR1="grep '^site'"
SITE_AGGR2="cut -f 2-"

CO_AGGR1="grep '^cooccurence'"
CO_AGGR2="cut -f 2-"

# Step 4: For all samples, take all coverage tuples for the sample and
#         from them calculate a normalization factor
NORMALIZE_AGGR="sort -k1,3"
NORMALIZE="python $SCR_DIR/normalize2.py"
SAMPLE_OUT=$INTERMEDIATE_DIR/samples
mkdir -p $SAMPLE_OUT
rm -f $SAMPLE_OUT/*
CHROM_SIZES=$PWD/chrom.sizes
cat $FASTA_IDX | cut -f -2 > $CHROM_SIZES
# Step 5: Collect all the norm factors together and write to file
# In Hadoop no partitioning or sorting (mapper only)
NORMALIZE_POST_AGGR="cat"
NORMALIZE_POST="python $SCR_DIR/normalize_post.py"

# Step 6: Walk over genome windows again (taking output from Step 2)
#         but this time, calculate per-position coverage vectors and
#         fit a linear model to each
WALK_FIT="python $SCR_DIR/walk_fit.py"

# Step 7: Given all the t-statistics, moderate them and emit moderated
#         t-stats
EBAYES_AGGR="cat"
EBAYES="python $SCR_DIR/ebayes.py"

# Step 8: Given all the moderated t-statistics, calculate the HMM
#         parameters to use in the next step
HMM_PARAMS_AGGR="cat"
HMM_PARAMS="python $SCR_DIR/hmm_params.py"

# Step 9: Given sorted bins of moderated t-statistics, and HMM
#         parameters, run the HMM
HMM_AGGR1="sort -n -k2,2"
HMM_AGGR2="sort -s -k1,1"
HMM="python $SCR_DIR/hmm.py"

# Step 10: Given the permutation outputs from the HMM step, this 
#          stores the coverage vectors for each permutation into separate files
PATH_AGGR1="sort -n -k3,3"
PATH_AGGR2="sort -s -k2,2"
PATH_AGGR3="sort -n -k1,1"
AGGR_PATH="python $SCR_DIR/aggr_path.py"
PERM_OUT=$INTERMEDIATE_DIR/permutations
mkdir -p $PERM_OUT

# Temporary files so we can form a DAG
WALK_IN_TMP=${TMPDIR}walk_in.tsv
HMM_IN_TMP=${TMPDIR}hmm_in.tsv

echo "Temporary file for walk_fit.py input is '$WALK_IN_TMP'" 1>&2
echo "Temporary file for hmm.py input is '$HMM_IN_TMP'" 1>&2

gzip -dc $INPUT_DIR/* \
	| $ALIGN_AGGR | $ALIGN \
		--ntasks=$NTASKS \
		--genomeLen=$GENOME_LEN \
		--bowtieExe $BOWTIE_EXE \
		--bowtieIdx=$BOWTIE_IDX \
		--readletLen $READLET_LEN \
		--readletIval $READLET_IVAL \
		--refseq=$GENOME \
 		--faidx=$FASTA_IDX \
		--splice-overlap=$SPLICE_OVERLAP \
		--exon-differentials \
		--exon-intervals \
		--verbose \
		--intron-partition-overlap=30 \
		-- -v 2 -m 1 -p 10 \
		| tee ${INTERMEDIATE_DIR}/align_out.tsv \
	| grep '^exon_diff' | cut -f 2- | $NORMALIZE_PRE_AGGR | tee $WALK_IN_TMP | $NORMALIZE_PRE \
		--ntasks=$NTASKS \
		--genomeLen=$GENOME_LEN \
	| $NORMALIZE_AGGR | tee ${INTERMEDIATE_DIR}/pre_normalize.tsv | $NORMALIZE \
		--percentile 0.75 \
		--out_dir $SAMPLE_OUT \
		--bigbed_exe $BIGBED_EXE \
		--chrom_sizes $CHROM_SIZES \
	| $NORMALIZE_POST_AGGR | $NORMALIZE_POST \
		--manifest $MANIFEST_FN > ${INTERMEDIATE_DIR}/norm.tsv

cp $WALK_IN_TMP ${INTERMEDIATE_DIR}/walk_in_input.tsv

cat ${INTERMEDIATE_DIR}/align_out.tsv \
	| grep '^intron' | $INTRON_AGGR2 | $INTRON_AGGR3 | $INTRON_AGGR4 | $INTRON \
		--ntasks=$NTASKS \
		--genomeLen=$GENOME_LEN \
		--refseq=$GENOME \
		--readletIval $READLET_IVAL \
		--readletLen $READLET_LEN  \
		--intron-partition-overlap=30 \
		--cluster-radius=$RADIUS \
		--per-span \
		--per-site \
		> ${INTERMEDIATE_DIR}/intron_out.tsv 

cat ${INTERMEDIATE_DIR}/intron_out.tsv \
    | grep '^site' | $SITE_AGGR2 | $SITE2BED > ${INTERMEDIATE_DIR}/splice_sites.bed

#cat ${INTERMEDIATE_DIR}/intron_out.tsv \
#    | grep '^cooccurence' | $CO_AGGR2 > ${INTERMEDIATE_DIR}/cooccurences.tab

cat ${INTERMEDIATE_DIR}/align_out.tsv \
    | grep '^intron' | $INTRON_AGGR2 | $INTRON_AGGR3 | $INTRON_AGGR4 | $INTRONS2BED > ${INTERMEDIATE_DIR}/flanks.bed

# cat ${INTERMEDIATE_DIR}/align_out.tsv \
#     | grep '^exon_ival' | $INTRON_AGGR2 | $INTRON_AGGR3 | $INTRON_AGGR4 | $EXONS2BED | $JUNCS2BED > ${INTERMEDIATE_DIR}/exons.bed

cat ${INTERMEDIATE_DIR}/align_out.tsv \
    | grep '^exon_ival' | cut -f 2- | sort -k1,3 | awk '{split($0,array,";")}{print array[1]"\t"$2"\t"$3"\t"$4}' > ${INTERMEDIATE_DIR}/exon.bed

#cat ${INTERMEDIATE_DIR}/align_out.tsv \
#    | grep '^exon' | $INTRON_AGGR2 | $INTRON_AGGR3 | $INTRON_AGGR4 > ${INTERMEDIATE_DIR}/align_out_exons.tsv
 
# cat $WALK_IN_TMP \  
# 	| tee ${INTERMEDIATE_DIR}/walk_fit_in.tsv | $WALK_FIT \
# 		--ntasks=$NTASKS \
# 		--genomeLen=$GENOME_LEN \
# 		--seed=777 \
# 		--permutations=$PERMUTATIONS \
# 		--permutations-out=${INTERMEDIATE_DIR}/permutations.tsv \
# 		--normals=${INTERMEDIATE_DIR}/norm.tsv \
# 	| $EBAYES_AGGR | tee ${INTERMEDIATE_DIR}/ebayes_in.tsv | $EBAYES \
# 		--ntasks=$NTASKS \
# 		--genomeLen=$GENOME_LEN \
# 		--hmm-overlap=$HMM_OVERLAP \
# 	| tee ${INTERMEDIATE_DIR}/hmm_in.tsv | $HMM_PARAMS_AGGR | $HMM_PARAMS \
# 		--null \
# 		--out ${INTERMEDIATE_DIR}/hmm_params.tsv 

# cat ${INTERMEDIATE_DIR}/hmm_in.tsv \
# 	| $HMM_AGGR1 | $HMM_AGGR2 | $HMM \
# 		--ntasks=$NTASKS \
# 		--genomeLen=$GENOME_LEN \
# 		--params ${INTERMEDIATE_DIR}/hmm_params.tsv \
# 		--hmm-overlap=$HMM_OVERLAP \
# 	| tee ${INTERMEDIATE_DIR}/hmm_out.tsv > hmm.out

# cat hmm.out | $PATH_AGGR1 | $PATH_AGGR2 | $PATH_AGGR3 | $AGGR_PATH \
#                 --out_dir $PERM_OUT \
#                 --bigbed_exe $BIGBED_EXE \
#                 --chrom_sizes $CHROM_SIZES \


# echo DONE 1>&2

# echo "Normalization file:" 1>&2
# cat ${INTERMEDIATE_DIR}/norm.tsv 1>&2

# echo "HMM parameter file:" 1>&2
# cat ${INTERMEDIATE_DIR}/hmm_params.tsv 1>&2

