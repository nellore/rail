#!/usr/bin/env bash
# This script generates some VERY small examples for the Rail-RNA manual
# It uses generate_bioreps.py with different parameters than those used for the Rail paper
# It also generates a dmel example without generate_bioreps.py
# $1: Flux executable
# $2: Drosophila gtf
# $3: Drosophila genome fasta dir
FLUX=$1
DM3GTF=$2
DM3GEN=$3
DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
cd $DIR
# Make "clean" Drosophila gtf; see http://sammeth.net/confluence/pages/viewpage.action?pageId=7013276
awk 'BEGIN{FS="\t";OFS="\t"}{split($NF,a," ");pfx="";s="";for(i=1;i<=length(a);i+=2){if(a[i]=="transcript_id"){pfx=a[i]" "a[i+1]}else{s=s" "a[i]" "a[i+1]}}if(pfx==""){print "[WARN] line "NR" without transcript_id!" > "/dev/stderr"}else{$NF=pfx""s;print$0} }' $2 > genes_clean.gtf
read -r -d '' var <<EOF
NB_MOLECULES\t5000000\n
LOAD_NONCODING\tYES\n
TSS_MEAN\t50\n
POLYA_SCALE\tNaN\n
POLYA_SHAPE\tNaN\n
FRAG_SUBSTRATE\tRNA\n
FRAG_METHOD\tUR\n
FRAG_UR_ETA\t350\n
FRAG_UR_D0\t1\n
RTRANSCRIPTION\tYES\n
RT_PRIMER\tRH\n
RT_LOSSLESS\tYES\n
RT_MIN\t500\n
RT_MAX\t5500\n
GC_MEAN\tNaN\n
PCR_PROBABILITY\t0.050000\n
FILTERING\tNO\n
READ_NUMBER\t50000\n
READ_LENGTH\t76\n
PAIRED_END\tYES\n
ERR_FILE\t76\n
FASTA\tYES\n
UNIQUE_IDS\tYES\n
GEN_DIR\t$DM3GEN\n
REF_FILE_NAME\tgenes_clean.gtf\n
SEED\t1
EOF
echo -e $var >dm3_example.par
$FLUX -x -l -s -p dm3_example.par
# Split Flux FASTQs
awk '(NR-1) % 8 < 4' dm3_example.fastq >dm3_example_1.fastq
awk '(NR-1) % 8 >= 4' dm3_example.fastq >dm3_example_2.fastq
rm dm3_example.fastq
read -r -d '' var <<EOF
http://verve.webfactional.com/dm3_example_1.fastq\t0\thttp://verve.webfactional.com/dm3_example_2.fastq\t0\tdm3_example-1-1
EOF
echo -e $var >hg19_example.manifest
# Grab the last two lines of the 112 manifest
tail -n 2 ../eval/GEUVADIS_112.manifest >src.manifest
cd ../eval
python generate_bioreps.py -c 20000 --single-end -l 76 -p 2 -o ../ex --manifest ../ex/src.manifest
# Rename files, removing redundant labels at the end
cd ../ex
mv NA07048_male_CEU_UU_6-1-2_sim.bed NA07048_male_CEU_UU.bed
mv NA07048_male_CEU_UU_6-1-2_sim.fastq NA07048_male_CEU_UU.fastq
mv NA19129_female_YRI_UU_6-1-1_sim.bed NA19129_female_YRI_UU.bed
mv NA19129_female_YRI_UU_6-1-1_sim.fastq NA19129_female_YRI_UU.fastq
rm src.manifest
rm *.lib
rm *.pro
rm *.log
read -r -d '' var <<EOF
http://verve.webfactional.com/NA07048_male_CEU_UU.fastq\t0\tNA07048_male_CEU_UU-1-1\nhttp://verve.webfactional.com/NA19129_female_YRI_UU.fastq\t0\tNA19129_female_YRI_UU_6-1-1
EOF
echo -e $var >hg19_example.manifest
# Make a bad FASTQ file to demonstrate how Rail spews errors and can be robust to them if the user wants
cat >bad.fastq <<EOF
@badrecord
ATACAGATGACAGATGACAGGGTAGAGACAAATAGACAGATGACGATGGACAGATGACAGATAGAACAGATAGAGA
+
IIIIIIIIIIIIII
@goodrecord
ATACAGATGACAGATGACAGGAGACAAATAGACAGATGACGATGGACAGATGACAGATAGAACAGATAGAGA
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
EOF
read -r -d '' var <<EOF
http://verve.webfactional.com/bad.fastq\t0\tbad-1-1
EOF
echo -e $var >bad.manifest
# Enter password
scp *.fastq *.bed *.manifest verve@verve.webfactional.com:/home/verve/webapps/burn1
