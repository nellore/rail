#!/usr/bin/env bash
# This script generates some VERY small human examples for the Rail-RNA manual
# It uses generate_bioreps.py with different parameters than those used for the Rail paper
DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
cd $DIR
mkdir -p small_human_examples
cd small_human_examples
# Grab the last two lines of the 112 manifest
tail -n 2 ../GEUVADIS_112.manifest >src.manifest
cd ..
python generate_bioreps.py -c 20000 --single-end -l 76 -p 2 -o ./small_human_examples --manifest ./small_human_examples/src.manifest
# Rename files, removing redundant labels at the end
mv NA07048_male_CEU_UU_6-1-2.bed NA07048_male_CEU_UU.bed
mv NA07048_male_CEU_UU_6-1-2.fastq NA07048_male_CEU_UU.fastq
mv NA19129_female_YRI_UU_6-1-1.bed NA19129_female_YRI_UU.bed
mv NA19129_female_YRI_UU_6-1-1.fastq NA19129_female_YRI_UU.fastq
rm src.manifest
read -r -d '' var <<EOF
NA07048_male_CEU_UU.fastq\t0\tNA07048_male_CEU_UU-1-1\nNA19129_female_YRI_UU.fastq\t0\tNA19129_female_YRI_UU_6-1-1
EOF
echo -e $var >example.manifest
# Enter password
scp *.fastq *.bed verve@verve.webfactional.com:/home/verve/webapps/burn1
read -r -d '' var <<EOF
http://verve.webfactional.com/NA07048_male_CEU_UU.fastq\t0\tNA07048_male_CEU_UU-1-1\nhttp://verve.webfactional.com/NA19129_female_YRI_UU.fastq\t0\tNA19129_female_YRI_UU_6-1-1
EOF
echo -e $var >web_example.manifest
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
bad.fastq\t0\tbad-1-1
EOF
echo -e $var >bad.manifest