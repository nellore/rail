#!/usr/bin/env bash
# This script generates some VERY small examples for the Rail-RNA manual
# It uses generate_bioreps.py with different parameters than those used for the Rail paper
# It also generates a dmel example without generate_bioreps.py
DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
cd $DIR
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
# Enter password
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
scp *.fastq *.bed verve@verve.webfactional.com:/home/verve/webapps/burn1
