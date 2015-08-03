#!/usr/bin/env bash
## Grabs results of experiment running Rail-RNA on 112 simulated datasets
## reflecting the coverage distributions of 112 GEUVADIS datasets from S3
## and computes performance measures on only those samples considered
## in run_all_small_data_sims_locally.sh. Also grabs transcript indexes
## and counts numbers of introns in each index with count_introns.py
## Requires AWS CLI
## $1: output directory of job WITH FILTER on S3
## $2: output directory of job WITHOUT FILTER on S3
## $3: directory in which to dump output
## $4: where to find Flux BEDs storing true read alignments from simulation
## $5: where to find Bowtie 1 index basename of genome
## Command we ran was sh grab_112_sim_results.sh s3://rail-eu-west-1/geuv112sim_v2.out s3://rail-eu-west-1/geuv112sim_v2.out.nofilter /dcl01/leek/data/railsims/112simresults /dcl01/leek/data/railsims /dcl01/leek/data/railsims/indexes_for_paper/genome

# Use the 20 samples that were randomly selected in create_single_sample_sim_commands.py
SAMPLES=( HG00115_male_GBR_UU_6-1-1 NA06984_male_CEU_UNIGE_1-1-1 HG00313_female_FIN_UNIGE_1-1-1 NA19095_female_YRI_LUMC_7-1-3 HG00117_male_GBR_LUMC_7-1-1 NA20768_female_TSI_HMGU_5-1-1 NA20582_female_TSI_ICMB_4-1-1 NA19130_male_YRI_HMGU_5-1-1 HG00139_male_GBR_LUMC_7-1-1 NA18486_male_YRI_LUMC_7-1-1 HG01790_female_GBR_MPIMG_3-1-1 NA12287_female_CEU_UNIGE_1-1-1 NA12287_female_CEU_MPIMG_3-1-2 HG00096_male_GBR_UNIGE_1-1-1 NA11831_male_CEU_LUMC_7-1-1 NA12874_male_CEU_UNIGE_1-1-1 HG00154_female_GBR_HMGU_5-1-1 NA07051_male_CEU_HMGU_5-1-1 NA12776_female_CEU_UU_6-1-1 NA20813_female_TSI_HMGU_5-1-1 )
# Minor idiosyncrasy in naming put the "sim" in the middle of the sample name
SAMPLENAMES=( HG00115_male_GBR_UU_6_sim-1-1 NA06984_male_CEU_UNIGE_1_sim-1-1 HG00313_female_FIN_UNIGE_1_sim-1-1 NA19095_female_YRI_LUMC_7_sim-1-3 HG00117_male_GBR_LUMC_7_sim-1-1 NA20768_female_TSI_HMGU_5_sim-1-1 NA20582_female_TSI_ICMB_4_sim-1-1 NA19130_male_YRI_HMGU_5_sim-1-1 HG00139_male_GBR_LUMC_7_sim-1-1 NA18486_male_YRI_LUMC_7_sim-1-1 HG01790_female_GBR_MPIMG_3_sim-1-1 NA12287_female_CEU_UNIGE_1_sim-1-1 NA12287_female_CEU_MPIMG_3_sim-1-2 HG00096_male_GBR_UNIGE_1_sim-1-1 NA11831_male_CEU_LUMC_7_sim-1-1 NA12874_male_CEU_UNIGE_1_sim-1-1 HG00154_female_GBR_HMGU_5_sim-1-1 NA07051_male_CEU_HMGU_5_sim-1-1 NA12776_female_CEU_UU_6_sim-1-1 NA20813_female_TSI_HMGU_5_sim-1-1 )

WITHFILTER=$1
WITHOUTFILTER=$2
OUTPUTDIR=$3
DATADIR=$4
BOWTIEIDX=$5

BOWTIE2INSPECT=/home/student/anellor1/raildotbio/bowtie2-2.2.5/bowtie2-inspect
RAILHOME=/home/student/anellor1/rail
PYTHON=/home/student/anellor1/raildotbio/pypy-2.5-linux_x86_64-portable/bin/pypy
SAMTOOLS=/home/student/anellor1/raildotbio/samtools-1.2/samtools
PERFORMANCE=perform

mkdir -p $OUTPUTDIR
cd $OUTPUTDIR
# Grab intron indexes with and without filter
mkdir -p withfilter
cd withfilter
aws s3 cp ${WITHFILTER}/cross_sample_results/isofrags.tar.gz ./
tar xvzf isofrags.tar.gz
$PYTHON $RAILHOME/eval/count_introns.py --bowtie1-idx $BOWTIEIDX --bowtie2-idx isofrags -t ${DATADIR} --bowtie2-inspect ${BOWTIE2INSPECT} >intron_count &
cd ..
mkdir -p withoutfilter
cd withoutfilter
aws s3 cp ${WITHOUTFILTER}/cross_sample_results/isofrags.tar.gz ./
tar xvzf isofrags.tar.gz
$PYTHON $RAILHOME/eval/count_introns.py --bowtie1-idx $BOWTIEIDX --bowtie2-idx isofrags -t ${DATADIR} --bowtie2-inspect ${BOWTIE2INSPECT} >intron_count &
cd ..
# Download/compute performances for all 20 samples with and without filter
for i in $(seq 0 19)
do
	cd withfilter
	mkdir -p ${SAMPLES[$i]}
	cd ${SAMPLES[$i]}
	for j in $(aws s3 ls ${WITHFILTER}/alignments/alignments.${SAMPLENAMES[$i]} | grep -v .bai | tr -s '[:blank:]' '\t' | cut -f4)
	do
		aws s3 cp ${WITHFILTER}/alignments/$j ./
	done
	(for k in *.bam; do $SAMTOOLS view $k; done | $PYTHON $RAILHOME/eval/spliced_read_recovery_performance.py -t ${DATADIR}/${SAMPLES[$i]}_sim.bed >$PERFORMANCE 2>${PERFORMANCE}_summary) &
	(for k in *.bam; do $SAMTOOLS view $k; done | $PYTHON $RAILHOME/eval/intron_recovery_performance.py -t ${DATADIR}/${SAMPLES[$i]}_sim.bed >${PERFORMANCE}_intron_recovery_summary) &
	(for k in *.bam; do $SAMTOOLS view $k; done | $PYTHON $RAILHOME/eval/mapping_accuracy.py -t ${DATADIR}/${SAMPLES[$i]}_sim.bed >${PERFORMANCE}_mapping_accuracy_summary) &
	(for k in *.bam; do $SAMTOOLS view $k; done | $PYTHON $RAILHOME/eval/mapping_accuracy.py -t ${DATADIR}/${SAMPLES[$i]}_sim.bed -c 0.1 >${PERFORMANCE}_mapping_accuracy_SC_summary) &
	cd ../..
	cd withoutfilter
	mkdir -p ${SAMPLES[$i]}
	cd ${SAMPLES[$i]}
	for j in $(aws s3 ls ${WITHOUTFILTER}/alignments/alignments.${SAMPLENAMES[$i]} | grep -v .bai | tr -s '[:blank:]' '\t' | cut -f4)
	do
		aws s3 cp ${WITHOUTFILTER}/alignments/$j ./
	done
	(for k in *.bam; do $SAMTOOLS view $k; done | $PYTHON $RAILHOME/eval/spliced_read_recovery_performance.py -t ${DATADIR}/${SAMPLES[$i]}_sim.bed >$PERFORMANCE 2>${PERFORMANCE}_summary) &
	(for k in *.bam; do $SAMTOOLS view $k; done | $PYTHON $RAILHOME/eval/intron_recovery_performance.py -t ${DATADIR}/${SAMPLES[$i]}_sim.bed >${PERFORMANCE}_intron_recovery_summary) &
	(for k in *.bam; do $SAMTOOLS view $k; done | $PYTHON $RAILHOME/eval/mapping_accuracy.py -t ${DATADIR}/${SAMPLES[$i]}_sim.bed >${PERFORMANCE}_mapping_accuracy_summary) &
	(for k in *.bam; do $SAMTOOLS view $k; done | $PYTHON $RAILHOME/eval/mapping_accuracy.py -t ${DATADIR}/${SAMPLES[$i]}_sim.bed -c 0.1 >${PERFORMANCE}_mapping_accuracy_SC_summary) &
	cd ../..
	wait
done
