#!/usr/bin/env bash
# $1: number of cores
# $2: output directory -- SPECIFY FULL PATH
# $3: where to find sample fastqs from generate_bioreps.py
# $4: sample name; this is the prefix of "_sim.fastq"
# $5: scratch directory; files are written here first, and relevant output is copied back to $2
# Ex: taskset -c 0,1,2,3 sh run_single_sample_tophat_sim.sh 4 ./myoutput NA11829_male_CEU_UU_6-1-1 /tmp
# See generate_bioreps.py for how sample data was generated.

# Specify number of parallel processes for each program
CORES=$1

# Specify FULL PATH to output directory
MAINOUTPUT=$2
mkdir -p ${MAINOUTPUT}

# Specify data directory; fastqs should be of the form [SAMPLE NAME]_sim.fastq; Flux beds should be
# of form [SAMPLE_NAME]_sim.bed
DATADIR=$3

# Specify sample name at command line
SAMPLE=$4

# Temp dir
SCRATCH=$5
mkdir -p ${SCRATCH}

## Specify locations of executables
# Used version 2.1.0 of TopHat; wrapped version 2.2.5 of Bowtie2 and version 1.1.1 of Bowtie
TOPHAT=/scratch0/langmead-fs1/shared/tophat-2.1.0.Linux_x86_64/tophat2
# Specify Python executable/loc of get_junctions.py; PyPy 2.5.0 was used
PYTHON=/home/anellor1/raildotbio/pypy-2.5-linux_x86_64-portable/bin/pypy
# Samtools v1.2 was used
SAMTOOLS=/home/anellor1/raildotbio/samtools-1.2/samtools

# Specify log filename for recording times
TIMELOG=${MAINOUTPUT}/tophat_times.log

## Specify locations of reference-related files
## See create_indexes.sh for index creation script
# Bowtie indexes
BOWTIE1IDX=/scratch0/langmead-fs1/indexes_for_paper/genome
BOWTIE2IDX=/scratch0/langmead-fs1/indexes_for_paper/genome

# Generic name of file measuring performance of a given alignment strategy
# Performance is computed with spliced_read_recovery_performance.py; refer to that file for details
PERFORMANCE=perform

## Specify location of annotation
# This is Gencode v12, which may be obtained at ftp://ftp.sanger.ac.uk/pub/gencode/release_12/gencode.v12.annotation.gtf.gz
ANNOTATION=/scratch0/langmead-fs1/geuvadis_sim/gencode.v12.annotation.gtf

# Flux outputs paired-end reads in one file; split files here
echo 'Splitting Flux FASTQs...'
awk '(NR-1) % 8 < 4' $DATADIR/${SAMPLE}_sim.fastq >${SCRATCH}/${SAMPLE}_sim_left.fastq
awk '(NR-1) % 8 >= 4' $DATADIR/${SAMPLE}_sim.fastq >${SCRATCH}/${SAMPLE}_sim_right.fastq

cd $SCRATCH
mkdir -p ${SAMPLE}
cd ${SAMPLE}
mkdir -p tophat

cd $MAINOUTPUT
mkdir -p ${SAMPLE}
SAMPLEOUTPUT=${MAINOUTPUT}/${SAMPLE}

# Run simulations
OUTPUT=$SCRATCH/${SAMPLE}
echo 'Running TopHat on sample '${SAMPLE}' with no annotation and in paired-end mode...'
echo '#'${SAMPLE}' TopHat noann paired' >>$TIMELOG
time ($TOPHAT -o $OUTPUT/tophat/noann_paired -p $CORES $BOWTIE2IDX ${SCRATCH}/${SAMPLE}_sim_left.fastq ${SCRATCH}/${SAMPLE}_sim_right.fastq 2>&1) 2>>$TIMELOG
echo 'Computing precision and recall...'
($SAMTOOLS view $OUTPUT/tophat/noann_paired/accepted_hits.bam | $PYTHON $RAILHOME/eval/spliced_read_recovery_performance.py -t $DATADIR/${SAMPLE}_sim.bed -g >$OUTPUT/tophat/noann_paired/$PERFORMANCE 2>$OUTPUT/tophat/noann_paired/${PERFORMANCE}_summary) &
($SAMTOOLS view $OUTPUT/tophat/noann_paired/accepted_hits.bam | $PYTHON $RAILHOME/eval/intron_recovery_performance.py -t $DATADIR/${SAMPLE}_sim.bed >$OUTPUT/tophat/noann_paired/${PERFORMANCE}_intron_recovery_summary) &
($SAMTOOLS view $OUTPUT/tophat/noann_paired/accepted_hits.bam | $PYTHON $RAILHOME/eval/mapping_accuracy.py -t $DATADIR/${SAMPLE}_sim.bed >${PERFORMANCE}_mapping_accuracy_summary) &
($SAMTOOLS view $OUTPUT/tophat/noann_paired/accepted_hits.bam | $PYTHON $RAILHOME/eval/mapping_accuracy.py -t $DATADIR/${SAMPLE}_sim.bed -c 0.1 >${PERFORMANCE}_mapping_accuracy_SC_summary) &
wait
echo 'Running TopHat on sample '${SAMPLE}' with annotation and in paired-end mode...'
echo '#'${SAMPLE}' TopHat ann paired' >>$TIMELOG
time ($TOPHAT -o $OUTPUT/tophat/ann_paired -G $ANNOTATION -p $CORES $BOWTIE2IDX ${SCRATCH}/${SAMPLE}_sim_left.fastq ${SCRATCH}/${SAMPLE}_sim_right.fastq 2>&1) 2>>$TIMELOG
echo 'Computing precision and recall...'
($SAMTOOLS view $OUTPUT/tophat/ann_paired/accepted_hits.bam | $PYTHON $RAILHOME/eval/spliced_read_recovery_performance.py -t $DATADIR/${SAMPLE}_sim.bed -g >$OUTPUT/tophat/ann_paired/$PERFORMANCE 2>$OUTPUT/tophat/ann_paired/${PERFORMANCE}_summary) &
($SAMTOOLS view $OUTPUT/tophat/ann_paired/accepted_hits.bam | $PYTHON $RAILHOME/eval/intron_recovery_performance.py -t $DATADIR/${SAMPLE}_sim.bed >$OUTPUT/tophat/ann_paired/${PERFORMANCE}_intron_recovery_summary) &
($SAMTOOLS view $OUTPUT/tophat/ann_paired/accepted_hits.bam | $PYTHON $RAILHOME/eval/mapping_accuracy.py -t $DATADIR/${SAMPLE}_sim.bed >${PERFORMANCE}_mapping_accuracy_summary) &
($SAMTOOLS view $OUTPUT/tophat/ann_paired/accepted_hits.bam | $PYTHON $RAILHOME/eval/mapping_accuracy.py -t $DATADIR/${SAMPLE}_sim.bed -c 0.1 >${PERFORMANCE}_mapping_accuracy_SC_summary) &
# Move TopHat results to final destination
mv $OUTPUT/tophat $SAMPLEOUTPUT