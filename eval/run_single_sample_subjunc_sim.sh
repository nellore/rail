#!/usr/bin/env bash
# $1: number of cores
# $2: output directory -- SPECIFY FULL PATH
# $3: where to find sample fastqs from generate_bioreps.py
# $4: sample name; this is the prefix of "_sim.fastq"
# $5: scratch directory; files are written here first, and relevant output is copied back to $2
# Ex: taskset -c 0,1,2,3 sh run_single_sample_subjunc_sim.sh 4 ./myoutput NA11829_male_CEU_UU_6-1-1 /tmp
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
# Use v1.4.6-p4 of Subread/Subjunc
SUBJUNC=/scratch0/langmead-fs1/shared/subread-1.4.6-p4-Linux-x86_64/bin/subjunc
# Used version 0.1.8 of Rail-RNA, but wrapped version 2.2.4 of Bowtie2 and version 1.1.1 of Bowtie
# Specify Python executable/loc of get_junctions.py; PyPy 2.5.0 was used
PYTHON=/home/anellor1/raildotbio/pypy-2.5-linux_x86_64-portable/bin/pypy
# Samtools v1.2 was used
SAMTOOLS=/home/anellor1/raildotbio/samtools-1.2/samtools

# Specify log filename for recording times
TIMELOG=${MAINOUTPUT}/small_data_times.log

## Specify locations of reference-related files
## See create_indexes.sh for index creation script
# Subjunc index
SUBJUNCIDX=/scratch0/langmead-fs1/indexes_for_paper/subreadgenome

# Generic name of file measuring performance of a given alignment strategy
# Performance is computed with spliced_read_recovery_performance.py; refer to that file for details
PERFORMANCE=perform

# Flux outputs paired-end reads in one file; split files here
echo 'Splitting Flux FASTQs...'
awk '(NR-1) % 8 < 4' $DATADIR/${SAMPLE}_sim.fastq >${SCRATCH}/${SAMPLE}_sim_left.fastq
awk '(NR-1) % 8 >= 4' $DATADIR/${SAMPLE}_sim.fastq >${SCRATCH}/${SAMPLE}_sim_right.fastq

RAILHOME=/scratch0/langmead-fs1/rail
cd $SCRATCH
mkdir -p ${SAMPLE}
cd ${SAMPLE}
mkdir -p subjunc

cd $MAINOUTPUT
mkdir -p ${SAMPLE}
SAMPLEOUTPUT=${MAINOUTPUT}/${SAMPLE}

# Run simulations
OUTPUT=$SCRATCH/${SAMPLE}
echo 'Running Subjunc on sample '${SAMPLE}' in paired-end mode...'
echo '#'${SAMPLE}' Subjunc' >>$TIMELOG
mkdir -p $OUTPUT/subjunc
cd $OUTPUT/subjunc
# Use Subjunc defaults
time (${SUBJUNC} -T $CORES -d 50 -D 600 -i ${SUBJUNCIDX} -r ${SCRATCH}/${SAMPLE}_sim_left.fastq -R ${SCRATCH}/${SAMPLE}_sim_right.fastq -o Aligned.out.sam) 2>>$TIMELOG
echo 'Computing precision and recall...'
(cat Aligned.out.sam | $PYTHON $RAILHOME/eval/spliced_read_recovery_performance.py -g -t $DATADIR/${SAMPLE}_sim.bed >$PERFORMANCE 2>${PERFORMANCE}_summary) &
(cat Aligned.out.sam | $PYTHON $RAILHOME/eval/intron_recovery_performance.py -t $DATADIR/${SAMPLE}_sim.bed >${PERFORMANCE}_intron_recovery_summary) &
(cat Aligned.out.sam | $PYTHON $RAILHOME/eval/mapping_accuracy.py -t $DATADIR/${SAMPLE}_sim.bed -g >${PERFORMANCE}_mapping_accuracy_summary) &
(cat Aligned.out.sam | $PYTHON $RAILHOME/eval/mapping_accuracy.py -t $DATADIR/${SAMPLE}_sim.bed -c 0.1 -g >${PERFORMANCE}_mapping_accuracy_SC_summary) &
wait
# Move Subjunc results to final destination
rm -rf ${SAMPLEOUTPUT}/subjunc
cp -r ${OUTPUT}/subjunc $SAMPLEOUTPUT
rm -rf ${OUTPUT}/subjunc
