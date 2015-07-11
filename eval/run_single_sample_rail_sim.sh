#!/usr/bin/env bash
# $1: number of cores
# $2: output directory -- SPECIFY FULL PATH
# $3: where to find sample fastqs from generate_bioreps.py
# $4: sample name; this is the prefix of "_sim.fastq"
# Ex: taskset -c 0,1,2,3 sh run_single_sample_rail_sim.sh 4 ./myoutput NA11829_male_CEU_UU_6-1-1 /tmp
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
# Used version 0.1.8 of Rail-RNA, but wrapped version 2.2.4 of Bowtie2 and version 1.1.1 of Bowtie
# Specify Python executable/loc of get_junctions.py; PyPy 2.5.0 was used
PYTHON=/home/anellor1/raildotbio/pypy-2.5-linux_x86_64-portable/bin/pypy
RAILHOME=/scratch0/langmead-fs1/rail
RAILRNA=rail-rna
# Samtools v1.2 was used
SAMTOOLS=/home/anellor1/raildotbio/samtools-1.2/samtools

# Specify log filename for recording times
TIMELOG=${MAINOUTPUT}/small_data_times.log

## Specify locations of reference-related files
## See create_indexes.sh for index creation script
# Bowtie indexes
BOWTIE1IDX=/scratch0/langmead-fs1/indexes_for_paper/genome
BOWTIE2IDX=/scratch0/langmead-fs1/indexes_for_paper/genome

# Generic name of file measuring performance of a given alignment strategy
# Performance is computed with spliced_read_recovery_performance.py; refer to that file for details
PERFORMANCE=perform

echo 'Splitting Flux FASTQs...'
awk '(NR-1) % 8 < 4' $DATADIR/${SAMPLE}_sim.fastq >${SCRATCH}/${SAMPLE}_sim_left.fastq
awk '(NR-1) % 8 >= 4' $DATADIR/${SAMPLE}_sim.fastq >${SCRATCH}/${SAMPLE}_sim_right.fastq

cd $SCRATCH
mkdir -p ${SAMPLE}
cd ${SAMPLE}
mkdir -p rail

cd $MAINOUTPUT
mkdir -p ${SAMPLE}
SAMPLEOUTPUT=${MAINOUTPUT}/${SAMPLE}

# Run simulations
OUTPUT=$SCRATCH/${SAMPLE}
echo 'Running Rail-RNA on sample '${SAMPLE}'...'
echo '#'${SAMPLE}' Rail-RNA' >>$TIMELOG
# Write manifest file
echo -e ${SCRATCH}/${SAMPLE}_sim_left.fastq'\t0\t'${SCRATCH}/${SAMPLE}_sim_right.fastq'\t0\t'${SAMPLE} >${SCRATCH}/${SAMPLE}.manifest
time ($RAILRNA go local -p $CORES -m ${SCRATCH}/${SAMPLE}.manifest -o $OUTPUT/rail --log $OUTPUT/rail.log -x $BOWTIE1IDX,$BOWTIE2IDX -f >/dev/null 2>&1) 2>>$TIMELOG
echo 'Computing precision and recall...'
(for i in $OUTPUT/rail/alignments/*.bam; do $SAMTOOLS view $i; done | $PYTHON $RAILHOME/eval/spliced_read_recovery_performance.py -t $DATADIR/${SAMPLE}_sim.bed >$OUTPUT/rail/$PERFORMANCE 2>$OUTPUT/rail/${PERFORMANCE}_summary) &
(for i in $OUTPUT/rail/alignments/*.bam; do $SAMTOOLS view $i; done | $PYTHON $RAILHOME/eval/intron_recovery_performance.py -t $DATADIR/${SAMPLE}_sim.bed >$OUTPUT/rail/${PERFORMANCE}_intron_recovery_summary) &
(for i in $OUTPUT/rail/alignments/*.bam; do $SAMTOOLS view $i; done | $PYTHON $RAILHOME/eval/mapping_accuracy.py -t $DATADIR/${SAMPLE}_sim.bed >${PERFORMANCE}_mapping_accuracy_summary) &
(for i in $OUTPUT/rail/alignments/*.bam; do $SAMTOOLS view $i; done | $PYTHON $RAILHOME/eval/mapping_accuracy.py -t $DATADIR/${SAMPLE}_sim.bed -c 0.1 >${PERFORMANCE}_mapping_accuracy_SC_summary) &
wait
# Move rail results to final destination
rm -rf ${SAMPLEOUTPUT}/rail
cp -r ${OUTPUT}/rail $SAMPLEOUTPUT
rm -rf ${OUTPUT}/rail