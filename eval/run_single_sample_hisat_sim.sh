#!/usr/bin/env bash
# $1: number of cores
# $2: output directory -- SPECIFY FULL PATH
# $3: where to find sample fastqs from generate_bioreps.py
# $4: sample name; this is the prefix of "_sim.fastq"
# $5: scratch directory; files are written here first, and relevant output is copied back to $2
# Ex: taskset -c 0,1,2,3 sh run_single_sample_hisat_sim.sh 4 ./myoutput NA11829_male_CEU_UU_6-1-1 /tmp
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
# Used version 0.1.6-beta of HISAT
HISAT=/scratch0/langmead-fs1/shared/hisat-0.1.6-beta/hisat
# Use HISAT's tool for extracting splice sites for its junction database
HISATSPLICE=/scratch0/langmead-fs1/shared/hisat-0.1.6-beta/extract_splice_sites.py
# Specify Python executable/loc of get_junctions.py; PyPy 2.5.0 was used
PYTHON=/home/anellor1/raildotbio/pypy-2.5-linux_x86_64-portable/bin/pypy
# Samtools v1.2 was used
SAMTOOLS=/home/anellor1/raildotbio/samtools-1.2/samtools

# Specify log filename for recording times
TIMELOG=${MAINOUTPUT}/hisat_times.log

## Specify locations of reference-related files
## See create_indexes.sh for index creation script
# HISAT index
HISATIDX=/scratch0/langmead-fs1/indexes_for_paper/hisatgenome

# Generic name of file measuring performance of a given alignment strategy
# Performance is computed with spliced_read_recovery_performance.py; refer to that file for details
PERFORMANCE=perform

## Specify location of annotation
# This is Gencode v12, which may be obtained at ftp://ftp.sanger.ac.uk/pub/gencode/release_12/gencode.v12.annotation.gtf.gz
ANNOTATION=/scratch0/langmead-fs1/geuvadis_sim/gencode.v12.annotation.gtf

## HISAT also uses a list of introns, but here we use a tool that came with HISAT to grab splice sites because coordinate system could be different
# Build HISAT junction index from GTF
HISATANNOTATION=$SCRATCH/junctions_for_hisat.txt

# Flux outputs paired-end reads in one file; split files here
echo 'Building annotation for HISAT from GTF...'
$PYTHON $HISATSPLICE $ANNOTATION >$HISATANNOTATION
echo 'Splitting Flux FASTQs...'
awk '(NR-1) % 8 < 4' $DATADIR/${SAMPLE}_sim.fastq >${SCRATCH}/${SAMPLE}_sim_left.fastq
awk '(NR-1) % 8 >= 4' $DATADIR/${SAMPLE}_sim.fastq >${SCRATCH}/${SAMPLE}_sim_right.fastq

cd $SCRATCH
mkdir -p ${SAMPLE}
cd ${SAMPLE}
mkdir -p hisat

cd $MAINOUTPUT
mkdir -p ${SAMPLE}
SAMPLEOUTPUT=${MAINOUTPUT}/${SAMPLE}

# Run simulations
OUTPUT=$SCRATCH/${SAMPLE}
echo 'Running HISAT on sample '${SAMPLE}' with no annotation and in paired-end mode...'
echo '#'${SAMPLE}' HISAT 1-pass noann paired' >>$TIMELOG
mkdir -p $OUTPUT/hisat/noann_paired_1pass
cd $OUTPUT/hisat/noann_paired_1pass
time ($HISAT -x $HISATIDX -1 ${SCRATCH}/${SAMPLE}_sim_left.fastq -2 ${SCRATCH}/${SAMPLE}_sim_right.fastq -p $CORES -S Aligned.out.sam --novel-splicesite-outfile novel_splice_sites.txt 2>&1) 2>>$TIMELOG
echo 'Computing precision and recall...'
(cat Aligned.out.sam | $PYTHON $RAILHOME/eval/spliced_read_recovery_performance.py -g -t $DATADIR/${SAMPLE}_sim.bed >$PERFORMANCE 2>${PERFORMANCE}_summary) &
(cat Aligned.out.sam | $PYTHON $RAILHOME/eval/intron_recovery_performance.py -t $DATADIR/${SAMPLE}_sim.bed >${PERFORMANCE}_intron_recovery_summary) &
(cat Aligned.out.sam | $PYTHON $RAILHOME/eval/mapping_accuracy.py -t $DATADIR/${SAMPLE}_sim.bed >${PERFORMANCE}_mapping_accuracy_summary) &
(cat Aligned.out.sam | $PYTHON $RAILHOME/eval/mapping_accuracy.py -t $DATADIR/${SAMPLE}_sim.bed -c 0.1 >${PERFORMANCE}_mapping_accuracy_SC_summary) &
wait
echo 'Running HISAT on sample '${SAMPLE}' with annotation and in paired-end mode...'
echo '#'${SAMPLE}' HISAT 1-pass ann paired' >>$TIMELOG
mkdir -p $OUTPUT/hisat/ann_paired_1pass
cd $OUTPUT/hisat/ann_paired_1pass
time ($HISAT -x $HISATIDX -1 ${SCRATCH}/${SAMPLE}_sim_left.fastq -2 ${SCRATCH}/${SAMPLE}_sim_right.fastq -p $CORES -S Aligned.out.sam --novel-splicesite-outfile novel_splice_sites.txt --novel-splicesite-infile $HISATANNOTATION 2>&1) 2>>$TIMELOG
echo 'Computing precision and recall...'
(cat Aligned.out.sam | $PYTHON $RAILHOME/eval/spliced_read_recovery_performance.py -g -t ${SCRATCH}/${SAMPLE}_sim.bed >$PERFORMANCE 2>${PERFORMANCE}_summary) &
(cat Aligned.out.sam | $PYTHON $RAILHOME/eval/intron_recovery_performance.py -t ${SCRATCH}/${SAMPLE}_sim.bed >${PERFORMANCE}_intron_recovery_summary) &
(cat Aligned.out.sam | $PYTHON $RAILHOME/eval/mapping_accuracy.py -t $DATADIR/${SAMPLE}_sim.bed >${PERFORMANCE}_mapping_accuracy_summary) &
(cat Aligned.out.sam | $PYTHON $RAILHOME/eval/mapping_accuracy.py -t $DATADIR/${SAMPLE}_sim.bed -c 0.1 >${PERFORMANCE}_mapping_accuracy_SC_summary) &
wait
echo 'Running second pass of HISAT on sample '${SAMPLE}' with no annotation and in paired-end mode...'
echo '#'${SAMPLE}' HISAT 2-pass noann paired' >>$TIMELOG
mkdir -p $OUTPUT/hisat/noann_paired_2pass
cd $OUTPUT/hisat/noann_paired_2pass
time ($HISAT -x $HISATIDX -1 ${SCRATCH}/${SAMPLE}_sim_left.fastq -2 ${SCRATCH}/${SAMPLE}_sim_right.fastq -p $CORES -S Aligned.out.sam --novel-splicesite-infile $OUTPUT/hisat/noann_paired_1pass/novel_splice_sites.txt 2>&1) 2>>$TIMELOG
echo 'Computing precision and recall...'
(cat Aligned.out.sam | $PYTHON $RAILHOME/eval/spliced_read_recovery_performance.py -g -t $DATADIR/${SAMPLE}_sim.bed >$PERFORMANCE 2>${PERFORMANCE}_summary) &
(cat Aligned.out.sam | $PYTHON $RAILHOME/eval/intron_recovery_performance.py -t $DATADIR/${SAMPLE}_sim.bed >${PERFORMANCE}_intron_recovery_summary) &
(cat Aligned.out.sam | $PYTHON $RAILHOME/eval/mapping_accuracy.py -t $DATADIR/${SAMPLE}_sim.bed >${PERFORMANCE}_mapping_accuracy_summary) &
(cat Aligned.out.sam | $PYTHON $RAILHOME/eval/mapping_accuracy.py -t $DATADIR/${SAMPLE}_sim.bed -c 0.1 >${PERFORMANCE}_mapping_accuracy_SC_summary) &
echo 'Running second pass of HISAT on sample '${SAMPLE}' with annotation and in paired-end mode...'
echo '#'${SAMPLE}' HISAT 2-pass ann paired' >>$TIMELOG
mkdir -p $OUTPUT/hisat/ann_paired_2pass
cd $OUTPUT/hisat/ann_paired_2pass
time ($HISAT -x $HISATIDX -1 ${SCRATCH}${SAMPLE}_sim_left.fastq -2 ${SCRATCH}/${SAMPLE}_sim_right.fastq -p $CORES -S Aligned.out.sam --novel-splicesite-infile $OUTPUT/hisat/ann_paired_1pass/novel_splice_sites.txt 2>&1) 2>>$TIMELOG
echo 'Computing precision and recall...'
(cat Aligned.out.sam | $PYTHON $RAILHOME/eval/spliced_read_recovery_performance.py -g -t ${SCRATCH}/${SAMPLE}_sim.bed >$PERFORMANCE 2>${PERFORMANCE}_summary) &
(cat Aligned.out.sam | $PYTHON $RAILHOME/eval/intron_recovery_performance.py -t ${SCRATCH}/${SAMPLE}_sim.bed >${PERFORMANCE}_intron_recovery_summary) &
(cat Aligned.out.sam | $PYTHON $RAILHOME/eval/mapping_accuracy.py -t $DATADIR/${SAMPLE}_sim.bed >${PERFORMANCE}_mapping_accuracy_summary) &
(cat Aligned.out.sam | $PYTHON $RAILHOME/eval/mapping_accuracy.py -t $DATADIR/${SAMPLE}_sim.bed -c 0.1 >${PERFORMANCE}_mapping_accuracy_SC_summary) &
wait
# Move all hisat results to final destination
mv $OUTPUT/hisat $SAMPLEOUTPUT