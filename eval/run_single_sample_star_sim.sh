#!/usr/bin/env bash
# $1: number of cores
# $2: output directory -- SPECIFY FULL PATH
# $3: where to find sample fastqs from generate_bioreps.py
# $4: sample name; this is the prefix of "_sim.fastq"
# $5: scratch directory; files are written here first, and relevant output is copied back to $2
# Ex: taskset -c 0,1,2,3 sh run_single_sample_star_sim.sh 4 ./myoutput NA11829_male_CEU_UU_6-1-1 /tmp
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
RAILHOME=/scratch0/langmead-fs1/rail

## Specify locations of executables
# Used version 2.4.2a of STAR
STAR=/scratch0/langmead-fs1/shared/STARv2.4.2a/STAR
# Specify Python executable/loc of get_junctions.py; PyPy 2.5.0 was used
PYTHON=/home/anellor1/raildotbio3/pypy-2.5-linux_x86_64-portable/bin/pypy
# Samtools v1.2 was used
SAMTOOLS=/home/anellor1/raildotbio/samtools-1.2/samtools

# Specify log filename for recording times
TIMELOG=${MAINOUTPUT}/star_times.log

## Specify locations of reference-related files
## See create_indexes.sh for index creation script
# STAR index
STARIDX=/scratch0/langmead-fs1/indexes_for_paper/star
# Overhang length for STAR; this should ideally be max read length - 1
OVERHANG=75

# Generic name of file measuring performance of a given alignment strategy
# Performance is computed with spliced_read_recovery_performance.py; refer to that file for details
PERFORMANCE=perform

## Specify location of annotation
# This is Gencode v12, which may be obtained at ftp://ftp.sanger.ac.uk/pub/gencode/release_12/gencode.v12.annotation.gtf.gz
ANNOTATION=/scratch0/langmead-fs1/geuvadis_sim/gencode.v12.annotation.gtf

## STAR requires its own annotation format that lists introns directly; conversion is done by get_junctions.py
# Build STAR junction index from GTF
STARANNOTATION=$SCRATCH/junctions_for_star.txt

echo 'Building annotation for STAR from GTF...'
cat $ANNOTATION | $PYTHON $RAILHOME/eval/get_junctions.py >$STARANNOTATION
# Flux outputs paired-end reads in one file; split files here
echo 'Splitting Flux FASTQs...'
awk '(NR-1) % 8 < 4' $DATADIR/${SAMPLE}_sim.fastq >${SCRATCH}/${SAMPLE}_sim_left.fastq
awk '(NR-1) % 8 >= 4' $DATADIR/${SAMPLE}_sim.fastq >${SCRATCH}/${SAMPLE}_sim_right.fastq

## Where Feb 2009 hg19 chromosome files are located; download them at http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz
# These are here because it's necessary to build a new STAR index for its 2-pass and annotation protocols
FADIR=/scratch0/langmead-fs1/shared/references/hg19/fasta

cd $SCRATCH
mkdir -p ${SAMPLE}
cd ${SAMPLE}
mkdir -p star
cd ..

# Create STAR index including annotation
echo 'Creating new STAR index including splice junctions...'
echo '#STAR index pre-1-pass ann' >>$TIMELOG
STARANNIDX=${SCRATCH}/${SAMPLE}/starannidx
mkdir -p $STARANNIDX
# Use --sjdbOverhang 75 because reads are 76 bases long! See p. 3 of STAR manual for details.
time ($STAR --runMode genomeGenerate --genomeDir $STARANNIDX --genomeFastaFiles $FADIR/chr{1..22}.fa $FADIR/chr{X,Y,M}.fa \
		--runThreadN $CORES --sjdbFileChrStartEnd $STARANNOTATION --sjdbOverhang $OVERHANG 2>&1) 2>>$TIMELOG

cd $MAINOUTPUT
mkdir -p ${SAMPLE}
SAMPLEOUTPUT=${MAINOUTPUT}/${SAMPLE}

# Run simulations
OUTPUT=$SCRATCH/${SAMPLE}
# STAR protocol for 2-pass execution w/ index construction is described on pp. 43-44 of the supplement of the RGASP spliced alignment paper
# (http://www.nature.com/nmeth/journal/v10/n12/extref/nmeth.2722-S1.pdf)
echo 'Running STAR on sample '${SAMPLE}' with no regenerated genome/no annotation and in paired-end mode...'
echo '#'${SAMPLE}' STAR 2-pass nogen noann paired' >>$TIMELOG
mkdir -p $OUTPUT/star/nogen_noann_paired_2pass
cd $OUTPUT/star/nogen_noann_paired_2pass
time ($STAR --genomeDir $STARIDX --readFilesIn ${SCRATCH}/${SAMPLE}_sim_left.fastq ${SCRATCH}/${SAMPLE}_sim_right.fastq --runThreadN $CORES --twopass1readsN -1 --sjdbOverhang $OVERHANG --twopassMode Basic >&1) 2>>$TIMELOG
echo 'Computing precision and recall...'
(cat Aligned.out.sam | $PYTHON $RAILHOME/eval/spliced_read_recovery_performance.py -g -t $DATADIR/${SAMPLE}_sim.bed >$PERFORMANCE 2>${PERFORMANCE}_summary) &
(cat Aligned.out.sam | $PYTHON $RAILHOME/eval/intron_recovery_performance.py -t $DATADIR/${SAMPLE}_sim.bed >${PERFORMANCE}_intron_recovery_summary) &
(cat Aligned.out.sam | $PYTHON $RAILHOME/eval/mapping_accuracy.py -t $DATADIR/${SAMPLE}_sim.bed -g >${PERFORMANCE}_mapping_accuracy_summary) &
(cat Aligned.out.sam | $PYTHON $RAILHOME/eval/mapping_accuracy.py -t $DATADIR/${SAMPLE}_sim.bed -c 0.1 -g >${PERFORMANCE}_mapping_accuracy_SC_summary) &
echo 'Running STAR on sample '${SAMPLE}' with no annotation and in paired-end mode...'
echo '#'${SAMPLE}' STAR 1-pass noann paired' >>$TIMELOG
mkdir -p $OUTPUT/star/noann_paired_1pass
cd $OUTPUT/star/noann_paired_1pass
time ($STAR --genomeDir $STARIDX --readFilesIn ${SCRATCH}/${SAMPLE}_sim_left.fastq ${SCRATCH}/${SAMPLE}_sim_right.fastq --runThreadN $CORES 2>&1) 2>>$TIMELOG
echo 'Computing precision and recall...'
(cat Aligned.out.sam | $PYTHON $RAILHOME/eval/spliced_read_recovery_performance.py -g -t $DATADIR/${SAMPLE}_sim.bed >$PERFORMANCE 2>${PERFORMANCE}_summary) &
(cat Aligned.out.sam | $PYTHON $RAILHOME/eval/intron_recovery_performance.py -t $DATADIR/${SAMPLE}_sim.bed >${PERFORMANCE}_intron_recovery_summary) &
(cat Aligned.out.sam | $PYTHON $RAILHOME/eval/mapping_accuracy.py -t $DATADIR/${SAMPLE}_sim.bed -g >${PERFORMANCE}_mapping_accuracy_summary) &
(cat Aligned.out.sam | $PYTHON $RAILHOME/eval/mapping_accuracy.py -t $DATADIR/${SAMPLE}_sim.bed -c 0.1 -g >${PERFORMANCE}_mapping_accuracy_SC_summary) &
echo 'Running STAR on sample '${SAMPLE}' with annotation and in paired-end mode...'
echo '#'${SAMPLE}' STAR 1-pass ann paired' >>$TIMELOG
mkdir -p $OUTPUT/star/ann_paired_1pass
cd $OUTPUT/star/ann_paired_1pass
time ($STAR --genomeDir $STARANNIDX --readFilesIn ${SCRATCH}/${SAMPLE}_sim_left.fastq ${SCRATCH}/${SAMPLE}_sim_right.fastq --runThreadN $CORES 2>&1) 2>>$TIMELOG
echo 'Computing precision and recall...'
(cat Aligned.out.sam | $PYTHON $RAILHOME/eval/spliced_read_recovery_performance.py -g -t $DATADIR/${SAMPLE}_sim.bed >$PERFORMANCE 2>${PERFORMANCE}_summary) &
(cat Aligned.out.sam | $PYTHON $RAILHOME/eval/intron_recovery_performance.py -t $DATADIR/${SAMPLE}_sim.bed >${PERFORMANCE}_intron_recovery_summary) &
(cat Aligned.out.sam | $PYTHON $RAILHOME/eval/mapping_accuracy.py -t $DATADIR/${SAMPLE}_sim.bed -g >${PERFORMANCE}_mapping_accuracy_summary) &
(cat Aligned.out.sam | $PYTHON $RAILHOME/eval/mapping_accuracy.py -t $DATADIR/${SAMPLE}_sim.bed -c 0.1 -g >${PERFORMANCE}_mapping_accuracy_SC_summary) &
echo 'Creating new STAR index for sample '${SAMPLE}' with annotation and in paired-end mode...'
echo '#'${SAMPLE}' STAR 1-pass ann paired index' >>$TIMELOG
STARIDXANNPAIRED=$OUTPUT/star/ann_paired_idx
mkdir -p $STARIDXANNPAIRED
time ($STAR --runMode genomeGenerate --genomeDir $STARIDXANNPAIRED --genomeFastaFiles $FADIR/chr{1..22}.fa $FADIR/chr{X,Y,M}.fa \
		--sjdbFileChrStartEnd $OUTPUT/star/ann_paired_1pass/SJ.out.tab --sjdbOverhang $OVERHANG --runThreadN $CORES 2>&1) 2>>$TIMELOG
echo 'Running second pass of STAR on sample '${SAMPLE}' with annotation and in paired-end mode...'
echo '#'${SAMPLE}' STAR 2-pass ann paired' >>$TIMELOG
mkdir -p $OUTPUT/star/ann_paired_2pass
cd $OUTPUT/star/ann_paired_2pass
time ($STAR --genomeDir $STARIDXANNPAIRED --readFilesIn ${SCRATCH}/${SAMPLE}_sim_left.fastq ${SCRATCH}/${SAMPLE}_sim_right.fastq --runThreadN $CORES 2>&1) 2>>$TIMELOG
echo 'Computing precision and recall...'
(cat Aligned.out.sam | $PYTHON $RAILHOME/eval/spliced_read_recovery_performance.py -g -t $DATADIR/${SAMPLE}_sim.bed >$PERFORMANCE 2>${PERFORMANCE}_summary) &
(cat Aligned.out.sam | $PYTHON $RAILHOME/eval/intron_recovery_performance.py -t $DATADIR/${SAMPLE}_sim.bed >${PERFORMANCE}_intron_recovery_summary) &
(cat Aligned.out.sam | $PYTHON $RAILHOME/eval/mapping_accuracy.py -t $DATADIR/${SAMPLE}_sim.bed -g >${PERFORMANCE}_mapping_accuracy_summary) &
(cat Aligned.out.sam | $PYTHON $RAILHOME/eval/mapping_accuracy.py -t $DATADIR/${SAMPLE}_sim.bed -c 0.1 -g >${PERFORMANCE}_mapping_accuracy_SC_summary) &
wait
# Move STAR results to final destination
mv $OUTPUT/star $SAMPLEOUTPUT
