#!/usr/bin/env bash

# Select two sample names for analysis. See generate_bioreps.py for how sample data was generated.
SAMPLE1=NA18861.1.M_120209_2
SAMPLE2=NA18508.1.M_111124_1

# Specify data directory; fastqs should be of the form [SAMPLE NAME]_sim.fastq
DATADIR=/scratch0/langmead-fs1/geuvadis_sim

## Specify locations of executables
# Use version 2.0.12 of TopHat; wrapped version 2.2.2 of Bowtie2
TOPHAT=/scratch0/langmead-fs1/shared/tophat-2.0.12.Linux_x86_64/tophat2
# Use version 2.3.0e of STAR
STAR=/scratch0/langmead-fs1/shared/STAR_2.3.0e.Linux_x86_64/STAR
# Use version 0.1.0 of Rail-RNA; wrapped version 2.2.2 of Bowtie2
RAILRNA=python\ /scratch0/langmead-fs1/rail/src

## Specify location of annotation
# This is Gencode v12, which may be obtained at ftp://ftp.sanger.ac.uk/pub/gencode/release_12/gencode.v12.annotation.gtf.gz
ANNOTATION=/scratch0/langmead-fs1/geuvadis_sim/gencode.v12.annotation.gtf

# Specify number of parallel processes for each program
CORES=32

# Specify output directory
OUTPUT=/scratch0/langmead-fs1/geuvadis_sim/local_out

# Specify log file for recording times
TIMELOG=/scratch0/langmead-fs1/geuvadis_sim/small_data_times.log

## Specify locations of reference-related files
## See create_indexes.sh for index creation script
# Bowtie indexes
BOWTIE1_IDX=/scratch0/langmead-fs1/indexes_for_paper/genome
BOWTIE2_IDX=/scratch0/langmead-fs1/indexes_for_paper/genome
# STAR index
STAR_REF=/scratch0/langmead-fs1/indexes_for_paper/star

# Flux outputs paired-end reads in one file; split files here
echo 'Splitting Flux FASTQs...'
for SAMPLE in {$SAMPLE1,$SAMPLE2}
do
	awk 'NR % 8 < 4' $DATADIR/$SAMPLE_sim.fastq > $DATADIR/$SAMPLE_sim_left.fastq
	awk 'NR % 8 >= 4' $DATADIR/$SAMPLE_sim.fastq > $DATADIR/$SAMPLE_sim_right.fastq
done

# Create output directories
mkdir -p $OUTPUT
mkdir -p $OUTPUT/tophat
mkdir -p $OUTPUT/star
mkdir -p $OUTPUT/rail

# Run simulations
for SAMPLE in {$SAMPLE1,$SAMPLE2}
do
	echo 'Running TopHat on sample $SAMPLE with no annotation and in single-end mode...'
	echo '#$SAMPLE TopHat noann single' >>$TIMELOG
	time ($TOPHAT -o $OUTPUT/tophat/noann_single -p $CORES $BOWTIE2_IDX $SAMPLE_sim.fastq) 2>>$TIMELOG
	echo 'Running TopHat on sample $SAMPLE with no annotation and in paired-end mode...'
	echo '#$SAMPLE TopHat noann paired' >>$TIMELOG
	time ($TOPHAT -o $OUTPUT/tophat/noann_paired -p $CORES $BOWTIE2_IDX $SAMPLE_sim_left.fastq $SAMPLE_sim_right.fastq) 2>>$TIMELOG
	echo 'Running TopHat on sample $SAMPLE with annotation and in single-end mode...'
	echo '#$SAMPLE TopHat ann single' >>$TIMELOG
	time ($TOPHAT -o $OUTPUT/tophat/ann_single -G $ANNOTATION -p $CORES $BOWTIE2_IDX $SAMPLE_sim.fastq) 2>>$TIMELOG
	echo 'Running TopHat on sample $SAMPLE with annotation and in paired-end mode...'
	echo '#$SAMPLE TopHat ann paired' >>$TIMELOG
	time ($TOPHAT -o $OUTPUT/tophat/noann_paired -G $ANNOTATION -p $CORES $BOWTIE2_IDX $SAMPLE_sim_left.fastq $SAMPLE_sim_right.fastq) 2>>$TIMELOG
done