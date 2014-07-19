#!/usr/bin/env bash
# THIS SCRIPT IS SUBSUMED BY run_all_small_data_sims_locally.sh
# It is included because Rail had to be rerun to time it
# RUN WITH TASKSET to limit CPU usage on multicore system for benchmarking!
# $1: number of cores
# $2: output directory -- SPECIFY FULL PATH
# Ex: taskset -c 0,1,2,3 sh run_all_small_data_sims_locally.sh 4 ./myoutput
# Select two sample names for analysis. See generate_bioreps.py for how sample data was generated.
SAMPLE1=NA18861.1.M_120209_2
SAMPLE2=NA18508.1.M_111124_1

# Specify data directory; fastqs should be of the form [SAMPLE NAME]_sim.fastq; Flux beds should be
# of form [SAMPLE_NAME]_sim.bed
DATADIR=/scratch0/langmead-fs1/geuvadis_sim

## Specify locations of executables
# Used version 0.1.0 of Rail-RNA; wrapped version 2.2.2 of Bowtie2
# Specify Python executable/loc of get_junctions.py; PyPy 2.2.1 was used
PYTHON=pypy
RAILHOME=/scratch0/langmead-fs1/rail
RAILRNA=$PYTHON\ $RAILHOME/src
# Samtools v0.1.19-44428cd was used (and wrapped by TopHat and Rail)
SAMTOOLS=samtools

# Specify number of parallel processes for each program
CORES=$1

# Specify FULL PATH to output directory
MAINOUTPUT=$2
mkdir -p $MAINOUTPUT

# Specify log filename for recording times
TIMELOG=$MAINOUTPUT/rail_times.log

## Specify locations of reference-related files
## See create_indexes.sh for index creation script
# Bowtie indexes
BOWTIE1IDX=/scratch0/langmead-fs1/indexes_for_paper/genome
BOWTIE2IDX=/scratch0/langmead-fs1/indexes_for_paper/genome

# Generic name of file measuring performance of a given alignment strategy
# Performance is computed with spliced_read_recovery_performance.py; refer to that file for details
PERFORMANCE=perform

cd $MAINOUTPUT
for SAMPLE in {$SAMPLE1,$SAMPLE2}
do
	mkdir -p ${SAMPLE}
	cd ${SAMPLE}
	mkdir -p rail
	cd ..
done

# Run simulations
for SAMPLE in {$SAMPLE1,$SAMPLE2}
do
	OUTPUT=$MAINOUTPUT/${SAMPLE}
	echo 'Running Rail-RNA on sample '${SAMPLE}'...'
	echo '#'${SAMPLE}' Rail-RNA' >>$TIMELOG
	# Write manifest file
	echo -e $DATADIR/${SAMPLE}_sim.fastq'\t0\t'${SAMPLE}'-0-0' >$MAINOUTPUT/${SAMPLE}.manifest
	time ($RAILRNA go local -p $CORES -m $MAINOUTPUT/${SAMPLE}.manifest -o $OUTPUT/rail2 --log $OUTPUT/rail2.log -1 $BOWTIE1IDX -2 $BOWTIE2IDX -f >/dev/null 2>&1) 2>>$TIMELOG
	echo 'Computing precision and recall...'
	$SAMTOOLS merge - $OUTPUT/rail/alignments/*.bam | $SAMTOOLS view - | $PYTHON $RAILHOME/eval/spliced_read_recovery_performance.py \
		-t $DATADIR/${SAMPLE}_sim.bed >$OUTPUT/rail2/$PERFORMANCE 2>$OUTPUT/rail2/${PERFORMANCE}_summary
done