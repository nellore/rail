#!/usr/bin/env bash
## Grabs results of experiment running Rail-RNA on 112 simulated datasets
## reflecting the coverage distributions of 112 GEUVADIS datasets from S3
## and computes performance measures on only those two samples considered
## in run_all_small_data_sims_locally.sh. Also grabs transcript indexes
## and counts numbers of introns in each index with count_introns.py
## Requires s3cmd, Bowtie 2 in PATH
## $1: output directory of job WITH FILTER on S3
## $2: output directory of job WITHOUT FILTER on S3
## $3: directory in which to dump output
## $4: where to find Flux BED storing true read alignments from simulation
## Command we ran was sh grab_112_sim_results.sh s3://rail-results/geuv112sim.out s3://rail-results/geuv112sim.out.nofilter /scratch0/langmead-fs1/geuvadis_sims_for_paper/fromemr /scratch0/langmead-fs1/geuvadis_sims_for_paper

SAMPLE1=NA19129_female_YRI_UU_6-1-1
SAMPLE2=NA07048_male_CEU_UU_6-1-2

DATADIR=$4

RAILHOME=/scratch0/langmead-fs1/rail
PYTHON=pypy
PERFORMANCE=perform

CWD=$(pwd)
mkdir -p $3
cd $3
mkdir -p withfilter
cd withfilter
for SAMPLE in {$SAMPLE1, $SAMPLE2}
do 
	s3cmd get $1/alignments/alignments.$SAMPLE*.bam
	s3cmd get $1/transcript_index/*
	tar xvzf *.tar.gz
	$PYTHON $RAILHOME/eval/count_introns.py --basename intron >intron_count
	(for i in $SAMPLE*.bam; do samtools view $i; done | $PYTHON $RAILHOME/eval/spliced_read_recovery_performance.py -t $DATADIR/${SAMPLE}_sim.bed >$PERFORMANCE 2>${PERFORMANCE}_summary) &
	(for i in $SAMPLE*.bam; do samtools view $i; done | $PYTHON $RAILHOME/eval/intron_recovery_performance.py -t $DATADIR/${SAMPLE}_sim.bed >${PERFORMANCE}_intron_recovery_summary) &
done
wait
cd ..
mkdir -p withoutfilter
cd withoutfilter
for SAMPLE in {$SAMPLE1, $SAMPLE2}
do 
	s3cmd get $1/alignments/alignments.$SAMPLE*.bam
	s3cmd get $1/transcript_index/*
	tar xvzf *.tar.gz
	$PYTHON $RAILHOME/eval/count_introns.py --basename intron >intron_count
	(for i in $SAMPLE*.bam; do samtools view $i; done | $PYTHON $RAILHOME/eval/spliced_read_recovery_performance.py -t $DATADIR/${SAMPLE}_sim.bed >$PERFORMANCE 2>${PERFORMANCE}_summary) &
	(for i in $SAMPLE*.bam; do samtools view $i; done | $PYTHON $RAILHOME/eval/intron_recovery_performance.py -t $DATADIR/${SAMPLE}_sim.bed >${PERFORMANCE}_intron_recovery_summary) &
done
cd $CWD