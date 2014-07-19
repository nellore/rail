#!/usr/bin/env bash
# This script is SUBSUMED BY run_all_small_data_sims_locally.sh; it was written to rerun spliced_read_recovery_performance.py
# after enhancements were made
# $1: output directory -- SPECIFY FULL PATH
# Select two sample names for analysis. Should match those in script run_all_small_data_sims_locally
SAMPLE1=NA18861.1.M_120209_2
SAMPLE2=NA18508.1.M_111124_1

# Specify Python executable/loc of get_junctions.py; PyPy 2.2.1 was used
PYTHON=pypy
# Samtools v0.1.19-44428cd was used
SAMTOOLS=samtools

# Specify FULL PATH to output directory
MAINOUTPUT=$1

# Generic name of file measuring performance of a given alignment strategy
# Performance is computed with spliced_read_recovery_performance.py; refer to that file for details
PERFORMANCE=perform

# Measure performances
for SAMPLE in {$SAMPLE1,$SAMPLE2}
do
	OUTPUT=$MAINOUTPUT/${SAMPLE}
	echo 'Computing precision and recall for TopHat on sample '${SAMPLE}' with no annotation and in single-end mode...'
	$SAMTOOLS view $OUTPUT/tophat/noann_single/accepted_hits.bam | $PYTHON $RAILHOME/eval/spliced_read_recovery_performance.py \
		-t $DATADIR/${SAMPLE}_sim.bed -g >$OUTPUT/tophat/noann_single/$PERFORMANCE 2>$OUTPUT/tophat/noann_single/${PERFORMANCE}_summary
	echo 'Computing precision and recall for TopHat on sample '${SAMPLE}' with no annotation and in paired-end mode...'
	$SAMTOOLS view $OUTPUT/tophat/noann_paired/accepted_hits.bam | $PYTHON $RAILHOME/eval/spliced_read_recovery_performance.py \
		-t $DATADIR/${SAMPLE}_sim.bed -g >$OUTPUT/tophat/noann_paired/$PERFORMANCE 2>$OUTPUT/tophat/noann_paired/${PERFORMANCE}_summary
	echo 'Computing precision and recall for TopHat on sample '${SAMPLE}' with annotation and in single-end mode...'
	$SAMTOOLS view $OUTPUT/tophat/ann_single/accepted_hits.bam | $PYTHON $RAILHOME/eval/spliced_read_recovery_performance.py \
		-t $DATADIR/${SAMPLE}_sim.bed -g >$OUTPUT/tophat/ann_single/$PERFORMANCE 2>$OUTPUT/tophat/ann_single/${PERFORMANCE}_summary
	echo 'Computing precision and recall for TopHat on sample '${SAMPLE}' with annotation and in paired-end mode...'
	$SAMTOOLS view $OUTPUT/tophat/ann_paired/accepted_hits.bam | $PYTHON $RAILHOME/eval/spliced_read_recovery_performance.py \
		-t $DATADIR/${SAMPLE}_sim.bed -g >$OUTPUT/tophat/ann_paired/$PERFORMANCE 2>$OUTPUT/tophat/ann_paired/${PERFORMANCE}_summary
	cd $OUTPUT/star/noann_single_1pass
	echo 'Computing precision and recall for STAR on sample '${SAMPLE}' with no annotation and in single-end mode......'
	cat Aligned.out.sam | $PYTHON $RAILHOME/eval/spliced_read_recovery_performance.py -g -t $DATADIR/${SAMPLE}_sim.bed >$PERFORMANCE 2>${PERFORMANCE}_summary
	cd $OUTPUT/star/noann_paired_1pass
	echo 'Computing precision and recall for STAR on sample '${SAMPLE}' with no annotation and in paired-end mode......'
	cat Aligned.out.sam | $PYTHON $RAILHOME/eval/spliced_read_recovery_performance.py -g -t $DATADIR/${SAMPLE}_sim.bed >$PERFORMANCE 2>${PERFORMANCE}_summary
	cd $OUTPUT/star/ann_single_1pass
	echo 'Computing precision and recall for STAR on sample '${SAMPLE}' with annotation and in single-end mode...'
	cat Aligned.out.sam | $PYTHON $RAILHOME/eval/spliced_read_recovery_performance.py -g -t $DATADIR/${SAMPLE}_sim.bed >$PERFORMANCE 2>${PERFORMANCE}_summary
	cd $OUTPUT/star/ann_paired_1pass
	echo 'Computing precision and recall for STAR on sample '${SAMPLE}' with annotation and in paired-end mode...'
	cat Aligned.out.sam | $PYTHON $RAILHOME/eval/spliced_read_recovery_performance.py -g -t $DATADIR/${SAMPLE}_sim.bed >$PERFORMANCE 2>${PERFORMANCE}_summary
	cd $OUTPUT/star/noann_single_2pass
	echo 'Computing precision and recall for second pass of STAR on sample '${SAMPLE}' with no annotation and in single-end mode...'
	cat Aligned.out.sam | $PYTHON $RAILHOME/eval/spliced_read_recovery_performance.py -g -t $DATADIR/${SAMPLE}_sim.bed >$PERFORMANCE 2>${PERFORMANCE}_summary
	cd $OUTPUT/star/noann_paired_2pass
	echo 'Computing precision and recall for second pass of STAR on sample '${SAMPLE}' with no annotation and in paired-end mode...'
	cat Aligned.out.sam | $PYTHON $RAILHOME/eval/spliced_read_recovery_performance.py -g -t $DATADIR/${SAMPLE}_sim.bed >$PERFORMANCE 2>${PERFORMANCE}_summary
	cd $OUTPUT/star/ann_single_2pass
	echo 'Computing precision and recall for second pass of STAR on sample '${SAMPLE}' with annotation and in single-end mode...'
	cat Aligned.out.sam | $PYTHON $RAILHOME/eval/spliced_read_recovery_performance.py -g -t $DATADIR/${SAMPLE}_sim.bed >$PERFORMANCE 2>${PERFORMANCE}_summary
	cd $OUTPUT/star/ann_paired_2pass
	echo 'Computing precision and recall for second pass of STAR on sample '${SAMPLE}' with annotation and in paired-end mode...'
	cat Aligned.out.sam | $PYTHON $RAILHOME/eval/spliced_read_recovery_performance.py -g -t $DATADIR/${SAMPLE}_sim.bed >$PERFORMANCE 2>${PERFORMANCE}_summary
	echo 'Computing precision and recall for Rail-RNA on sample '${SAMPLE}'...'
	$SAMTOOLS merge - $OUTPUT/rail/alignments/*.bam | $SAMTOOLS view - | $PYTHON $RAILHOME/eval/spliced_read_recovery_performance.py \
		-t $DATADIR/${SAMPLE}_sim.bed >$OUTPUT/rail/$PERFORMANCE 2>$OUTPUT/rail/${PERFORMANCE}_summary
done