#!/usr/bin/env bash
# This script computes performances for the two samples chosen in run_all_small_data_sims_locally.sh
# after they were analyzed together with the other 18 bioreplicates on AWS (see submit_big_data_sims.sh
# and the manifest file geuvadis_sim.manifest).

SAMPLE1=YRI-7-0 # this is sample NA18861.1.M_120209_2 
SAMPLE2=YRI-0-0 # this is sample NA18508.1.M_111124_1

# Specify Python executable/loc of get_junctions.py; PyPy 2.2.1 was used
PYTHON=pypy
# Samtools v0.1.19-44428cd was used
SAMTOOLS=samtools

# Specify FULL PATH to output directory
MAINOUTPUT=$1

# Generic name of file measuring performance of a given alignment strategy
# Performance is computed with spliced_read_recovery_performance.py; refer to that file for details
PERFORMANCE=perform

# Specify data directory; fastqs should be of the form [SAMPLE NAME]_sim.fastq; Flux beds should be
# of form [SAMPLE_NAME]_sim.bed
DATADIR=/scratch0/langmead-fs1/geuvadis_sim
# Where the S3 BAMs were downloaded
BAMDIR=/scratch0/langmead-fs1/geuvadis_sim/rail_combined_biorep_results

RAILHOME=/scratch0/langmead-fs1/rail

echo 'Computing precision and recall for Rail-RNA on sample NA18861.1.M_120209_2'
(for i in $BAMDIR/*YRI-7-0*.bam; do $SAMTOOLS view $i; done) | $PYTHON $RAILHOME/eval/spliced_read_recovery_performance.py \
	-t $DATADIR/NA18861.1.M_120209_2_sim.bed >$MAINOUTPUT/$PERFORMANCE_NA18861.1.M_120209_2 2>$MAINOUTPUT/${PERFORMANCE}_NA18861.1.M_120209_2_summary &
echo 'Computing precision and recall for Rail-RNA on sample NA18508.1.M_111124_1'
(for i in $BAMDIR/*YRI-0-0*.bam; do $SAMTOOLS view $i; done) | $PYTHON $RAILHOME/eval/spliced_read_recovery_performance.py \
	-t $DATADIR/NA18508.1.M_111124_1_sim.bed >$MAINOUTPUT/$PERFORMANCE_NA18508.1.M_111124_1 2>$MAINOUTPUT/${PERFORMANCE}_NA18508.1.M_111124_1_summary &