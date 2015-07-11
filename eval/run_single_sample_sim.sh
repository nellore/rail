#!/usr/bin/env bash
# $1: number of cores
# $2: output directory -- SPECIFY FULL PATH
# $3: where to find sample fastqs from generate_bioreps.py
# $4: sample name; this is the prefix of "_sim.fastq"
# $5: scratch directory; files are written here first, and relevant output is copied back to $2
# Ex: taskset -c 0,1,2,3 sh run_all_small_data_sims.sh 4 ./myoutput NA11829_male_CEU_UU_6-1-1 /tmp
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
# Used version 2.4.2a of STAR
STAR=/scratch0/langmead-fs1/shared/STARv2.4.2a/STAR
# Used version 0.1.6-beta of HISAT
HISAT=/scratch0/langmead-fs1/shared/hisat-0.1.6-beta/hisat
# Use HISAT's tool for extracting splice sites for its junction database
HISATSPLICE=/scratch0/langmead-fs1/shared/hisat-0.1.6-beta/extract_splice_sites.py
# Use v1.4.6-p4 of Subread/Subjunc
SUBJUNC=/scratch0/langmead-fs1/shared/subread-1.4.6-p4-Linux-x86_64/bin/subjunc
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
# STAR index
STARIDX=/scratch0/langmead-fs1/indexes_for_paper/star
# HISAT index
HISATIDX=/scratch0/langmead-fs1/indexes_for_paper/hisatgenome
# Subjunc index
SUBJUNCIDX=/scratch0/langmead-fs1/indexes_for_paper/subreadgenome
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

## HISAT also uses a list of introns, but here we use a tool that came with HISAT to grab splice sites because coordinate system could be different
# Build HISAT junction index from GTF
HISATANNOTATION=$SCRATCH/junctions_for_hisat.txt

# Flux outputs paired-end reads in one file; split files here
echo 'Building annotation for HISAT from GTF...'
$PYTHON $HISATSPLICE $ANNOTATION >$HISATANNOTATION
echo 'Building annotation for STAR from GTF...'
cat $ANNOTATION | $PYTHON $RAILHOME/eval/get_junctions.py >$STARANNOTATION
echo 'Splitting Flux FASTQs...'
awk '(NR-1) % 8 < 4' $DATADIR/${SAMPLE}_sim.fastq >${SCRATCH}/${SAMPLE}_sim_left.fastq
awk '(NR-1) % 8 >= 4' $DATADIR/${SAMPLE}_sim.fastq >${SCRATCH}/${SAMPLE}_sim_right.fastq

## Where Feb 2009 hg19 chromosome files are located; download them at http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz
# These are here because it's necessary to build a new STAR index for its 2-pass and annotation protocols
FADIR=/scratch0/langmead-fs1/shared/references/hg19/fasta

# Create STAR index including annotation
echo 'Creating new STAR index including splice junctions...'
echo '#STAR index pre-1-pass ann' >>$TIMELOG
STARANNIDX=${SCRATCH}/starannidx
mkdir -p $STARANNIDX
# Use --sjdbOverhang 75 because reads are 76 bases long! See p. 3 of STAR manual for details.
time ($STAR --runMode genomeGenerate --genomeDir $STARANNIDX --genomeFastaFiles $FADIR/chr{1..22}.fa $FADIR/chr{X,Y,M}.fa \
		--runThreadN $CORES --sjdbFileChrStartEnd $STARANNOTATION --sjdbOverhang $OVERHANG 2>&1) 2>>$TIMELOG

cd $SCRATCH
mkdir -p ${SAMPLE}
cd ${SAMPLE}
mkdir -p hisat
mkdir -p tophat
mkdir -p star
mkdir -p rail
mkdir -p subjunc
cd ..

cd $MAINOUTPUT
mkdir -p ${SAMPLE}
$SAMPLEOUTPUT=${MAINOUTPUT}/${SAMPLE}

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
(cat Aligned.out.sam | $PYTHON $RAILHOME/eval/mapping_accuracy.py -t $DATADIR/${SAMPLE}_sim.bed >${PERFORMANCE}_mapping_accuracy_summary) &
(cat Aligned.out.sam | $PYTHON $RAILHOME/eval/mapping_accuracy.py -t $DATADIR/${SAMPLE}_sim.bed -c 0.1 >${PERFORMANCE}_mapping_accuracy_SC_summary) &
wait
# Move Subjunc results to final destination
mv $OUTPUT/subjunc $SAMPLEOUTPUT
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
mv $OUTPUT/rail $SAMPLEOUTPUT
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
# STAR protocol for 2-pass execution w/ index construction is described on pp. 43-44 of the supplement of the RGASP spliced alignment paper
# (http://www.nature.com/nmeth/journal/v10/n12/extref/nmeth.2722-S1.pdf)
echo 'Running STAR on sample '${SAMPLE}' with no annotation and in paired-end mode...'
echo '#'${SAMPLE}' STAR 1-pass noann paired' >>$TIMELOG
mkdir -p $OUTPUT/star/noann_paired_1pass
cd $OUTPUT/star/noann_paired_1pass
time ($STAR --genomeDir $STARIDX --readFilesIn ${SCRATCH}/${SAMPLE}_sim_left.fastq ${SCRATCH}/${SAMPLE}_sim_right.fastq --runThreadN $CORES 2>&1) 2>>$TIMELOG
echo 'Computing precision and recall...'
(cat Aligned.out.sam | $PYTHON $RAILHOME/eval/spliced_read_recovery_performance.py -g -t $DATADIR/${SAMPLE}_sim.bed >$PERFORMANCE 2>${PERFORMANCE}_summary) &
(cat Aligned.out.sam | $PYTHON $RAILHOME/eval/intron_recovery_performance.py -t $DATADIR/${SAMPLE}_sim.bed >${PERFORMANCE}_intron_recovery_summary) &
(cat Aligned.out.sam | $PYTHON $RAILHOME/eval/mapping_accuracy.py -t $DATADIR/${SAMPLE}_sim.bed >${PERFORMANCE}_mapping_accuracy_summary) &
(cat Aligned.out.sam | $PYTHON $RAILHOME/eval/mapping_accuracy.py -t $DATADIR/${SAMPLE}_sim.bed -c 0.1 >${PERFORMANCE}_mapping_accuracy_SC_summary) &
wait
echo 'Running STAR on sample '${SAMPLE}' with annotation and in paired-end mode...'
echo '#'${SAMPLE}' STAR 1-pass ann paired' >>$TIMELOG
mkdir -p $OUTPUT/star/ann_paired_1pass
cd $OUTPUT/star/ann_paired_1pass
time ($STAR --genomeDir $STARANNIDX --readFilesIn ${SCRATCH}/${SAMPLE}_sim_left.fastq ${SCRATCH}/${SAMPLE}_sim_right.fastq --runThreadN $CORES 2>&1) 2>>$TIMELOG
echo 'Computing precision and recall...'
(cat Aligned.out.sam | $PYTHON $RAILHOME/eval/spliced_read_recovery_performance.py -g -t $DATADIR/${SAMPLE}_sim.bed >$PERFORMANCE 2>${PERFORMANCE}_summary) &
(cat Aligned.out.sam | $PYTHON $RAILHOME/eval/intron_recovery_performance.py -t $DATADIR/${SAMPLE}_sim.bed >${PERFORMANCE}_intron_recovery_summary) &
(cat Aligned.out.sam | $PYTHON $RAILHOME/eval/mapping_accuracy.py -t $DATADIR/${SAMPLE}_sim.bed >${PERFORMANCE}_mapping_accuracy_summary) &
(cat Aligned.out.sam | $PYTHON $RAILHOME/eval/mapping_accuracy.py -t $DATADIR/${SAMPLE}_sim.bed -c 0.1 >${PERFORMANCE}_mapping_accuracy_SC_summary) &
wait
echo 'Creating new STAR index for sample '${SAMPLE}' with no annotation and in paired-end mode...'
echo '#'${SAMPLE}' STAR 1-pass noann paired index' >>$TIMELOG
STARIDXNOANNPAIRED=$OUTPUT/star/noann_paired_idx
mkdir -p $STARIDXNOANNPAIRED
time ($STAR --runMode genomeGenerate --genomeDir $STARIDXNOANNPAIRED --genomeFastaFiles $FADIR/chr{1..22}.fa $FADIR/chr{X,Y,M}.fa \
		--sjdbFileChrStartEnd $OUTPUT/star/noann_paired_1pass/SJ.out.tab --sjdbOverhang $OVERHANG --runThreadN $CORES 2>&1) 2>>$TIMELOG
echo 'Creating new STAR index for sample '${SAMPLE}' with annotation and in paired-end mode...'
echo '#'${SAMPLE}' STAR 1-pass ann paired index' >>$TIMELOG
STARIDXANNPAIRED=$OUTPUT/star/ann_paired_idx
mkdir -p $STARIDXANNPAIRED
time ($STAR --runMode genomeGenerate --genomeDir $STARIDXANNPAIRED --genomeFastaFiles $FADIR/chr{1..22}.fa $FADIR/chr{X,Y,M}.fa \
		--sjdbFileChrStartEnd $OUTPUT/star/ann_paired_1pass/SJ.out.tab --sjdbOverhang $OVERHANG --runThreadN $CORES 2>&1) 2>>$TIMELOG
echo 'Running second pass of STAR on sample '${SAMPLE}' with no annotation and in paired-end mode...'
echo '#'${SAMPLE}' STAR 2-pass noann paired' >>$TIMELOG
mkdir -p $OUTPUT/star/noann_paired_2pass
cd $OUTPUT/star/noann_paired_2pass
time ($STAR --genomeDir $STARIDXNOANNPAIRED --readFilesIn ${SCRATCH}/${SAMPLE}_sim_left.fastq ${SCRATCH}/${SAMPLE}_sim_right.fastq --runThreadN $CORES 2>&1) 2>>$TIMELOG
echo 'Computing precision and recall...'
(cat Aligned.out.sam | $PYTHON $RAILHOME/eval/spliced_read_recovery_performance.py -g -t $DATADIR/${SAMPLE}_sim.bed >$PERFORMANCE 2>${PERFORMANCE}_summary) &
(cat Aligned.out.sam | $PYTHON $RAILHOME/eval/intron_recovery_performance.py -t $DATADIR/${SAMPLE}_sim.bed >${PERFORMANCE}_intron_recovery_summary) &
(cat Aligned.out.sam | $PYTHON $RAILHOME/eval/mapping_accuracy.py -t $DATADIR/${SAMPLE}_sim.bed >${PERFORMANCE}_mapping_accuracy_summary) &
(cat Aligned.out.sam | $PYTHON $RAILHOME/eval/mapping_accuracy.py -t $DATADIR/${SAMPLE}_sim.bed -c 0.1 >${PERFORMANCE}_mapping_accuracy_SC_summary) &
wait
echo 'Running second pass of STAR on sample '${SAMPLE}' with annotation and in paired-end mode...'
echo '#'${SAMPLE}' STAR 2-pass ann paired' >>$TIMELOG
mkdir -p $OUTPUT/star/ann_paired_2pass
cd $OUTPUT/star/ann_paired_2pass
time ($STAR --genomeDir $STARIDXANNPAIRED --readFilesIn ${SCRATCH}/${SAMPLE}_sim_left.fastq ${SCRATCH}/${SAMPLE}_sim_right.fastq --runThreadN $CORES 2>&1) 2>>$TIMELOG
echo 'Computing precision and recall...'
(cat Aligned.out.sam | $PYTHON $RAILHOME/eval/spliced_read_recovery_performance.py -g -t $DATADIR/${SAMPLE}_sim.bed >$PERFORMANCE 2>${PERFORMANCE}_summary) &
(cat Aligned.out.sam | $PYTHON $RAILHOME/eval/intron_recovery_performance.py -t $DATADIR/${SAMPLE}_sim.bed >${PERFORMANCE}_intron_recovery_summary) &
(cat Aligned.out.sam | $PYTHON $RAILHOME/eval/mapping_accuracy.py -t $DATADIR/${SAMPLE}_sim.bed >${PERFORMANCE}_mapping_accuracy_summary) &
(cat Aligned.out.sam | $PYTHON $RAILHOME/eval/mapping_accuracy.py -t $DATADIR/${SAMPLE}_sim.bed -c 0.1 >${PERFORMANCE}_mapping_accuracy_SC_summary) &
wait
echo 'Running STAR on sample '${SAMPLE}' with no regenerated genome/no annotation and in paired-end mode...'
echo '#'${SAMPLE}' STAR 2-pass nogen noann paired' >>$TIMELOG
mkdir -p $OUTPUT/star/nogen_noann_paired_2pass
cd $OUTPUT/star/nogen_noann_paired_2pass
time ($STAR --genomeDir $STARIDX --readFilesIn ${SCRATCH}/${SAMPLE} ${SCRATCH}/${SAMPLE} --runThreadN $CORES --twopass1readsN 50000000 --sjdbOverhang $OVERHANG >&1) 2>>$TIMELOG
echo 'Computing precision and recall...'
(cat Aligned.out.sam | $PYTHON $RAILHOME/eval/spliced_read_recovery_performance.py -g -t $DATADIR/${SAMPLE}_sim.bed >$PERFORMANCE 2>${PERFORMANCE}_summary) &
(cat Aligned.out.sam | $PYTHON $RAILHOME/eval/intron_recovery_performance.py -t $DATADIR/${SAMPLE}_sim.bed >${PERFORMANCE}_intron_recovery_summary) &
(cat Aligned.out.sam | $PYTHON $RAILHOME/eval/mapping_accuracy.py -t $DATADIR/${SAMPLE}_sim.bed >${PERFORMANCE}_mapping_accuracy_summary) &
(cat Aligned.out.sam | $PYTHON $RAILHOME/eval/mapping_accuracy.py -t $DATADIR/${SAMPLE}_sim.bed -c 0.1 >${PERFORMANCE}_mapping_accuracy_SC_summary) &
wait
# Move STAR results to final destination
mv $OUTPUT/star $SAMPLEOUTPUT
