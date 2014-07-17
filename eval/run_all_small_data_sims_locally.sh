#!/usr/bin/env bash
# RUN WITH TASKSET to limit CPU usage on multicore system for benchmarking!
# Ex: taskset -c 0,1,2,3 sh run_all_small_data_sims_locally.sh
# Select two sample names for analysis. See generate_bioreps.py for how sample data was generated.
SAMPLE1=NA18861.1.M_120209_2
SAMPLE2=NA18508.1.M_111124_1

# Specify data directory; fastqs should be of the form [SAMPLE NAME]_sim.fastq
DATADIR=/scratch0/langmead-fs1/geuvadis_sim

## Specify locations of executables
# Used version 2.0.12 of TopHat; wrapped version 2.2.2 of Bowtie2
TOPHAT=/scratch0/langmead-fs1/shared/tophat-2.0.12.Linux_x86_64/tophat2
# Used version 2.3.0e of STAR
STAR=/scratch0/langmead-fs1/shared/STAR_2.3.0e.Linux_x86_64/STAR
# Used version 0.1.0 of Rail-RNA; wrapped version 2.2.2 of Bowtie2
# Specify Python executable/loc of get_junctions.py; PyPy 2.2.1 was used
PYTHON=pypy
RAILHOME=/scratch0/langmead-fs1/rail
RAILRNA=$PYTHON\ $RAILHOME/src
# Samtools v0.1.19-44428cd was used
SAMTOOLS=samtools

# Specify number of parallel processes for each program
CORES=4

# Specify FULL PATH to output directory
MAINOUTPUT=/scratch0/langmead-fs1/geuvadis_sim/4cpus
mkdir -p $MAINOUTPUT

# Specify log filename for recording times
TIMELOG=$MAINOUTPUT/small_data_times.log

## Specify locations of reference-related files
## See create_indexes.sh for index creation script
# Bowtie indexes
BOWTIE1IDX=/scratch0/langmead-fs1/indexes_for_paper/genome
BOWTIE2IDX=/scratch0/langmead-fs1/indexes_for_paper/genome
# STAR index
STARIDX=/scratch0/langmead-fs1/indexes_for_paper/star
# Overhang length for STAR; this should ideally be max read length - 1
OVERHANG=75

# Flux outputs paired-end reads in one file; split files here
echo 'Splitting Flux FASTQs...'
for SAMPLE in {$SAMPLE1,$SAMPLE2}
do
	awk '(NR-1) % 8 < 4' $DATADIR/${SAMPLE}_sim.fastq >$DATADIR/${SAMPLE}_sim_left.fastq
	awk '(NR-1) % 8 >= 4' $DATADIR/${SAMPLE}_sim.fastq >$DATADIR/${SAMPLE}_sim_right.fastq
done

## Specify location of annotation
# This is Gencode v12, which may be obtained at ftp://ftp.sanger.ac.uk/pub/gencode/release_12/gencode.v12.annotation.gtf.gz
ANNOTATION=/scratch0/langmead-fs1/geuvadis_sim/gencode.v12.annotation.gtf

## STAR requires its own annotation format that lists introns directly; conversion is done by get_junctions.py
# Build STAR junction index from GTF
STARANNOTATION=$MAINOUTPUT/junctions_for_star.txt
echo 'Building annotation for STAR from GTF...'
cat $ANNOTATION | $PYTHON $RAILHOME/eval/get_junctions.py >$STARANNOTATION

## Where Feb 2009 hg19 chromosome files are located; download them at http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz
# These are here because it's necessary to build a new STAR index for its 2-pass and annotation protocols
FADIR=/scratch0/langmead-fs1/shared/references/hg19/fasta

# Create STAR index including annotation
echo 'Creating new STAR index including splice junctions...'
echo '#STAR index pre-1-pass ann' >>$TIMELOG
STARANNIDX=$MAINOUTPUT/starannidx
mkdir -p $STARANNIDX
# Use --sjdbOverhang 75 because reads are 76 bases long! See p. 3 of STAR manual for details.
time ($STAR --runMode genomeGenerate --genomeDir $STARANNIDX --genomeFastaFiles $FADIR/chr{1..22}.fa $FADIR/chr{X,Y,M}.fa \
		--runThreadN $CORES --sjdbFileChrStartEnd $STARANNOTATION --sjdbOverhang $OVERHANG 2>&1) 2>>$TIMELOG

cd $MAINOUTPUT
for SAMPLE in {$SAMPLE1,$SAMPLE2}
do
	mkdir -p ${SAMPLE}
	cd ${SAMPLE}
	mkdir -p tophat
	mkdir -p star
	mkdir -p rail
	cd ..
done

# Run simulations
for SAMPLE in {$SAMPLE1,$SAMPLE2}
do
	OUTPUT=$MAINOUTPUT/${SAMPLE}
	echo 'Running TopHat on sample '${SAMPLE}' with no annotation and in single-end mode...'
	echo '#'${SAMPLE}' TopHat noann single' >>$TIMELOG
	time ($TOPHAT -o $OUTPUT/tophat/noann_single -p $CORES $BOWTIE2IDX $DATADIR/${SAMPLE}_sim.fastq 2>&1) 2>>$TIMELOG
	echo 'Running TopHat on sample '${SAMPLE}' with no annotation and in paired-end mode...'
	echo '#'${SAMPLE}' TopHat noann paired' >>$TIMELOG
	time ($TOPHAT -o $OUTPUT/tophat/noann_paired -p $CORES $BOWTIE2IDX $DATADIR/${SAMPLE}_sim_left.fastq $DATADIR/${SAMPLE}_sim_right.fastq 2>&1) 2>>$TIMELOG
	echo 'Running TopHat on sample '${SAMPLE}' with annotation and in single-end mode...'
	echo '#'${SAMPLE}' TopHat ann single' >>$TIMELOG
	time ($TOPHAT -o $OUTPUT/tophat/ann_single -G $ANNOTATION -p $CORES $BOWTIE2IDX $DATADIR/${SAMPLE}_sim.fastq 2>&1) 2>>$TIMELOG
	echo 'Running TopHat on sample '${SAMPLE}' with annotation and in paired-end mode...'
	echo '#'${SAMPLE}' TopHat ann paired' >>$TIMELOG
	time ($TOPHAT -o $OUTPUT/tophat/ann_paired -G $ANNOTATION -p $CORES $BOWTIE2IDX $DATADIR/${SAMPLE}_sim_left.fastq $DATADIR/${SAMPLE}_sim_right.fastq 2>&1) 2>>$TIMELOG
	# STAR protocols for execution are described on pp. 43-44 of the supplement of the RGASP spliced alignment paper
	# (http://www.nature.com/nmeth/journal/v10/n12/extref/nmeth.2722-S1.pdf)
	echo 'Running STAR on sample '${SAMPLE}' with no annotation and in single-end mode...'
	echo '#'${SAMPLE}' STAR 1-pass noann single' >>$TIMELOG
	mkdir -p $OUTPUT/star/noann_single_1pass
	cd $OUTPUT/star/noann_single_1pass
	time ($STAR --genomeDir $STARIDX --readFilesIn $DATADIR/${SAMPLE}_sim.fastq --runThreadN $CORES 2>&1) 2>>$TIMELOG
	echo 'Running STAR on sample '${SAMPLE}' with no annotation and in paired-end mode...'
	echo '#'${SAMPLE}' STAR 1-pass noann paired' >>$TIMELOG
	mkdir -p $OUTPUT/star/noann_paired_1pass
	cd $OUTPUT/star/noann_paired_1pass
	time ($STAR --genomeDir $STARIDX --readFilesIn $DATADIR/${SAMPLE}_sim_left.fastq $DATADIR/${SAMPLE}_sim_right.fastq --runThreadN $CORES 2>&1) 2>>$TIMELOG
	echo 'Running STAR on sample '${SAMPLE}' with annotation and in single-end mode...'
	echo '#'${SAMPLE}' STAR 1-pass ann single' >>$TIMELOG
	mkdir -p $OUTPUT/star/ann_single_1pass
	cd $OUTPUT/star/ann_single_1pass
	time ($STAR --genomeDir $STARANNIDX --readFilesIn $DATADIR/${SAMPLE}_sim.fastq --runThreadN $CORES 2>&1) 2>>$TIMELOG
	echo 'Running STAR on sample '${SAMPLE}' with annotation and in paired-end mode...'
	echo '#'${SAMPLE}' STAR 1-pass ann paired' >>$TIMELOG
	mkdir -p $OUTPUT/star/ann_paired_1pass
	cd $OUTPUT/star/ann_paired_1pass
	time ($STAR --genomeDir $STARANNIDX --readFilesIn $DATADIR/${SAMPLE}_sim_left.fastq $DATADIR/${SAMPLE}_sim_right.fastq --runThreadN $CORES 2>&1) 2>>$TIMELOG
	echo 'Creating new STAR index for sample '${SAMPLE}' with no annotation and in single-end mode...'
	echo '#'${SAMPLE}' STAR 1-pass noann single index' >>$TIMELOG
	STARIDXNOANNSINGLE=$OUTPUT/star/noann_single_idx
	mkdir -p $STARIDXNOANNSINGLE
	time ($STAR --runMode genomeGenerate --genomeDir $STARIDXNOANNSINGLE --genomeFastaFiles $FADIR/chr{1..22}.fa $FADIR/chr{X,Y,M}.fa \
			--sjdbFileChrStartEnd $OUTPUT/star/noann_single_1pass/SJ.out.tab --sjdbOverhang $OVERHANG --runThreadN $CORES 2>&1) 2>>$TIMELOG
	echo 'Creating new STAR index for sample '${SAMPLE}' with no annotation and in paired-end mode...'
	echo '#'${SAMPLE}' STAR 1-pass noann paired index' >>$TIMELOG
	STARIDXNOANNPAIRED=$OUTPUT/star/noann_paired_idx
	mkdir -p $STARIDXNOANNPAIRED
	time ($STAR --runMode genomeGenerate --genomeDir $STARIDXNOANNPAIRED --genomeFastaFiles $FADIR/chr{1..22}.fa $FADIR/chr{X,Y,M}.fa \
			--sjdbFileChrStartEnd $OUTPUT/star/noann_paired_1pass/SJ.out.tab --sjdbOverhang $OVERHANG --runThreadN $CORES 2>&1) 2>>$TIMELOG
	echo 'Creating new STAR index for sample '${SAMPLE}' with annotation and in single-end mode...'
	echo '#'${SAMPLE}' STAR 1-pass ann single index' >>$TIMELOG
	STARIDXANNSINGLE=$OUTPUT/star/ann_single_idx
	mkdir -p $STARIDXANNSINGLE
	time ($STAR --runMode genomeGenerate --genomeDir $STARIDXANNSINGLE --genomeFastaFiles $FADIR/chr{1..22}.fa $FADIR/chr{X,Y,M}.fa \
			--sjdbFileChrStartEnd $OUTPUT/star/ann_single_1pass/SJ.out.tab --sjdbOverhang $OVERHANG --runThreadN $CORES 2>&1) 2>>$TIMELOG
	echo 'Creating new STAR index for sample '${SAMPLE}' with annotation and in paired-end mode...'
	echo '#'${SAMPLE}' STAR 1-pass ann paired index' >>$TIMELOG
	STARIDXANNPAIRED=$OUTPUT/star/ann_paired_idx
	mkdir -p $STARIDXANNPAIRED
	time ($STAR --runMode genomeGenerate --genomeDir $STARIDXANNPAIRED --genomeFastaFiles $FADIR/chr{1..22}.fa $FADIR/chr{X,Y,M}.fa \
			--sjdbFileChrStartEnd $OUTPUT/star/ann_paired_1pass/SJ.out.tab --sjdbOverhang $OVERHANG --runThreadN $CORES 2>&1) 2>>$TIMELOG
	echo 'Running second pass of STAR on sample '${SAMPLE}' with no annotation and in single-end mode...'
	echo '#'${SAMPLE}' STAR 2-pass noann single' >>$TIMELOG
	mkdir -p $OUTPUT/star/noann_single_2pass
	cd $OUTPUT/star/noann_single_2pass
	time ($STAR --genomeDir $STARIDXNOANNSINGLE --readFilesIn $DATADIR/${SAMPLE}_sim.fastq --runThreadN $CORES 2>&1) 2>>$TIMELOG
	echo 'Running second pass of STAR on sample '${SAMPLE}' with no annotation and in paired-end mode...'
	echo '#'${SAMPLE}' STAR 2-pass noann paired' >>$TIMELOG
	mkdir -p $OUTPUT/star/noann_paired_2pass
	cd $OUTPUT/star/noann_paired_2pass
	time ($STAR --genomeDir $STARIDXNOANNPAIRED --readFilesIn $DATADIR/${SAMPLE}_sim_left.fastq $DATADIR/${SAMPLE}_sim_right.fastq --runThreadN $CORES 2>&1) 2>>$TIMELOG
	echo 'Running second pass of STAR on sample '${SAMPLE}' with annotation and in single-end mode...'
	echo '#'${SAMPLE}' STAR 2-pass ann single' >>$TIMELOG
	mkdir -p $OUTPUT/star/ann_single_2pass
	cd $OUTPUT/star/ann_single_2pass
	time ($STAR --genomeDir $STARIDXANNSINGLE --readFilesIn $DATADIR/${SAMPLE}_sim.fastq --runThreadN $CORES 2>&1) 2>>$TIMELOG
	echo 'Running second pass of STAR on sample '${SAMPLE}' with annotation and in paired-end mode...'
	echo '#'${SAMPLE}' STAR 2-pass ann paired' >>$TIMELOG
	mkdir -p $OUTPUT/star/ann_paired_2pass
	cd $OUTPUT/star/ann_paired_2pass
	time ($STAR --genomeDir $STARIDXANNPAIRED --readFilesIn $DATADIR/${SAMPLE}_sim_left.fastq $DATADIR/${SAMPLE}_sim_right.fastq --runThreadN $CORES 2>&1) 2>>$TIMELOG
	echo 'Running Rail-RNA on sample '${SAMPLE}'...'
	echo '#'${SAMPLE}' Rail-RNA' >>$TIMELOG
	# Write manifest file
	echo -e $DATADIR/${SAMPLE}_sim.fastq'\t0\t'${SAMPLE}'-0-0' >$MAINOUTPUT/${SAMPLE}.manifest
	time ($RAILRNA go local -m $MAINOUTPUT/${SAMPLE}.manifest -o $OUTPUT/rail --log $OUTPUT/rail.log -1 $BOWTIE1IDX -2 $BOWTIE2IDX -f >/dev/null 2>&1) 2>>$TIMELOG
done