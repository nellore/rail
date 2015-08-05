#!/usr/bin/env bash
# $1: number of cores
# $2: output directory -- SPECIFY FULL PATH
# $3: where to find sample fastqs
# $4: scratch directory; files are written here first, and relevant output is copied back to $2
# Do perform our timing simulation, we chose the GEUVADIS sample {ERR205018_1.fastq.gz, ERR205018_2.fastq.gz} at random and ran:
# taskset -c 1,2,3,4,5,6,7,8 sh run_single_sample_sim.sh 8 /scratch0/langmead-fs1/geuvadis_sims_for_paper_v2/8core_timing /scratch0/langmead-fs1/data/big_public_datasets/geuvadis /scratch2/langmead-fs1/tmp_timing
# taskset -c 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16 sh run_single_sample_sim.sh 16 /scratch0/langmead-fs1/geuvadis_sims_for_paper_v2/16core_timing /scratch0/langmead-fs1/data/big_public_datasets/geuvadis /scratch2/langmead-fs1/tmp_timing
# sh run_single_sample_sim.sh 32 /scratch0/langmead-fs1/geuvadis_sims_for_paper_v2/32core_timing /scratch0/langmead-fs1/data/big_public_datasets/geuvadis /scratch2/langmead-fs1/tmp_timing
# on a machine with 32 cores.

# Specify number of parallel processes for each program
CORES=$1

# Specify FULL PATH to output directory
MAINOUTPUT=$2
mkdir -p ${MAINOUTPUT}

# Specify data directory; fastqs should be of the form [SAMPLE NAME]_sim.fastq; Flux beds should be
# of form [SAMPLE_NAME]_sim.bed
DATADIR=$3

# Temp dir
SCRATCH=$4
SAMPLE=ERR205018
mkdir -p ${SCRATCH}
cd ${SCRATCH}
echo 'Writing FASTQs...'
gzip -cd $DATADIR/ERR205018_1.fastq.gz >${SCRATCH}/${SAMPLE}_sim_left.fastq
gzip -cd $DATADIR/ERR205018_2.fastq.gz >${SCRATCH}/${SAMPLE}_sim_right.fastq

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
SAMPLEOUTPUT=${MAINOUTPUT}/${SAMPLE}

# Run simulations
OUTPUT=$SCRATCH/${SAMPLE}
echo 'Running Subjunc on sample '${SAMPLE}' in paired-end mode...'
echo '#'${SAMPLE}' Subjunc' >>$TIMELOG
mkdir -p $OUTPUT/subjunc
cd $OUTPUT/subjunc
# Use Subjunc defaults
time (${SUBJUNC} -T $CORES -d 50 -D 600 -i ${SUBJUNCIDX} -r ${SCRATCH}/${SAMPLE}_sim_left.fastq -R ${SCRATCH}/${SAMPLE}_sim_right.fastq -o Aligned.out.sam) 2>>$TIMELOG
# Move Subjunc results to final destination
rm -rf ${SAMPLEOUTPUT}/subjunc
cp -r ${OUTPUT}/subjunc $SAMPLEOUTPUT
rm -rf ${OUTPUT}/subjunc
echo 'Running HISAT on sample '${SAMPLE}' with no annotation and in paired-end mode...'
echo '#'${SAMPLE}' HISAT 1-pass noann paired' >>$TIMELOG
mkdir -p $OUTPUT/hisat/noann_paired_1pass
cd $OUTPUT/hisat/noann_paired_1pass
time ($HISAT -x $HISATIDX -1 ${SCRATCH}/${SAMPLE}_sim_left.fastq -2 ${SCRATCH}/${SAMPLE}_sim_right.fastq -p $CORES -S Aligned.out.sam --novel-splicesite-outfile novel_splice_sites.txt 2>&1) 2>>$TIMELOG
echo 'Running HISAT on sample '${SAMPLE}' with annotation and in paired-end mode...'
echo '#'${SAMPLE}' HISAT 1-pass ann paired' >>$TIMELOG
mkdir -p $OUTPUT/hisat/ann_paired_1pass
cd $OUTPUT/hisat/ann_paired_1pass
time ($HISAT -x $HISATIDX -1 ${SCRATCH}/${SAMPLE}_sim_left.fastq -2 ${SCRATCH}/${SAMPLE}_sim_right.fastq -p $CORES -S Aligned.out.sam --novel-splicesite-outfile novel_splice_sites.txt --novel-splicesite-infile $HISATANNOTATION 2>&1) 2>>$TIMELOG
echo 'Running second pass of HISAT on sample '${SAMPLE}' with no annotation and in paired-end mode...'
echo '#'${SAMPLE}' HISAT 2-pass noann paired' >>$TIMELOG
mkdir -p $OUTPUT/hisat/noann_paired_2pass
cd $OUTPUT/hisat/noann_paired_2pass
time ($HISAT -x $HISATIDX -1 ${SCRATCH}/${SAMPLE}_sim_left.fastq -2 ${SCRATCH}/${SAMPLE}_sim_right.fastq -p $CORES -S Aligned.out.sam --novel-splicesite-infile $OUTPUT/hisat/noann_paired_1pass/novel_splice_sites.txt 2>&1) 2>>$TIMELOG
echo 'Running second pass of HISAT on sample '${SAMPLE}' with annotation and in paired-end mode...'
echo '#'${SAMPLE}' HISAT 2-pass ann paired' >>$TIMELOG
mkdir -p $OUTPUT/hisat/ann_paired_2pass
cd $OUTPUT/hisat/ann_paired_2pass
time ($HISAT -x $HISATIDX -1 ${SCRATCH}/${SAMPLE}_sim_left.fastq -2 ${SCRATCH}/${SAMPLE}_sim_right.fastq -p $CORES -S Aligned.out.sam --novel-splicesite-infile $OUTPUT/hisat/ann_paired_1pass/novel_splice_sites.txt 2>&1) 2>>$TIMELOG
# Move all hisat results to final destination
rm -rf ${SAMPLEOUTPUT}/hisat
cp -r ${OUTPUT}/hisat $SAMPLEOUTPUT
rm -rf ${OUTPUT}/hisat
echo 'Running Rail-RNA on sample '${SAMPLE}', outputting only BAM...'
echo '#'${SAMPLE}' Rail-RNA only bam' >>$TIMELOG
# Write manifest file
echo -e ${SCRATCH}/${SAMPLE}_sim_left.fastq'\t0\t'${SCRATCH}/${SAMPLE}_sim_right.fastq'\t0\t'${SAMPLE}'-1-1' >${SCRATCH}/${SAMPLE}.manifest
time ($RAILRNA go local -p $CORES -m ${SCRATCH}/${SAMPLE}.manifest -o $OUTPUT/railb --log $OUTPUT/railb.log -x $BOWTIE1IDX,$BOWTIE2IDX -f -d bam >/dev/null 2>&1) 2>>$TIMELOG
# Move railb results to final destination
rm -rf ${SAMPLEOUTPUT}/railb
cp -r ${OUTPUT}/railb $SAMPLEOUTPUT
rm -rf ${OUTPUT}/railb
echo 'Running Rail-RNA on sample '${SAMPLE}' with all default outputs...'
echo '#'${SAMPLE}' Rail-RNA all default outputs' >>$TIMELOG
# Write manifest file
echo -e ${SCRATCH}/${SAMPLE}_sim_left.fastq'\t0\t'${SCRATCH}/${SAMPLE}_sim_right.fastq'\t0\t'${SAMPLE}'-1-1' >${SCRATCH}/${SAMPLE}.manifest
time ($RAILRNA go local -p $CORES -m ${SCRATCH}/${SAMPLE}.manifest -o $OUTPUT/raild --log $OUTPUT/raild.log -x $BOWTIE1IDX,$BOWTIE2IDX -f >/dev/null 2>&1) 2>>$TIMELOG
# Move raild results to final destination
rm -rf ${SAMPLEOUTPUT}/raild
cp -r ${OUTPUT}/raild $SAMPLEOUTPUT
rm -rf ${OUTPUT}/raild
echo 'Running TopHat on sample '${SAMPLE}' with no annotation and in paired-end mode...'
echo '#'${SAMPLE}' TopHat noann paired' >>$TIMELOG
time ($TOPHAT -o $OUTPUT/tophat/noann_paired -p $CORES $BOWTIE2IDX ${SCRATCH}/${SAMPLE}_sim_left.fastq ${SCRATCH}/${SAMPLE}_sim_right.fastq 2>&1) 2>>$TIMELOG
echo 'Running TopHat on sample '${SAMPLE}' with annotation and in paired-end mode...'
echo '#'${SAMPLE}' TopHat ann paired' >>$TIMELOG
time ($TOPHAT -o $OUTPUT/tophat/ann_paired -G $ANNOTATION -p $CORES $BOWTIE2IDX ${SCRATCH}/${SAMPLE}_sim_left.fastq ${SCRATCH}/${SAMPLE}_sim_right.fastq 2>&1) 2>>$TIMELOG
# Move TopHat results to final destination
rm -rf ${SAMPLEOUTPUT}/tophat
cp -r ${OUTPUT}/tophat $SAMPLEOUTPUT
rm -rf ${OUTPUT}/tophat
# STAR protocol for 2-pass execution w/ index construction is described on pp. 43-44 of the supplement of the RGASP spliced alignment paper
# (http://www.nature.com/nmeth/journal/v10/n12/extref/nmeth.2722-S1.pdf)
echo 'Running STAR on sample '${SAMPLE}' with no regenerated genome/no annotation and in paired-end mode...'
echo '#'${SAMPLE}' STAR 2-pass nogen noann paired' >>$TIMELOG
mkdir -p $OUTPUT/star/nogen_noann_paired_2pass
cd $OUTPUT/star/nogen_noann_paired_2pass
time ($STAR --genomeDir $STARIDX --readFilesIn ${SCRATCH}/${SAMPLE}_sim_left.fastq ${SCRATCH}/${SAMPLE}_sim_right.fastq --runThreadN $CORES --twopass1readsN -1 --sjdbOverhang $OVERHANG --twopassMode Basic >&1) 2>>$TIMELOG
echo 'Running STAR on sample '${SAMPLE}' with no annotation and in paired-end mode...'
echo '#'${SAMPLE}' STAR 1-pass noann paired' >>$TIMELOG
mkdir -p $OUTPUT/star/noann_paired_1pass
cd $OUTPUT/star/noann_paired_1pass
time ($STAR --genomeDir $STARIDX --readFilesIn ${SCRATCH}/${SAMPLE}_sim_left.fastq ${SCRATCH}/${SAMPLE}_sim_right.fastq --runThreadN $CORES 2>&1) 2>>$TIMELOG
echo 'Running STAR on sample '${SAMPLE}' with annotation and in paired-end mode...'
echo '#'${SAMPLE}' STAR 1-pass ann paired' >>$TIMELOG
mkdir -p $OUTPUT/star/ann_paired_1pass
cd $OUTPUT/star/ann_paired_1pass
time ($STAR --genomeDir $STARANNIDX --readFilesIn ${SCRATCH}/${SAMPLE}_sim_left.fastq ${SCRATCH}/${SAMPLE}_sim_right.fastq --runThreadN $CORES 2>&1) 2>>$TIMELOG
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
