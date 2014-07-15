#!/usr/bin/env bash
## Create indexes; use hg19 chr{1..22} + chr{X,Y,M}
REF_DIR=/scratch0/langmead-fs1/indexes_for_paper
# Where Feb 2009 hg19 chromosome files are located; download them at http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz
FA_DIR=/scratch0/langmead-fs1/shared/references/hg19/fasta

# Where to find Bowtie build executables; versions 1.0.0 and 2.2.2 are used here
BOWTIE_BUILD=/scratch0/langmead-fs1/shared/bowtie-1.0.0/bowtie-build
BOWTIE2_BUILD=/scratch0/langmead-fs1/shared/bowtie2-2.2.2/bowtie2-build
# Where to find STAR; version 2.3.0e is used here
STAR=/scratch0/langmead-fs1/shared/STAR_2.3.0e.Linux_x86_64/STAR

# Number of cores
CORES=32

mkdir -p $REF_DIR
cd $FA_DIR
cat chr{1..22}.fa chr{X,Y,M}.fa >$REF_DIR/hg19.fa
cd $REF_DIR
$STAR --runMode genomeGenerate --genomeDir $REF_DIR/star --genomeFastaFiles hg19.fa --runThreadN $CORES
$BOWTIE_BUILD hg19.fa genome || { echo 'Problem encountered building Bowtie 1 index'; exit 1; }
$BOWTIE2_BUILD hg19.fa genome || { echo 'Problem encountered building Bowtie 2 index'; exit 1; }