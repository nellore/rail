#!/usr/bin/env bash
## Create indexes; use hg19 chr{1..22} + chr{X,Y,M}
REFDIR=/scratch0/langmead-fs1/indexes_for_paper
# Where Feb 2009 hg19 chromosome files are located; download them at http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz
FADIR=/scratch0/langmead-fs1/shared/references/hg19/fasta

# Where to find Bowtie build executables; versions 1.0.0 and 2.2.2 are used here
BOWTIEBUILD=/scratch0/langmead-fs1/shared/bowtie-1.0.0/bowtie-build
BOWTIE2BUILD=/scratch0/langmead-fs1/shared/bowtie2-2.2.2/bowtie2-build
# Where to find STAR; version 2.3.0e is used here
STAR=/scratch0/langmead-fs1/shared/STAR_2.3.0e.Linux_x86_64/STAR

# Number of cores
CORES=32

mkdir -p $REFDIR
cd $FADIR
cat chr{1..22}.fa chr{X,Y,M}.fa >$REFDIR/genome.fa
cd $REFDIR
mkdir -p $REFDIR/star
$STAR --runMode genomeGenerate --genomeDir $REFDIR/star --genomeFastaFiles genome.fa --runThreadN $CORES
$BOWTIEBUILD genome.fa genome || { echo 'Problem encountered building Bowtie 1 index'; exit 1; }
$BOWTIE2BUILD genome.fa genome || { echo 'Problem encountered building Bowtie 2 index'; exit 1; }