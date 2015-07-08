#!/usr/bin/env bash
## Create indexes; use hg19 chr{1..22} + chr{X,Y,M}
REFDIR=/scratch0/langmead-fs1/indexes_for_paper
# Where Feb 2009 hg19 chromosome files are located; download them at http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz
FADIR=/scratch0/langmead-fs1/shared/references/hg19/fasta

# Where to find Bowtie build executables; versions 1.1.1 and 2.2.4 are used here
BOWTIEBUILD=/scratch0/langmead-fs1/shared/bowtie-1.1.1/bowtie-build
BOWTIE2BUILD=/scratch0/langmead-fs1/shared/bowtie2-2.2.4/bowtie2-build
# Where to find STAR; version 2.4.0j is used here
STAR=/scratch0/langmead-fs1/shared/STAR-STAR_2.4.0j/bin/Linux_x86_64_static
# Where to find HISAT build executable; version 0.1.1-beta is used
HISATBUILD=/scratch0/langmead-fs1/shared/hisat-0.1.1-beta/hisat-build
#Where to find Subread build executable; version 1.4.6-p4 is used
SUBREADBUILD=/scratch0/langmead-fs1/shared/subread-1.4.6-p4-Linux-x86_64/bin/subread-buildindex

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
$HISATBUILD genome.fa hisatgenome || { echo 'Problem encountered building HISAT index'; exit 1; }
$SUBREADBUILD -o subreadgenome genome.fa || { echo 'Problem encountered building Subread index'; exit 1; }
