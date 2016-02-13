#!/usr/bin/env bash
# Runs timing simulations on an i2.8xlarge EC2 machine

# Update
sudo yum update

# Set up partition
sudo mkfs.ext4 /dev/xvdi
sudo mount -t ext4 /dev/xvdi /mnt

# Download all required software
cd /mnt
mkdir aligners
sudo chown ec2-user aligners
cd aligners
# STAR v2.4.2a
wget https://github.com/alexdobin/STAR/blob/ff38b7836c132bcd91ca2cd3a23e10d59efc5382/bin/Linux_x86_64_static/STAR?raw=true
# HISAT v0.1.6-beta
wget http://www.ccb.jhu.edu/software/hisat/downloads/hisat-0.1.6-beta-Linux_x86_64.zip
unzip hisat-0.1.6-beta-Linux_x86_64.zip
# Subred 
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

# Get data
cd /mnt
sudo mkdir -p geuvadis_data
sudo chown ec2-user geuvadis_data
cd geuvadis_data
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR188/ERR18846$i/ERR188460_1.fastq.gz;
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR188/ERR18846$i/ERR188460_2.fastq.gz;
for i in 0..9; do wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR188/ERR18846$i/ERR18846$i_2.fastq.gz; done;
cd ..
# Install Rail; used v0.2.2
sudo yum install gcc
sudo yum install ncurses-devel ncurses
sudo yum install zlib-devel
(INSTALLER=/var/tmp/$(cat /dev/urandom | env LC_CTYPE=C tr -cd 'a-f0-9' | head -c 32); curl http://verve.webfactional.com/rail -o $INSTALLER; python $INSTALLER -m || true; rm -f $INSTALLER)
sudo mkdir -p rail_sim
sudo chown ec2-user rail_sim
cd rail_sim
TAB="$(printf '\t')"
cat >rail_sim.manifest <<EOF
/mnt/geuvadis_data/ERR188460_1.fastq.gz${TAB}0${TAB}/mnt/geuvadis_data/ERR188460_2.fastq.gz${TAB}0${TAB}ERR188460
EOF
sudo mkdir -p /mnt/tmp
sudo chown ec2-user /mnt/tmp
(time env TMPDIR=/mnt/tmp rail-rna go local -x /mnt/indexes/genome -m rail_sim.manifest -p 32) 2>railtime.txt

