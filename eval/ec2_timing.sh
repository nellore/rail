#!/usr/bin/env bash
# Runs timing simulations on an i2.8xlarge EC2 machine
# Be sure to create all indexes with create_indexes.sh first and put them in /mnt/indexes

# Set up partition
sudo mkfs.ext4 /dev/xvdi
sudo mount -t ext4 /dev/xvdi /mnt

# Get data
cd /mnt
sudo mkdir -p geuvadis_data
sudo chown ec2-user geuvadis_data
cd geuvadis_data
for i in 0..9; do wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR188/ERR18846$i/ERR18846$i_1.fastq.gz; done;
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
/mnt/geuvadis_data/ERR188461_1.fastq.gz${TAB}0${TAB}/mnt/geuvadis_data/ERR188461_2.fastq.gz${TAB}0${TAB}ERR188461
/mnt/geuvadis_data/ERR188462_1.fastq.gz${TAB}0${TAB}/mnt/geuvadis_data/ERR188462_2.fastq.gz${TAB}0${TAB}ERR188462
/mnt/geuvadis_data/ERR188463_1.fastq.gz${TAB}0${TAB}/mnt/geuvadis_data/ERR188463_2.fastq.gz${TAB}0${TAB}ERR188463
/mnt/geuvadis_data/ERR188464_1.fastq.gz${TAB}0${TAB}/mnt/geuvadis_data/ERR188464_2.fastq.gz${TAB}0${TAB}ERR188464
/mnt/geuvadis_data/ERR188465_1.fastq.gz${TAB}0${TAB}/mnt/geuvadis_data/ERR188465_2.fastq.gz${TAB}0${TAB}ERR188465
/mnt/geuvadis_data/ERR188466_1.fastq.gz${TAB}0${TAB}/mnt/geuvadis_data/ERR188466_2.fastq.gz${TAB}0${TAB}ERR188466
/mnt/geuvadis_data/ERR188467_1.fastq.gz${TAB}0${TAB}/mnt/geuvadis_data/ERR188467_2.fastq.gz${TAB}0${TAB}ERR188467
/mnt/geuvadis_data/ERR188468_1.fastq.gz${TAB}0${TAB}/mnt/geuvadis_data/ERR188468_2.fastq.gz${TAB}0${TAB}ERR188468
/mnt/geuvadis_data/ERR188469_1.fastq.gz${TAB}0${TAB}/mnt/geuvadis_data/ERR188469_2.fastq.gz${TAB}0${TAB}ERR188469
EOF
sudo mkdir -p /mnt/tmp
sudo chown ec2-user /mnt/tmp
(time env TMPDIR=/mnt/tmp rail-rna go local -x /mnt/indexes/genome -m rail_sim.manifest -p 32) 2>railtime.txt

