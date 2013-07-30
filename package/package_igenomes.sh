#!/bin/sh

#
# package_igenomes.sh
#
# Run from just above the database subdirectory (ususally named UCSC/NCBI/etc)
# to create a reference archive for Tornado
#

# Get assembly name and database name from the enclosing directories
pw=`pwd`
genome=`basename $pw`
tmp=`dirname $pw`
db=`basename $tmp`

mkdir -p .archive/index .archive/gtf .archive/fasta

cp Sequence/BowtieIndex/*.ebwt .archive/index
cp Annotation/Genes/genes.gtf .archive/gtf
cp Sequence/WholeGenomeFasta/genome.fa* .archive/fasta

cd .archive

tar -zcvf ${genome}_${db}.tar.gz .
mv ${genome}_${db}.tar.gz ..
cd ..

echo "DONE -- you can 'rm -rf .archive'"
