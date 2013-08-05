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

# Put index, gene annotations, and fasta+index in separate archives, since we
# might not always need all three
for i in index gtf fasta ; do
	tar -zcvf ${genome}_${db}.$i.tar.gz $i
	mv ${genome}_${db}.$i.idx.tar.gz ..
done

cd ..

echo "DONE -- you can 'rm -rf .archive'"
