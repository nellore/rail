#!/bin/sh

# e.g., for drosophila: sh igenomes_setup.sh Drosophila_melanogaster UCSC dm3

if [ -z "$IGENOMES_DIR" ] ; then
	if [ -z "$RAIL_HOME" ] ; then
		echo "Neither IGENOMES_DIR nor RAIL_HOME set"
		exit 10
	fi
	IGENOMES_DIR=$RAIL_HOME/igenomes
fi

if [ ! -d "$IGENOMES_DIR" ] ; then
	echo "No such iGenomes directory as $IGENOMES_DIR"
fi

SPECIES=$1

if [ -z "$SPECIES" ] ; then
	echo "No species specified; specify as first arg"
	exit 1
fi

DB=$2

if [ -z "$DB" ] ; then
	echo "No database specified; specify as second arg"
	exit 2
fi

ASM=$3

if [ -z "$ASM" ] ; then
	echo "No assmebly name specified; specify as third arg"
	exit 3
fi

FULL="$IGENOMES_DIR/$SPECIES/$DB/$ASM"
if [ ! -d "$FULL" ] ; then
	echo "No such iGenomes directory as $FULL"
fi

rm -f Chromosomes genes.fixed.gtf

ln -s -f $FULL/Sequence/Chromosomes Chromosomes

python gtf_flux_fix.py \
	--gtf $FULL/Annotation/Genes/genes.gtf \
	--fasta $FULL/Sequence/Chromosomes/*.fa \
	> genes.fixed.gtf
