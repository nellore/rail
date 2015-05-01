#!/usr/bin/env bash
CWD=$(pwd)
cd $1
mkdir -p intermediate
for i in *junctions*.bed; do (python $CWD/grab_introns.py $i >intermediate/$i.intermediate); done
python write_table.py $CWD/intermediate | gzip >GEUVADIS_introns_v4.04.2015.tsv.gz
gzip -cd GEUVADIS_introns_v4.04.2015.tsv.gz | awk '{print $1}' | tail -n +2 | awk -F ';' '{print substr($1, 1, length($1)-1) "\t" $2 "\t" $3 "\t" substr($1, length($1), length($1)) }' >intron_table.tsv