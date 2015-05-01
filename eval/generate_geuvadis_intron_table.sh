#!/usr/bin/env bash
CWD=$(pwd)
cd $1
mkdir -p intermediate
for i in *junctions*.bed; do (python $CWD/grab_introns.py $i >intermediate/$i.intermediate); done
python write_table.py $CWD/intermediate | gzip >geuvadis_intron_table.tsv.gz