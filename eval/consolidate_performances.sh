#!/usr/bin/env bash
# Consolidates performances from all sims, local and elastic, in one directory
# $1: directory in which sims were performed
# $2: directory in which to put performances
# We ran sh consolidate_performances.sh /scratch0/langmead-fs1/geuvadis_sims_for_paper /scratch0/langmead-fs1/geuvadis_sims_for_paper/performances
cd $1
mkdir -p performances; for i in $(find . -name perform\*); do cp $i performances/${i//\//_}; done
cd performances
for i in .*; do mv $i ${i//._/}; done