#!/usr/bin/env bash
# Consolidates performances from sims for all in a single directory "performances" in the same directory as $1
# $1: directory in which sims were performed
# We ran sh consolidate_performances.sh /scratch0/langmead-fs1/geuvadis_sims_for_paper_v2/8core
# sh consolidate_performances.sh /dcl01/leek/data/railsims/112simresults
# sh consolidate_performances.sh /dcl01/leek/data/railsims/8core
# Then we copied only those files in each "performances" directory the with the sample name NA20768_female_TSI_HMGU_5-1-1 to the local directory
# /Users/eterna/results_v2, which is inside the SetDirectory of the Mathematica 10 notebook performance.nb
# You will have to change the argument of SetDirectory there so it matches the directory to which you copy the performance files.
cd $1
mkdir -p performances; for i in $(find . -name perform\*); do cp $i performances/${i//\//_}; done
cd performances
for i in .*; do mv $i ${i//._/}; done