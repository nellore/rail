#!/usr/bin/env bash
# Consolidates performances from all sims, local and elastic, in one directory
# $1: directory in which sims were performed
cd $1
mkdir performances; for i in $(find . -name \*perform\*); do cp $i performances/; done
cd performances
for i in .*; do mv $i ${i//._/}; done