#!/bin/sh

## Usage
# sh run-all.sh railGEU

# Define variables
EXPERIMENT=$1

mkdir -p ${EXPERIMENT}/CoverageInfo
mkdir -p ${EXPERIMENT}/regionMatrix

sh step1-fullCoverage.sh ${EXPERIMENT}

sh step6-regionMatrix.sh ${EXPERIMENT}
sh step6b-regionMatrix-merge.sh ${EXPERIMENT}
