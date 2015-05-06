#!/bin/bash
#$ -cwd
#$ -m e
#$ -l mem_free=2G,h_vmem=20G,h_fsize=50G
#$ -N derfinder_analysis
#$ -pe local 12

echo '**** Job starts ****'
date

## Run analysis on ERs
module load R/3.1.x
Rscript derfinder_analysis.R

echo '**** Job ends ****'
date
