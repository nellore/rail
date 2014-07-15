#!/usr/bin/env bash

# Select two sample names for analysis
SAMPLE1=NA18861.1.M_120209_2
SAMPLE2=NA18508.1.M_111124_1

# Specify data directory; fastqs should be of the form [SAMPLE NAME]_sim.fastq
DATADIR=/scratch0/langmead-fs1/geuvadis_sim

# Flux outputs paired-end reads in one file; split them here
awk 

# Specify locations of executables
# Use version 2.0.12 of TopHat
TOPHAT=/scratch0/langmead-fs1/shared/tophat-2.0.12.Linux_x86_64/tophat2
# Use version 2.3.0e of STAR
STAR=/scratch0/langmead-fs1/shared/STAR_2.3.0e.Linux_x86_64/STAR
# Use version 0.1.0 of Rail-RNA
RAILRNA=python\ /scratch0/langmead-fs1/rail/src

