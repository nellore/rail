#!/usr/bin/env bash
# Script run on Johns Hopkins cluster to generate PBS scripts for running simulations
# May have to be adjusted by the user to suit other clusters
# STAR/subjunc simulations were not run on cluster because of memory requirement; they were run separately on a high-memory node
DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
cd $DIR
mkdir -p pbs_scripts
cd pbs_scripts
i=0; IFS=$'\n'; for j in $(python ../create_single_sample_sim_commands.py --aligners tophat rail); do python -c "print '#!/bin/bash\n#PBS -q batch\n#PBS -l walltime=72:00:00\n#PBS -l nodes=1:ppn=8\n#PBS -j n\n#PBS -m abe\n#PBS -M anellore@gmail.com\n$j'" >pbs_$i; let i=i+1; done