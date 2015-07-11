#!/usr/bin/env bash
# Script run on Johns Hopkins cluster to generate PBS scripts for running simulations
# May have to be adjusted by the user to suit other clusters
DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
cd $DIR
mkdir -p pbs_scripts
cd pbs_scripts
for i in $(python ../create_single_sample_sim_commands.py); do cat <(echo -e "
#!/bin/tcsh
#PBS -q batch
#PBS -l walltime=48:00
#PBS -l nodes=1:ppn=8
#PBS -j n
") <(echo $i); done