#!/usr/bin/env bash
# This script was run to generate all commands for tophat, subjunc, star, hisat, and rail
DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
python create_single_sample_sim_commands.py --aligners tophat subjunc star hisat rail --num-processes 8 --scratch /scratch2/langmead-fs1/tmp