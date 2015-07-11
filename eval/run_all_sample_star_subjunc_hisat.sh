#!/usr/bin/env bash
# This script was run to generate all commands for hisat, star, and subjunc
DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
python create_single_sample_sim_commands.py --aligners subjunc star hisat --num-processes 8 --scratch /scratch2/langmead-fs1/tmp