#!/usr/bin/env bash
# This script was run to perform STAR and Subjunc simulations on a node with lots of memory
DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
for i in $(python create_single_sample_sim_commands.py --aligners subjunc star --num-processes 8 --scratch /scratch0/langmead-fs1/tmp); do eval $i; done