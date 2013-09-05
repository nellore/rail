#!/bin/bash
FLUX=flux-simulator

#./flux-simulator/bin/flux-simulator -p parameters/sailfish.par
SIM_OUT=simulation_output
python run_sim.py \
    --output-dir=$SIM_OUT \
    --flux-path=$FLUX \
    --gtf-file=$PWD/genes.fixed.gtf \
    --chromosomes=$PWD/Chromosomes \
    --num-reads=15000000 \
    --num-molecules=5000000 \
    --pcr-probability=1 \
    --num-samples=2
    
