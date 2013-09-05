#!/bin/bash
FLUX=flux_sim/flux-simulator

#./flux-simulator/bin/flux-simulator -p parameters/sailfish.par
SIM_OUT=simulation_output
python run_sim.py \
    --output-dir=$SIM_OUT \
    --flux-path=$FLUX \
    --gtf-file=$PWD/genes.fixed.gtf \
    --chromosomes=Chromosomes \
    --num-samples=1

