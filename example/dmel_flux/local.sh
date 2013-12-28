#!/bin/sh

[ -z "$RAIL_HOME" ] && echo "Set RAIL_HOME" && exit 1
[ -z "$IGENOMES_HOME" ] && echo "Set IGENOMES_HOME" && exit 1

python $RAIL_HOME/src/driver/rail-rna.py \
    --local \
    --start-with-align --no-differential \
    --manifest dmel_flux.manifest \
    --input preprocessed_reads \
    --intermediate intermediate \
    --output local_out \
    --reference $IGENOMES_HOME/Drosophila_melanogaster/UCSC/dm3 \
    --igenomes \
    --stop-after-align \
    --num-processes=6 \
    --keep-all \
    $*
