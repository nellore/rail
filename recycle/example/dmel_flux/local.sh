#!/bin/sh

[ -z "$RAIL_HOME" ] && echo "Set RAIL_HOME" && exit 1
[ -z "$IGENOMES_HOME" ] && echo "Set IGENOMES_HOME" && exit 1

python absoluteize.py < dmel_flux.manifest > dmel_flux.abs.manifest

python $RAIL_HOME/src/driver/rail-rna.py \
	--local \
        --no-differential \
        --manifest dmel_flux.abs.manifest \
	--input preprocess  \
	--intermediate intermediate \
	--output local_out \
        --reference $IGENOMES_HOME/Drosophila_melanogaster/UCSC/dm3 \
        --igenomes \
        --min-cap-query-size 8 \
        --readlet-interval 4 \
        --readlet-length 25 \
        --num-processes 25 \
	--cap-search-window-size 1000 \
        --output-bam-by-chromosome \
        --keep-all \
	$*
