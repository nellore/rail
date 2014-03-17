#!/bin/sh

[ -z "$RAIL_HOME" ] && echo "Set RAIL_HOME" && exit 1
[ -z "$IGENOMES_HOME" ] && echo "Set IGENOMES_HOME" && exit 1

python absoluteize.py < dmel_flux.manifest > dmel_flux.abs.manifest

python $RAIL_HOME/src/driver/rail-rna.py \
	--local \
	--no-differential \
        --manifest dmel_flux.abs.manifest \
	--input preprocessed_reads \
	--intermediate intermediate \
	--output local_out \
        --reference $IGENOMES_HOME/Drosophila_melanogaster/UCSC/dm3 \
        --igenomes \
        --readlet-interval 3 \
        --num-processes 25 \
	--cap-search-window-size 0 \
        --keep-all \
	$*
