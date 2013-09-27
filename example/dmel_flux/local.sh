#!/bin/sh

[ -z "$TORNADO_HOME" ] && echo "Set TORNADO_HOME" && exit 1
[ -z "$IGENOMES_HOME" ] && echo "Set IGENOMES_HOME" && exit 1

python $TORNADO_HOME/src/driver/tornado.py \
	--local \
	--start-with-align --no-differential \
	--manifest dmel_flux.manifest \
	--input preprocessed_reads \
	--intermediate intermediate \
	--output local_out \
	--reference $IGENOMES_HOME/Drosophila_melanogaster/UCSC/dm3 \
	--igenomes \
	--num-processes=6 \
	--keep-all \
	$*
