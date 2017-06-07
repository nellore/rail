#!/bin/sh

[ "$1" = "bowtie" -o "$1" = "bowtie2" ] || { echo "Must specify bowtie or bowtie2 as 1st arg" && exit 1; }
[ -z "$2" ] && echo "Must specify branch as 2nd arg" && exit 1

docker run -t -i --rm -v `pwd`:/io \
    phusion/holy-build-box-64:latest \
    /hbb_exe_gc_hardened/activate-exec bash \
    /io/static_tbb.bash $1 $2
