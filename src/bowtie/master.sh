#!/bin/sh

set -e

sh hbb_build.sh bowtie2 rail_load_balancing ff9a75028d4c8e7ee5ba3b5b42ccc3aefa9acf90
nm=`ls -1 bowtie2-*.zip`
[ -z "${nm}" -o ! -f "${nm}" ] && echo "Failed" && exit 1
unzip ${nm}
nm=`echo ${nm} | sed 's/\.zip$//'`
nm_short=`echo ${nm} | sed 's/[a-zA-Z]*$//'`
mv "${nm}" "${nm_short}"
nm="${nm_short}"
rm -rf ${nm}/example ${nm}/doc ${nm}/scripts
rm -f ${nm}/*-debug
rm -f *.zip
zip -r "${nm}.zip" "${nm}"

sh hbb_build.sh bowtie rel1.2.1 f053cb6521c4020a736578430e53fb3b00990c10
nm="bowtie-bin.zip"
[ -z "${nm}" -o ! -f "${nm}" ] && echo "Failed" && exit 1
unzip "${nm}"
nm=`find . -name 'bowtie-*' -type d | sed 's/^\.\///'`
rm -rf ${nm}/SeqAn-1.1 ${nm}/doc ${nm}/scripts ${nm}/genomes ${nm}/indexes ${nm}/reads
rm -f ${nm}/*-debug
rm -f bowtie-bin.zip
zip -r "${nm}.zip" "${nm}"
