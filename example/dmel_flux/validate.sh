UTIL=$TORNADO_HOME/src/util
CHECK=$TORNADO_HOME/src/check

BEDSITES=intermediate/splice_sites.bed
GENOME=$IGENOMES_HOME/$SPECIES/$DB/$ASM/Sequence/WholeGenomeFasta/genome.fa
ANNOTATIONS=$IGENOMES_HOME/$SPECIES/$DB/$ASM/Annotation/Genes/genes.gtf
echo *.bed
rm fly.sites
cat *.bed | python $UTIL/junc2site.py | sort -k2,2 | sort -k1,1 > fly.sites

SIM_READS=fly.sites
VALIDATE_OUT=validation_output
mkdir $VALIDATE_OUT

echo $ANNOTATIONS
python $CHECK/validate2.py \
        --refseq=$GENOME \
	--bed-file=$BEDSITES \
        --gtf=$ANNOTATIONS \
        --actual-sites=$SIM_READS \
        --out-dir=$VALIDATE_OUT \
	$*

VALIDATE_OUT=validation_output
#False positives track
FALSE_POSITIVES=$VALIDATE_OUT/false_positives
FALSE_NEGATIVES=$VALIDATE_OUT/false_negatives
ANNOTATED_SITES=$VALIDATE_OUT/annotated
DETECTED_SITES=$VALIDATE_OUT/../intermediate/splice_sites
EXONS=$VALIDATE_OUT/../intermediate/exon
GENES=$IGENOMES_DIR/$SPECIES/$DB/$ASM/Annotation/Genes/genes.gtf

#Sort and index all of the files
igvtools sort $FALSE_POSITIVES.bed $FALSE_POSITIVES'_sort.bed'
igvtools sort $FALSE_NEGATIVES.bed $FALSE_NEGATIVES'_sort.bed'
igvtools sort $ANNOTATED_SITES.bed $ANNOTATED_SITES'_sort.bed'
igvtools sort $DETECTED_SITES.bed $DETECTED_SITES'_sort.bed'
igvtools sort $EXONS.bed $EXONS'_sort.bed'
igvtools index $FALSE_POSITIVES'_sort.bed'
igvtools index $FALSE_NEGATIVES'_sort.bed'
igvtools index $ANNOTATED_SITES'_sort.bed'
igvtools index $DETECTED_SITES'_sort.bed'
igvtools index $EXONS'_sort.bed'

#Create batch script for the IGV to run
echo "new" > batch.txt
echo "load $FALSE_POSITIVES""_sort.bed">> batch.txt
echo "load $FALSE_NEGATIVES""_sort.bed">> batch.txt
echo "load $ANNOTATED_SITES""_sort.bed">> batch.txt
echo "load $DETECTED_SITES""_sort.bed">> batch.txt
echo "load $EXONS""_sort.bed">> batch.txt
echo "load $GENES">> batch.txt
IGV_HOME=$TORNADO_HOME/tools/igv/IGV
echo "prefix=`dirname $(readlink $0 || echo $0)`">$IGV_HOME/igv.sh
echo "exec java -Xmx2000m -Dproduction=true -Dapple.laf.useScreenMenuBar=true -Djava.net.preferIPv4Stack=true -jar "$IGV_HOME/$prefix"/igv.jar -b "$PWD/batch.txt" "$@" &">>$IGV_HOME/igv.sh

sh $IGV_HOME/igv.sh
