#First gather all of the tracks
IGV_TRACKS=igv_tracks
mkdir $IGV_TRACKS

IGENOME=/media/91356d9c-b219-4f38-a6a3-ff802758da9c/home/morton/Documents/Bioinformatics/JHU/Drosophila_melanogaster/UCSC/dm3
VALIDATE_OUT=/media/91356d9c-b219-4f38-a6a3-ff802758da9c/home/morton/Documents/Bioinformatics/JHU/tornado/tools/flux_sim/small_test/validation_output
#False positives track
FALSE_POSITIVES=$VALIDATE_OUT/false_positives
FALSE_NEGATIVES=$VALIDATE_OUT/false_negatives
ANNOTATED_SITES=$VALIDATE_OUT/annotated
DETECTED_SITES=$VALIDATE_OUT/../intermediate/splice_sites
GENES=$IGENOME/Annotation/Genes/genes.gtf
IGV_HOME=/media/91356d9c-b219-4f38-a6a3-ff802758da9c/home/morton/Documents/Bioinformatics/JHU/tornado/tools/igv/IGV

#Sort and index all of the files
igvtools sort $FALSE_POSITIVES.bed $FALSE_POSITIVES'_sort.bed'
igvtools sort $FALSE_NEGATIVES.bed $FALSE_NEGATIVES'_sort.bed'
igvtools sort $ANNOTATED_SITES.bed $ANNOTATED_SITES'_sort.bed'
igvtools sort $DETECTED_SITES.bed $DETECTED_SITES'_sort.bed'
igvtools index $FALSE_POSITIVES'_sort.bed'
igvtools index $FALSE_NEGATIVES'_sort.bed'
igvtools index $ANNOTATED_SITES'_sort.bed'
igvtools index $DETECTED_SITES'_sort.bed'

#Create batch script for the IGV to run
echo "new" > batch.txt
echo "load $FALSE_POSITIVES""_sort.bed">> batch.txt
echo "load $FALSE_NEGATIVES""_sort.bed">> batch.txt
echo "load $ANNOTATED_SITES""_sort.bed">> batch.txt
echo "load $DETECTED_SITES""_sort.bed">> batch.txt
echo "load $GENES">> batch.txt
sh $IGV_HOME/igv.sh