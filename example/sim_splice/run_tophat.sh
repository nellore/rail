BOWTIE_INDEX=../drosophila/Drosophila_melanogaster/UCSC/dm3/Sequence/BowtieIndex/genome
OUT_DIR=tophat_out
FASTA_SRC=/damsl/projects/myrna2/software/tornado/src/fasta
UTIL=/damsl/projects/myrna2/software/tornado/src/util
B=*.tab
for tab in $TAB
do
    file=`basename $tab`
    cat $tab | python $FASTA_SRC/tab2fastq.py > $file.fastq
done


tophat2  $BOWTIE_INDEX *.fastq --bowtie1 -o $OUT_DIR
cat $OUT_DIR/junctions.bed | python $UTIL/junc2site.py > $OUT_DIR/tophat_sites.bed
