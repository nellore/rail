CHECK=../src/check


BEDSITES=$SIM_HOME/intermediate/splice_sites.bed
TOPHAT_SITES=$SIM_HOME/tophat_out/tophat_sites.bed
#IGENOME=/media/91356d9c-b219-4f38-a6a3-ff802758da9c/home/morton/Documents/Bioinformatics/JHU/Drosophila_melanogaster/UCSC/dm3
GENOME=$IGENOMES_DIR/$SPECIES/$DB/$ASM/Sequence/WholeGenomeFasta/genome.fa
ANNOTATIONS=$IGENOMES_DIR/$SPECIES/$DB/$ASM/Annotation/Genes/genes.gtf
echo $SIM_HOME/*.lib
cat $SIM_HOME/*.lib > $SIM_HOME/fly.lib
LIB=$SIM_HOME/fly.lib
VALIDATE_OUT=$SIM_HOME/validation_output
mkdir $VALIDATE_OUT

# Don't re-run run.sh unless it seems to be needed
if [ ! -f $BEDSITES -o ! ] ; then
	sh run.sh
fi
echo $ANNOTATIONS
python $CHECK/validate2.py \
        --refseq=$GENOME \
	--bed-file=$BEDSITES \
        --gtf=$ANNOTATIONS \
        --lib-file=$LIB \
        --out-dir=$VALIDATE_OUT \
	$*
