
SIM=../../src/simulate
ANNOTATIONS=../drosophila/Drosophila_melanogaster/UCSC/dm3/Annotation/Genes/genes.gtf
GENOME=../drosophila/Drosophila_melanogaster/UCSC/dm3/Sequence/WholeGenomeFasta/genome.fa
DRIVERS=$PWD/../../drivers
FASTAIDX=../drosophila/Drosophila_melanogaster/UCSC/dm3/Sequence/WholeGenomeFasta/genome.fa.fai
ANNOT_SITES=fly.genes
python $SIM/splice_sim.py --fasta=$GENOME --gtf=$ANNOTATIONS --output-prefix=$PWD/fly --num-xscripts=10 --num-nucs=1000000 --alternative-spliced=0 --stranded=1 --readmm_rate=0 --snp_rate=0 --indel_rate=0 --variants_file=$DRIVERS/variants.txt --chrsizes=$FASTAIDX --output-annotations=$ANNOT_SITES

rm fly.genes.gtf
#IFS=$'\n';
for line in `cat $ANNOT_SITES`
do
    echo $line
    cat $ANNOTATIONS | grep $line | grep exon >> fly.genes.gtf
done