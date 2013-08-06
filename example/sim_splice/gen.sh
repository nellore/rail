
SIM=../../src/simulate
ANNOTATIONS=../drosophila/Drosophila_melanogaster/UCSC/dm3/Annotation/Genes/genes.gtf
GENOME=../drosophila/Drosophila_melanogaster/UCSC/dm3/Sequence/WholeGenomeFasta/genome.fa
DRIVERS=/media/jamie/3cd05818-b444-4592-81e3-e1fa9121c2c9/Documents/summer2013/workspace/tornado/drivers
FASTAIDX=../drosophila/Drosophila_melanogaster/UCSC/dm3/Sequence/WholeGenomeFasta/genome.fa.fai
ANNOT_SITES=fly.genes.gtg
python $SIM/splice_sim.py --fasta=$GENOME --gtf=$ANNOTATIONS --output-prefix=$PWD/fly --num-xscripts=100 --num-nucs=10000000 --alternative-spliced=0 --stranded=1 --readmm_rate=0 --snp_rate=0 --indel_rate=0 --variants_file=$DRIVERS/variants.txt --chrsizes=$FASTAIDX --output-annotations=$ANNOT_SITES