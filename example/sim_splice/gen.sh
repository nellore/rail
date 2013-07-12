
SIM=../../src/simulate/
ANNOTATIONS=../drosophila/Drosophila_melanogaster/UCSC/dm3/Annotation/Genes/genes.gtf
GENOME=../drosophila/Drosophila_melanogaster/UCSC/dm3/Sequence/WholeGenomeFasta/genome.fa

python $SIM/splice_sim.py --fasta=$GENOME --gtf=$ANNOTATIONS --output-prefix=fly --num-xscripts=2 --num-nucs=10000000