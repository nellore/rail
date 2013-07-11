
SIM=/damsl/projects/myrna2/software/tornado/src/simulate
ANNOTATIONS=/damsl/projects/myrna2/software/tornado/example/drosophila/Drosophila_melanogaster/UCSC/dm3/Annotation/Genes/genes.gtf
GENOME=/damsl/projects/myrna2/software/tornado/example/drosophila/Drosophila_melanogaster/UCSC/dm3/Sequence/WholeGenomeFasta/genome.fa

python $SIM/splice_sim.py --fasta=$GENOME --gtf=$ANNOTATIONS --output-prefix=fly --num-xscripts=2 --num-nucs=10000000