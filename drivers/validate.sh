CHECK=../src/check
SIM_SPLICE=../example/sim_splice
XSCRIPTS=$SIM_SPLICE/fly.xscripts
BEDSITES=$SIM_SPLICE/intermediate/splice_sites.bed
SITES=$SIM_SPLICE/fly.sites
GENOME=../example/drosophila/Drosophila_melanogaster/UCSC/dm3/Sequence/WholeGenomeFasta/genome.fa
python $CHECK/validate.py --xscripts-file=$XSCRIPTS --bed-file=$BEDSITES --radius=10 --refseq=$GENOME --sites-file=$SITES
