CHECK=/media/jamie/3cd05818-b444-4592-81e3-e1fa9121c2c9/Documents/summer2013/workspace/tornado/src/check
SIM_SPLICE=/media/jamie/3cd05818-b444-4592-81e3-e1fa9121c2c9/Documents/summer2013/workspace/tornado/example/sim_splice
XSCRIPTS=$SIM_SPLICE/fly.xscripts
BEDSITES=intermediate/splice_sites.bed
GENOME=/media/jamie/3cd05818-b444-4592-81e3-e1fa9121c2c9/Documents/summer2013/workspace/tornado/example/drosophila/Drosophila_melanogaster/UCSC/dm3/Sequence/WholeGenomeFasta/genome.fa
python $CHECK/validate.py --xscripts-file=$XSCRIPTS --bed-file=$BEDSITES --radius=10 --refseq=$GENOME
