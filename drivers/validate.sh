CHECK=../src/check
SIM_SPLICE=../example/sim_splice
XSCRIPTS=$SIM_SPLICE/fly.xscripts
BEDSITES=$SIM_SPLICE/intermediate/splice_sites.bed
SITES=$SIM_SPLICE/fly.sites
GENOME=../example/drosophila/Drosophila_melanogaster/UCSC/dm3/Sequence/WholeGenomeFasta/genome.fa
FLANKS=$SIM_SPLICE/intermediate/flanks.tab
cat $SIM_SPLICE/intermediate/align_out.tsv | grep intron > $FLANKS   #really just the intron.py input

python $CHECK/validate.py --xscripts-file=$XSCRIPTS --bed-file=$BEDSITES --radius=10 --refseq=$GENOME --sites-file=$SITES --flank-seqs=$FLANKS
