CHECK=../src/check
SIM_SPLICE=../example/sim_splice
#SIM_SPLICE=../example/splice_short
XSCRIPTS=$SIM_SPLICE/fly.xscripts

BEDSITES=$SIM_SPLICE/intermediate/splice_sites.bed
TOPHAT_SITES=$SIM_SPLICE/tophat_out/tophat_sites.bed
SITES=$SIM_SPLICE/fly.sites
GENOME=../example/drosophila/Drosophila_melanogaster/UCSC/dm3/Sequence/WholeGenomeFasta/genome.fa
FLANKS=$SIM_SPLICE/intermediate/flanks.tab
COV=$SIM_SPLICE/fly.cov

cat $SIM_SPLICE/intermediate/align_out.tsv | grep intron > $FLANKS   #really just the intron.py input

# echo "Myrna Aligner validation"
sh run.sh
python $CHECK/validate.py --xscripts-file=$XSCRIPTS --bed-file=$BEDSITES --radius=10 --refseq=$GENOME --sites-file=$SITES --flank-seqs=$FLANKS --coverage-file=$COV --window-radius=100
#python $CHECK/validate.py --xscripts-file=$XSCRIPTS --bed-file=$BEDSITES --radius=10 --refseq=$GENOME --sites-file=$SITES --flank-seqs=$FLANKS --coverage-file=$COV --window-radius=100 --false-negatives  | sed -e '/False/,/Region/!d' | paste -d " " - - | uniq > myrna_out

# python $CHECK/validate.py --xscripts-file=$XSCRIPTS --bed-file=$BEDSITES --radius=10 --refseq=$GENOME --sites-file=$SITES --flank-seqs=$FLANKS --coverage-file=$COV --window-radius=100 --false-negatives  > myrna_false_negatives

# python $CHECK/validate.py \
#     --xscripts-file=$XSCRIPTS \
#     --bed-file=$BEDSITES \
#     --radius=10 \
#     --refseq=$GENOME \
#     --sites-file=$SITES \
#     --flank-seqs=$FLANKS \
#     --coverage-file=$COV \
#     --window-radius=100 \
#     --false-positives  \
#     | sed -e '/False/,/Region/!d' | paste -d " " - - | sort | uniq | cut -d' ' -f 4 \
#     | perl -ne 'm/(\S+):(\d+)-(\d+)/; print "(refID=='\''$1'\'' or ( abs(st-$2)<=radius and abs(en-$3)<=radius)) or\n" '

# echo "Tophat validation"
# python $CHECK/validate.py --xscripts-file=$XSCRIPTS --bed-file=$TOPHAT_SITES --radius=10 --refseq=$GENOME --sites-file=$SITES --flank-seqs=$FLANKS --coverage-file=$COV --window-radius=100 --false-negatives | sed -e '/False/,/Region/!d' | paste -d " " - - | uniq > tophat_out

#echo $MYRNA_OUT
#echo $TOPHAT_OUT
#diff -y <(echo $MYRNA_OUT) <(echo $TOPHAT_OUT)
#diff -y myrna_out tophat_out
