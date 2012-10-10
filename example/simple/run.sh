SCR_DIR=../../src/rnawesome

cat *.tab \
	| python $SCR_DIR/align.py \
		--bowtieArgs '-v 2 -m 1' \
		--bowtieExe $HOME/software/bowtie-0.12.8/bowtie \
		--bowtieIdx=$HOME/software/bowtie-0.12.8/indexes/e_coli \
		--readletLen 20 \
		--readletIval 2 \
		--manifest simple.manifest \
		| python $SCR_DIR/splice.py \
		--ntasks=10 \
		--genomeLen=1000 \
		--manifest simple.manifest \
	| sort -n -k2,2 | sort -s -k1,1 \
	| python $SCR_DIR/merge.py \
		--manifest simple.manifest \
	| sort -k1,1 \
	| python $SCR_DIR/walk.py \
		--manifest simple.manifest \
		--ntasks=10 \
		--genomeLen=1000
