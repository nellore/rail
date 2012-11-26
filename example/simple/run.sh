SCR_DIR=../../src/rnawesome

# Step 1: Readletize input reads and use Bowtie to align the readlets 
ALIGN_AGGR="cat"
ALIGN="python $SCR_DIR/align.py"

# Step 2: 
# (actually, should be merged with Step 1 in a Hadoop pipeline)
SPLICE_AGGR="cat"
SPLICE="python $SCR_DIR/splice.py"

# Step 3: Take all the intervals in a given genome window for all
#         samples and make a coverage table for the interval
MERGE_AGGR1="sort -n -k2,2"
MERGE_AGGR2="sort -s -k1,1"
MERGE="python $SCR_DIR/merge.py"

# Step 4: Walk over genome windows and...
WALK_AGGR="sort -k1,1"
WALK="python $SCR_DIR/walk.py"

# Step 5: Normalize
NORMALIZE_AGGR="sort -k1,1"
NORMALIZE="python $SCR_DIR/normalize.py"

cat *.tab \
	| $ALIGN_AGGR | $ALIGN \
		--bowtieArgs '-v 2 -m 1' \
		--bowtieExe $HOME/software/bowtie-0.12.8/bowtie \
		--bowtieIdx=$HOME/software/bowtie-0.12.8/indexes/e_coli \
		--readletLen 20 \
		--readletIval 2 \
		--manifest simple.manifest \
	| $SPLICE_AGGR | $SPLICE \
		--ntasks=10 \
		--genomeLen=1000 \
		--manifest simple.manifest \
	| $MERGE_AGGR1 | $MERGE_AGGR2 | $MERGE \
		--manifest simple.manifest \
	| $WALK_AGGR | $WALK \
		--manifest simple.manifest \
		--ntasks=10 \
		--genomeLen=1000 \
		--columnize \
	| $NORMALIZE_AGGR | $NORMALIZE \
		--percentile 0.75
