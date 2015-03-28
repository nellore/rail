#!/usr/bin/env bash
# Uses awk to turn "perform" data output by spliced_read_recovery_performance.py
# -- as embedded in the script run_all_small_data_sims_locally -- into two matrices:
# 1) rows = true coverages of introns; columns = (true coverage, # relevant, # retrieved, # relevant and retrieved) instances;
#	here, an instance is an alignment overlapping introns whose minimally covered intron has coverage = row coverage
# 2) rows = retrieved coverages of introns; columns = (retrieved coverage # relevant, # retrieved, # relevant and retrieved) instances;
#	here, an instance is an alignment overlapping introns whose minimally covered intron has coverage = row coverage.
# May break if too much memory is used. If that happens, remove the &s below and run sequentially, or "wait" in chunks.

TRUETHRESHOLD=100 # max number of rows in matrix 1
RETRIEVEDTHRESHOLD=100 # max number of rows in matrix 2

INPUT=/scratch0/langmead-fs1/geuvadis_sims_for_paper/performances # where to find output of harvest_spliced_read_performance.py
OUTPUT=$INPUT/harvest

mkdir -p $OUTPUT

for perform in $INPUT/*perform
do
	awk '{ 
			for (i = 1; i <= '$TRUETHRESHOLD'; i++) {
				if ($5 == i) {
					rel[i] += $1;
					ret[i] += $2;
					if ($1 == $2) {
						relret[i] += $1
					}
				}
			}
		} END {
			for (i = 1; i <='$TRUETHRESHOLD'; i++) {
				printf "%d\t%d\t%d\t%d\n",i,rel[i],ret[i],relret[i]
			}
		}' $perform >$OUTPUT/$(basename "$perform")_true &
	awk '{ 
			for (i = 1; i <= '$RETRIEVEDTHRESHOLD'; i++) {
				if ($8 >= i) {
					rel[i] += $1;
					ret[i] += $2;
					if ($1 == $2) {
						relret[i] += $1
					}
				}
			}
		} END {
			for (i = 1; i <='$RETRIEVEDTHRESHOLD'; i++) {
				printf "%d\t%d\t%d\t%d\n",i,rel[i],ret[i],relret[i]
			}
		}' $perform >$OUTPUT/$(basename "$perform")_retrieved &
done
wait