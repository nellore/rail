#!gmake

all: dooplicity-unit-tests rna-unit-tests

.PHONY: dooplicity-unit-tests
dooplicity-unit-tests:
	make -C src/dooplicity tests

# NOTE: requires bowtie tools in PATH
.PHONY: rna-unit-tests
rna-unit-tests:
	make -C src/rna/steps tests
