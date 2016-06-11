#!gmake

all: dooplicity-unit-tests rna-tests

.PHONY: dooplicity-unit-tests
dooplicity-unit-tests:
	make -C src/dooplicity tests

.PHONY: rna-unit-tests
rna-tests:
	make -C src/rna/steps tests
