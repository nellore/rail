#!gmake

.PHONY: tests

tests:
	for i in `grep 'import unittest' *.py | sed 's/:.*//'` ; do python $$i ; done
