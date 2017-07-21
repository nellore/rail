#!gmake

.PHONY: tests

tests:
	grep 'import unittest' *.py  | grep -v 'hadoop_runner.py' | sed 's/:.*//' | xargs -I % sh -c "echo %; python %;"
	python hadoop_runner.py --test
