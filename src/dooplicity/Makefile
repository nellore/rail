#!gmake

.PHONY: tests

tests:
	grep 'import unittest' *.py | sed 's/:.*//' | xargs -I % sh -c "echo %; python %;"
