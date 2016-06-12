#!gmake

.PHONY: tests

tests:
	ls *.py | xargs -I % sh -c "echo %; python % --test;"
