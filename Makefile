#!gmake

.PHONY: cmodules
cmodules: src/alignment/_nw.so

src/alignment/_nw.so: src/alignment/nw_wrap.c
	cd src/alignment && python setup.py build_ext --inplace

src/alignment/nw_wrap.c: src/alignment/nw.i
	cd src/alignment && swig -python nw.i

.PHONY: clean
clean:
	rm -f src/alignment/nw_wrap.c \
	      src/alignment/nw.py \
	      src/alignment/_nw.so
	rm -rf src/alignment/build
