#!gmake

.PHONY: cmodules
cmodules: alignment/_nw.so

alignment/_nw.so: alignment/nw_wrap.c
	cd alignment && python setup.py build_ext --inplace

alignment/nw_wrap.c: alignment/nw.i
	cd alignment && swig -python nw.i

.PHONY: clean
clean:
	rm -f alignment/nw_wrap.c \
	      alignment/nw.py \
	      alignment/_nw.so
	rm -rf alignment/build
