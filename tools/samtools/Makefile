NM=samtools
VER=0.1.19
AR=$(NM)-$(VER).tar.bz2
URL=http://sourceforge.net/projects/samtools/files/$(NM)/$(VER)/$(AR)/download

VERSION: samtools/samtools
	-$< 2> $@

samtools/samtools: samtools/Makefile
	cd samtools && $(MAKE)

samtools/Makefile: $(AR)
	tar xvfj $(AR)
	mv $(NM)-$(VER) samtools
	touch samtools/*

$(AR):
	wget -O $(AR) $(URL)
