URL_DIR=http://tophat.cbcb.umd.edu/downloads
VER=2.0.9
NORM_NM=tophat
NM=$(NORM_NM)-$(VER)
PKG=$(NM).tar.gz
URL=$(URL_DIR)/$(PKG)
SAMTOOLS_DIR=../samtools

VERSION: $(NORM_NM)/src/sam_juncs
	-$< 2>&1 | tee $@

$(NORM_NM)/src/sam_juncs: $(NORM_NM)/Makefile
	make -C $(NORM_NM)

$(NORM_NM)/Makefile: $(NORM_NM) $(SAMTOOLS_DIR)/samtools/libbam.a
	mkdir -p $(NORM_NM)/libbam/lib $(NORM_NM)/libbam/include/bam
	cp $(SAMTOOLS_DIR)/samtools/libbam.a $(NORM_NM)/libbam/lib
	cp $(SAMTOOLS_DIR)/samtools/*.h $(NORM_NM)/libbam/include/bam
	cd $(NORM_NM) && ./configure --with-bam=`pwd`/libbam

$(NORM_NM): $(PKG)
	tar zxvf $<
	mv $(NM) $@

$(PKG):
	wget $(URL)

.PHONY: clean
clean:
	rm -rf VERSION $(NM) $(PKG) $(NORM_NM)
