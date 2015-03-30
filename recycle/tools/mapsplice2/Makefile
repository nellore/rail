URL_DIR=http://protocols.netlab.uky.edu/~zeng
VER=v2.1.6
NORM_NM=MapSplice
NM=$(NORM_NM)-$(VER)
PKG=$(NM).zip
URL=$(URL_DIR)/$(PKG)

#VERSION:

$(NORM_NM)/install-root/bin/tophat: $(NORM_NM)/Makefile
	make -C $(NORM_NM)
	make -C $(NORM_NM) install

$(NORM_NM)/Makefile: $(NORM_NM) $(SAMTOOLS_DIR)/samtools/libbam.a
	mkdir -p $(NORM_NM)/libbam/lib $(NORM_NM)/libbam/include/bam
	cp $(SAMTOOLS_DIR)/samtools/libbam.a $(NORM_NM)/libbam/lib
	cp $(SAMTOOLS_DIR)/samtools/*.h $(NORM_NM)/libbam/include/bam
	cd $(NORM_NM) && ./configure --prefix=`pwd`/install-root --with-bam=`pwd`/libbam

$(NORM_NM): $(PKG)
	unzip $<
	mv $(NM) $@

$(PKG):
	wget $(URL)

.PHONY: clean
clean:
	rm -rf VERSION $(NM) $(PKG) $(NORM_NM)
