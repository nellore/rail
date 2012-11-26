URL_DIR=http://sammeth.net/artifactory/barna/barna/barna.simulator
VER=1.2
NORM_NM=flux-simulator
NM=$(NORM_NM)-$(VER)
PKG=$(NM).tgz
URL=$(URL_DIR)/$(VER)/$(PKG)

VERSION: $(NORM_NM)/bin/flux-simulator
	$< 2>/dev/null | tee $@

$(NORM_NM)/bin/flux-simulator: $(NORM_NM)
	touch $@

$(NORM_NM): $(NM)
	mv $< $@

$(NM): $(PKG)
	tar zxvf $<

$(PKG):
	wget $(URL)
