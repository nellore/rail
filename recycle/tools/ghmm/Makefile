VER=0.9-rc1
NM=ghmm
NM_VER=ghmm-$(VER)
ARC=$(NM_VER).tar.gz
URL=http://sourceforge.net/projects/ghmm/files/ghmm/ghmm%20$(VER)/$(ARC)/download

# When configuring, use --with-rng=bsd to avoid warning messages at Python import time 

INSTALLED: $(NM)
	cd $(NM) && ./configure --with-rng=bsd
	sudo make -C $(NM) install
	touch $@

$(NM): $(ARC)
	tar zxvf $(ARC)
	mv $(NM_VER) $@

$(ARC):
	wget -O $@ $(URL)
