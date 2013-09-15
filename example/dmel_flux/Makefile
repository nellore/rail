ifeq (,$(TORNADO_HOME))
$(error "Set TORNADO_HOME")
endif

ifeq (,$(IGENOMES_HOME))
$(error "Set IGENOMES_HOME")
endif

FLUX_DIR=$(TORNADO_HOME)/tools/flux_sim
IGV_DIR=$(TORNADO_HOME)/tools/igv
OUTPUT_DIR=.
NUM_READS=100000
NUM_MOLECULES=500000

all: $(OUTPUT_DIR)/preprocessed_reads $(TORNADO_HOME)/tools/igv/IGV/igv.jar

$(OUTPUT_DIR)/preprocessed_reads: $(OUTPUT_DIR)/dmel_flux.manifest $(OUTPUT_DIR)/dmel_flux-0.fastq
	mkdir -p $@
	cat $< | python $(TORNADO_HOME)/src/rnawesome/preprocess.py --gzip-output --push=$@

# Run flux, many in parallel so we have essentially tech reps
$(OUTPUT_DIR)/dmel_flux-0.fastq: $(OUTPUT_DIR)/Chromosomes $(OUTPUT_DIR)/genes.fixed.gtf $(FLUX_DIR)/VERSION
	python $(TORNADO_HOME)/src/simulate/flux.py \
		--name=dmel_flux \
		--output-dir=$(OUTPUT_DIR) \
		--flux-path=$(FLUX_DIR)/flux-simulator \
		--gtf-file=$(OUTPUT_DIR)/genes.fixed.gtf \
		--chromosomes=$(OUTPUT_DIR)/Chromosomes \
		--num-reads=$(NUM_READS) \
		--num-molecules=$(NUM_MOLECULES) \
		--num-samples=6 \
		--num-processes=6

# Make the link to the directory with 1 FASTA file per chromosome
$(OUTPUT_DIR)/Chromosomes: $(IGENOMES_HOME)/Drosophila_melanogaster
	cd $(OUTPUT_DIR) && ln -s -f $(IGENOMES_HOME)/Drosophila_melanogaster/UCSC/dm3/Sequence/Chromosomes Chromosomes

# Make sure flux fixed gtf file is in place
$(OUTPUT_DIR)/genes.fixed.gtf: $(IGENOMES_HOME)/Drosophila_melanogaster $(FLUX_DIR)/VERSION
	python $(FLUX_DIR)/gtf_flux_fix.py \
		--fasta $(IGENOMES_HOME)/Drosophila_melanogaster/UCSC/dm3/Sequence/Chromosomes/*.fa \
		--gtf $(IGENOMES_HOME)/Drosophila_melanogaster/UCSC/dm3/Annotation/Genes/genes.gtf \
		> $(OUTPUT_DIR)/genes.unsorted.fixed.gtf
	$(FLUX_DIR)/flux-simulator/bin/flux-simulator \
		-t sortGTF \
		-i $(OUTPUT_DIR)/genes.unsorted.fixed.gtf \
		-o $(OUTPUT_DIR)/genes.fixed.gtf

# Make sure flux is installed
$(FLUX_DIR)/VERSION:
	$(MAKE) -C $(FLUX_DIR)

# Make sure IGV is installed
$(IGV_DIR)/IGV/igv.jar:
	cd $(IGV_DIR) && sh download.sh

# Make sure D. melanogaster archive from iGenomes is present
IG_URL=ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com
$(IGENOMES_HOME)/Drosophila_melanogaster:
	mkdir -p $(IGENOMES_HOME)
	cd $(IGENOMES_HOME) && wget $(IG_URL)/Drosophila_melanogaster/UCSC/dm3/Drosophila_melanogaster_UCSC_dm3.tar.gz
	cd $(IGENOMES_HOME) && tar zxvf Drosophila_melanogaster_UCSC_dm3.tar.gz

.PHONY: clean
clean:
	rm -f *.bed *.fastq *.lib *.pro *.par *.gtf Chromosomes *.manifest chrom.sizes
	rm -rf preprocessed_reads intermediate local_out
