ifeq (,$(RAIL_HOME))
$(error "Set RAIL_HOME")
endif

ifeq (,$(IGENOMES_HOME))
$(error "Set IGENOMES_HOME")
endif

ifeq (,$(TOPHAT_HOME))
$(error "Set TOPHAT_HOME")
endif

FLUX_DIR=$(RAIL_HOME)/tools/flux_sim
IGV_DIR=$(RAIL_HOME)/tools/igv
OUTPUT_DIR=.
NREADS_PER_SAMP=1000000
NMOLECULES_PER_SAMP=500000
NSAMP=1
TOPHAT_EXE=$(TOPHAT_HOME)/tophat2

all: $(OUTPUT_DIR)/preprocessed_reads $(RAIL_HOME)/tools/igv/IGV/igv.jar

tophat_out/accepted_hits.bam: $(OUTPUT_DIR)/dmel_flux-0.fastq $(TOPHAT_EXE)
	$(TOPHAT_EXE) $(IGENOMES_HOME)/Drosophila_melanogaster/UCSC/dm3/Sequence/Bowtie2Index/genome dmel_flux-*.fastq

$(TOPHAT_EXE):
	make -C $(TOPHAT_DIR)

$(OUTPUT_DIR)/preprocessed_reads: $(OUTPUT_DIR)/dmel_flux.manifest $(OUTPUT_DIR)/dmel_flux-0.fastq
	mkdir -p $@
	cat $< | python $(RAIL_HOME)/src/rail-rna/preprocess.py --gzip-output --push=$@

# Run flux, many in parallel so we have essentially tech reps
$(OUTPUT_DIR)/dmel_flux-0.fastq: $(OUTPUT_DIR)/Chromosomes $(OUTPUT_DIR)/genes.fixed.gtf $(FLUX_DIR)/VERSION
	python $(RAIL_HOME)/src/simulate/flux.py \
		--name=dmel_flux \
		--output-dir=$(OUTPUT_DIR) \
		--flux-path=$(FLUX_DIR)/flux-simulator \
		--gtf-file=$(OUTPUT_DIR)/genes.fixed.gtf \
		--chromosomes=$(OUTPUT_DIR)/Chromosomes \
		--num-reads=$(NREADS_PER_SAMP) \
		--num-molecules=$(NMOLECULES_PER_SAMP) \
		--num-samples=$(NSAMP) \
		--num-processes=$(NSAMP)

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
	rm -f *.bed *.fastq *.lib *.pro *.par *.gtf Chromosomes chrom.sizes
	rm -rf preprocessed_reads intermediate local_out
