#       AUTHOR: Martin Bens, bensmartin@gmail.com
# ORGANIZATION: FLI Jena

FRAMA_VERSION := 20160626_development

SHELL := /bin/bash

.SECONDEXPANSION:

.DELETE_ON_ERROR:

.PHONY: all check assembly primary-processing blast blast-summary annotation \
	vector fasta-conversion catalog-summary functional

# IMPORT PARAMETERS {{{

ifdef FILE
include $(FILE)
else
$(error Use "-FILE" parameter to specify configuration file)
endif

ifndef OUTPUT_DIR
$(error OUTPUT_DIR not defined!)
endif

# }}}
# VARIABLES & FUNCTIONS {{{
# -----------------------------------------------------------------------------

# warning: does not follow symbolic links but directory doesn't have to exist
OUTPUT_DIR := $(abspath $(OUTPUT_DIR))

REF_TRANSCRIPTOME := $(realpath $(REF_TRANSCRIPTOME))

# directories
LOG_DIR      := $(OUTPUT_DIR)/logs
TABLE_DIR    := $(OUTPUT_DIR)/tables
TRANS_DIR    := $(OUTPUT_DIR)/transcripts
TRINITY_DIR  := $(OUTPUT_DIR)/trinity
BLAST_DIR    := $(OUTPUT_DIR)/blast
BBH_DIR      := $(OUTPUT_DIR)/bbh
ASSEMBLY_DIR := $(OUTPUT_DIR)/assembly
REF_DIR      := $(OUTPUT_DIR)/reference

# assembly
RAW_ASSEMBLY    := $(ASSEMBLY_DIR)/Trinity
PRIPRO_ASSEMBLY := $(ASSEMBLY_DIR)/Trinity_preprocessed
CDHIT_ASSEMBLY  := $(ASSEMBLY_DIR)/cdhit/Trinity.cdhit
TGICL_ASSEMBLY  := $(ASSEMBLY_DIR)/tgicl/Trinity.tgicl

# reference
MASKED_REF     := $(REF_DIR)/ref_transcriptome.masked.fa
MASKED_REF_CDS := $(REF_DIR)/ref_transcriptome_cds.masked.fa

# final output
FINAL_RAW_ASSEMBLY  := $(OUTPUT_DIR)/assembly
FINAL_ASSEMBLY      := $(OUTPUT_DIR)/assembly_pripro

FINAL_CDS      := $(OUTPUT_DIR)/transcriptome_CDS
FINAL_MRNA     := $(OUTPUT_DIR)/transcriptome_mRNA
FINAL_GENBANK  := $(OUTPUT_DIR)/transcriptome

FINAL_SUMMARY := $(OUTPUT_DIR)/annotation
FINAL_KEGG    := $(OUTPUT_DIR)/annotation_KEGG
FINAL_GO      := $(OUTPUT_DIR)/annotation_GO


# we need space
NOTHING :=
space := $(NOTHING) $(NOTHING)
define newline


endef

# stdout
DATE    = `date +'%y-%m-%d %H:%M:%S'`
P_INFO  = $(DATE) [INF]
P_LOG   = $(DATE) [LOG]    $(LOG_DIR)
P_OUT   = $(DATE) [OUT]   
P_ERR   = $(DATE) [ERR]

# funcitons
NDEF_ANY = $(filter undefined,$(foreach v,$(1),$(origin $(v))))
DEF_ANY  = $(filter-out undefined,$(foreach v,$(1),$(origin $(v))))
FILTER   = $(foreach v,$(2),$(if $(findstring $(1),$(v)),$(v),))

# extend path {{{

# get PATH_ variables
MYPATHS := $(call FILTER,PATH_,$(.VARIABLES))
MYPATHS := $(strip $(MYPATHS))

# appends values to PATH
EPATHS := $(foreach v,$(MYPATHS),$(realpath $($(v))))
EPATHS := $(strip $(EPATHS))
EPATHS := $(subst $(space),:,$(EPATHS)):$(PATH)
export PATH := $(EPATHS)

# }}}

# commands
COMMAND_SUM_BLAST :=
ifeq ($(NCBI_BLAST),1)
COMMAND_SUM_BLAST := perl genblasta2table.pl
else
COMMAND_SUM_BLAST := perl wublast_average.pl
endif


COMMAND_TRINITY := $(shell export PATH="$(EPATHS)" && { command -v Trinity || command -v Trinity.pl; })
COMMAND_CDHIT   := $(shell export PATH="$(EPATHS)" && { command -v cd-hit-est || command -v cdhit-est; })

# determine trinity version
PERL := perl -ne 'chomp; s/.+[vr](.+)$$/$$1/g; s/[_.]//g; print length'
TRINITY_VERSION := $(shell $(COMMAND_TRINITY) --version | grep version | $(PERL))

# }}}
# DEFAULT PARAMETER {{{
# -----------------------------------------------------------------------------
OPT_R := --slave --vanilla

TMP_DIR ?= $(OUTPUT_DIR)/tmp

.INTERMEDIATE: $(wildcard $(TMP_DIR)/*)

OPT_TRINITY                  += --no_cleanup --CPU $(OPT_CPUS)
OPT_REPEAT_REF_TRANSCRIPTOME += -xsmall -par $(OPT_CPUS)
OPT_REPEAT_ASSEMBLY          += -xsmall -par $(OPT_CPUS)
OPT_TGICL                    += -c $(OPT_CPUS)
OPT_CD_HIT_EST               += -T $(OPT_CPUS)

# default: RepeatMasker {{{
NCBI_BLAST ?= 1

#ifeq ($(NCBI_BLAST),1)
#OPT_REPEAT_REF_TRANSCRIPTOME += -engine ncbi
#OPT_REPEAT_ASSEMBLY          += -engine ncbi
#else
#OPT_REPEAT_REF_TRANSCRIPTOME += -engine wublast
#OPT_REPEAT_ASSEMBLY          += -engine wublast
#endif
# }}}
# default: Blast {{{

ENGINE := sge
ifndef OPT_SGE
ENGINE  := parallel
OPT_SGE := $(OPT_CPUS)
endif

BLAST_PARAMETER := -j $(OPT_SGE) -format tabular -engine $(ENGINE)

ifeq ($(NCBI_BLAST),1)
OPT_BLAST        := $(shell perl parameter_blast.pl "$(OPT_GENBLAST)" "$(OPT_BLAST)")
BLAST_PARAMETER  += -path bin=$(dir $(shell export PATH="$(EPATHS)" && command -v genblasta))/ -flavor genblasta -parameter '$(OPT_BLAST)'
else
OPT_BLAST        += -wordmask=seg lcmask -topcomboN 3 -cpus 1
BLAST_PARAMETER  += -path bin=$(dir $(shell export PATH="$(EPATHS)" && command -v blastn))/ -flavor wu -parameter '$(OPT_BLAST)'
endif
# }}}
# default: annotation.pl {{{

PATH_GENSCAN_MAT := $(abspath $(PATH_GENSCAN_MAT))
ANNOTATION_PARAMETER := $(OPT_PREDICTCDS) -genscan-matrix $(PATH_GENSCAN_MAT) -cpus $(OPT_CPUS)

# OPTIONAL: homolog groups
ifneq ($(call DEF_ANY, ORTHOLOG_TABLE ORTHOLOG_CDS),)
ORTHOLOG_TABLE := $(realpath $(ORTHOLOG_TABLE))
ORTHOLOG_CDS   := $(realpath $(ORTHOLOG_CDS))
ANNOTATION_PARAMETER += -ortholog-table $(ORTHOLOG_TABLE) -ortholog-cds $(ORTHOLOG_CDS)
endif

# OPTIONAL: need former contig name of cluster
ifneq (,$(findstring tgicl,$(ASSEMBLY_PREPROCESS)))
ANNOTATION_PARAMETER += -tgicl $(TGICL_ASSEMBLY).fa.clstr
endif

# }}}
# default: SBH {{{

OPT_SBH ?= -identity 99.9 -coverage 99.9

# adds '=' if missing
$(eval OPT_SBH := $(shell echo $(OPT_SBH) | perl -pe 's/(-.+?)[\s=](.+?)/$$1=$$2/g'))
# }}}
# primary processing {{{
ifdef ASSEMBLY_PREPROCESS
ASSEMBLY_PREPROCESS := $(addprefix complete.,$(ASSEMBLY_PREPROCESS))
ASSEMBLY_PREPROCESS := $(addprefix $(ASSEMBLY_DIR)/,$(ASSEMBLY_PREPROCESS))
endif
# }}}
# sort order {{{
# evalue (asc) > highscore (desc) > target coverage (desc) > total score (desc) > query coverage (desc) > query name (asc) > target name (asc)
SORT_ORDER=-k4,4g -k6,6nr -k16,16nr -k5,5nr -k15,15nr
# }}}
# GENERIC RULES {{{
# -----------------------------------------------------------------------------

%/.c :
	mkdir -p $(@D)
	touch $@

%.fa.fai: %.fa
	samtools faidx $<

# }}}
# }}}
# READ DATA (FASTQ and FASTA only) {{{
# -----------------------------------------------------------------------------

READ_DIR := $(realpath $(READ_DIR))

#  Trinity accepts only one specific format at once
READ_EXTENSIONS := *.fq *.fastq *.fa
READ_EXTENSIONS := $(addprefix $(READ_DIR)/,$(READ_EXTENSIONS))
READ_EXTENSIONS := $(READ_EXTENSIONS) $(addsuffix .gz,$(READ_EXTENSIONS))
READS := $(wildcard $(READ_EXTENSIONS))

# creates symbolic links to read data
READ_LINKS := $(subst $(READ_DIR)/,$(OUTPUT_DIR)/reads/,$(READS))
$(READ_LINKS) : $$(subst $(OUTPUT_DIR)/reads/,$(READ_DIR)/,$$@) | $(OUTPUT_DIR)/reads/.c
	ln -f -s $< $@

# }}}
# PIPELINE {{{
# -----------------------------------------------------------------------------

all: info check assembly primary-processing blast blast-summary annotation vector \
	fasta-conversion catalog-summary functional
	@printf "$(P_INFO) Pipeline complete.\n"
	@printf "\n"

info: | $(OUTPUT_DIR)/.c
	@printf "\n"
	@printf "$(P_INFO) Starting FRAMA...\n"
	@printf "$(P_OUT) $(OUTPUT_DIR)/parameters.cfg [parameter file]\n"
	@printf "$(P_OUT) $(OUTPUT_DIR) [output directory]\n\n"
	@printf "#\n# FRAMA Version: $(FRAMA_VERSION)\n#\n" > $(OUTPUT_DIR)/parameters.cfg
	@cat $(FILE) >> $(OUTPUT_DIR)/parameters.cfg

check: $(FILE) | $(LOG_DIR)/.c $(TABLE_DIR)/.c $(TRANS_DIR)/.c \
	$(TRINITY_DIR)/.c $(BLAST_DIR)/.c $(BBH_DIR)/.c $(ASSEMBLY_DIR)/.c \
	$(REF_DIR)/.c 
	@printf "$(P_INFO) Checking required software and modules\n"
	@printf "$(P_LOG)/software_check.log\n"
	perl init.pl $< > $(LOG_DIR)/software_check.log
	@printf "\n"

assembly: $(RAW_ASSEMBLY).fa
	@printf "$(P_INFO) Assembly complete.\n"

primary-processing: $(PRIPRO_ASSEMBLY).fa
	@printf "$(P_INFO) Primary processing complete.\n\n"

blast : $(BLAST_DIR)/raw_CDS_trin2ref.csv $(BLAST_DIR)/raw_CDS_ref2trin.csv \
	$(BLAST_DIR)/raw_mRNA_ref2trin.csv $(BLAST_DIR)/raw_mRNA_trin2ref.csv
	@printf "$(P_INFO) BLAST complete.\n\n"

blast-summary : $(BLAST_DIR)/avg_CDS_trin2ref.csv $(BLAST_DIR)/avg_CDS_ref2trin.csv \
	$(BLAST_DIR)/avg_mRNA_ref2trin.csv $(BLAST_DIR)/avg_mRNA_trin2ref.csv
	@printf "$(P_INFO) Blast preprocessing complete.\n\n"

annotation : $(ASSEMBLY_DIR)/annotated.gbk
	@printf "$(P_INFO) Annotation complete\n\n"

vector: $(FINAL_GENBANK).gbk
	@printf "$(P_INFO) Contamination check compete.\n"
	@printf "$(P_OUT) $<\n\n"

fasta-conversion: $(FINAL_CDS).fa $(FINAL_MRNA).fa
	@printf "$(P_INFO) Output files: Transcript Catalog.\n"
	@printf "$(P_OUT) $(word 1,$^)\n"
	@printf "$(P_OUT) $(word 2,$^)\n\n"

catalog-summary : $(FINAL_SUMMARY).pdf
	@printf "$(P_INFO) Statistics complete.\n\n"

functional : $(FINAL_GO).pdf $(FINAL_KEGG).pdf
	@printf "$(P_INFO) Functional annotation.\n\n"

# }}}
# MASKING {{{
# -----------------------------------------------------------------------------

# - mask repeats only if path to RepeatMasker is available (create symlinks otherwise)
# - performance boost: concatenating sequences for RepeatMasker
# - repeat library specified in configuration file

%.long.fa : %.fa
	perl fasta_concatenate.pl -file $< -out $(basename $<).long

%.masked.fa : %.long.fa.masked
	perl fasta_unconcat.pl -file $< -index $(basename $(basename $<)).index -o $@

# Repeat mask assembly
$(PRIPRO_ASSEMBLY).long.fa.masked : $(PRIPRO_ASSEMBLY).long.fa
	@printf "$(P_INFO) Repeat masking assembly.\n"
ifeq ($(REPEAT),0)
	ln -f -s $< $@
else
	@printf "$(P_LOG)/rm_assembly.log\n\n"
	cd $(ASSEMBLY_DIR) && RepeatMasker $< $(OPT_REPEAT_ASSEMBLY) &> $(LOG_DIR)/rm_assembly.log
endif

# Repeat mask reference
$(REF_DIR)/ref_transcriptome.long.fa.masked : $(REF_DIR)/ref_transcriptome.long.fa
	@printf "$(P_INFO) Repeat masking reference.\n"
ifdef REF_TRANSCRIPTOME_FASTA_MASKED
	ln -f -s $(REF_TRANSCRIPTOME_FASTA_MASKED) $@
else ifeq ($(REPEAT),0)
	ln -f -s $< $@
else
	@printf "$(P_LOG)/repeat_masking_reference.log\n\n"
	cd $(REF_DIR) && RepeatMasker $< $(OPT_REPEAT_REF_TRANSCRIPTOME) &> $(LOG_DIR)/repeat_masking_reference.log
endif

# Repeat mask reference CDS
$(MASKED_REF_CDS) : $(REF_TRANSCRIPTOME) $(MASKED_REF)
ifndef REF_TRANSCRIPTOME_FASTA_CDS_MASKED
	perl genbank_CDS.pl -genbank $< -input $(word 2,$^) -output $@
else
	ln -f -s $(REF_TRANSCRIPTOME_FASTA_CDS_MASKED) $@
endif

# }}}
# Genbank to fasta {{{
$(REF_DIR)/ref_transcriptome.fa : $(REF_TRANSCRIPTOME)
ifndef REF_TRANSCRIPTOME_FASTA
	perl genbank2fasta.pl -i $< -o $@ -feature CDS
else
	ln -f -s $(REF_TRANSCRIPTOME_FASTA) $@
endif

$(REF_DIR)/ref_transcriptome_cds.fa : $(REF_TRANSCRIPTOME)
ifndef REF_TRANSCRIPTOME_FASTA_CDS
	perl genbank2fasta.pl -i $< -o $@ -feature CDS -feature-sequence
else
	ln -f -s $(REF_TRANSCRIPTOME_FASTA_CDS) $@
endif
# }}}
# TRANSCRIPTOME ASSEMBLY {{{
# -----------------------------------------------------------------------------

# - Trinity handles precomputed results by itself
# - JVM decides how many garbage collection threads to use => we limit
#   GC-threads to one and start half of allowed process in parallel

$(TRINITY_DIR)/Trinity.fasta : $(READ_LINKS)
	@printf "$(P_INFO) Starting Trinity.\n"

	$(eval options := --output $(dir $@))

ifeq ($(OPT_READTYPE),pe)
	$(eval read_R1 := $(strip $(FILTER $^,"_R1")))
	$(eval read_R2 := $(strip $(FILTER $^,"_R2")))
	$(eval options += --left $(read_R1) --right $(read_R2))
else
	$(eval options += --single $^)
endif

ifneq (,$(findstring Trinity.pl,$(COMMAND_TRINITY)))
	@printf "$(P_LOG)/inchworm.log\n"
	ulimit -s unlimited; $(COMMAND_TRINITY) $(options) --no_run_chrysalis $(OPT_TRINITY) &> $(LOG_DIR)/inchworm.log

	@printf "$(P_LOG)/chrysalis.log\n"
	ulimit -s unlimited; $(COMMAND_TRINITY) $(options) --no_run_quantifygraph $(OPT_TRINITY) $(OPT_BUTTERFLY) &> $(LOG_DIR)/chrysalis.log

	@printf "$(P_LOG)/quanitfyGraph.log\n"
	ulimit -s unlimited; cat $(dir $@)/chrysalis/quantifyGraph_commands | parallel --max-procs $(OPT_CPUS) --joblog $(LOG_DIR)/quantifyGraph.log

	Rscript parallel_checkLog.R $(LOG_DIR)/quantifyGraph.log

	@printf "$(P_LOG)/butterfly.log\n"
	ulimit -s unlimited; cat $(dir $@)/chrysalis/butterfly_commands | parallel --max-procs $(shell echo '$(OPT_CPUS)/2' | bc) --joblog $(LOG_DIR)/butterfly.log

	Rscript parallel_checkLog.R $(LOG_DIR)/butterfly.log
	ulimit -s unlimited; find $(dir $@)/chrysalis -name "*allProbPaths.fasta" -exec cat {} + > $@
else
	@printf "$(P_LOG)/Trinity.log\n\n"
	ulimit -s unlimited; $(COMMAND_TRINITY) $(options) $(OPT_TRINITY) > $(LOG_DIR)/Trinity.log
endif

ifneq (,$(findstring genome_guided_bam,$(OPT_TRINITY)))
	mv $(TRINITY_DIR)/Trinity-GG.fasta $@
endif

$(RAW_ASSEMBLY).fa : $(TRINITY_DIR)/Trinity.fasta | $(ASSEMBLY_DIR)/.c
ifeq (,$(findstring Trinity.pl,$(COMMAND_TRINITY)))
	perl -pe 's/^>(.+?)\|(.+?)\s.+/>$$2.$$1/' $< > $@
else
	cp $< $@
endif

# }}}
# PRIMARY PROCESSING {{{
# -----------------------------------------------------------------------------

TMP_ASSEMBLY := $(TMP_DIR)/Trinity_tmp

$(FINAL_ASSEMBLY).fa : $(PRIPRO_ASSEMBLY).masked.fa
	cp $< $@

$(ASSEMBLY_DIR)/complete.primary_processing : $(ASSEMBLY_DIR)/complete.tmp $(ASSEMBLY_PREPROCESS) $(PRIPRO_ASSEMBLY).fa
	touch $@

$(PRIPRO_ASSEMBLY).fa: $(RAW_ASSEMBLY).fa $(ASSEMBLY_DIR)/complete.tmp $(ASSEMBLY_PREPROCESS)
	if [ -a $@.index ]; then rm $@.index; fi
	cp $(TMP_ASSEMBLY).fa $@

# Create temporary file which all preprocessing steps will use
$(ASSEMBLY_DIR)/complete.tmp : $(RAW_ASSEMBLY).fa | $(TMP_DIR)/.c
	cp $(RAW_ASSEMBLY).fa $(TMP_ASSEMBLY).fa
	touch $@

# CD-HIT {{{

$(ASSEMBLY_DIR)/complete.cd-hit : $(CDHIT_ASSEMBLY).fa
	@printf "$(P_INFO) Clustering sequences complete.\n\n"
	cp $< $(TMP_ASSEMBLY).fa
	touch $@

$(CDHIT_ASSEMBLY).fa : $(ASSEMBLY_DIR)/complete.tmp | $(LOG_DIR)/.c $(dir $(CDHIT_ASSEMBLY))/.c
	@printf "\n"
	@printf "$(P_INFO) Clustering sequences: CDHIT-EST\n"
	@printf "$(P_LOG)/cdhit.log\n\n"
	$(COMMAND_CDHIT) -i $(TMP_ASSEMBLY).fa -o $@ $(OPT_CD_HIT_EST) > $(LOG_DIR)/cdhit.log

# }}}
# TGICL {{{

$(ASSEMBLY_DIR)/complete.tgicl : $(TGICL_ASSEMBLY).fa
	@printf "$(P_INFO) Meta assembly complete.\n\n"
	cp $< $(TMP_ASSEMBLY).fa
	touch $@

$(TGICL_ASSEMBLY).fa : | $(TMP_ASSEMBLY).fa $(dir $(TGICL_ASSEMBLY))/.c
	@printf "\n";
	@printf "$(P_INFO) Meta assembly: TGICL\n"
	@printf "$(P_LOG)/tgicl.log\n"

	$(eval TMPFILE := $(ASSEMBLY_DIR)/tgicl/Trinity.tmp.fa)
	cp $(TMP_ASSEMBLY).fa $(ASSEMBLY_DIR)/tgicl/tmp
	cd $(ASSEMBLY_DIR)/tgicl; tgicl tmp &> $(LOG_DIR)/tgicl.log
	#$(OPT_TGICL);
	bash tgicl.sh $(ASSEMBLY_DIR)/tgicl/tmp $@

# }}}
# }}}
# GENE SYMBOL ASSIGNMENT {{{
# -----------------------------------------------------------------------------

$(TMP_DIR)/best_hit_fusion.csv : $(TMP_DIR)/bbh_sbh.csv $(BBH_DIR)/fusion.csv
	cut -f1 $(word 2,$^) \
		| grep -F -f - -w $(word 1,$^) \
		| perl -ane '$$F[4] .= ",part0"; print join "\t", @F; print "\n"' \
		> $@

$(TMP_DIR)/xx.csv : $(TMP_DIR)/bbh_sbh.csv $(BBH_DIR)/fusion.csv $(TMP_DIR)/best_hit_fusion.csv
	{ cut -f1 $(TMP_DIR)/best_hit_fusion.csv | grep -F -v -f - -w $(word 1,$^); cat $(word 2,$^) $(word 3,$^); } > $@

# final annotation (contig => refseq)
$(TABLE_DIR)/annotation.csv : $(TMP_DIR)/xx.csv $(TABLE_DIR)/reference.csv
	perl repCol.pl -table1 $< -table2 $(word 2,$^) -add -column 2 -tkey 1 -tvalue 5 > $@

# BBH
$(BBH_DIR)/bbh.csv : $(BLAST_DIR)/best_CDS_trin2ref.csv $(BLAST_DIR)/best_CDS_ref2trin.csv | $(BBH_DIR)/.c
	@printf "$(P_INFO) Looking for best bidirectional hits.\n"
	perl BBH.pl -a $(word 1,$^) -b $(word 2,$^) -o $@

# SBH: query,target,strandQ,strandT,type,start,end
$(BBH_DIR)/sbh.csv : $(BLAST_DIR)/best_CDS_trin2ref_SBH.csv
	perl -ane 'print "$$F[0]\t$$F[1]\t$$F[7]\t1\tsbh\t$$F[8]\t$$F[9]\n";' $< > $@

# SBH: blast
$(BLAST_DIR)/best_CDS_trin2ref_SBH.csv : $(BLAST_DIR)/best_CDS_trin2ref.csv $(BBH_DIR)/bbh.csv
	@printf "$(P_INFO) Looking for single best hits\n"
	cat $< | perl -sane 'print if ($$F[15] > $$coverage && $$F[16] > $$identity);' -- $(OPT_SBH) | sort -k2,2V $(SORT_ORDER) -k1,1V | perl SBH.pl -annotated $(BBH_DIR)/bbh.csv > $@

# FUSION: best unknown gene
$(BBH_DIR)/fusion.csv : $(TMP_DIR)/fusion_full_report_sym.csv $(TMP_DIR)/bbh_sbh_sym.csv
	grep -F -v "^#" $< | grep -F -w -v -f $(word 2,$^) | perl -ane 'print "$$F[2]\t$$F[0]\t$$F[3]\t1\tsbh,part1\t$$F[4]\t$$F[5]\n";' | perl first_occurence.pl > $@

# FUSION full report
$(TMP_DIR)/fusion_full_report.csv : $(BLAST_DIR)/avg_CDS_ref2trin.csv $(TMP_DIR)/bbh_sbh.csv $(BLAST_DIR)/avg_CDS_trin2ref.csv $(REF_DIR)/ref_transcriptome_cds.fa.fai
	@printf "$(P_INFO) Looking for misassembled contigs\n"
	sort -k2,2V $(SORT_ORDER) -k1,1V $< | cut -f1,2,8,11,12,15,17 | perl blast_filter_fusion.pl -bbh $(word 2,$^) -hit-rev $(word 3,$^) -lq $(word 4,$^) $(OPT_FUSION) > $@

$(TMP_DIR)/bbh_sbh.csv : $(BBH_DIR)/bbh.csv $(BBH_DIR)/sbh.csv
	cat $^ > $@

$(TMP_DIR)/bbh_sbh_sym.csv : $(TMP_DIR)/bbh_sbh.csv $(TABLE_DIR)/reference.csv
	perl repCol.pl -table1 $< -table2 $(word 2,$^) -column 2 -tkey 1 -tvalue 5 | cut -f2 | sort -u > $@

$(TMP_DIR)/fusion_full_report_sym.csv : $(TMP_DIR)/fusion_full_report.csv $(TABLE_DIR)/reference.csv
	perl repCol.pl -table1 $< -table2 $(word 2,$^) -column 1 -tkey 1 -tvalue 5 -add > $@

# }}}
# BLAST {{{
# -----------------------------------------------------------------------------

# BLAST mRNA
$(BLAST_DIR)/raw_mRNA_ref2trin.csv : $(PRIPRO_ASSEMBLY).masked.fa $(MASKED_REF) | $(BLAST_DIR)/.c
	@printf "$(P_INFO) BLAST\n"
	@printf "$(P_INFO) query=reference database=assembly\n"
	@printf "$(P_LOG)/blast_mRNA2trin.log\n"
	@printf "$(P_OUT) $@\n\n"
	perl parallel_blast.pl -q $(word 2,$^) -db $< -o $@ $(BLAST_PARAMETER) > /dev/null

$(BLAST_DIR)/raw_mRNA_trin2ref.csv : $(MASKED_REF) $(PRIPRO_ASSEMBLY).masked.fa | $(BLAST_DIR)/.c
	@printf "$(P_INFO) BLAST\n"
	@printf "$(P_INFO) query=assembly database=reference\n"
	@printf "$(P_LOG)/blast_trin2mRNA.log\n"
	@printf "$(P_OUT) $@\n\n"
	perl parallel_blast.pl -q $(word 2,$^) -db $< -o $@ $(BLAST_PARAMETER) > /dev/null 2> $(LOG_DIR)/blast_trin2mRNA.log

# BLAST CDS
$(BLAST_DIR)/raw_CDS_ref2trin.csv : $(PRIPRO_ASSEMBLY).masked.fa $(MASKED_REF_CDS) | $(BLAST_DIR)/.c
	@printf "$(P_INFO) BLAST\n"
	@printf "$(P_INFO) query=reference database=assembly\n"
	@printf "$(P_LOG)/blast_cds2trin.log\n"
	@printf "$(P_OUT) $@\n\n"
	perl parallel_blast.pl -q $(word 2,$^) -db $< -o $@ $(BLAST_PARAMETER) > /dev/null 2> $(LOG_DIR)/blast_cds2trin.log

$(BLAST_DIR)/raw_CDS_trin2ref.csv : $(MASKED_REF_CDS) $(PRIPRO_ASSEMBLY).masked.fa | $(BLAST_DIR)/.c 
	@printf "$(P_INFO) BLAST\n"
	@printf "$(P_INFO) query=assembly database=reference\n"
	@printf "$(P_LOG)/blast_trin2cds.log\n"
	@printf "$(P_OUT) $@\n\n"
	perl parallel_blast.pl -q $(word 2,$^) -db $< -o $@ $(BLAST_PARAMETER) > /dev/null 2> $(LOG_DIR)/blast_cds2trin.log

# BLAST mRNA
$(BLAST_DIR)/avg_mRNA_ref2trin.csv :$(BLAST_DIR)/raw_mRNA_ref2trin.csv $(PRIPRO_ASSEMBLY).masked.fa.fai $(REF_DIR)/ref_transcriptome.fa.fai
	@printf "$(P_INFO) query=reference database=assembly\n"
	@printf "$(P_LOG)/blast_$(notdir $@).log\n"
	@printf "$(P_OUT) $@\n\n"
	$(COMMAND_SUM_BLAST) -input $< -lq $(word 3,$^) -lt $(word 2,$^) > $@ 2> $(LOG_DIR)/blast_$(notdir $@).log

$(BLAST_DIR)/avg_mRNA_trin2ref.csv : $(BLAST_DIR)/raw_mRNA_trin2ref.csv $(REF_DIR)/ref_transcriptome.fa.fai $(PRIPRO_ASSEMBLY).masked.fa.fai
	@printf "$(P_INFO) query=assembly database=reference\n"
	@printf "$(P_LOG)/blast_$(notdir $@).log\n"
	@printf "$(P_OUT) $@\n\n"
	$(COMMAND_SUM_BLAST) -input $< -lq $(word 3,$^) -lt $(word 2,$^) > $@ 2> $(LOG_DIR)/blast_$(notdir $@).log

# BLAST CDS
$(BLAST_DIR)/avg_CDS_ref2trin.csv : $(BLAST_DIR)/raw_CDS_ref2trin.csv $(PRIPRO_ASSEMBLY).masked.fa.fai $(REF_DIR)/ref_transcriptome_cds.fa.fai
	@printf "$(P_INFO) query=reference database=assembly\n"
	@printf "$(P_LOG)/blast_$(notdir $@).log\n"
	@printf "$(P_OUT) $@\n\n"
	$(COMMAND_SUM_BLAST) -input $< -lq $(word 3,$^) -lt $(word 2,$^) > $@ 2> $(LOG_DIR)/blast_$(notdir $@).log

$(BLAST_DIR)/avg_CDS_trin2ref.csv : $(BLAST_DIR)/raw_CDS_trin2ref.csv $(REF_DIR)/ref_transcriptome_cds.fa.fai $(PRIPRO_ASSEMBLY).masked.fa.fai
	@printf "$(P_INFO) query=assembly  database=reference\n"
	@printf "$(P_LOG)/blast_$(notdir $@).log\n"
	@printf "$(P_OUT) $@\n\n"
	$(COMMAND_SUM_BLAST) -input $< -lq $(word 3,$^) -lt $(word 2,$^) > $@ 2> $(LOG_DIR)/blast_$(notdir $@).log

# BEST
$(BLAST_DIR)/sorted_%.csv : $$(subst sorted,avg,$$@)
	sort $(SORT_ORDER) -k1,1V -k2,2V -o $@ $<

$(BLAST_DIR)/best_%.csv : $$(subst best,sorted,$$@)
	if [ ! -s $< ]; then  printf "$(P_ERR) Blast result empty:\n$(P_ERR)    $<\n"; exit 1; fi
	perl first_occurence.pl -i $< -o $@

# }}}
# FIND SCAFFOLDING CANDIDATES {{{
# -----------------------------------------------------------------------------

$(TMP_DIR)/all_contigs.txt : $(TABLE_DIR)/annotation.csv
	cut -f1 $< > $@

$(TMP_DIR)/all_references.txt : $(TABLE_DIR)/annotation.csv
	cut -f2 $< > $@

# remove already assigned queries and targets from blast results
$(TMP_DIR)/avg_mRNA_unassigned_contigs_trin2ref.csv : $(BLAST_DIR)/avg_mRNA_trin2ref.csv $(TMP_DIR)/all_contigs.txt $(TMP_DIR)/all_references.txt
	perl blast_filterIDs.pl -blast $< -ids $(word 3,$^) -output $@.tmp -type extract -seq subject
	perl blast_filterIDs.pl -blast $@.tmp -ids $(word 2,$^) -output $@ -type remove -seq query

# filter contigs according to specified options
$(TMP_DIR)/avg_mRNA_fragments_trin2ref.csv : $(TMP_DIR)/avg_mRNA_unassigned_contigs_trin2ref.csv $(PRIPRO_ASSEMBLY).fa.fai
	perl blast_filter.pl -blast $< -query-index $(word 2,$^) $(OPT_FRAGMENTS) | sort $(SORT_ORDER) -k1,1V -k2,2V -o $@

# best hit only
$(BLAST_DIR)/scaffolding_candidates_trin2ref.csv : $(TMP_DIR)/avg_mRNA_fragments_trin2ref.csv
	perl first_occurence.pl -i $< | sort -k2,2V > $@

# }}}
# CDS + MRNA BOUNDARY {{{
# -----------------------------------------------------------------------------

$(BLAST_DIR)/avg_CDS_ref2trin_sym.csv : $(BLAST_DIR)/avg_CDS_ref2trin.csv $(TABLE_DIR)/reference.csv
	perl repCol.pl -table1 $< -table2 $(word 2,$^) -column 1 -tkey 1 -tvalue 5 -add > $@

$(OUTPUT_DIR)/annotation_sym.csv : $(TABLE_DIR)/annotation.csv $(TABLE_DIR)/reference.csv
	perl repCol.pl -table1 $< -table2 $(word 2,$^) -column 2 -tkey 1 -tvalue 5 -add > $@

$(ASSEMBLY_DIR)/annotated.gbk : $(PRIPRO_ASSEMBLY).fa $(REF_TRANSCRIPTOME) \
	$(TABLE_DIR)/annotation.csv $(TRANS_DIR)/.c $(TRINITY_DIR)/dir2comp.index \
	$(BLAST_DIR)/avg_CDS_ref2trin_sym.csv \
	$(BLAST_DIR)/scaffolding_candidates_trin2ref.csv $(PATH_GENSCAN_MAT)
	@printf "\n$(P_INFO) Annotation\n"
	@printf "$(P_LOG)/annotation.log\n"
	@printf "$(P_LOG)/annotation.err\n\n"
	perl annotation.pl \
		-new  \
		-taxon-assemble $(SPECIES_ASS_ID) \
		-taxon-ortholog $(SPECIES_REF_ID) \
		-output-dir  $(TRANS_DIR) \
		-output-file $@ \
		-trinity     $(PRIPRO_ASSEMBLY).fa \
		-reference   $(REF_TRANSCRIPTOME) \
		-assignment  $(TABLE_DIR)/annotation.csv \
		-read-index  $(TRINITY_DIR)/dir2comp.index \
		-scaffolding $(BLAST_DIR)/scaffolding_candidates_trin2ref.csv \
		-blast 		 $(BLAST_DIR)/avg_CDS_ref2trin_sym.csv \
		$(ANNOTATION_PARAMETER) \
			> $(LOG_DIR)/annotation.log \
			2> $(LOG_DIR)/annotation.err

ifeq (3,$(TRINITY_VERSION))
# at least since Version 2.1
ifneq (,$(findstring genome_guided_bam,$(OPT_TRINITY)))
# genome guided
$(TRINITY_DIR)/dir2comp.csv : $(TRINITY_DIR)/Trinity.fasta
	perl find_ggreadlist.pl $(TRINITY_DIR) $(TRINITY_DIR)/Trinity.fasta | perl -pe 's/TRINITY_GG_(\d+)_.+?\s/c$$1\t/g' | sort -k1,1V  | first_occurence.pl > $@
else
# de novo
$(TRINITY_DIR)/dir2comp.csv : $(TRINITY_DIR)/Trinity.fasta $(TRINITY_DIR)/partitioned_reads.files.list
	cat $(TRINITY_DIR)/partitioned_reads.files.list | perl -ne 'chomp; $$a = $$_; $$a =~ s/^.+\/(c\d+).trinity.reads.fa$$/$$1/g; print "$$a\t$$_\n"' > $@
endif
else
# Versions: 2013 - 2014
$(TRINITY_DIR)/dir2comp.csv : $(TRINITY_DIR)/Trinity.fasta
	find $(TRINITY_DIR) -name "*.reads.tmp" | perl -ne 'chomp; $$a = $$_; $$a =~ s/^.+\/(c\d+).reads.tmp$$/$$1/g; print "$$a\t$$_\n"' > $@
endif

$(TRINITY_DIR)/dir2comp.index : $(TRINITY_DIR)/dir2comp.csv
	perl store_hash.pl -input $< -c1 1 -c2 2 -noarray -output $@

# }}}
# Vector contamination (or simple include in annotation.pl?) {{{
# -----------------------------------------------------------------------------

$(BLAST_DIR)/vecscreen.csv : $(FINAL_CDS).fa
	@printf "$(P_INFO) Vector screening\n"
	perl blast_VecScreen.pl -sequences $< -outstem $(basename $@) -cpu $(OPT_CPUS) -database $(PATH_BLASTDB) -ncbi $(NCBI_BLAST)

$(FINAL_GENBANK).gbk : $(ASSEMBLY_DIR)/annotated.gbk $(BLAST_DIR)/vecscreen.csv
	if [ ! -s $< ]; then  printf "$(P_ERR) Annotation failed:\n$(P_ERR)    $<\n"; exit 1; fi
	perl contamination_feature.pl -input $< -output $@ -contamination $(lastword $^)

# }}}
# FASTA CONVERSION {{{
# -----------------------------------------------------------------------------

$(FINAL_CDS).fa : $(ASSEMBLY_DIR)/annotated.gbk
	perl genbank2fasta.pl -i $< -o $@ -feature CDS -feature-sequence

$(FINAL_MRNA).fa : $(ASSEMBLY_DIR)/annotated.gbk
	perl genbank2fasta.pl -i $< -o $@ -feature mRNA -feature-sequence

# }}}
# SUMMARIZE {{{
# -----------------------------------------------------------------------------

$(FINAL_SUMMARY).pdf : $(FINAL_SUMMARY).csv $(TABLE_DIR)/reference.csv $(PRIPRO_ASSEMBLY).fa.fai | $(TABLE_DIR)/.c
	@printf "$(P_INFO) Summarize results.\n"
	@printf "$(P_LOG)/summary.log\n"
	@printf "$(P_OUT) $@\n\n"
	Rscript $(OPT_R) summarize.R $< $(TABLE_DIR)/reference.csv $(PRIPRO_ASSEMBLY).fa.fai $(TABLE_DIR) $@ &> $(LOG_DIR)/summary.log

$(FINAL_SUMMARY).csv : $(FINAL_GENBANK).gbk | $(TABLE_DIR)/.c
	@printf "$(P_INFO) Create summary table.\n"
	@printf "$(P_OUT) $@\n\n"
	perl genbank2table.pl -genbank $< -assembly -o $@

$(TABLE_DIR)/reference.csv : $(REF_TRANSCRIPTOME) | $(TABLE_DIR)/.c
	perl genbank2table.pl -genbank $< -o $@

$(FINAL_CDS).csv : $(FINAL_GENBANK).gbk | $(TABLE_DIR)/.c
	@printf "$(P_INFO) Extract CDS coordinates.\n"
	@printf "$(P_OUT) $@\n\n"
	perl genbank2table.pl -genbank $< -cds -o $@

# }}}
# FUNCTIONAL ANNOTATION {{{
# -----------------------------------------------------------------------------

$(FINAL_GO).pdf : $(FINAL_SUMMARY).csv
	@printf "$(P_INFO) GO annotation.\n"
	@printf "$(P_LOG)/GO.log\n\n"
	Rscript $(OPT_R) GO.R $< $(OPT_ANNOTATION) $(TABLE_DIR) $@ &> $(LOG_DIR)/GO.log

$(FINAL_KEGG).pdf : $(FINAL_SUMMARY).csv
	@printf "$(P_INFO) KEGG annotation.\n"
	@printf "$(P_LOG)/KEGG.log\n\n"
	Rscript $(OPT_R) KEGG.R $< $(OPT_ANNOTATION) $(TABLE_DIR) $@ &> $(LOG_DIR)/KEGG.log

# }}}

