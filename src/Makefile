#       AUTHOR: Martin Bens, bensmartin@gmail.com
# ORGANIZATION: FLI Jena

.SECONDEXPANSION:

ifdef FILE
include $(FILE)
else
$(error Use "-FILE" parameter to specify configuration file)
endif

export PATH := $(shell echo $(PATH_BAMTOOLS):$(PATH_BOWTIE):$(PATH_CD_HIT_EST):$(PATH_EMBOSS):$(PATH_GENSCAN):$(PATH_MAFFT):$(PATH_PERL):$(PATH_REPEATMASKER):$(PATH_SAMTOOLS):$(PATH_TGICL):$(PATH_TRINITY):$(PATH_WUBLAST):$(PATH_XDFORMAT):$(PATH_RSCRIPT):$(PATH) | perl -pe "s/:{2,}/:/g; s/^://;")

ifndef OUTPUT_DIR
	$(error OUTPUT_DIR not defined!)
endif


# -----------------------------------------------------------------------------
# DIRECTORIES
# -----------------------------------------------------------------------------

LOG_DIR := $(OUTPUT_DIR)/logs
TABLE_DIR := $(OUTPUT_DIR)/tables
TRANS_DIR := $(OUTPUT_DIR)/transcripts
TRINITY_DIR := $(OUTPUT_DIR)/trinity

BLAST_DIR := $(OUTPUT_DIR)/blast
BLASTDB_DIR := $(OUTPUT_DIR)/db
BBH_DIR := $(OUTPUT_DIR)/bbh
ASSEMBLY_DIR := $(OUTPUT_DIR)/assembly
REF_DIR := $(OUTPUT_DIR)/reference
MASK_DIR := $(OUTPUT_DIR)/masking

DATE = `date +'%y-%m-%d %H:%M:%S'`:

ifndef DB_VEC
DB_VEC := UniVec_Core
endif

# default: parallel, switch to GNU parallel if OPT_MAX_SGE is not set
ENGINE := sge

ifndef OPT_MAX_SGE
ENGINE := parallel
OPT_MAX_SGE := $(OPT_CPUS)
endif

# use output directory for temporary files
ifndef TMP_DIR
TMP_DIR := $(OUTPUT_DIR)/tmp
endif

ANNOTATION_CUSTOM_PARAMETER := 
ifdef ORTHOLOG_TABLE
	ANNOTATION_CUSTOM_PARAMETER := -ortholog-table $(ORTHOLOG_TABLE) -ortholog-cds $(ORTHOLOG_CDS)
endif

ifneq (,$(findstring tgicl,$(ASSEMBLY_PREPROCESS)))
	ANNOTATION_CUSTOM_PARAMETER += -tgicl $(ASSEMBLY_DIR)/Trinity.tgicl.fa.clstr
endif
$(info $(ANNOTATION_CUSTOM_PARAMETER))

# -----------------------------------------------------------------------------
# DEFAULT PARAMETER
# -----------------------------------------------------------------------------

# Trinity
OPT_TRINITY += --no_cleanup --CPU $(OPT_CPUS)

# Blast
OPT_WUBLAST_BLASTN += -wordmask=seg lcmask -topcomboN 3 -cpus 1

# RepeatMasker
OPT_REPEAT_REF_TRANSCRIPTOME += -xsmall -par $(OPT_CPUS)
OPT_REPEAT_ASSEMBLY += -xsmall -par $(OPT_CPUS)

# TGICL
OPT_TGICL += -c $(OPT_CPUS)

# CD-HIT-EST
OPT_CD_HIT_EST:= -T $(OPT_CPUS)

ifdef ASSEMBLY_PREPROCESS
	ASSEMBLY_PREPROCESS := $(addprefix $(LOG_DIR)/,$(ASSEMBLY_PREPROCESS))
	ASSEMBLY_PREPROCESS := $(addsuffix .finished,$(ASSEMBLY_PREPROCESS))
endif
ASSEMBLY_PREPROCESS := $(LOG_DIR)/trinity-tmp.finished $(ASSEMBLY_PREPROCESS)

# defines space and newline variables
NOTHING :=
space := $(NOTHING) $(NOTHING)
define newline


endef

# -----------------------------------------------------------------------------
# READ DATA (FASTQ and FASTA only)
# -----------------------------------------------------------------------------

# at the moment, trinity accepts only one type of file at a time
READS := $(wildcard $(READ_DIR)/*.fq)
READS += $(wildcard $(READ_DIR)/*.fastq)
READS += $(wildcard $(READ_DIR)/*.fq)
READS += $(wildcard $(READ_DIR)/*.fa)
READS += $(wildcard $(READ_DIR)/*.fastq.gz)
READS += $(wildcard $(READ_DIR)/*.fq.gz)
READS += $(wildcard $(READ_DIR)/*.fq.gz)
READS += $(wildcard $(READ_DIR)/*.fa.gz)

# we don't want to copy reads
READ_LINKS := $(subst $(READ_DIR)/,$(OUTPUT_DIR)/reads/,$(READS))

# -----------------------------------------------------------------------------
# Pipeline
# -----------------------------------------------------------------------------

.PHONY : all check trinity preprocess-assembly blast blast-sum annotation vector fasta-conversion pipeline-stats functional

all: check trinity preprocess-assembly blast blast-sum annotation vector fasta-conversion pipeline-stats functional
	@echo $(DATE) Pipeline finished.

check: $(FILE) $(LOG_DIR)/.c $(REF_DIR)/.c
	@echo $(DATE) Check required software and modules
	@echo -e "$(DATE)\tLog: $(LOG_DIR)/software_check.log"
	$(PATH_PERL)perl init.pl $< 1> $(LOG_DIR)/software_check.log

trinity: $(TRINITY_DIR)/Trinity.fa
	@echo $(DATE) Trinity finished.

preprocess-assembly: $(LOG_DIR)/preprocessing.finished $(OUTPUT_DIR)/assembly.fa
	@echo $(DATE) Primary processing finished.

blast : $(BLAST_DIR)/raw_CDS_trin2ref.csv.csv $(BLAST_DIR)/raw_CDS_ref2trin.csv $(BLAST_DIR)/raw_mRNA_ref2trin.csv $(BLAST_DIR)/raw_mRNA_trin2ref.csv
	@echo $(DATE) Blast finished.

blast-sum : $(BLAST_DIR)/avg_CDS_trin2ref.csv $(BLAST_DIR)/avg_CDS_ref2trin.csv $(BLAST_DIR)/avg_mRNA_ref2trin.csv $(BLAST_DIR)/avg_mRNA_trin2ref.csv
	@echo $(DATE) Blast preprocessing finished.

annotation : $(ASSEMBLY_DIR)/annotated.gbk
	@echo $(DATE) Annotation finished.

vector: $(OUTPUT_DIR)/transcript_catalogue.gbk
	@echo $(DATE) Annotation of possible contamination finished.

fasta-conversion: $(OUTPUT_DIR)/sequences-cds.fa $(OUTPUT_DIR)/sequences-mRNA.fa
	@echo $(DATE) Fasta conversion finished.

pipeline-stats : $(OUTPUT_DIR)/summary.pdf
	@echo $(DATE) Statistics finished.

functional : $(OUTPUT_DIR)/summary_go.pdf $(OUTPUT_DIR)/summary_kegg.pdf
	@echo $(DATE) Functional annotation.


# -----------------------------------------------------------------------------
# GENERAL RULES
# -----------------------------------------------------------------------------

# Create directory
$(OUTPUT_DIR)/%/.c :
	mkdir -p $(@D)
	touch $@

# Index fasta
%.fa.fai: %.fa 
	samtools faidx $<

# -----------------------------------------------------------------------------
# MASKING
# -----------------------------------------------------------------------------

# - mask repeats only if path to RepeatMasker is available (create symlinks otherwise)
# - switched to NCBI Blast due to performance reason (and failing masking
#   runs with wublast)
# - performance boost: concatenating sequences for RepeatMasker
# - repeat library specified in configuration file

$(ASSEMBLY_DIR)/Trinity_preprocessed.fa.masked : $(ASSEMBLY_DIR)/Trinity_preprocessed.fa $(MASK_DIR)/.c $(TMP_DIR)/.c
ifeq ($(REPEAT),0) 
	ln -f -s $< $@
else
	@echo $(DATE) Repeat masking assembly.
	@echo -e "$(DATE)\tLog: $(LOG_DIR)/repeat_masking_assembly.log"
	$(PATH_PERL)perl fasta_concatenate.pl -file $< -out $(ASSEMBLY_DIR)/Trinity_preprocessed_long
	cd $(OUTPUT_DIR) && $(PATH_REPEATMASKER)RepeatMasker $(ASSEMBLY_DIR)/Trinity_preprocessed_long.fa $(OPT_REPEAT_ASSEMBLY) -dir $(MASK_DIR) &> $(LOG_DIR)/repeat_masking_assembly.log
	$(PATH_PERL)perl fasta_unconcat.pl -file $(MASK_DIR)/Trinity_preprocessed_long.fa.masked -index $(ASSEMBLY_DIR)/Trinity_preprocessed_long.index -o $@
endif 

$(OUTPUT_DIR)/assembly.fa : $(ASSEMBLY_DIR)/Trinity_preprocessed.fa.masked
	cp $< $@

# -----------------------------------------------------------------------------
# BLAST DBs
# -----------------------------------------------------------------------------

$(BLASTDB_DIR)/ref_transcriptome.fa.xnd : $(REF_DIR)/ref_transcriptome.fa $(BLASTDB_DIR)/.c
	$(PATH_XDFORMAT)xdformat -n -o $(basename $@) $< &> $(basename $@).log

$(BLASTDB_DIR)/ref_transcriptome_cds.fa.xnd : $(REF_DIR)/ref_transcriptome_cds.fa $(BLASTDB_DIR)/.c
	$(PATH_XDFORMAT)xdformat -n -o $(basename $@) $< &> $(basename $@).log

$(BLASTDB_DIR)/ref_transcriptome.fa.masked.xnd : $(REF_DIR)/ref_transcriptome.fa.masked $(BLASTDB_DIR)/.c
	$(PATH_XDFORMAT)xdformat -n -o $(basename $@) $< &> $(basename $@).log

$(BLASTDB_DIR)/ref_transcriptome_cds.fa.masked.xnd : $(REF_DIR)/ref_transcriptome_cds.fa.masked $(BLASTDB_DIR)/.c
	$(PATH_XDFORMAT)xdformat -n -o $(basename $@)  $< &> $(basename $@).log

$(BLASTDB_DIR)/Trinity_preprocessed.fa.xnd : $(ASSEMBLY_DIR)/Trinity_preprocessed.fa $(BLASTDB_DIR)/.c
	$(PATH_XDFORMAT)xdformat -n -o $(basename $@) $< &> $(basename $@).log

$(BLASTDB_DIR)/Trinity_preprocessed.fa.masked.xnd : $(ASSEMBLY_DIR)/Trinity_preprocessed.fa.masked $(BLASTDB_DIR)/.c
	$(PATH_XDFORMAT)xdformat -n -o $(basename $@) $< &> $(basename $@).log

# -----------------------------------------------------------------------------
# REFERENCE FASTA FILES
# -----------------------------------------------------------------------------

$(REF_DIR)/ref_transcriptome.fa : $(REF_TRANSCRIPTOME)
ifndef REF_TRANSCRIPTOME_FASTA
	$(PATH_PERL)perl genbank2fasta.pl -i $< -o $@ -feature CDS
else
	ln -f -s $(REF_TRANSCRIPTOME_FASTA) $@
endif

$(REF_DIR)/ref_transcriptome_cds.fa : $(REF_TRANSCRIPTOME)
ifndef REF_TRANSCRIPTOME_FASTA_CDS
	$(PATH_PERL)perl genbank2fasta.pl -i $< -o $@ -feature CDS -feature-sequence
else
	ln -f -s $(REF_TRANSCRIPTOME_FASTA_CDS) $@
endif

$(REF_DIR)/ref_transcriptome.fa.masked : $(REF_DIR)/ref_transcriptome.fa $(MASK_DIR)/.c $(TMP_DIR)/.c
ifdef REF_TRANSCRIPTOME_FASTA_MASKED
	ln -f -s $(REF_TRANSCRIPTOME_FASTA_MASKED) $@
else ifeq ($(REPEAT),0)
	ln -f -s $< $@
else
	@echo $(DATE) Repeat masking reference.
	@echo -e "$(DATE)\tLog: $(LOG_DIR)/repeat_masking_reference.log"
	$(PATH_PERL)perl fasta_concatenate.pl -file $< -out $(TMP_DIR)/$(basename $(notdir $<))
	cd $(OUTPUT_DIR) && $(PATH_REPEATMASKER)RepeatMasker $(TMP_DIR)/$(notdir $<) -dir $(MASK_DIR) $(OPT_REPEAT_REF_TRANSCRIPTOME) &> $(LOG_DIR)/repeat_masking_reference.log
	$(PATH_PERL)perl fasta_unconcat.pl -file $(MASK_DIR)/$(notdir $@) -index $(TMP_DIR)/$(basename $(notdir $<)).index > $@
endif

$(REF_DIR)/ref_transcriptome_cds.fa.masked : $(REF_TRANSCRIPTOME) $(REF_DIR)/ref_transcriptome.fa.masked 
ifndef REF_TRANSCRIPTOME_FASTA_CDS_MASKED
	$(PATH_PERL)perl genbank_CDS.pl -genbank $< -input $(word 2,$^) -output $@
else
	ln -f -s $(REF_TRANSCRIPTOME_FASTA_CDS_MASKED) $@
endif

# -----------------------------------------------------------------------------
# READS
# -----------------------------------------------------------------------------

# creates symbolic links to read data
$(READ_LINKS) : $$(subst $(OUTPUT_DIR)/reads/,$(READ_DIR)/,$$@) $(OUTPUT_DIR)/reads/.c
	ln -f -s $< $@

# -----------------------------------------------------------------------------
# TRANSCRIPTOME ASSEMBLY
# -----------------------------------------------------------------------------

# - Trinity handles precomputed results by itself
# - JVM decides how many garbage collection threads to use => we limit
#   GC-threads to one and start half of allowed process in parallel

$(TRINITY_DIR)/Trinity.fa : $(TRINITY_DIR)/.c $(READ_LINKS)
	@echo $(DATE) Starting Trinity.
	@echo -e "$(DATE)\tLog: $(LOG_DIR)/inchworm.log"

	$(eval read_R1 := $(subst $(newline),$(space),$(shell for i in $(filter-out $<,$^); do echo $$i | grep "_R1"; done)))
	$(eval read_R2 := $(subst $(newline),$(space),$(shell for i in $(filter-out $<,$^); do echo $$i | grep "_R2"; done)))
	$(eval TRINITY_DIR := $(dir $@))
ifeq ($(OPT_READTYPE),pe)
	ulimit -s unlimited; $(PATH_TRINITY)Trinity.pl --left $(read_R1) --right $(read_R2) --output $(TRINITY_DIR) --no_run_chrysalis $(OPT_TRINITY) &> $(LOG_DIR)/inchworm.log 
	@echo -e "$(DATE)\tLog: $(LOG_DIR)/chrysalis.log"
	ulimit -s unlimited; $(PATH_TRINITY)Trinity.pl --left $(read_R1) --right $(read_R2) --output $(TRINITY_DIR) --no_run_quantifygraph $(OPT_TRINITY) $(OPT_BUTTERFLY)  &> $(LOG_DIR)/chrysalis.log 
else
	ulimit -s unlimited; $(PATH_TRINITY)Trinity.pl --single $(filter-out $<,$^) --output $(TRINITY_DIR) --no_run_chrysalis $(OPT_TRINITY) &> $(LOG_DIR)/inchworm.log 
	@echo -e "$(DATE)\tLog: $(LOG_DIR)/chrysalis.log"
	ulimit -s unlimited; $(PATH_TRINITY)Trinity.pl --single $(filter-out $<,$^) --output $(TRINITY_DIR) --no_run_quantifygraph $(OPT_TRINITY) $(OPT_BUTTERFLY)  &> $(LOG_DIR)/chrysalis.log 
endif
	@echo -e "$(DATE)\tLog: $(LOG_DIR)/quanitfyGraph.log"
	ulimit -s unlimited; cat $(TRINITY_DIR)chrysalis/quantifyGraph_commands | parallel --max-procs $(OPT_CPUS) --joblog $(LOG_DIR)/quantifyGraph.log
	$(PATH_RSCRIPT)Rscript parallel_checkLog.R $(LOG_DIR)/quantifyGraph.log
	@echo -e "$(DATE)\tLog: $(LOG_DIR)/butterfly.log"
	ulimit -s unlimited; cat $(TRINITY_DIR)/chrysalis/butterfly_commands | parallel --max-procs $(shell echo '$(OPT_CPUS)/2' | bc) --joblog $(LOG_DIR)/butterfly.log
	$(PATH_RSCRIPT)Rscript parallel_checkLog.R $(LOG_DIR)/butterfly.log
	ulimit -s unlimited; find $(TRINITY_DIR)/chrysalis -name "*allProbPaths.fasta" -exec cat {} + > $@

$(ASSEMBLY_DIR)/Trinity.fa : $(TRINITY_DIR)/Trinity.fa $(ASSEMBLY_DIR)/.c
	cp $< $@

# -----------------------------------------------------------------------------
# TRANSCRIPTOME PREPROCESSING
# -----------------------------------------------------------------------------

$(LOG_DIR)/preprocessing.finished : $(ASSEMBLY_PREPROCESS)
	touch $@

$(ASSEMBLY_DIR)/Trinity_preprocessed.fa: 
	if [ -a $@.index ]; then rm $@.index; fi
	cp $(ASSEMBLY_DIR)/Trinity.tmp.fa $@


# Create temporary file which all preprocessing steps will use
$(LOG_DIR)/trinity-tmp.finished : $(ASSEMBLY_DIR)/Trinity.fa
	cp $(ASSEMBLY_DIR)/Trinity.fa $(ASSEMBLY_DIR)/Trinity.tmp.fa
	touch $@

# CD-HIT
$(LOG_DIR)/cd-hit.finished : $(ASSEMBLY_DIR)/Trinity.cdhit.fa
	cp $< $(ASSEMBLY_DIR)/Trinity.tmp.fa
	touch $@

$(ASSEMBLY_DIR)/Trinity.cdhit.fa : $(ASSEMBLY_DIR)/Trinity.tmp.fa $(LOG_DIR)/.c
	@echo $(DATE) Clustering sequences: CDHIT-EST
	@echo -e "$(DATE)\tLog: $(LOG_DIR)/cdhit.log"
	$(PATH_CD_HIT_EST)cd-hit-est -i $< -o $@ $(OPT_CD_HIT_EST) &> $(LOG_DIR)/cdhit.log

# TGICL
$(LOG_DIR)/tgicl.finished : $(ASSEMBLY_DIR)/Trinity.tgicl.fa
	cp $< $(ASSEMBLY_DIR)/Trinity.tmp.fa
	touch $@

$(ASSEMBLY_DIR)/Trinity.tgicl.fa : $(ASSEMBLY_DIR)/Trinity.tmp.fa $(ASSEMBLY_DIR)/tgicl/.c
	@echo $(DATE) Meta assembly: TGICL
	$(eval TMPFILE := $(ASSEMBLY_DIR)/tgicl/Trinity.tmp.fa)
	cd $(ASSEMBLY_DIR)/tgicl; $(PATH_TGICL)tgicl $(ASSEMBLY_DIR)/Trinity.tmp.fa $(OPT_TGICL); 
	tgicl.sh $< $(ASSEMBLY_DIR)/tgicl/Trinity.tmp.fa $(ASSEMBLY_DIR)/tgicl $@
	
# -----------------------------------------------------------------------------
# ASSEMBLY ANNOTATION
# -----------------------------------------------------------------------------

# final annotation (contig => refseq)
$(TABLE_DIR)/annotation.csv : $(BBH_DIR)/bbh_sbh.csv $(BBH_DIR)/fusion.csv $(TABLE_DIR)/reference.csv
	cut -f1 $(word 2,$^) | grep -F -f - -w $(word 1,$^) | $(PATH_PERL)perl -ane '$$F[4] .= ",part0"; print join "\t", @F; print "\n"' > $(BBH_DIR)/best_hit_fusion.csv
	{ cut -f1 $(BBH_DIR)/best_hit_fusion.csv | grep -F -v -f - -w $(word 1,$^); cat $(word 2,$^) $(BBH_DIR)/best_hit_fusion.csv; } > $@.tmp
	$(PATH_PERL)perl repCol.pl -table1 $@.tmp -table2 $(word 3,$^) -add -column 2 -tkey 1 -tvalue 5 > $@
	rm $@.tmp
	
# BBH
$(BBH_DIR)/bbh.csv : $(BLAST_DIR)/best_CDS_trin2ref.csv $(BLAST_DIR)/best_CDS_ref2trin.csv $(BBH_DIR)/.c
	@echo $(DATE) Looking for BBHs
	$(PATH_PERL)perl BBH.pl -a $(word 1,$^) -b $(word 2,$^) -o $@
	
# SBH: query,target,strandQ,strandT,type,start,end
$(BBH_DIR)/sbh.csv : $(BLAST_DIR)/best_CDS_trin2ref_SBH.csv
	$(PATH_PERL)perl -ane 'print "$$F[0]\t$$F[1]\t$$F[7]\t1\tsbh\t$$F[8]\t$$F[9]\n";' $< > $@
	
# SBH: blast 
$(BLAST_DIR)/best_CDS_trin2ref_SBH.csv : $(BLAST_DIR)/best_CDS_trin2ref.csv $(BBH_DIR)/bbh.csv 
	@echo $(DATE) Looking for SBHs
	cat $< | $(PATH_PERL)perl -sane 'print if ($$F[15] > $$coverage && $$F[16] > $$identity);' -- $(OPT_SBH) | sort -k2,2V $(SORT_ORDER) -k1,1V | $(PATH_PERL)perl SBH.pl -annotated $(BBH_DIR)/bbh.csv > $@

# FUSION: best unknown gene
$(BBH_DIR)/fusion.csv : $(BBH_DIR)/fusion_full_report_sym.csv $(BBH_DIR)/bbh_sbh_sym.csv
	grep -F -v "^#" $< | grep -F -w -v -f $(word 2,$^) | $(PATH_PERL)perl -ane 'print "$$F[2]\t$$F[0]\t$$F[3]\t1\tsbh,part1\t$$F[4]\t$$F[5]\n";' | $(PATH_PERL)perl first_occurence.pl > $@

# FUSION full report
$(BBH_DIR)/fusion_full_report.csv : $(BLAST_DIR)/avg_CDS_ref2trin.csv $(BBH_DIR)/bbh_sbh.csv $(BLAST_DIR)/avg_CDS_trin2ref.csv $(REF_DIR)/ref_transcriptome_cds.fa.fai
	@echo $(DATE) Looking for misassembled contigs
	sort -k2,2V $(SORT_ORDER) -k1,1V $< | cut -f1,2,8,11,12,15,17 | $(PATH_PERL)perl blast_filter_fusion.pl -bbh $(word 2,$^) -hit-rev $(word 3,$^) -lq $(word 4,$^) $(OPT_FUSION) > $@

$(BBH_DIR)/bbh_sbh.csv : $(BBH_DIR)/bbh.csv $(BBH_DIR)/sbh.csv
	cat $^ > $@

$(BBH_DIR)/bbh_sbh_sym.csv : $(BBH_DIR)/bbh_sbh.csv $(TABLE_DIR)/reference.csv
	$(PATH_PERL)perl repCol.pl -table1 $< -table2 $(word 2,$^) -column 2 -tkey 1 -tvalue 5 | cut -f2 | sort -u > $@

$(BBH_DIR)/fusion_full_report_sym.csv : $(BBH_DIR)/fusion_full_report.csv $(TABLE_DIR)/reference.csv
	$(PATH_PERL)perl repCol.pl -table1 $< -table2 $(word 2,$^) -column 1 -tkey 1 -tvalue 5 -add > $@

# -----------------------------------------------------------------------------
#  BLAST
# -----------------------------------------------------------------------------

BLAST_PARAMETER=-j $(OPT_MAX_SGE) -format tabular -flavor wu -parameter '$(OPT_WUBLAST_BLASTN)' -engine $(ENGINE) 

# evalue (asc) > highscore (desc) > target coverage (desc) > total score (desc) > query coverage (desc) > query name (asc) > target name (asc)
SORT_ORDER=-k4,4g -k6,6nr -k16,16nr -k5,5nr -k15,15nr

# MAPPING

# BLAST mRNA
$(BLAST_DIR)/raw_mRNA_ref2trin.csv : $(BLASTDB_DIR)/Trinity_preprocessed.fa.masked.xnd $(REF_DIR)/ref_transcriptome.fa.masked $(BLAST_DIR)/.c $(ASSEMBLY_DIR)/Trinity_preprocessed.fa.masked
	@echo $(DATE) blasting mRNA Q:reference DB:assembly
	@echo -e "$(DATE)\tLog: $(LOG_DIR)/blast_mRNA2trin.log"
	$(PATH_PERL)perl sge.pl -t blastn -q $(word 2,$^) -db $(basename $<) -o $@ $(BLAST_PARAMETER) > /dev/null 2> $(LOG_DIR)/blast_mRNA2trin.log

$(BLAST_DIR)/raw_mRNA_trin2ref.csv : $(BLASTDB_DIR)/ref_transcriptome.fa.masked.xnd $(ASSEMBLY_DIR)/Trinity_preprocessed.fa.masked $(BLAST_DIR)/.c $(REF_DIR)/ref_transcriptome.fa.masked
	@echo $(DATE) blasting mRNA Q:assembly DB:reference
	@echo -e "$(DATE)\tLog: $(LOG_DIR)/blast_trin2mRNA.log"
	$(PATH_PERL)perl sge.pl -t blastn -q $(word 2,$^) -db $(basename $<) -o $@ $(BLAST_PARAMETER) > /dev/null 2> $(LOG_DIR)/blast_trin2mRNA.log

# BLAST CDS
$(BLAST_DIR)/raw_CDS_ref2trin.csv : $(BLASTDB_DIR)/Trinity_preprocessed.fa.masked.xnd $(REF_DIR)/ref_transcriptome_cds.fa.masked $(BLAST_DIR)/.c $(ASSEMBLY_DIR)/Trinity_preprocessed.fa.masked
	@echo $(DATE) blasting CDS Q:reference DB:assembly
	@echo -e "$(DATE)\tLog: $(LOG_DIR)/blast_cds2trin.log"
	$(PATH_PERL)perl sge.pl -t blastn -q $(word 2,$^) -db $(basename $<) -o $@ $(BLAST_PARAMETER) > /dev/null 2> $(LOG_DIR)/blast_cds2trin.log

$(BLAST_DIR)/raw_CDS_trin2ref.csv.csv : $(BLASTDB_DIR)/ref_transcriptome_cds.fa.masked.xnd $(ASSEMBLY_DIR)/Trinity_preprocessed.fa.masked $(BLAST_DIR)/.c $(REF_DIR)/ref_transcriptome_cds.fa.masked
	@echo $(DATE) blasting CDS Q:assembly DB:reference
	@echo -e "$(DATE)\tLog: $(LOG_DIR)/blast_trin2cds.log"
	$(PATH_PERL)perl sge.pl -t blastn -q $(word 2,$^) -db $(basename $<) -o $@ $(BLAST_PARAMETER) > /dev/null 2> $(LOG_DIR)/blast_trin2cds.log

# SUMMARIZE BLAST

# BLAST mRNA
$(BLAST_DIR)/avg_mRNA_ref2trin.csv :$(BLAST_DIR)/raw_mRNA_ref2trin.csv $(ASSEMBLY_DIR)/Trinity_preprocessed.fa.fai $(REF_DIR)/ref_transcriptome.fa.fai 
	@echo $(DATE) average HSPs mRNA Q:reference DB:assembly
	@echo -e "$(DATE)\tLog: $(LOG_DIR)/blast_$(notdir $@).log"
	$(PATH_PERL)perl wublast_average.pl -blast $< -lq $(word 3,$^) -lt $(word 2,$^) > $@ 2> $(LOG_DIR)/blast_$(notdir $@).log

$(BLAST_DIR)/avg_mRNA_trin2ref.csv : $(BLAST_DIR)/raw_mRNA_trin2ref.csv $(REF_DIR)/ref_transcriptome.fa.fai $(ASSEMBLY_DIR)/Trinity_preprocessed.fa.fai
	@echo $(DATE) average HSPs mRNA Q:assembly DB:reference
	@echo -e "$(DATE)\tLog: $(LOG_DIR)/blast_$(notdir $@).log"
	$(PATH_PERL)perl wublast_average.pl -blast $< -lq $(word 3,$^) -lt $(word 2,$^) > $@ 2> $(LOG_DIR)/blast_$(notdir $@).log

# BLAST CDS
$(BLAST_DIR)/avg_CDS_ref2trin.csv : $(BLAST_DIR)/raw_CDS_ref2trin.csv $(ASSEMBLY_DIR)/Trinity_preprocessed.fa.fai $(REF_DIR)/ref_transcriptome_cds.fa.fai
	@echo $(DATE) average HSPs CDS Q:reference DB:assembly
	@echo -e "$(DATE)\tLog: $(LOG_DIR)/blast_$(notdir $@).log"
	$(PATH_PERL)perl wublast_average.pl -blast $< -lq $(word 3,$^) -lt $(word 2,$^) > $@ 2> $(LOG_DIR)/blast_$(notdir $@).log

$(BLAST_DIR)/avg_CDS_trin2ref.csv : $(BLAST_DIR)/raw_CDS_trin2ref.csv.csv $(REF_DIR)/ref_transcriptome_cds.fa.fai $(ASSEMBLY_DIR)/Trinity_preprocessed.fa.fai
	@echo $(DATE) average HSPs CDS Q:assembly  DB:reference
	@echo -e "$(DATE)\tLog: $(LOG_DIR)/blast_$(notdir $@).log"
	$(PATH_PERL)perl wublast_average.pl -blast $< -lq $(word 3,$^) -lt $(word 2,$^) > $@ 2> $(LOG_DIR)/blast_$(notdir $@).log

# BEST

$(BLAST_DIR)/sorted_%.csv : $$(subst sorted,avg,$$@)
	sort $(SORT_ORDER) -k1,1V -k2,2V -o $@ $<

$(BLAST_DIR)/best_%.csv : $$(subst best,sorted,$$@) 
	$(PATH_PERL)perl first_occurence.pl -i $< -o $@

# -----------------------------------------------------------------------------
# TRANSCRIPTOME SCAFFOLDING
# -----------------------------------------------------------------------------

$(BBH_DIR)/all_contigs.txt : $(TABLE_DIR)/annotation.csv
	cut -f1 $< > $@

$(BBH_DIR)/all_references.txt : $(TABLE_DIR)/annotation.csv
	cut -f2 $< > $@

# remove already assigned queries and targets from blast results
$(BLAST_DIR)/avg_mRNA_unassigned_contigs_trin2ref.csv : $(BLAST_DIR)/avg_mRNA_trin2ref.csv $(BBH_DIR)/all_contigs.txt $(BBH_DIR)/all_references.txt
	@echo $(DATE) Filter blast for assigned contigs
	$(PATH_PERL)perl blast_filterIDs.pl -blast $< -ids $(word 3,$^) -output $@.tmp -type extract -seq subject
	$(PATH_PERL)perl blast_filterIDs.pl -blast $@.tmp -ids $(word 2,$^) -output $@ -type remove -seq query 

# filter contigs according to specified options
$(BLAST_DIR)/avg_mRNA_fragments_trin2ref.csv : $(BLAST_DIR)/avg_mRNA_unassigned_contigs_trin2ref.csv $(ASSEMBLY_DIR)/Trinity_preprocessed.fa.fai
	@echo $(DATE) Filter blast based on fragment criteria
	$(PATH_PERL)perl blast_filter.pl -blast $< -query-index $(word 2,$^) $(OPT_FRAGMENTS) | sort $(SORT_ORDER) -k1,1V -k2,2V -o $@

# best hit only
$(BLAST_DIR)/best_mRNA_fragments_trin2ref.csv : $(BLAST_DIR)/avg_mRNA_fragments_trin2ref.csv
	$(PATH_PERL)perl first_occurence.pl -i $< | sort -k2,2V > $@

# -----------------------------------------------------------------------------
# Annotation
# -----------------------------------------------------------------------------

$(BLAST_DIR)/avg_CDS_ref2trin_sym.csv : $(BLAST_DIR)/avg_CDS_ref2trin.csv $(TABLE_DIR)/reference.csv
	$(PATH_PERL)repCol.pl -table1 $< -table2 $(word 2,$^) -column 1 -tkey 1 -tvalue 5 -add > $@

$(OUTPUT_DIR)/annotation_sym.csv : $(TABLE_DIR)/annotation.csv $(TABLE_DIR)/reference.csv
	$(PATH_PERL)repCol.pl -table1 $< -table2 $(word 2,$^) -column 2 -tkey 1 -tvalue 5 -add > $@

$(ASSEMBLY_DIR)/annotated.gbk : $(ASSEMBLY_DIR)/Trinity_preprocessed.fa $(REF_TRANSCRIPTOME) $(TABLE_DIR)/annotation.csv $(TRANS_DIR)/.c $(TRINITY_DIR)/dir2comp.index $(BLAST_DIR)/avg_CDS_ref2trin_sym.csv $(BLAST_DIR)/best_mRNA_fragments_trin2ref.csv $(ORTHOLOG_TABLE) $(ORTHOLOG_FASTA)
	@echo "$(DATE) Annotation"
	@echo -e "$(DATE)\tLog: $(LOG_DIR)/annotation.log"
	@echo -e "$(DATE)\tError: $(LOG_DIR)/annotation.err"
	$(PATH_PERL)perl annotation.pl -new $(ANNOTATION_CUSTOM_PARAMETER) -taxon-assemble $(SPECIES_ASS_ID) -taxon-ortholog $(SPECIES_REF_ID) -trinity $< -reference $(word 2,$^) -assignment $(word 3,$^) -output-dir $(dir $(word 4,$^)) -read-index $(word 5,$^) -blast $(word 6,$^) -cpus $(OPT_CPUS) -scaffolding $(word 7,$^) -output-file $@ -genscan-matrix $(PATH_GENSCAN_MAT) $(OPT_PREDICTCDS) > $(LOG_DIR)/annotation.log 2> $(LOG_DIR)/annotation.err

# index path to read files
$(TRINITY_DIR)/dir2comp.csv : $(TRINITY_DIR)/.c
	@echo $(DATE) Indexing trinity readfiles.
	$(PATH_PERL)perl trinity_index_reads.pl -dir $(dir $<) -output $@

$(TRINITY_DIR)/dir2comp.index : $(TRINITY_DIR)/dir2comp.csv
	$(PATH_PERL)perl store_hash.pl -input $< -c1 1 -c2 2 -noarray -output $@

# -----------------------------------------------------------------------------
# Vector contamination (or simple include in annotation.pl?)
# -----------------------------------------------------------------------------

$(BLAST_DIR)/vecscreen.csv : $(OUTPUT_DIR)/sequences-cds.fa
	@echo $(DATE) Vector screening
	$(PATH_PERL)perl blast_VecScreen.pl -sequences $< -outstem $(basename $@) -cpu $(OPT_CPUS) -database $(DB_VEC)

$(OUTPUT_DIR)/transcript_catalogue.gbk : $(ASSEMBLY_DIR)/annotated.gbk $(BLAST_DIR)/vecscreen.csv
	$(PATH_PERL)perl contamination_feature.pl -input $< -output $@ -contamination $(lastword $^)

# -----------------------------------------------------------------------------
# FASTA Conversion
# -----------------------------------------------------------------------------

$(OUTPUT_DIR)/sequences-cds.fa : $(ASSEMBLY_DIR)/annotated.gbk
	@echo $(DATE) Extract CDS.
	$(PATH_PERL)perl genbank2fasta.pl -i $< -o $@ -feature CDS -feature-sequence

$(OUTPUT_DIR)/sequences-mRNA.fa : $(ASSEMBLY_DIR)/annotated.gbk
	@echo $(DATE) Extract mRNA.
	$(PATH_PERL)perl genbank2fasta.pl -i $< -o $@ -feature mRNA -feature-sequence

$(OUTPUT_DIR)/sequences-source.fa : $(ASSEMBLY_DIR)/annotated.gbk
	@echo $(DATE) Extract contigs/scaffolds.
	$(PATH_PERL)perl genbank2fasta.pl -i $< -o $@

$(OUTPUT_DIR)/sequences-prot.fa : $(ASSEMBLY_DIR)/annotated.gbk
	@echo $(DATE) Extract proteins.
	$(PATH_PERL)perl genbank2fasta.pl -i $< -o $@ -protein

# -----------------------------------------------------------------------------
# SUMMARIZE
# -----------------------------------------------------------------------------

$(OUTPUT_DIR)/summary.pdf : $(TABLE_DIR)/summary.csv $(TABLE_DIR)/reference.csv $(ASSEMBLY_DIR)/Trinity_preprocessed.fa.fai $(TABLE_DIR)/.c
	@echo $(DATE) Summarize results.
	@echo -e "$(DATE)\tLog: $(LOG_DIR)/summary.log"
	$(PATH_RSCRIPT)Rscript --vanilla summarize.R $(TABLE_DIR)/summary.csv $(TABLE_DIR)/reference.csv $(ASSEMBLY_DIR)/Trinity_preprocessed.fa.fai $(TABLE_DIR) $(OUTPUT_DIR)/summary.pdf > $(LOG_DIR)/summary.log

$(TABLE_DIR)/summary.csv : $(OUTPUT_DIR)/transcript_catalogue.gbk $(TABLE_DIR)/.c
	@echo $(DATE) Create summary table.
	$(PATH_PERL)perl genbank2table.pl -genbank $< -assembly -o $@

$(TABLE_DIR)/reference.csv : $(REF_TRANSCRIPTOME) $(TABLE_DIR)/.c
	$(PATH_PERL)perl genbank2table.pl -genbank $< -o $@

$(TABLE_DIR)/cds.csv : $(OUTPUT_DIR)/transcript_catalogue.gbk  $(TABLE_DIR)/.c
	@echo $(DATE) Extract CDS coordinates.
	$(PATH_PERL)perl genbank2table.pl -genbank $< -cds -o $@

# -----------------------------------------------------------------------------
# FUNCTIONAL ANNOTATION
# -----------------------------------------------------------------------------

$(OUTPUT_DIR)/summary_go.pdf : $(TABLE_DIR)/summary.csv 
	@echo $(DATE) GO annotation.
	@echo -e "$(DATE)\tLog: $(LOG_DIR)/GO.log"
	$(PATH_RSCRIPT)Rscript --vanilla GO.R $< $(OPT_ANNOTATION) $(TABLE_DIR) $@ 2> $(LOG_DIR)/GO.log

$(OUTPUT_DIR)/summary_kegg.pdf : $(TABLE_DIR)/summary.csv 
	@echo $(DATE) KEGG annotation.
	@echo -e "$(DATE)\tLog: $(LOG_DIR)/KEGG.log"
	$(PATH_RSCRIPT)Rscript --vanilla KEGG.R $< $(OPT_ANNOTATION) $(TABLE_DIR) $@ 2> $(LOG_DIR)/KEGG.log

# -----------------------------------------------------------------------------
# CLEANING
# -----------------------------------------------------------------------------

full-cleanup :
	@echo $(DATE) This might take a while.
	@rm -rf $(TRANS_DIR)
	@rm -rf $(BBH_DIR)
	@rm -rf $(LOG_DIR)
	@rm -rf $(TRINITY_DIR)
	@rm -rf $(BLAST_DIR)
	@rm -rf $(ASSEMBLY_DIR)
	@rm -rf $(BLASTDB_DIR)
	@rm -rf $(REF_DIR)
	@rm -rf $(MASK_DIR)
	@rm -rf $(OUTPUT_DIR)/reads
	@rm -rf $(LOG_DIR)

cleanup :
	@echo $(DATE) This might take a while.
	@rm -rf $(BBH_DIR)
	@rm -rf $(LOG_DIR)
	@rm -rf $(BLAST_DIR)
	@rm -rf $(ASSEMBLY_DIR)
	@rm -rf $(BLASTDB_DIR)
	@rm -rf $(REF_DIR)
	@rm -rf $(MASK_DIR)
	@rm -rf $(LOG_DIR)

.INTERMEDIATE: $(TMP_DIR)


