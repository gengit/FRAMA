# get functional annotation

# description
#   summary.csv
#   annotation packages
#   output directory for tables
#   output file

sink("/dev/null")

arguments = commandArgs(T)

require(GO.db,quietly = T, warn.conflicts = F)
require(plyr,quietly = T, warn.conflicts = F)
require(ggplot2,quietly = T, warn.conflicts = F)
require(arguments[2], character.only = T,quietly = T, warn.conflicts = F)

packages = ls(paste0("package:", arguments[2]))

sum.table = read.table(
  arguments[1],
  header = T,
  sep ="\t",
  colClasses = c("character",     #assembly
                 "character",     #hit_type
                 "character",     #ortholog
                 "character",     #symbol
                 "integer",       #length
                 "integer",       #cds_length
                 "integer",       #num_fragments
                 "character",     #fragments
                 "character",     #cds_status
                 "character",     #cds_prediction_tool
                 "integer",       #clipped
                 "character"      #contamination
  ),
  na.string = "NA"
)[, c("assembly", "ortholog")]

######################################################################
# annotate each contig with GO ID if any
######################################################################

# geneid, go_id, evidence, ontology
go.frame = toTable(get(packages[grep("GO$",packages)]))
# geneid, accession
refseq.frame = toTable(get(packages[grep("REFSEQ$",packages)]))
# go_id, Term, Ontology, Definition, Synonym
go.desc = toTable(GOTERM)
go.desc = go.desc[,c(-1, -4, -6, -7)]
go.desc = go.desc[!duplicated(go.desc$go_id), ]
# add geneid
sum.table = merge(refseq.frame, sum.table, by.x ="accession", by.y = "ortholog", all = F)
# add go term
sum.table = merge(go.frame, sum.table, by.x ="gene_id", by.y = "gene_id", all = F)
if (nrow(sum.table) == 0) {
    # there was nothing wrong, but no annotated genes have been assembled
    file.create(arguments[4])
    sink(NULL)
    print("No genes assembled with annotated GO Term.")
    quit()
    
}
# so no annotated genes found...
# add description
sum.table = merge(go.desc, sum.table, by.x ="go_id", by.y = "go_id", all = F)
write.table(sum.table, file.path(arguments[3], "gene_ontology.csv"), row.names = F, quote = F, sep = "\t")

######################################################################
# number of genes in each path (for assembly and annotation)
######################################################################

# genes per GO in annotation
annotation.genes.per.goid = data.frame(table(go.frame$go_id))
colnames(annotation.genes.per.goid) = c("go_id", "genes")

tmp = sum.table[,c("go_id", "gene_id")]
tmp = tmp[!duplicated(tmp), ]
assembly.genes.per.goid =  data.frame(table(tmp$go_id))
colnames(assembly.genes.per.goid) = c("go_id", "genes")

genes.per.goid = merge(annotation.genes.per.goid, assembly.genes.per.goid, all.x = T, by = "go_id")
colnames(genes.per.goid) = c("go_id", "annotation", "assembly")
genes.per.goid$proportion = round((genes.per.goid$assembly / genes.per.goid$annotation) * 100, 2)


# add term name
genes.per.goid = merge(genes.per.goid, go.desc, all.x=T, all.y=F)
write.table(genes.per.goid[,c("go_id", "annotation", "assembly", "proportion", "Term")], file.path(arguments[3], "gene_ontology_genes_per_path.csv"), row.names = F, quote = F, sep = "\t")

######################################################################
# total GOIDs and total genes (in assembly and annotation)
######################################################################

ontologies = c("BP","MF","CC")
names(ontologies) = ontologies

annotation.summary = ldply(ontologies, function(x) {
  subtable = subset(go.frame, Ontology == x)
  if (nrow(subtable) > 0) {
      c( length(unique(as.character(subtable$go_id))), length(unique(as.character(subtable$gene_id))))
  } else {
      c ( 0,0)
  }
})
colnames(annotation.summary) = c("ontology", "go_ids", "gene_ids")

assembly.summary = ldply(ontologies, function(x) {
  subtable = subset(sum.table, Ontology == x)
  if (nrow(subtable) > 0) {
      c( length(unique(as.character(subtable$go_id))), length(unique(as.character(subtable$gene_id))))
  } else {
      c ( 0,0)
  }
  #c(length(unique(as.character(subtable$go_id))),
  #length(unique(as.character(subtable$gene_id))))
})
colnames(assembly.summary) = c("ontology", "go_ids", "gene_ids")

df = merge(assembly.summary, annotation.summary, by = "ontology")
colnames(df) = c("ontology", "assembly.goids", "assembly.genes", "annotation.goids", "annotation.genes")
df$proportion_goids = round((df$assembly.goids / df$annotation.goids) * 100, 2)
df$proportion_genes = round((df$assembly.genes / df$annotation.genes) * 100, 2)
write.table(df, file.path(arguments[3], "gene_ontology_genes_per_ontology.csv"), sep = "\t", quote = F, row.names = F)

######################################################################
# Histogram with covered genes for each go term
######################################################################
pdf(file = arguments[4], width = 20, height = 8)
print(ggplot(genes.per.goid, aes(x = proportion)) +
geom_bar(binwidth = 5) +
theme_bw() +
xlab("percentage of genes found") +
ylab("number of GO Terms") +
scale_x_continuous(breaks = seq(0,100, by = 5)) +
ggtitle("GO Terms"))
dev.off()
