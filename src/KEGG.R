# get functional annotation

# description
#   summary.csv
#   annotation packages
#   output directory for tables
#   output file

sink("/dev/null")

arguments = commandArgs(T)

require(plyr,quietly = T, warn.conflicts = F)
require(ggplot2,quietly = T, warn.conflicts = F)
require(KEGG.db,quietly = T, warn.conflicts = F)
require(arguments[2], character.only = T,quietly = T, warn.conflicts = F)

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

packages = ls(paste("package:", arguments[2], sep = ""))

paths = toTable(get(packages[grep("PATH$",packages)]))
refseqs = toTable(get(packages[grep("REFSEQ$",packages)]))
path.name = toTable(KEGGPATHID2NAME)

sum.table = merge(refseqs, sum.table, by.x ="accession", by.y = "ortholog", all = F)
sum.table = merge(paths, sum.table, by.x ="gene_id", by.y = "gene_id", all = F)
if (nrow(sum.table) == 0) {
    # there was nothing wrong, but no annotated genes have been assembled
    file.create(arguments[4])
    sink(NULL)
    print("No genes assembled with annotated KEGG Pathway.")
    quit()
    
}
sum.table = merge(path.name, sum.table)

write.table(sum.table,file.path(arguments[3], "kegg.csv"), row.names = F, col.names =T, quote = F, sep = "\t")

# genes per GO in annotation
annotation.genes.per.pathid = data.frame(table(paths$path_id))
colnames(annotation.genes.per.pathid) = c("path_id", "genes")

tmp = sum.table[,c("path_id", "gene_id")]
tmp = tmp[!duplicated(tmp), ]
assembly.genes.per.pathid =  data.frame(table(tmp$path_id))
colnames(assembly.genes.per.pathid) = c("path_id", "genes")

genes.per.pathid = merge(annotation.genes.per.pathid, assembly.genes.per.pathid, all.x = T, by = "path_id")
colnames(genes.per.pathid) = c("path_id", "annotation", "assembly")
genes.per.pathid$proportion = round((genes.per.pathid$assembly / genes.per.pathid$annotation) * 100, 2)

# add term name
genes.per.pathid = merge(genes.per.pathid, path.name, all.x=T, all.y=F)
write.table(genes.per.pathid, file.path(arguments[3], "kegg_genes_per_path.csv"), row.names = F, quote = F, sep = "\t")

assembly.paths = length(unique(as.character(sum.table$path_id)))
assembly.genes = length(unique(as.character(sum.table$gene_id)))
annotation.paths = length(unique(as.character(paths$path_id)))
annotation.genes = length(unique(as.character(paths$gene_id)))

df = data.frame(assembly.paths = assembly.paths, assembly.genes = assembly.genes, annotation.paths = annotation.paths, annotation.genes = annotation.genes)
df$proportion.paths = round((df$assembly.paths / df$annotation.paths ) * 100, 2)
df$proportion.genes = round((df$assembly.genes / df$annotation.genes ) * 100, 2)

write.table(df, file.path(arguments[3], "kegg_covered.csv"), sep = "\t", row.names = F, quote =F)

######################################################################
# Histogram with covered genes for each go term
######################################################################
pdf(file = arguments[4], width = 20, height = 8)
print(ggplot(genes.per.pathid, aes(x = proportion)) +
geom_bar(binwidth = 5) +
theme_bw() +
xlab("percentage of genes found") +
ylab("number of KEGG pathway") +
scale_x_continuous(breaks = seq(0,100, by = 5)) +
ggtitle("KEGG Pathways coverage"))
dev.off()



