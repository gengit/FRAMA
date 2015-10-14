#!/usr/bin/env Rscript --slave
#
# Description:
#
# Arguemnts:
#   summary.csv
#   reference.csv
#   Trinity_preprocessed.fa.fai
#   output directory
#   output pdf

# We could move this to init.pl, but as summary stuff is not an essential part,
# check for installed packages here.

# run as quietly as possible
#require(plyr, quietly = T, warn.conflicts = F)
#require(reshape, quietly = T, warn.conflicts = F)
#require(ggplot2, quietly = T, warn.conflicts = F)

packages = c("ggplot2", "plyr", "reshape", "gridExtra", "grid")
for (i in packages) {
    if(i %in% installed.packages() == FALSE) {
        # lib="../lib/R/packages?"
        print("Trying to install packages automatically")
        install.packages(i) 
    } else {
        require(i, quietly = T, warn.conflicts = F, character = T)
    }
}

sink("/dev/null")

N50  <-  function(data, percentage = 50) {
  percentage = percentage/100;
  total.bases = sum(data)
  sorted.contig.lengths = sort(data, dec= T)
  cumsum.contig.lengths = cumsum(sorted.contig.lengths)
  n50.index = which(cumsum.contig.lengths > (total.bases * percentage))[1]
  n50 = sorted.contig.lengths[n50.index]
  n50
}

assembly.stat <- function(names,lengths) {
  n50 = N50(lengths)

  numbers = as.character(c(
                 sum(lengths),
                 length(lengths), 
                 n50,
                 length(which(lengths >= n50)),
                 round(mean(lengths),3), 
                 round(sd(lengths),3),
                 round(sqrt(sum(lengths^2)/length(lengths)), 3),
                 min(lengths), 
                 max(lengths), 
                 sum(lengths < 300), 
                 sum(lengths > 1000), 
                 sum(lengths > 10000)
  ))
  labels = c(
             "total length",
             "# total contigs", 
             "N50",
             "# contigs >= N50",
             "mean length",
             "sd length",
             "RMS",
             "min length",
             "max length",
             "# contigs < 300bp",
             "# contigs > 1000bp",
             "# contigs > 10000bp"
  )

  df = data.frame(description = labels, number = numbers, stringsAsFactors = F)

  if (!is.null(names)) {
      numbers = length(unique(laply(as.character(names), function(x) { strsplit(x, '.', fixed = T)[[1]][1] })))
  } else {
      numbers = c(numbers, "NA")
  }
  df = rbind(df, c("graphs", numbers))
  df

}

arguments = commandArgs(T)

options(width=300)
uniqL <- function(x) {
  length(unique(as.character(x)))
}

# WORLD ORDER
order_prediction = c( "alignment", "genscan", "longest_orf")

order_cds = c(
 "full_length",
 "full_length;intron",
 "partial_3prime",
 "partial_3prime;intron",
 "partial_5prime",
 "partial_5prime;intron",
 "partial",
 "partial;intron",
 "full_length;pseudogene",
 "partial_3prime;pseudogene",
 "partial_5prime;pseudogene",
 "partial;pseudogene"
)

order_cds_short = c(
 "full_length(;intron|)",
 "3prime(;intron|)",
 "5prime(;intron|)",
 "partial(;intron|)$",
 "pseudogene",
 "intron"
 )

order_hit_type = c( "BBH", "SBH")


# transcript catalogue
#
# EXAMPLE:
#
#
#             assembly hit_type  ortholog symbol length cds_length num_fragments fragments     cds_status cds_prediction_tool clipped contamination
# 1 c10.graph_c0_seq10      BBH NM_005329   HAS3   1881       1734             0      <NA>    full_length           alignment       0          <NA>
# 2  c11.graph_c0_seq1      BBH NM_003218  TERF1    828        786             0      <NA> partial_3prime           alignment       0          <NA>
# 3   c9.graph_c0_seq1      BBH NM_005328   HAS2   4623       1659             0      <NA>    full_length           alignment       0          <NA>
#
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
)

# extra column with trinity ids
sum.table$original_name = gsub("_PART\\d", "", sum.table$assembly, perl = T)


# reference transcriptome
#
# EXAMPLE:
#
#   accession nucl_version length cds_length symbol geneid      ccds   protein prot_version     organism
# 1 NM_003218            3   2900       1260  TERF1   7013  CCDS6210 NP_003209            2 Homo sapiens
# 2 NM_005328            2   3275       1659   HAS2   3037  CCDS6335 NP_005319            1 Homo sapiens
# 3 NM_005329            2   4220       1662   HAS3   3038 CCDS10871 NP_005320            2 Homo sapiens
# 4 NM_001523            2   2116       1737   HAS1   3036 CCDS12838 NP_001514            2 Homo sapiens
#
reference.table = read.table(
  arguments[2],
  header = T,
  sep ="\t",
  colClasses = c("character", #accession
                 "integer",   #nucl_version
                 "integer",   #length
                 "integer",   #cds_length
                 "character", #symbol
                 "character",   #geneid
                 "character", #ccds
                 "character", #protein
                 "integer",   #prot_version
                 "character"  #organism
 ),
 na.strings = "NA"
)[, c("accession", "length", "cds_length", "geneid", "protein")]

colnames(reference.table) = c("accession", "ortholog_length", "ortholog_cds_length", "geneid", "protein")

# processed assembly
#
# EXAMPLE:
# 
#             contig length
# 1 c2.graph_c0_seq1    236
# 2 c3.graph_c0_seq1    221
# 3 c5.graph_c0_seq1   1107
# 4 c5.graph_c0_seq2    310
# 5 c6.graph_c0_seq1    620
# 6 c6.graph_c0_seq2    522
#
contig.table = read.table(
  arguments[3],
  sep = "\t", 
  header = F 
)[,c(1,2)]
colnames(contig.table) = c("contig", "contig_length")

sum.table = merge(reference.table, sum.table, all.y = T, by.x = "accession", by.y = "ortholog")
sum.table = merge(contig.table, sum.table, all.y = T, by.x = "contig", by.y = "original_name")

# adding "cds_reconstruction" column
sum.table$cds_reconstruction = round(((sum.table$cds_length / sum.table$ortholog_cds_length) * 100), 3)
# possible that CDS is longer than orthologs length
sum.table$cds_reconstruction[sum.table$cds_reconstruction > 100] = 100.000 

# best CDS based on order_cds
sum.table = ldply(order_cds, function(x) { 
    sub.t = sum.table[which(sum.table$cds_status == x),] 
    sub.t[order(sub.t$cds_length, decreasing = T), ]
})
sum.table.bestCDS = sum.table[!duplicated(as.character(sum.table$symbol)), ]
stopifnot(nrow(sum.table.bestCDS) == length(unique(sum.table$symbol))) # sanity

# SAVE SUMMARY/SUMMARY.BESY
write.table(sum.table, file.path(arguments[4], "summary_detailed.csv"), sep = "\t", quote = F, row.names = F)
write.table(sum.table.bestCDS, file.path(arguments[4], "summary_best_detailed.csv"), sep = "\t", quote = F, row.names = F)

######################################################################
# Transcript catalogue overview
######################################################################

final.df = data.frame(description = rep(NA, 25), transcripts = rep(NA, 25), genes = rep(NA, 25), stringsAsFactors = F)

final.df[1,] = (c("No. of transcripts/genes",nrow(sum.table), nrow(sum.table.bestCDS)))

# SCAFFOLDS
index = which(sum.table$num_fragments > 0)
if (length(index) > 0) {
    final.df[2, ] = c("No. of scaffolds", sum(sum.table$num_fragments > 0), sum(sum.table.bestCDS$num_fragments > 0))
} else {
    final.df[2, ] = c("No. of scaffolds", 0, 0)
}

# longest/shortest transcript contig (not scaffolded)
final.df[3,] = c("Max contig length",
               sum.table[order(sum.table$contig_length, sum.table$symbol, decreasing = T), ][1, c("contig_length")],
               sum.table.bestCDS[order(sum.table.bestCDS$contig_length, sum.table.bestCDS$symbol, decreasing = T), ][1, c("contig_length")]
)
final.df[4,] = c("Min contig length",
               sum.table[order(sum.table$contig_length, sum.table$symbol, decreasing = F), ][1, c("contig_length")],
               sum.table.bestCDS[order(sum.table.bestCDS$contig_length, sum.table.bestCDS$symbol, decreasing = F), ][1, c("contig_length")]
)
final.df[5,] = c("Max scaffold length",
               sum.table[order(sum.table$length, sum.table$symbol, decreasing = T), ][1, c("length")],
               sum.table.bestCDS[order(sum.table.bestCDS$length, sum.table.bestCDS$symbol, decreasing = T), ][1, c("length")]
)
final.df[6,] = c("Min scaffold length",
               sum.table[order(sum.table$length, sum.table$symbol, decreasing = F), ][1, c("length")],
               sum.table.bestCDS[order(sum.table.bestCDS$length, sum.table.bestCDS$symbol, decreasing = F), ][1, c("length")]
)

final.df[7, ] = c("CDS reconstructed > 90%" ,
               sum(sum.table$cds_reconstruction > 90),
               sum(sum.table.bestCDS$cds_reconstruction > 90)
)
final.df[8, ] = c("CDS reconstructed > 99%" ,
               sum(sum.table$cds_reconstruction > 99),
               sum(sum.table.bestCDS$cds_reconstruction > 99)
)

# completeness
a = ldply(order_cds, function(possibleValue) { c(sum(sum.table$cds_status == possibleValue), sum(sum.table.bestCDS$cds_status == possibleValue)) })
colnames(a) = c("transcripts", "genes")

a$description = order_cds

final.df[9:20, ] = a[, c("description", "transcripts", "genes")]

# prediction tools
a = ldply(order_prediction, function(possibleValue) { c(sum(sum.table$cds_prediction_tool == possibleValue), sum(sum.table.bestCDS$cds_prediction == possibleValue)) })
colnames(a) = c("transcripts", "genes")
a$description = order_prediction

final.df[21:23, ] = a[, c("description", "transcripts", "genes")]


# fusion
final.df[24,] = c("fusion contigs", sum(duplicated(as.character(sum.table$contig))), sum(duplicated(as.character(sum.table.bestCDS$contig))))

# clipping
final.df[25,] = c("clipped contigs", sum(sum.table$clipped > 0),sum(sum.table.bestCDS$clipped > 0))

write.table(final.df, file.path(arguments[4], "overview.csv"), sep = "\t", quote = F, row.names = F)

######################################################################
# Assembly statistic
######################################################################

a = assembly.stat(as.character(sum.table$contig),sum.table$length)
b = assembly.stat(as.character(sum.table$ortholog),sum.table$ortholog_length)
c = assembly.stat(as.character(contig.table$contig),contig.table$contig_length)

final.df1 = data.frame(description = a$description, total.assembly = c$number, transcript.catalogue  = a$number, orthologous.transcripts = b$number)
write.table(final.df1, file.path(arguments[4], "transcriptome_statistic.csv"), sep = "\t", quote = F, row.names = F)

######################################################################
######################################################################
# Graphics
######################################################################
######################################################################

pdf(file = arguments[5], width = 20, height = 8)

final.df$transcripts = as.numeric(final.df$transcripts)
final.df$genes = as.numeric(final.df$genes)
#final.df1$total.assembly = as.numeric(final.df1$total.assembly)
#final.df1$transcript.catalogue = as.numeric(final.df1$transcript.catalogue)
#final.df1$orthologous.transcripts = as.numeric(final.df1$orthologous.transcripts)

tG1<-tableGrob(
              format(final.df, big.mark = ","),
              #core.just="right",
              #col.just="right",
              #show.rownames = F,
              #h.even.alpha = 0,
              #gpar.rowtext = gpar(col="black", equal.width = TRUE, show.vlines = F, show.hlines = F, separator="grey")
)
tG2<-tableGrob(
              format(final.df1, big.mark = ","),
              #core.just="right",
              #col.just="right",
              #show.rownames = F,
              #h.even.alpha = 0,
              #gpar.rowtext = gpar(col="black", equal.width = TRUE, show.vlines = F, show.hlines = F, separator="grey")
)

grid.arrange(tG2, tG1, ncol = 2)

######################################################################
# Pie chart with hit type
######################################################################
par(mfrow=c(1,2))
comp = laply(order_hit_type, function(possibleValue) { sum(sum.table$hit_type == possibleValue) })
names(comp) = order_hit_type
labels = paste(names(comp), paste0("(", comp, ")"), sep = "\n")
a <- pie(comp, labels = labels, main = "hit type (transcript)")

comp = laply(order_hit_type, function(possibleValue) { sum(sum.table.bestCDS$hit_type == possibleValue) })
names(comp) = order_hit_type
labels = paste(names(comp), paste0("(", comp, ")"), sep = "\n")
b <- pie(comp, labels = labels, main = "hit type (best transcript per gene)")

######################################################################
# Pie chart with CDS completeness
######################################################################
par(mfrow=c(1,2))

comp = laply(order_cds, function(possibleValue) { sum(sum.table$cds_status == possibleValue) })
names(comp) = order_cds
completeness = c(sum(comp[c(1,2,9)]), sum(comp[c(4,5,10)]), sum(comp[c(7,8,11)]), sum(comp[c(9,10,12)]))
names(completeness) = c("full_length", "partial 3'", "partial 5'", "partial")
labels = paste(names(completeness), paste0("(", completeness, ")"), sep = "\n")
pie(completeness, labels = labels, main = "CDS completeness (transcript)")

comp = laply(order_cds, function(possibleValue) { sum(sum.table.bestCDS$cds_status == possibleValue) })
names(comp) = order_cds
completeness = c(sum(comp[c(1,2,9)]), sum(comp[c(4,5,10)]), sum(comp[c(7,8,11)]), sum(comp[c(9,10,12)]))
names(completeness) = c("full_length", "partial 3'", "partial 5'", "partial")
labels = paste(names(completeness), paste0("(", completeness, ")"), sep = "\n")
pie(completeness, labels = labels, main = "CDS completeness (best transcript per gene)")

######################################################################
# Pie chart with CDS prediction tool
######################################################################

par(mfrow=c(1,2))
comp = laply(order_prediction, function(possibleValue) { sum(as.character(sum.table$cds_prediction_tool) == possibleValue) })
names(comp) = order_prediction
labels = paste(names(comp), paste0("(", comp, ")"), sep = "\n")
pie(comp, labels = labels, main = "CDS inference method (transcript)")

comp = laply(order_prediction, function(possibleValue) { sum(sum.table.bestCDS$cds_prediction_tool == possibleValue) })
names(comp) = order_prediction
labels = paste(names(comp), paste0("(", comp, ")"), sep = "\n")
pie(comp, labels = labels, main = "CDS inference method (best transcript per gene)")


######################################################################
# Histogram contig length
######################################################################

print(ggplot(sum.table, aes(x=contig_length)) +
geom_bar(stat = "bin", binwidth = 100) +
scale_x_continuous(breaks = round(seq(200, max(sum.table$contig_length), by =300),1)) +
theme_bw() +
xlab("length in bp") +
ylab("No. of contigs") +
ggtitle("contig length histogram (transcript)") +
theme(axis.title.y=element_text(vjust=1.5)) + 
theme(axis.text.x=element_text(angle = 90, vjust=0.5))
)

print(ggplot(sum.table.bestCDS, aes(x=contig_length)) +
geom_bar(stat = "bin", binwidth = 100) +
scale_x_continuous(breaks = round(seq(200, max(sum.table$contig_length), by =300),1)) +
theme_bw() +
xlab("length in bp") +
ylab("No. of contigs") +
ggtitle("contig length histogram (best transcript per gene)") +
theme(axis.title.y=element_text(vjust=1.5)) + 
theme(axis.text.x=element_text(angle = 90, vjust=0.5))
)

######################################################################
# Histogram CDS length
######################################################################

print(ggplot(sum.table, aes(x=cds_length)) +
geom_bar(stat = "bin", binwidth = 100) +
scale_x_continuous(breaks = round(seq(200, max(sum.table$contig_length), by =300),1)) +
theme_bw() +
xlab("CDS length") +
ylab("No. of contigs") +
ggtitle("CDS length histogram (transcript)") +
theme(axis.title.y=element_text(vjust=1.5)) + 
theme(axis.text.x=element_text(angle = 90, vjust=0.5))
)

print(ggplot(sum.table.bestCDS, aes(x=cds_length)) +
geom_bar(stat = "bin", binwidth = 100) +
scale_x_continuous(breaks = round(seq(200, max(sum.table$contig_length), by =300),1)) +
theme_bw() +
xlab("CDS length") +
ylab("No. of contigs") +
ggtitle("CDS length histogram (best transcript per gene)") +
theme(axis.title.y=element_text(vjust=1.5)) + 
theme(axis.text.x=element_text(angle = 90, vjust=0.5))
)


######################################################################
# Histogram CDS reconstruction
######################################################################

print(ggplot(sum.table, aes(x=cds_reconstruction)) +
stat_bin(aes(y=(..count../sum(..count..))*100), binwidth = 1)  +
stat_bin(aes(y=..count../sum(..count..)), binwidth = 1)  +
scale_fill_brewer(palette="Spectral") +
scale_x_continuous(breaks = seq(0, 100, by = 5)) +
scale_y_continuous(breaks = seq(0, 100, by = 5)) +
theme_bw() +
theme(axis.line = element_line(colour = "black")) +
ylab("proportion of genes (%)") +
xlab("reconstructed CDS length in comparison to reference transcript (%)") +
ggtitle("CDS reconstruction (transcript)")
)

print(ggplot(sum.table.bestCDS, aes(x=cds_reconstruction)) +
      stat_bin(aes(y=(..count../sum(..count..))*100), binwidth = 1)  +
      stat_bin(aes(y=..count../sum(..count..)), binwidth = 1)  +
      scale_fill_brewer(palette="Spectral") +
      scale_x_continuous(breaks = seq(0, 100, by = 5)) +
      scale_y_continuous(breaks = seq(0, 100, by = 5)) +
      theme_bw() +
      theme(axis.line = element_line(colour = "black")) +
      ylab("proportion of genes (%)") +
      xlab("reconstructed CDS length in comparison to reference transcript (%)") + 
      ggtitle("CDS reconstruction (best transcript per gene)")
      )

######################################################################
# Histogram length improvement by scaffolding
######################################################################

scaffold.table = sum.table[which(sum.table$num_fragments > 0),]
if (nrow(scaffold.table) > 0) {
    scaffold.table$length_improvement = (scaffold.table$contig_length / scaffold.table$length) * 100
    scaffold.table$length_improvement = round((scaffold.table$length/scaffold.table$contig_length), 1)

    print(ggplot(scaffold.table, aes(x=length_improvement)) +
          geom_bar(stat = "bin", binwidth = 0.1) +
          scale_fill_brewer(palette="Spectral") +
          scale_x_continuous(breaks = round(seq(1, max(scaffold.table$length_improvement), by = 0.2),1)) +
          theme_bw() +
          xlab("length increase") +
          ylab("No. of contigs") +
          ggtitle("Scaffolding (transcript)") +
          theme(axis.title.y=element_text(vjust=1.5)) + 
          theme(axis.text.x=element_text(angle = 90, vjust=0.5))
          )
}

scaffold.table.bestCDS = sum.table.bestCDS[which(sum.table.bestCDS$num_fragments > 0),]
if (nrow(scaffold.table.bestCDS) > 0) {
    scaffold.table.bestCDS$length_improvement = (scaffold.table.bestCDS$contig_length / scaffold.table.bestCDS$length) * 100
    scaffold.table.bestCDS$length_improvement = round((scaffold.table.bestCDS$length/scaffold.table.bestCDS$contig_length), 1)

    print(ggplot(scaffold.table.bestCDS, aes(x=length_improvement)) +
          geom_bar(stat = "bin", binwidth = 0.1) +
          scale_fill_brewer(palette="Spectral") +
          scale_x_continuous(breaks = round(seq(1, max(scaffold.table.bestCDS$length_improvement), by = 0.2),1)) +
          theme_bw() +
          xlab("length increase") +
          ylab("No. of contigs") +
          ggtitle("Scaffolding (best transcript per gene)") +
          theme(axis.title.y=element_text(vjust=1.5)) + 
          theme(axis.text.x=element_text(angle = 90, vjust=0.5))
          )
}

######################################################################
# Histogram length clipping
######################################################################

sum.table.clipped = sum.table[sum.table$clipped > 0, ]
if (nrow(sum.table.clipped) > 0) {
    print(ggplot(sum.table.clipped, aes(x=clipped)) +           
          geom_bar(stat = "bin", binwidth = 50) +               
          scale_x_continuous(breaks = seq(0, max(sum.table.clipped$clipped), by = 200)) +
          theme_bw() +                                          
          theme(axis.line = element_line(colour = "black")) +
          theme(axis.text.x=element_text(angle=90, vjust = 0.5)) +
          xlab("number of clipped bases") +                     
          ylab("number of contigs") +                           
          ggtitle("3' UTR clipping (transcript)")               
          )
}

sum.table.bestCDS.clipped = sum.table.bestCDS[sum.table.bestCDS$clipped > 0, ]
if (nrow(sum.table.bestCDS.clipped) > 0) {
    print(ggplot(sum.table.bestCDS.clipped, aes(x=clipped)) +           
          geom_bar(stat = "bin", binwidth = 50) +               
          scale_x_continuous(breaks = seq(0, max(sum.table.clipped$clipped), by = 200)) +
          theme_bw() +                                          
          theme(axis.line = element_line(colour = "black")) +
          theme(axis.text.x=element_text(angle=90, vjust = 0.5)) +
          xlab("number of clipped bases") +                     
          ylab("number of contigs") +                           
          ggtitle("3' UTR clipping (best transcript per gene)")               
          )
}

dev.off()

