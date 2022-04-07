####call library####
library(DESeq2)
library(tidyverse)
library(RColorBrewer)
library(pheatmap)
library(tximport)
library(jsonlite)
library(ggplot2)
library(ggrepel)
library(readr)

####List all directories containing data####
samples <- list.files(full.names = T, pattern="quant$")

####Obtain a vector of all filenames including the path####
files <- file.path(samples, "quant.sf")

####call library####
library(DESeq2)
library(tidyverse)
library(RColorBrewer)
library(pheatmap)
library(tximport)
library(jsonlite)
library(ggplot2)
library(ggrepel)
library(readr)

####List all directories containing data####
samples <- list.files(full.names = T, pattern="quant$")

####Obtain a vector of all filenames including the path####
files <- file.path(samples, "quant.sf")

####Since all quant files have the same name it is useful to have names for each element####
names(files) <- c("invitro_rep1", "invitro_rep2","invitro_rep3",
                  "invivo_37C_rep1", "invivo_37C_rep2", "nvivo_37C_rep3",
                  "invivo_30C_rep1", "invivo_30C_rep2", "invivo_30C_rep3",
                  "invivo_47C_rep1", "invivo_47C_rep2", "invivo_47C_rep3",
                  "invivo_hs_rep1", "invivo_hs_rep2", "invivo_hs_rep3")

####Load the annotation table for Eco tRNA####
tx.exp <- tximport(files, type = "salmon", txOut = TRUE)
head(tx.exp$counts)

tx2gene <- data.frame(
  TXNAME = rownames(tx.exp$counts),
  GENEID = sapply(strsplit(rownames(tx.exp$counts), '\\.'), '[', 1)
)
head(tx2gene)

gene.exp <- summarizeToGene(tx.exp, tx2gene, countsFromAbundance = "scaledTPM")
head(gene.exp$counts)

write.csv(gene.exp$counts, file = "quant_Data_salmon.csv")
write.csv(gene.exp$abundance, file = "quant_abundance_salmon.csv")

table <- as.data.frame(gene.exp$counts)
write.table(table, file = "gene_exp_counts.txt", col.names = T, row.names = T, sep = "\t")

