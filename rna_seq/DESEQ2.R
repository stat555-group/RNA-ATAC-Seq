####call library####
library(DESeq2)
library(tidyverse)
library(RColorBrewer)
library(tximport)
library(jsonlite)
library(ggplot2)
library(readr)
library(ggrepel)
library(extrafont)

####read salmon output####
data <- read.csv("quant_Data_salmon_37C47C.csv", sep = ",", header = T, row.names = 1)
data2 <- read.csv("quant_Data_salmon_30CHS.csv", sep = ",", header = T, row.names = 1)
data3 <- read.csv("quant_Data_salmon_47Chs.csv", sep = ",", header = T, row.names = 1)

####Assign experimental variables####
dt <- as.matrix(data) 
dt1 <- dt[apply(dt,1,sum)>6,]
group = data.frame(con = factor(c("37C", "37C",
                                  "37C", "47C", "47C",
                                  "47C")))

dt2 <- as.matrix(data2)
dt3 <- dt2[apply(dt2,1,sum)>6,]
group2 = data.frame(con2 = factor(c("30C", "30C",
                                    "30C", "HS", "HS",
                                    "HS")))

dt4 <- as.matrix(data3)
dt5 <- dt4[apply(dt2,1,sum)>6,]
group3 = data.frame(con3 = factor(c("47C","47C","47C","HS","HS","HS")))

####Run DEseq2####
dds <- DESeqDataSetFromMatrix(countData = round(dt1), colData = group, design = ~con)
dds <- DESeq(dds)
res <- results(dds)
res$gene_id <- row.names(res)
write.table(res,file="DESeq2_result_37C47C.txt",sep="\t",row.names=F,col.names=T,quote=F)

dds2 <- DESeqDataSetFromMatrix(countData = round(dt3), colData = group2, design = ~con2)
dds2 <- DESeq(dds2)
res2 <- results(dds2)
res2$gene_id <- row.names(res2)
write.table(res2,file="DESeq2_result_30CHS.txt",sep="\t",row.names=F,col.names=T,quote=F)

dds3 <- DESeqDataSetFromMatrix(countData = round(dt5), colData = group3, design = ~con3)
dds3 <- DESeq(dds3)
res3 <- results(dds3)
res3$gene_id <- row.names(res3)
write.table(res3,file="DESeq2_result_47CHS.txt",sep="\t",row.names=F,col.names=T,quote=F)


####ggplot for 37C47C####
DESeq_37C47C <- read.csv("DESeq2_result_37C47C.txt", sep = "\t", header = T)

# The significantly differentially expressed genes are the ones found in the upper-left and upper-right corners.
# Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2FoldChange respectively positive or negative)
# add a column of NAs
DESeq_37C47C$Expression <- "NO"
# if log2Foldchange > 1 and pvalue < 0.05, set as "UP" 
DESeq_37C47C$Expression[DESeq_37C47C$log2FoldChange > 1 & DESeq_37C47C$pvalue < 0.05] <- "UP"
# if log2Foldchange < -1 and pvalue < 0.05, set as "DOWN"
DESeq_37C47C$Expression[DESeq_37C47C$log2FoldChange < -1 & DESeq_37C47C$pvalue < 0.05] <- "DOWN"
# add a column of NAs
DESeq_37C47C$label <- NA
# Create a new column "label" to DESeq_37C47, that will contain the name of genes differentially expressed (NA in case they are not)
DESeq_37C47C$label[DESeq_37C47C$Expression != "NO"] <- DESeq_37C47C$gene_id[DESeq_37C47C$Expression != "NO"]

#volcano plot
volcanoplot <- ggplot(DESeq_37C47C, aes(x=log2FoldChange, y=-log10(pvalue), col = Expression, label = label))+
  geom_point()+
  theme_classic()+
  theme(text = element_text(size = 20, family = "Arial", color = "black"))+
  geom_vline(xintercept=c(-1, 1), col="grey", linetype = "dashed") +
  geom_hline(yintercept=-log10(0.05), col="grey", linetype = "dashed")+
  scale_color_manual(values = c("blue", "grey", "red"))+
  xlim(-4, 4)+
  geom_text_repel(aes(fontface = "bold"))+
  theme(legend.position = "top")+
  labs(x=expression(log[2](Fold~change)), y=expression(-log[10](p-value)))
volcanoplot

####ggplot for 30CHS####
DESeq_30CHS <- read.csv("DESeq2_result_30CHS.txt", sep = "\t", header = T)

# The significantly differentially expressed genes are the ones found in the upper-left and upper-right corners.
# Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2FoldChange respectively positive or negative)
# add a column of NAs
DESeq_30CHS$Expression <- "NO"
# if log2Foldchange > 1 and pvalue < 0.05, set as "UP" 
DESeq_30CHS$Expression[DESeq_30CHS$log2FoldChange > 1 & DESeq_30CHS$pvalue < 0.05] <- "UP"
# if log2Foldchange < -1 and pvalue < 0.05, set as "DOWN"
DESeq_30CHS$Expression[DESeq_30CHS$log2FoldChange < -1 & DESeq_30CHS$pvalue < 0.05] <- "DOWN"
# add a column of NAs
DESeq_30CHS$label <- NA
# Create a new column "label" to DESeq_37C47, that will contain the name of genes differentially expressed (NA in case they are not)
DESeq_30CHS$label[DESeq_30CHS$Expression != "NO"] <- DESeq_30CHS$gene_id[DESeq_30CHS$Expression != "NO"]

#volcano plot
volcanoplot2 <- ggplot(DESeq_30CHS, aes(x=log2FoldChange, y=-log10(pvalue), col = Expression, label = label))+
  geom_point()+
  theme_classic()+
  theme(text = element_text(size = 20, family = "Arial", color = "black"))+
  geom_vline(xintercept=c(-1, 1), col="grey", linetype = "dashed") +
  geom_hline(yintercept=-log10(0.05), col="grey", linetype = "dashed")+
  scale_color_manual(values = c("blue", "grey", "red"))+
  xlim(-4, 4)+
  geom_text_repel(aes(fontface = "bold"))+
  theme(legend.position = "top")+
  labs(x=expression(log[2](Fold~change)), y=expression(-log[10](p-value)))
volcanoplot2

####ggplot for 47CHS####
DESeq_47CHS <- read.csv("DESeq2_result_47CHS.txt", sep = "\t", header = T)

# The significantly differentially expressed genes are the ones found in the upper-left and upper-right corners.
# Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2FoldChange respectively positive or negative)
# add a column of NAs
DESeq_47CHS$Expression <- "NO"
# if log2Foldchange > 1 and pvalue < 0.05, set as "UP" 
DESeq_47CHS$Expression[DESeq_47CHS$log2FoldChange > 1 & DESeq_47CHS$pvalue < 0.05] <- "UP"
# if log2Foldchange < -1 and pvalue < 0.05, set as "DOWN"
DESeq_47CHS$Expression[DESeq_47CHS$log2FoldChange < -1 & DESeq_47CHS$pvalue < 0.05] <- "DOWN"
# add a column of NAs
DESeq_47CHS$label <- NA
# Create a new column "label" to DESeq_47CHS, that will contain the name of genes differentially expressed (NA in case they are not)
DESeq_47CHS$label[DESeq_47CHS$Expression != "NO"] <- DESeq_47CHS$gene_id[DESeq_47CHS$Expression != "NO"]

#volcano plot
volcanoplot3 <- ggplot(DESeq_47CHS, aes(x=log2FoldChange, y=-log10(pvalue), col = Expression, label = label))+
  geom_point()+
  theme_classic()+
  theme(text = element_text(size = 20, family = "Arial", color = "black"))+
  geom_vline(xintercept=c(-1, 1), col="grey", linetype = "dashed") +
  geom_hline(yintercept=-log10(0.05), col="grey", linetype = "dashed")+
  scale_color_manual(values = c("blue", "grey", "red"))+
  xlim(-4, 4)+
  geom_text_repel(aes(fontface = "bold"))+
  theme(legend.position = "top")+
  labs(x=expression(log[2](Fold~change)), y=expression(-log[10](p-value)))
volcanoplot3




#### ggsave ####
library(svglite)
#ggsave(filename = 'DGE_37C_37C.png', path = getwd(), plot = volcanoplot, scale = 2.5, width = 5, height = 3.5, units = "cm", dpi = 300)
#ggsave(filename = 'DGE_30CHS.png', path = getwd(), plot = volcanoplot2, scale = 2.5, width = 5, height = 3.5, units = "cm", dpi = 300)
ggsave(filename = 'DGE_47CHS.png', path = getwd(), plot = volcanoplot3, scale = 2.5, width = 5, height = 3.5, units = "cm", dpi = 300)
#ggsave(filename = 'DGE_37C_37C.svg', path = getwd(), plot = volcanoplot, scale = 2.5, width = 5, height =3.5, units = "cm", dpi = 300)
#ggsave(filename = 'DGE_30CHS.svg', path = getwd(), plot = volcanoplot2, scale = 2.5, width = 5, height = 3.5, units = "cm", dpi = 300)



