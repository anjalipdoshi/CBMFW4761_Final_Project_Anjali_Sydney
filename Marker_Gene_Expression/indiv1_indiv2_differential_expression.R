## This code runs edgeR to analyze differential expression of the 10 marker genes in two individuals. 

library('edgeR')
library('dplyr')
library('ggplot2')
library('ggrepel')
library('qvalue')
library("gridExtra")
library("ggthemes")
library(RColorBrewer)

setwd("/Users/sydneyblattman/Dropbox/ncbi_data/Project_Code/Marker_Gene_Expression/")

data1 <- read.delim("SRR769418_normalized_matrix.txt",sep="\t")
data2 <- read.delim("SRR769429_normalized_matrix.txt",sep="\t")
rownames(data1) <- data1$ID
rownames(data2) <- data2$ID
data1 <- subset(data1, select = -c(ID,phylum))
data2 <- subset(data2, select = -c(ID,phylum))

data1_t = as.data.frame(t(data1))
sub_data1_t = as.data.frame(data1_t$meta_mOTU_v2_6557)
rownames(sub_data1_t) = rownames(data1_t)
colnames(sub_data1_t) = 'meta_mOTU_v2_6557_SRR769418'

data2_t = as.data.frame(t(data2))
sub_data1_t$meta_mOTU_v2_6557_SRR769429 = data2_t$meta_mOTU_v2_6557

data = sub_data1_t
samples <- c('indiv_1','indiv_2')
group <- c('indiv_1','indiv_2')

Counts <- data
dgList <- DGEList(counts=Counts, genes=rownames(Counts), group=group)
## dispersion estimated because there are no replicates
et <- exactTest(dgList, pair=c('indiv_1','indiv_2'), dispersion = 0.5^2)
result <- cbind(bnum=rownames(et$table), et$table )

result$genes <- rownames(result)
volc <- ggplot(result, aes(logFC, -log10(PValue),color=genes)) + #volcanoplot with log2Foldchange versus pvalue
  geom_point(alpha=1,size=6) + #add points colored by significance
  theme_classic() + scale_color_brewer(palette='Spectral') + 
  ggtitle('Individual 1 vs Individual 2 meta OTU 6557 Expression')
volc + xlab('log(fold change)') + ylab('-log10(p-value)')
