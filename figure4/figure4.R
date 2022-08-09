#pooja.singh09@gmail.com
#code for figure 4A in Singh et al GWAS manuscript
# figure 4B is made in https://string-db.org

library(dplyr)
library(tidyr)
library(stringr)
library(data.table)
library(pheatmap)
library(reshape2)
options(stringsAsFactors = FALSE)
library(dendextend)
library(ggplot2)
setwd(".")



###  RHAB I vs R

data1 <- read.table("Rh_RvsS_fdr0.05.txt.locus.FREQ", header=T, sep=" ")
dim(data1)
rownames(data1)[2] <- "ERF1-like"

cormat <- cor(t(data1), method="spearman", use = "pairwise.complete.obs") 
r2 <- (cormat)^2		



pdf('Rh_RvsS_sig_snps_r2.pdf')
#pheatmap(r2, legend=T, legend_labels="r2", color = c("white", "lightyellow", "orange", "red"),breaks=c(0, 0.3, 0.6, 0.8, 1.0))
pheatmap(r2, legend=T)
dev.off()

svg('Rh_RvsS_sig_snps_r2.svg')
pheatmap(r2, legend=T)
dev.off()



hclust <- hclust(dist(r2), method = "complete")
pdf("Rh_RvsS_sig_snps_r2.denogram.pdf", w=10, h=5)
par(mar = c(2,2,2,10))
as.dendrogram(hclust) %>%
  plot(horiz = TRUE, cex=0.5)
dev.off()


