#pooja.singh09@gmail.com
#may2020
# this script calcualtes the average LD per gene using the LD calculation fitting the Hill and Weir model using all SNPs in significant GWAS genes

input1 ="/data/home/pooja/GWAS/SNC/SNC-varscan_all_bedfiles_SNP_maf_RD-recalculated_R-column-pairs_I-R-test_CMH-test-results.txt.padjust.sig.anno.ranked"
input2 = "/data/home/pooja/GWAS/SNC/LD_decay/IR/SNC_IR_all_SNP_spearmans_cor_r2.txt"
input3 = "/data/home/pooja/GWAS/SNC/SNC-varscan_all_bedfiles_SNP_maf_RD-recalculated_R-column-pairs_I-T-test_CMH-test-results.txt.padjust.sig.anno.ranked"
input4 = "/data/home/pooja/GWAS/SNC/LD_decay/IT/SNC_IT_all_SNP_spearmans_cor_r2.txt"
input5 = "/data/home/pooja/GWAS/SNC/SNC-varscan_all_bedfiles_SNP_maf_RD-recalculated_R-column-pairs_T-R-test_CMH-test-results.txt.padjust.sig.anno.ranked"
input6 = "/data/home/pooja/GWAS/SNC/LD_decay/TR/SNC_TR_all_SNP_spearmans_cor_r2.txt"

#####load libs########

library(stringr)
library(dplyr)


##################### I vs R
setwd("/data/home/pooja/GWAS/SNC/LD_decay")

df1 <-read.table(input1,header=T)
df1$scaff <- stringr::str_split_fixed(df1$snp, "-", 2) [,1]
df2 <- df1[,c(3,4)]

dfcor <-read.table(input2,header=T)
dfcor$scaff <- stringr::str_split_fixed(dfcor$snp1, "-", 2) [,1]


annotated_dfcor <- merge(dfcor, df2, by = "scaff") ## annotate LD estimates

list <- data.frame(unique(annotated_dfcor$gene))
colnames(list) <- c("genes")

out<- matrix(nrow=nrow(list),ncol=6)
colnames(out) <- c("gene", "avrg_spearmans_r2", "avrg_distance", "avrg_fpoints_predicted", "avrg_residual", "hdecay_distance")

for (i in 2:nrow(list)){
	
	subset <- subset(annotated_dfcor, gene == list[i,1])
	out[i,1] <- as.character(list[i,1])
	out[i,2] <- mean(subset$spearmans_r2)
	out[i,3] <- mean(subset$distance)
	out[i,4] <- mean(subset$fpoints_predicted)
	out[i,5] <- mean(subset$residual)
	h.decay <- (max(subset$spearmans_r2, na.rm =T))/2
	out[i,6] <- subset$distance[which.min(abs(subset$spearmans_r2-h.decay))]
	print(h.decay)
	print(out[i,6])
	
	
}

out1 <- na.omit(out)

snps_in_genes <- data.frame(table(df1$gene))[-1,] # count the number of signifcant hits per gene, for snps that his genes obviously
colnames(snps_in_genes) <- c("gene", "number_of_sig_SNPs")

gene_avrg_snps_ld <- merge(out1, snps_in_genes, by="gene")
write.table(gene_avrg_snps_ld, "SNC_IvR_avrg_rho_fpoints_residual_numberofSNPs.txt", row.names=F, quote=F, sep="\t")



##################### I vs T
setwd("/data/home/pooja/GWAS/SNC/LD_decay")

df1 <-read.table(input3,header=T)
df1$scaff <- stringr::str_split_fixed(df1$snp, "-", 2) [,1]
df2 <- df1[,c(3,4)]

dfcor <-read.table(input4,header=T)
dfcor$scaff <- stringr::str_split_fixed(dfcor$snp1, "-", 2) [,1]


annotated_dfcor <- merge(dfcor, df2, by = "scaff") ## annotate LD estimates

list <- data.frame(unique(annotated_dfcor$gene))
colnames(list) <- c("genes")

out<- matrix(nrow=nrow(list),ncol=6)
colnames(out) <- c("gene", "avrg_spearmans_r2", "avrg_distance", "avrg_fpoints_predicted", "avrg_residual", "hdecay_distance")

for (i in 2:nrow(list)){
	
	subset <- subset(annotated_dfcor, gene == list[i,1])
	out[i,1] <- as.character(list[i,1])
	out[i,2] <- mean(subset$spearmans_r2)
	out[i,3] <- mean(subset$distance)
	out[i,4] <- mean(subset$fpoints_predicted)
	out[i,5] <- mean(subset$residual)
	h.decay <- (max(subset$fpoints_predicted))/2
        out[i,6] <- subset$distance[which.min(abs(subset$fpoints_predicted-h.decay))]
	
}

out1 <- na.omit(out)

snps_in_genes <- data.frame(table(df1$gene))[-1,] # count the number of signifcant hits per gene, for snps that his genes obviously
colnames(snps_in_genes) <- c("gene", "number_of_sig_SNPs")

gene_avrg_snps_ld <- merge(out1, snps_in_genes, by="gene")
write.table(gene_avrg_snps_ld, "SNC_IvT_avrg_rho_fpoints_residual_numberofSNPs.txt", row.names=F, quote=F, sep="\t")


##################### T vs R
setwd("/data/home/pooja/GWAS/SNC/LD_decay")

df1 <-read.table(input5,header=T)
df1$scaff <- stringr::str_split_fixed(df1$snp, "-", 2) [,1]
df2 <- df1[,c(3,4)]

dfcor <-read.table(input6,header=T)
dfcor$scaff <- stringr::str_split_fixed(dfcor$snp1, "-", 2) [,1]


annotated_dfcor <- merge(dfcor, df2, by = "scaff") ## annotate LD estimates

list <- data.frame(unique(annotated_dfcor$gene))
colnames(list) <- c("genes")

out<- matrix(nrow=nrow(list),ncol=6)
colnames(out) <- c("gene", "avrg_spearmans_r2", "avrg_distance", "avrg_fpoints_predicted", "avrg_residual", "hdecay_distance")

for (i in 2:nrow(list)){
	
	subset <- subset(annotated_dfcor, gene == list[i,1])
	out[i,1] <- as.character(list[i,1])
	out[i,2] <- mean(subset$spearmans_r2)
	out[i,3] <- mean(subset$distance)
	out[i,4] <- mean(subset$fpoints_predicted)
	out[i,5] <- mean(subset$residual)
	h.decay <- (max(subset$fpoints_predicted))/2
        out[i,6] <- subset$distance[which.min(abs(subset$fpoints_predicted-h.decay))]
	
}

out1 <- na.omit(out)

snps_in_genes <- data.frame(table(df1$gene))[-1,] # count the number of signifcant hits per gene, for snps that his genes obviously
colnames(snps_in_genes) <- c("gene", "number_of_sig_SNPs")

gene_avrg_snps_ld <- merge(out1, snps_in_genes, by="gene")
write.table(gene_avrg_snps_ld, "SNC_TvR_avrg_rho_fpoints_residual_numberofSNPs.txt", row.names=F, quote=F, sep="\t")


############### final plot ##############

ir <- read.table("SNC_IvR_avrg_rho_fpoints_residual_numberofSNPs.txt", header=T)
it <- read.table("SNC_IvT_avrg_rho_fpoints_residual_numberofSNPs.txt", header=T)
tr <- read.table("SNC_TvR_avrg_rho_fpoints_residual_numberofSNPs.txt", header=T)

pdf("SNC_SNPs_vs_LD.pdf")
par(mfrow=(c(3,2)))
plot(ir$number_of_sig_SNPs,ir$avrg_residual, main="IvR no. of SNPs in gene vs \n average LD residual of SNPs in each gene", pch = 16, col="orange", xlab="No. of sig SNPs in gene", ylab="residual LD (r2)" )
abline(lm(ir$avrg_residual ~ ir$number_of_sig_SNPs))

plot(ir$number_of_sig_SNPs,ir$avrg_fpoints_predicted, main="IvR no. of SNPs in gene vs \n average predicted LD of SNPs in each gene", pch = 16, col="orange", xlab="No. of sig SNPs in gene", ylab = "LD (r2)")
abline(lm(ir$avrg_fpoints_predicted ~ ir$number_of_sig_SNPs))

plot(it$number_of_sig_SNPs,it$avrg_residual, main="IvT no. of SNPs in gene vs \n average LD residual of SNPs in each gene", pch = 16, col="green", xlab="No. of sig SNPs in gene",ylab="residual LD (r2)")
abline(lm(ir$avrg_residual ~ ir$number_of_sig_SNPs))

plot(it$number_of_sig_SNPs,it$avrg_fpoints_predicted, main="IvT no. of SNPs in gene vs \n average predicted LD of SNPs in each gene", pch = 16, col="green", xlab="No. of sig SNPs in gene",ylab = "LD (r2)")
abline(lm(ir$avrg_fpoints_predicted ~ ir$number_of_sig_SNPs))


plot(tr$number_of_sig_SNPs,tr$avrg_residual, main="TvR no. of SNPs in gene vs \n average LD residual of SNPs in each gene", pch = 16, col="blue", xlab="No. of sig SNPs in gene",ylab="residual LD (r2)")
abline(lm(tr$avrg_residual ~ tr$number_of_sig_SNPs))

plot(tr$number_of_sig_SNPs,tr$avrg_fpoints_predicted, main="TvR no. of SNPs in gene vs \n average predicted LD of SNPs in each gene", pch = 16, col="blue", xlab="No. of sig SNPs in gene",ylab = "LD (r2)")
abline(lm(tr$avrg_fpoints_predicted ~ tr$number_of_sig_SNPs))

dev.off()
