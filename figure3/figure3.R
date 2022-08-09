

library(ggplot2)
library(reshape2)
library(ggplot2)
library(ggridges)
library(data.table)
library(dplyr)
library(cowplot)



ir <- read.table("SNC_IR_sig_SNP_spearmans_cor_r2.txt", header=T)
ir = ir[order(ir$distance),] 


it <- read.table("SNC_IT_sig_SNP_spearmans_cor_r2.txt", header=T)
it = it[order(it$distance),] 

tr <- read.table("SNC_TR_sig_SNP_spearmans_cor_r2.txt", header=T)
tr = tr[order(tr$distance),] 

all <- read.table("SNC_IR_all_SNP_spearmans_cor_r2.txt", header=T)
all = all[order(all$distance),] 


ir$info <- "IvsR"
it$info <- "IvsT"
tr$info <- "TvsR"
all$info <- "background"

df <- rbind(ir, tr, it, all)

p <- ggplot(df, aes(x=info, y=residual)) + 
  geom_violin(trim=F)


A <- p + geom_boxplot(width=0.1, outlier.shape = NA, col="red") + xlab("") + ylab(bquote('LD'~(r^2) ~ 'residual')) + theme_bw()



ir2 <- read.table("SNC_IvR_avrg_rho_fpoints_residual_numberofSNPs.txt", header=T)
it2 <- read.table("SNC_IvT_avrg_rho_fpoints_residual_numberofSNPs.txt", header=T)
tr2 <- read.table("SNC_TvR_avrg_rho_fpoints_residual_numberofSNPs.txt", header=T)

cor.test(it2$number_of_sig_SNPs, it2$avrg_residual, alternative = "greater")
cor.test(tr2$number_of_sig_SNPs, tr2$avrg_residual, alternative = "greater")
cor.test(ir2$number_of_sig_SNPs, ir2$avrg_residual, alternative = "greater")


B <- ggplot(data=ir2, aes(x=number_of_sig_SNPs,y=avrg_residual)) + 
	geom_point(color='orange')+
	geom_smooth(method='lm', formula= y~x, col="black", se=F)+
	annotate("text", label = "cor = 0.04 \n p.value = 0.18",x=9, y=0.4, vjust = 0)+
	xlab("No. of significant IvsR SNPs in gene") + ylab(bquote('average residual LD '~(r^2)~'across gene')) +
	theme_bw()

D <- ggplot(data=tr2, aes(x=number_of_sig_SNPs,y=avrg_residual)) + 
	geom_point(color='lightblue')+
	geom_smooth(method='lm', formula= y~x, col="black", se=F)+
	annotate("text", label = "cor = 0.12 \n p.value = 0.005",x=5, y=0.5)+
	xlab("No. of significant TvsR SNPs in gene") + ylab(bquote('average residual LD '~(r^2)~'across gene')) +
	theme_bw()

C <- ggplot(data=it2, aes(x=number_of_sig_SNPs,y=avrg_residual)) + 
	geom_point(color='lightgreen')+
	geom_smooth(method='lm', formula= y~x, col="black", se=F)+
	annotate("text", label = "cor = 0.0003 \n p.value = 0.5",x=2, y=0.3)+
	xlab("No. of significant IvsT SNPs in gene") + ylab(bquote('average residual LD '~(r^2)~'across gene')) +
	theme_bw()


pdf("fig3.pdf")
plot_grid(A, B, C, D, labels = c('A', 'B', 'C', 'D'), label_size = 12, nrow=2, ncol=2)
dev.off()

svg("fig3.svg")
plot_grid(A, B, C, D, labels = c('A', 'B', 'C', 'D'), label_size = 12, nrow=2, ncol=2)
dev.off
