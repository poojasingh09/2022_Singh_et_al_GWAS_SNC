#pooja.singh09@gmail.com
#code for figure 1 in Singh et al GWAS manuscript
#the output svg file from this script was edited in inkscape

library("SuperExactTest")
library(plyr)
ir <- read.table("SNC-varscan_all_bedfiles_SNP_maf_RD-recalculated_R-column-pairs_I-R-test_CMH-test-results.txt.padjust.sig", header=F)
it <- read.table("SNC-varscan_all_bedfiles_SNP_maf_RD-recalculated_R-column-pairs_I-T-test_CMH-test-results.txt.padjust.sig", header=F)
tr <- read.table("SNC-varscan_all_bedfiles_SNP_maf_RD-recalculated_R-column-pairs_T-R-test_CMH-test-results.txt.padjust.sig", header=F)
rhab <- read.table("Rh_data-varscan_all_bedfiles_SNP_maf_RD-recalculated_CMH-test-results.txt.padjust.sig", header=F)


list <- list(ir$V1, tr$V1, it$V1, rhab$V1)
names(list) <- c("SNC I vs R", "SNC T vs R", "SNC I vs T", "Rh I vs R")
total = as.numeric(summary(list)[1,1]) + as.numeric(summary(list)[2,1]) + as.numeric(summary(list)[3,1]) + as.numeric(summary(list)[4,1])
res=supertest(list, n=total)

svg("overlap_plots_SNC.svg")
plot(res, Layout="landscape", sort.by="size",  degree=2:4, margin=c(0.5,5,1,2),ylab="significant SNPs (p.adjust < 0.05)")
dev.off()


#test for enrichment of overlap (this is what I want)
n_A = length(list$`SNC I vs R`);n_B = length(list$`SNC T vs R`); n_C = 672965; n_A_B = 84
phyper(84-1, n_A, 672904, n_B, lower.tail = F, log.p = FALSE)
E = (n_A*n_B)/672904
print(E)

n_A = length(list$`SNC I vs R`);n_B = length(list$`SNC I vs T`); n_C = 672965; n_A_B = 3
phyper(3-1, n_A, 672904, n_B, lower.tail = F, log.p = FALSE)
E = (n_A*n_B)/672904
print(E)

## test for depletion
n_A = length(list$`I vs R`);n_B = length(list$`T vs R`); n_C = 672965; n_A_B = 84
phyper(84, n_A, 672904, n_B, lower.tail = T, log.p = FALSE)


n_A = length(list$`I vs R`);n_B = length(list$`I vs T`); n_C = 672965; n_A_B = 3
phyper(3, n_A, 672904, n_B, lower.tail = T, log.p = FALSE)



