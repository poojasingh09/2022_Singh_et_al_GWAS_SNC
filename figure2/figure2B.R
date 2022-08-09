##pooja.singh09@gmail.com
##may 2021
##count number of significant SNPs per gene for SNC manuscript Singh et al


library("ggplot2")
library("ggprism")
library(cowplot)

ir <- read.table("SNC-varscan_all_bedfiles_SNP_maf_RD-recalculated_R-column-pairs_I-R-test_CMH-test-results.txt.padjust.sig.anno.ranked.txt", header=T)
ir2 <- ir[,c(1,3)]
ir_count <- data.frame(table(ir2$gene))[-1,]
ir_count$group <- "SNC IvsR"


tr <- read.table("SNC-varscan_all_bedfiles_SNP_maf_RD-recalculated_R-column-pairs_T-R-test_CMH-test-results.txt.padjust.sig.anno.ranked.txt", header=T)
tr2 <- tr[,c(1,3)]
tr_count <- data.frame(table(tr2$gene))[-1,]
tr_count$group <- "SNC TvsR"

it <- read.table("SNC-varscan_all_bedfiles_SNP_maf_RD-recalculated_R-column-pairs_I-T-test_CMH-test-results.txt.padjust.sig.anno.ranked.txt", header=T)
it2 <- it[,c(1,3)]
it_count <- data.frame(table(it2$gene))[-1,]
it_count$group <- "SNC IvsT"


df <- rbind(ir_count, tr_count, it_count)


svg("figure2B_snpspergene.svg")
p <- ggplot(df, aes(x = Freq, fill = group)) + theme_prism()                     # Draw overlaying histogram
p+geom_histogram(position = "identity", alpha = 0.3, bins = 50) + xlab("No. of SNPs per gene") + ylab("No. of genes")+ theme(legend.title = element_blank())
dev.off()



# Histogram by group in ggplot2


svg("figure2B_snpspergene_v2.svg")
ggplot(df, aes(x = Freq, fill = group, colour = group)) + 
  geom_histogram(position = "dodge", alpha = 0.3, bins = 16) + xlab("No. of SNPs per gene") + ylab("No. of genes")+ theme(legend.title = element_blank()) + theme_prism()
dev.off()

# Histogram by group in ggplot2 with inset



dfa <- df[df$Freq <= 3,]
dfb <- df[df$Freq > 3,]

group.colors <- c("SNC IvsR" = "orange", "SNC TvsR" = "skyblue", "SNC IvsT" ="green")

main.plot <- ggplot(dfa, aes(x = Freq, fill = group, colour = group)) + 
  geom_histogram(position = "dodge", alpha = 0.3, binwidth=0.5, bins=3, colour="black") + xlab("No. of SNPs per gene") + ylab("No. of genes")+ theme(legend.title = element_blank()) + theme_prism() + scale_fill_manual(values = group.colors)

inset.plot <- ggplot(dfb, aes(x = Freq, fill = group, colour = group)) + 
  geom_histogram(position = "dodge", alpha = 0.3, bins = 13, show.legend = FALSE, colour="black") + xlab("No. of SNPs per gene") + ylab("No. of genes") + theme_prism() +  scale_fill_manual(values = group.colors) +theme(axis.text=element_text(size=8), axis.title=element_text(size=8))

plot.with.inset <-
  ggdraw() +
  draw_plot(main.plot) +
  draw_plot(inset.plot, x = .4, y = 0.6, width = .4, height = .4)
plot.with.inset


svg("figure2B_snpspergene_v3.svg")
plot.with.inset
dev.off()


png("figure2B_snpspergene_v3.png")
plot.with.inset
dev.off()



## density plot
plot(density(it_count$Freq), col="lightgreen")
lines(density(ir_count$Freq), col="orange")
lines(density(tr_count$Freq), col="skyblue")


