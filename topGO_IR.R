#pooja.singh09@gmail.com
#feb2021
#SNC GO enrichment

### I vs R

# read in gene2gomappings and first format out GO annno

library(topGO)
library(tidyverse)
library(stringr)


##BP
godata <- read.table("/data/projects/annotation/brandon_EnTAP/DF/Psme.1_0-DF.merged.final-entap_results.BP", header=F)

go <- godata %>%
 group_by(V1) %>%
 summarize(text = str_c(V2, collapse = ", "))
godf <- as.data.frame(go)
colnames(godf) <- NULL
write.table(godf, "/data/projects/annotation/brandon_EnTAP/DF/Psme.1_0-DF.merged.final-entap_results.BP.formatted", row.names=F, quote=F, sep="\t")

geneID2GO<-readMappings("/data/projects/annotation/brandon_EnTAP/DF/Psme.1_0-DF.merged.final-entap_results.BP.formatted", IDsep=", ", sep="\t")

geneNames <- unique(names(geneID2GO))
str(head(geneID2GO))
length(geneNames)

# read in significant hits from GWAS

sig<-read.table("/data/home/pooja/GWAS/SNC/annotate/SNC_IvsR_fdr0.05.bed.anno.genes4GO",header=F)
dim(sig)
geneList <- factor(as.integer(geneNames %in% as.character(sig$V1)))
names(geneList) <- geneNames
length(geneList)


GOdata_BP <- new("topGOdata", ontology = "BP", allGenes = geneList,annot = annFUN.gene2GO, gene2GO = geneID2GO) 
test.stat <- new("weightCount", testStatistic = GOFisherTest, name = "Fisher test", sigRatio = "ratio")
resultWeight <- getSigGroups(GOdata_BP,test.stat)
allRes_BP<-GenTable(GOdata_BP, weight=resultWeight, orderBy="weight",topNodes=500)
write.table(allRes_BP, file="/data/home/pooja/GWAS/SNC/annotate/SNC_IvsR_fdr0.05.bed.anno.genes4GO.BPenrichment", quote=FALSE, row.names=FALSE, sep="\t")




##CC
godata <- read.table("/data/projects/annotation/brandon_EnTAP/DF/Psme.1_0-DF.merged.final-entap_results.CC", header=F)

go <- godata %>%
 group_by(V1) %>%
 summarize(text = str_c(V2, collapse = ", "))
godf <- as.data.frame(go)
colnames(godf) <- NULL
write.table(godf, "/data/projects/annotation/brandon_EnTAP/DF/Psme.1_0-DF.merged.final-entap_results.CC.formatted", row.names=F, quote=F, sep="\t")

geneID2GO<-readMappings("/data/projects/annotation/brandon_EnTAP/DF/Psme.1_0-DF.merged.final-entap_results.CC.formatted", IDsep=", ", sep="\t")

geneNames <- unique(names(geneID2GO))
str(head(geneID2GO))
length(geneNames)


GOdata_CC <- new("topGOdata", ontology = "CC", allGenes = geneList,annot = annFUN.gene2GO, gene2GO = geneID2GO)
test.stat <- new("weightCount", testStatistic = GOFisherTest, name = "Fisher test", sigRatio = "ratio")
resultWeight <- getSigGroups(GOdata_CC,test.stat)
allRes_CC<-GenTable(GOdata_CC, weight=resultWeight, orderBy="weight",topNodes=500)
write.table(allRes_CC, file="/data/home/pooja/GWAS/SNC/annotate/SNC_IvsR_fdr0.05.bed.anno.genes4GO.CCenrichment", quote=FALSE, row.names=FALSE, sep="\t")


#MF
godata <- read.table("/data/projects/annotation/brandon_EnTAP/DF/Psme.1_0-DF.merged.final-entap_results.MF", header=F)

go <- godata %>%
 group_by(V1) %>%
 summarize(text = str_c(V2, collapse = ", "))
godf <- as.data.frame(go)
colnames(godf) <- NULL
write.table(godf, "/data/projects/annotation/brandon_EnTAP/DF/Psme.1_0-DF.merged.final-entap_results.MF.formatted", row.names=F, quote=F, sep="\t")

geneID2GO<-readMappings("/data/projects/annotation/brandon_EnTAP/DF/Psme.1_0-DF.merged.final-entap_results.MF.formatted", IDsep=", ", sep="\t")

geneNames <- unique(names(geneID2GO))
str(head(geneID2GO))
length(geneNames)



GOdata_MF <- new("topGOdata", ontology = "MF", allGenes = geneList,annot = annFUN.gene2GO, gene2GO = geneID2GO)
test.stat <- new("weightCount", testStatistic = GOFisherTest, name = "Fisher test", sigRatio = "ratio")
resultWeight <- getSigGroups(GOdata_MF,test.stat)
allRes_MF<-GenTable(GOdata_MF, weight=resultWeight, orderBy="weight",topNodes=500)
write.table(allRes_MF, file="/data/home/pooja/GWAS/SNC/annotate/SNC_IvsR_fdr0.05.bed.anno.genes4GO.MFenrichment", quote=FALSE, row.names=FALSE, sep="\t")
