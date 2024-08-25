library(DESeq2)
#sigs <- na.omit(res.crypt)
library(dplyr)
#sigs <- filter(sigs, padj < 0.05)
library(clusterProfiler)
library(org.Mm.eg.db)
library(AnnotationDbi)

wgenes <- readRDS("txi.kallisto.tsv$counts_crypt_ensid.RDS")
g <- rownames(read.csv("h_sig.csv", row.names = 1))

wgenes <- read.csv("D:/WGCNA/ex_c_map.csv", row.names = 2)
wgenes$X <- NULL
totalg <- read.csv("EdgeRActive_filterbyExpr.csv")[,2]

ballcomsumer <- rownames(wgenes)

table(wgenes$unlist.dyanmicColors.) # blue     brown     green       red turquoise    yellow 
# Or X if reading from BiCor
wids <- rownames(filter(wgenes, x == "darkgreen"))
# real <- as.data.frame(totalg)
# rownames(real) <- totalg
DGE <- read.csv("D:/DGE/CryptLimmaDGE.csv")
DGE <- readRDS("ensCounts.RDS")
sigs <- filter(DGE, adj.P.Val < 0.05)
# bigboy <- union(ballcomsumer, totalg)

library(biomaRt)

mart <- useEnsembl("ensembl","mmusculus_gene_ensembl", version = 110, host="https://asia.ensembl.org")

Martgenes <- getBM(attributes=c("ensembl_gene_id", "external_gene_name"),
                   filters = 'ensembl_gene_id', values = ballcomsumer, mart = mart, bmHeader = TRUE)

universeSymbols <- Martgenes$`Gene name`

# universeSymbols <- readRDS("univserseSymbols.RDS")

library(org.Mm.eg.db)

genes_to_test <- filter(sigs, logFC > 0.5)$X
genes_to_test <- wids
GO_results <- enrichGO(gene = g, OrgDb = "org.Mm.eg.db", keyType = "ENSEMBL", ont = "ALL", universe = rownames(wgenes), pvalueCutoff = 0.01)
head(GO_results$BgRatio)
write.csv(as.data.frame(GO_results), "Cryptupwithproperuniverse_0_01.csv")
fit <- plot(barplot(GO_results, showCategory = 15))

library(ggplot2)

library(dplyr)

df <- as.data.frame(GO_results)
attach(df)
df$padj <- -log(df$p.adjust)
GO_results_gg <- filter(df, p.adjust < 0.01)[order(padj),]


ggplot(GO_results_gg[1:15,], aes(reorder(Description, Count), Count)) +
  geom_col(aes(fill=-log10(p.adjust))) +
  coord_flip() +
  labs(x="Ontology", y="Hits in DGE",
       title="Gene Ontology Analysis of DGEs (adj p < 0.01)") + 
  theme_minimal()

png("Cryptupwithproperuniverse_0_01.png", res = 250, width = 2800, height = 3600)
print(fit)
dev.off()

genes_to_test <- filter(sigs, logFC < 0.5)$X
GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "ALL")
write.csv(as.data.frame(GO_results), "CryptGODown.csv")
fit <- plot(barplot(GO_results, showCategory = 15))

png("CryptGODown.png", res = 250, width = 2800, height = 3600)
print(fit)
dev.off()

# COLON

DGE <- read.csv("D:/DGE/ColonLimmaDGE.csv")
sigs <- filter(DGE, adj.P.Val < 0.05)

genes_to_test <- filter(sigs, logFC > 0.5)$X
GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "ALL")
write.csv(as.data.frame(GO_results), "ColonGOUp.csv")
fit <- plot(barplot(GO_results, showCategory = 15))

png("ColonGOUp.png", res = 250, width = 2800, height = 3600)
print(fit)
dev.off()

genes_to_test <- filter(sigs, logFC < 0.5)$X
GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "ALL")
write.csv(as.data.frame(GO_results), "ColonGODown.csv")
fit <- plot(barplot(GO_results, showCategory = 15))

png("ColonGODown.png", res = 250, width = 2800, height = 3600)
print(fit)
dev.off()

# VILLI 

DGE <- read.csv("D:/DGE/VilliLimmaDGE.csv")
sigs <- filter(DGE, adj.P.Val < 0.05)

genes_to_test <- filter(sigs, logFC > 0.5)$X
GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "ALL")
write.csv(as.data.frame(GO_results), "VilliGOUp.csv")
fit <- plot(barplot(GO_results, showCategory = 15))

png("VilliGOUp.png", res = 250, width = 2800, height = 3600)
print(fit)
dev.off()

genes_to_test <- filter(sigs, logFC < 0.5)$X
GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "ALL")
write.csv(as.data.frame(GO_results), "VilliGODown.csv")
fit <- plot(barplot(GO_results, showCategory = 15))

png("VilliGODown.png", res = 250, width = 2800, height = 3600)
print(fit)
dev.off()

filepaths <- c("C:/Users/notaj/OneDrive/Desktop/EigenModule/brown _list.rds", "C:/Users/notaj/OneDrive/Desktop/EigenModule/darkmagenta _list.rds", "C:/Users/notaj/OneDrive/Desktop/EigenModule/darkturquoise _list.rds", "C:/Users/notaj/OneDrive/Desktop/EigenModule/floralwhite _list.rds", "C:/Users/notaj/OneDrive/Desktop/EigenModule/green _list.rds", "C:/Users/notaj/OneDrive/Desktop/EigenModule/ivory _list.rds", "C:/Users/notaj/OneDrive/Desktop/EigenModule/lightyellow _list.rds", "C:/Users/notaj/OneDrive/Desktop/EigenModule/turquoise _list.rds")
brown <- readRDS(filepaths[1]) #did
darkmagenta <- readRDS(filepaths[2]) #NULL ontology #No real pathways #did
darkturquoise <- readRDS(filepaths[3]) #did
floralwhite <- readRDS(filepaths[4]) #NULL ontology #No real pathways #did
green <- readRDS(filepaths[5]) #did #Wnt
ivory <- readRDS(filepaths[6]) #NULL ontology #No real pathways #did
lightyellow <- readRDS(filepaths[7]) #did
turquoise <- readRDS(filepaths[8]) #division #did

#Enrichr

GO_results <- enrichGO(gene = green, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "ALL")
write.csv(as.data.frame(GO_results), "GO_resultsGreen.png.csv")
fit <- plot(barplot(GO_results, showCategory = 15))

png("Green.png", res = 250, width = 2800, height = 3600)
print(fit)
dev.off()
