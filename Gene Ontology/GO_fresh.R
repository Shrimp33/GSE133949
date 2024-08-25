library(DESeq2)
library(dplyr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(AnnotationDbi)

h_sig_crypt <- filter(read.csv("h_sig.csv", row.names = 1), adj.P.Val < 0.01)
rownames(h_sig_crypt) <- h_sig_crypt$symbol
h_sig_villi <- filter(read.csv("h_sig_villi.csv", row.names = 1), adj.P.Val < 0.01)
h_sig_colon <- filter(read.csv("h_sig_colon.csv", row.names = 1), adj.P.Val < 0.01)

crypt_backgroud <- rownames(readRDS("Pass_cutoff_crypt_genes.RDS")$counts)
colon_backgroud <- rownames(readRDS("Colon_counts_pass_cutoff.rds")$counts)
villi_backgroud <- rownames(readRDS("Villi_counts_pass_cuttoff.RDS")$counts)

# mapping <- read.csv("ensmbl_to_symbol110.csv", row.names = 1)
# universeSym <- readRDS("activezFPKM.RDS")
# sym <- mapping[unlist(universeSym),]

GO_results_crypt <- enrichGO(gene = rownames(h_sig_crypt), OrgDb = "org.Mm.eg.db", keyType = "SYMBOL",
                       ont = "ALL", universe = crypt_backgroud, pvalueCutoff = 0.01)

fit_crypt <- plot(barplot(GO_results_crypt, showCategory = 15))

write.csv(GO_results_crypt, "GO_results_crypt.csv")

library(ggplot2)

ggplot(head(GO_results_crypt, 15), aes(reorder(Description, Count), Count)) +
  geom_col(aes(fill=-log10(p.adjust))) +
  coord_flip() +
  labs(x="Ontology", y="Count",
       title="Gene Ontology Analysis of DGEs (adj p < 0.01)") + 
  theme_minimal()

ggsave("GO_crypt.png", width = 10, height = 7.5)

GO_results_villi <- enrichGO(gene = rownames(h_sig_villi), OrgDb = "org.Mm.eg.db", keyType = "SYMBOL",
                       ont = "ALL", universe = villi_backgroud, pvalueCutoff = 0.01)

write.csv(GO_results_villi, "GO_results_villi.csv")

fit_villi <- plot(barplot(GO_results_villi, showCategory = 15))

ggplot(GO_results_villi[1:15,], aes(reorder(Description, Count), Count)) +
  geom_col(aes(fill=-log10(p.adjust))) +
  coord_flip() +
  labs(x="Ontology", y="Count",
       title="Gene Ontology Analysis of DGEs (adj p < 0.01)") + 
  theme_minimal()

ggsave("GO_villi.png", width = 10, height = 7.5)

GO_results_colon <- enrichGO(gene = rownames(h_sig_colon), OrgDb = "org.Mm.eg.db", keyType = "SYMBOL",
                       ont = "ALL", universe = colon_backgroud, pvalueCutoff = 0.01)

write.csv(GO_results_colon, "GO_results_colon.csv")

fit_colon <- plot(barplot(GO_results_colon, showCategory = 15))

library(scales)
library(stringr)

ggplot(head(GO_results_colon,15), aes(reorder(Description, Count), Count)) +
  geom_col(aes(fill=-log10(p.adjust))) +
  coord_flip() +
  labs(x="Ontology", y="Count",
       title="Gene Ontology Analysis of DGEs (adj p < 0.01)") + 
  theme_minimal() +
  scale_x_discrete(labels = function(Description) str_wrap(Description, width = 40))

ggsave("GO_colon.png", width = 10, height = 7.5)
