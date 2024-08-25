# library(tximportData)
library(tximport)
library(rhdf5)
library(GenomicFeatures)
library(tidyr)

dir <- list.dirs("F:/KYS/Kallisto quant on collection 766 and collection 805_ Abundances (tabular)")
n <- list.files("F:/KYS/Kallisto quant on collection 766 and collection 805_ Abundances (tabular)")
files <- file.path(dir,n)
files<-files[grepl("tabular", files)]

files<-files[1:13]
n<-n[1:13]

names(files) <- files

# txdb <- makeTxDbFromGFF(file="D:/W/P/China/KYS/Mus_musculus.GRCm39.110.chr.gff3", format="gff3", dataSource="ENSMBL")
# saveDb(x=txdb, file = "genomic.TxDb")
txdb <- loadDb("genomic.TxDb")
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME") 

txi.kallisto.tsv  <- tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreAfterBar = TRUE, ignoreTxVersion=TRUE)
head(txi.kallisto.tsv$counts)
# txi.kallisto.tsv$abundance

# library(org.Mm.eg.db)
# db <- org.Mm.eg.db
# 
# counts <- txi.kallisto.tsv$counts
# gene_name <- mapIds(db, keys=rownames(counts), column="SYMBOL", keytype="ENSEMBL")
# counts <- as.data.frame(counts)
# counts[,'gene_name'] <- gene_name
# 
# # Assuming counts$gene_name is your vector
# gene_name_vector <- counts$gene_name  # Replace with your actual vector
# 
# # Find duplicated elements
# duplicated_elements <- gene_name_vector[duplicated(gene_name_vector)]
# 
# # Extract unique duplicated elements
# unique_duplicated_elements <- unique(duplicated_elements)
# 
# 
# # Sum rows by gene_name
# summed_counts <- counts %>%
#   group_by(gene_name) %>%
#   summarise_all(sum)
# 
# count <- as.data.frame(summed_counts)
# 
# library(ape)
# ensmb <- read.gff("D:/W/P/China/KYS/Mus_musculus.GRCm39.110.chr.gff3", na.strings = c(".", "?"), GFF3 = TRUE)
# 
# library(stringr)
# 
# ensgenes <- NULL
# # Extract specific fields using string manipulation
# ensgenes <- str_extract(ensmb$attributes, "gene:([^;]+)")
# ensgenes<-as.data.frame(ensgenes)
# ensgenes$gene_name <- str_extract(ensmb$attributes, "Name=([^;]+)")
# #ensmb$biotype <- str_extract(ensmb$attributes, "biotype=([^;]+)")
# #ensmb$description <- str_extract(ensmb$attributes, "description=([^;]+)")
# #ensmb$logic_name <- str_extract(ensmb$attributes, "logic_name=([^;]+)")
# #ensmb$version <- str_extract(ensmb$attributes, "version=([^;]+)")
# 
# # Remove prefixes
# ensgenes$ensgenes <- sub("gene:", "", ensgenes$ensgenes)
# ensgenes$gene_name <- sub("Name=", "", ensgenes$gene_name)
# ensgenes <- ensgenes %>% dplyr::filter(!grepl("-", gene_name))
# ensgenes <- drop_na(ensgenes)
# 
# # Print the extracted information
# head(ensmb[, c("gene_id", "gene_name", "biotype", "description", "logic_name", "version")])
# 
# # Assuming counts$gene_name is your vector
# gene_name_vector <- ensgenes$gene_name  # Replace with your actual vector
# 
# # Find duplicated elements
# duplicated_elements <- gene_name_vector[duplicated(gene_name_vector)]
# 
# # Extract unique duplicated elements
# unique_duplicated_elements <- unique(duplicated_elements)
# 
# counts <- txi.kallisto.tsv$counts
# gene_name <- mapIds(db, keys=rownames(counts), column="SYMBOL", keytype="ENSEMBL")
# counts <- as.data.frame(counts)
# counts[,'gene_name'] <- gene_name

library(biomaRt)
library(annotables)

counts <- txi.kallisto.tsv$counts
# input list of Ensembl ID's
ensembl.ids <- rownames(counts)


# method 1: biomaRt
listEnsembl(mirror = "useast")
ensembl <- useEnsembl(biomart = "genes", mirror = "useast", host="https://useast.ensembl.org")
datasets <- listDatasets(ensembl)

listEnsemblArchives()
listMarts(host="https://jul2023.archive.ensembl.org")
ensembl.con <- useMart("ensembl_mart_110", dataset = 'mmusculus_gene_ensembl', version=110, host="https://jul2023.archive.ensembl.org")

attr <- listAttributes(ensembl.con)
filters <- listFilters(ensembl.con)

Martgenes <- getBM(attributes = c('ensembl_gene_id','external_gene_name'),
      filters = "ensembl_gene_id",
      values = ensembl.ids,
      mart = ensembl.con)

counts <- as.data.frame(counts)
counts$Symbol <- Martgenes$external_gene_name

library(dplyr)

# Assuming counts$gene_name is your vector
gene_name_vector <- counts$Symbol  # Replace with your actual vector

# Find duplicated elements
duplicated_elements <- gene_name_vector[duplicated(gene_name_vector)]

# Extract unique duplicated elements
unique_duplicated_elements <- unique(duplicated_elements)


# Sum rows by gene_name
summed_counts <- counts %>%
  group_by(Symbol) %>%
  summarise_all(sum)

count <- as.data.frame(summed_counts)


library(edgeR)

counts <- count

rownames(counts) <- count$Symbol
counts$Symbol <- NULL
colnames(counts) <- c("c1","t1","c2","c3","c4","c5","c6","t2","t3","t4","t5","t6", "t7")

d0 <- DGEList(counts)
d0 <- calcNormFactors(d0)

cutoff <- 1
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,] 
dim(d) # number of genes left $ Dropped from 56479 to 17261 genes

counts <- d
n <- colnames(counts) # names
n <- substr(n, 1, 1)
design <- model.matrix(~0+n)
mm <- design
y <- voom(d, mm, plot = T)

plotMDS(d)

dim(mm)
dim(d)

#fit <- lmFit(y, design)
# head(coef(fit))
contrasts.matrix <- makeContrasts(nt - nc, levels=design)
# fit2 <- contrasts.fit(fit, contrasts.matrix)
# fit2 <- eBayes(fit2)

# Inshallah
# topTable(fit2)
# decideTests(fit2)


# Bomb if bad
aw <- arrayWeights(y,design)
fit <- lmFit(y, design,
             weights = aw)
contrasts.matrix <- makeContrasts(nt - nc, levels=design)
fit2 <- contrasts.fit(fit, contrasts.matrix)
fit2 <- eBayes(fit2)

topTable(fit2)
decideTests(fit2)
# End of outlier bombing

full_results <- topTable(fit2, number=Inf)

full_results$Symbol <- rownames(full_results)

write.csv(full_results, "full_results_villi.csv")

h_sig <- filter(full_results, adj.P.Val < 0.05)
write.csv(h_sig, "h_sig_villi.csv")

library(ggplot2)
library(ggrepel)
p_cutoff <- 0.05
fc_cutoff <- 0
topN <- 20

library(dplyr)

# full_results <- read.csv("full_results_crypt.csv", row.names = 1)

sd <- filter(full_results, adj.P.Val < 0.05)
u <- filter(sd, logFC > 0)
d <- filter(sd, logFC < 0)

u <- u[order(u$adj.P.Val, decreasing = FALSE), ]
d <- d[order(d$adj.P.Val, decreasing = FALSE), ]

stori <- c(u$Symbol[1:15], d$Symbol[1:10], c("Cyp24a1"))

# full_results %>%
#   mutate(Significant = adj.P.Val < p_cutoff, abs(logFC) > fc_cutoff) %>%
#   mutate(Up = logFC > fc_cutoff & adj.P.Val < p_cutoff) %>%
#   mutate(ColorCategory = case_when(
#     Up ~ "Upregulated",
#     abs(logFC) < fc_cutoff | adj.P.Val > p_cutoff ~ "Not Significant",
#     TRUE ~ "Downregulated"
#   )) %>%
#   mutate(Rank = 1:n(), Label = ifelse(((Rank < topN) & adj.P.Val < p_cutoff), symbol,"")) %>%
#   ggplot(aes(x = logFC, y = -log10(adj.P.Val), color = ColorCategory, label = Label)) +
#   geom_point() +
#   geom_text_repel(col = "black", box.padding = 0.5, max.overlaps = Inf, size=5) +  # Increase max.overlaps
#   scale_color_manual(values = c("blue", "black", "red")) +  # Define colors
#   labs(x = "log2FC", y = "-log10(adj.P.Val)") +   # Add labels to axes
#  theme(text=element_text(size=10))

# dev.new(width = 800, height = 800)

library(ggbreak) 
library(patchwork)

# write.csv(full_results, "full_results_crypt.csv")

# stori <- stori[!stori %in% c("Cyp24a1", "Epha4", "D130043K22Rik", "Rnf149", "Ifngr1", "5033404E19Rik")]


stori <- stori[!stori %in% c("Myo19")]

# Myo
# Pdlim
# Mmgt


p <- full_results %>%
  mutate(Significant = adj.P.Val < p_cutoff, abs(logFC) > fc_cutoff) %>%
  mutate(Up = logFC > fc_cutoff & adj.P.Val < p_cutoff) %>%
  mutate(ColorCategory = case_when(
    Up ~ "Upregulated",
    abs(logFC) < fc_cutoff | adj.P.Val > p_cutoff ~ "Not Significant",
    TRUE ~ "Downregulated"
  )) %>%
  mutate(Rank = 1:n(), Label = ifelse(Symbol %in% stori, Symbol,"")) %>%
  ggplot(aes(x = logFC, y = -log10(adj.P.Val), color = ColorCategory, label = Label)) +
  geom_point() +
  geom_text_repel(aes(size = ifelse(Label %in% c("Msi1", "Mex3a", "Lgr5", "Smoc2", "Tert"), 24/.pt, 14/.pt)),
                  col = "black", box.padding = 1, max.overlaps = Inf, min.segment.length = unit(0, 'lines'), force = 1) +  
  geom_segment(aes(x = logFC, y = -log10(adj.P.Val), xend = logFC, yend = -log10(adj.P.Val)), color = "grey50", alpha = 0.5) + # Add lines connecting labels to points
  scale_color_manual(values = c("blue", "black", "red")) +  
  scale_size_identity() +  
  labs(x = "log2FC", y = "-log10(adj.P.Val)") +   
  theme(text=element_text(size=14))

# plot <- p + scale_x_break(c(6, 12.5), scales=0.15, expand = expansion(add=0.4), ticklabels = c(12.5, 13))

plot <- p + scale_x_break(c(6, 11), scales=0.2, expand = expansion(add=0.4), ticklabels = c(11, 11.5, 12, 12.5))

ggsave("full_results_villi_scaled_manual.png", plot, width = 9.38, height = 7.5, units = "in", dpi = 300)

num_sig <- nrow(filter(full_results, adj.P.Val < 0.05))
num_up <- nrow(filter((filter(full_results, adj.P.Val < 0.05)), logFC > 0))
num_down <- nrow(filter((filter(full_results, adj.P.Val < 0.05)), logFC < 0))

num_sig
num_up
num_down
