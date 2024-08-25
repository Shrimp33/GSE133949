# library(edgeR)

# # counts <- as.matrix(readRDS("txi.kallisto.tsv$counts.RDS"))
# counts <- readRDS("CryptCounts.RDS")
# rownames(counts) <- count$id
# counts$Symbol <- NULL
# counts$id <- NULL
# colnames(counts) <- c("c1","t1","c2","c3","c4","c5","c6","t2","t3","t4","t5","t6")

# d0 <- DGEList(counts)
# d0 <- calcNormFactors(d0)

# # library(biomaRt)
# # 
# # ensembl = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
# # gene_lengths <- getBM(attributes = c("ensembl_gene_id", "gene_biotype", "transcript_length"),
# #                       filters = "ensembl_gene_id",
# #                       values = counts$id,
# #                       mart = ensembl)
# # saveRDS(gene_lengths, "gene_lengths_for_kallisto.RDS")
# # 
# # gene_lengths_unique <- aggregate(transcript_length ~ ensembl_gene_id, gene_lengths, mean) # Take the average of isoforms
# # saveRDS(gene_lengths_unique, "gene_lengths_for_kallisto_avg.RDS")
# # 
# # avg_gn_ln <- mean(gene_lengths_unique$transcript_length)
# # 
# # cpmC <- cpm(d0, log=FALSE)

# # merged_data <- merge(gene_lengths_unique, cpmC, by.x = "ensembl_gene_id", by.y = "row.names")
# # merged_data$transcript_length <- merged_data$transcript_length / avg_gn_ln
# # 
# # tmm <- merged_data[,-(1:2)]
# # rownames(tmm) <- merged_data$ensembl_gene_id
# # 
# # tmm_norm <- tmm / merged_data[,"transcript_length"]


# # Question of what threshhold to drop?
# #counts #merged_data #rpk #norm.counts.edger #norm.counts.rpk_edger #tpm #norm.counts.deseq

# # 
# # 
# # cutoff <- 3
# # # drop <- which(apply(norm.counts.edger, 1, max) < cutoff)
# # 
# # # Calculate median count for each gene
# # median_counts <- apply(norm.counts.rpk_edger, 1, mean)
# # 
# # # Identify genes with median counts less than cutoff
# # drop <- which(median_counts < cutoff)
# # 
# # d <- d0[-drop,]
# # 
# # dim(d) # number of genes left $ Dropped from 56478 to 14910 genes

# ags <- intersect(rownames(d0), unlist(activeGenes))

# d <- d0[ags,]
# dim(d)

# counts <- d
# n <- colnames(counts) # names
# n <- substr(n, 1, 1)
# design <- model.matrix(~0+n)
# mm <- design
# y <- voom(d, mm, plot = T)

# plotMDS(norm.counts.deseq)

# dim(mm)
# dim(d)

# #fit <- lmFit(y, design)
# # head(coef(fit))
# contrasts.matrix <- makeContrasts(nt - nc, levels=design)
# # fit2 <- contrasts.fit(fit, contrasts.matrix)
# # fit2 <- eBayes(fit2)

# # Inshallah
# # topTable(fit2)
# # decideTests(fit2)


# # Bomb if bad
# aw <- arrayWeights(y,design)
# fit <- lmFit(y, design,
#              weights = aw)
# # contrasts <- makeContrasts(nt - nc, levels=design)
# fit2 <- contrasts.fit(fit, contrasts.matrix)
# fit2 <- eBayes(fit2)
# # End of outlier bombing

# full_results <- topTable(fit2, number=Inf)
# counts <- readRDS("CryptCounts.RDS")

# enssym <- counts[, c("Symbol", "id")]
# full_results <- merge(full_results, enssym, by.x='row.names', by.y="id", all=FALSE)

# library(tximportData)
# library(tximport)
# library(rhdf5)
# library(GenomicFeatures)
# library(tidyr)

# dir <- list.dirs("D:/W/P/China/KYS/Kallisto quant on collection 766 and collection 805_ Abundances (tabular)")
# n <- list.files("D:/W/P/China/KYS/Kallisto quant on collection 766 and collection 805_ Abundances (tabular)")
# files <- file.path(dir,n)
# files<-files[grepl("tabular", files)]

# files<-files[14:25]
# n<-n[14:25]

# names(files) <- files

# txdb <- makeTxDbFromGFF(file="D:/W/P/China/KYS/Mus_musculus.GRCm39.110.chr.gff3", format="gff3", dataSource="ENSMBL")
# saveDb(x=txdb, file = "genomic.TxDb")
# k <- keys(txdb, keytype = "TXNAME")
# tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME") 

# txi.kallisto.tsv  <- tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreAfterBar = TRUE, ignoreTxVersion=TRUE)
# head(txi.kallisto.tsv$counts)
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

# library(biomaRt)
# library(annotables)

# counts <- txi.kallisto.tsv$counts
# # input list of Ensembl ID's
# ensembl.ids <- rownames(counts)


# # method 1: biomaRt
# listEnsembl()
# ensembl <- useEnsembl(biomart = "genes")
# datasets <- listDatasets(ensembl)

# ensembl.con <- useMart("ensembl", dataset = 'mmusculus_gene_ensembl')

# attr <- listAttributes(ensembl.con)
# filters <- listFilters(ensembl.con)

# Martgenes <- getBM(attributes = c('ensembl_gene_id','external_gene_name'),
#                    filters = "ensembl_gene_id",
#                    values = ensembl.ids,
#                    mart = ensembl.con)

# counts <- as.data.frame(counts)
# counts$Symbol <- Martgenes$external_gene_name

# library(dplyr)

# # Assuming counts$gene_name is your vector
# gene_name_vector <- counts$Symbol  # Replace with your actual vector

# # Find duplicated elements
# duplicated_elements <- gene_name_vector[duplicated(gene_name_vector)]

# # Extract unique duplicated elements
# unique_duplicated_elements <- unique(duplicated_elements)


# # Sum rows by gene_name
# summed_counts <- counts %>%
#   group_by(Symbol) %>%
#   summarise_all(sum)

# count <- as.data.frame(summed_counts)


library(edgeR)

counts <- count

rownames(counts) <- count$Symbol
counts$Symbol <- NULL
counts$id <- NULL
colnames(counts) <- c("c1","t1","c2","c3","c4","c5","c6","t2","t3","t4","t5","t6")

d0 <- DGEList(counts)
d0 <- calcNormFactors(d0)

cutoff <- 1
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,] 
dim(d)

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
# contrasts <- makeContrasts(nt - nc, levels=design)
fit2 <- contrasts.fit(fit, contrasts.matrix)
fit2 <- eBayes(fit2)
# End of outlier bombing

full_results <- topTable(fit2, number=Inf)
full_results$Symbol <- rownames(full_results)

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

stori <- c(u$Symbol[1:15], d$Symbol[1:10], c("Mex3a", "Msi1", "Lgr5", "Smoc2", "Tert"))

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

stori <- stori[!stori %in% c("Cyp24a1", "Epha4", "D130043K22Rik", "Rnf149", "Ifngr1", "5033404E19Rik")]

#Cyp
# 19Rik
# D130043K22Rik
# Epha4
# Ifngr1

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
                  col = "black", box.padding = 1.5, max.overlaps = Inf, min.segment.length = unit(0, 'lines')) +  
  geom_segment(aes(x = logFC, y = -log10(adj.P.Val), xend = logFC, yend = -log10(adj.P.Val)), color = "grey50", alpha = 0.5) + # Add lines connecting labels to points
  scale_color_manual(values = c("blue", "black", "red")) +  
  scale_size_identity() +  
  labs(x = "log2FC", y = "-log10(adj.P.Val)") +   
  theme(text=element_text(size=14))

plot <- p + scale_x_break(c(6, 12.5), scales=0.15, expand = expansion(add=0.4), ticklabels = c(12.5, 13))

ggsave("full_results_crypt_scaled_manual.png", plot, width = 9.38, height = 7.5, units = "in", dpi = 300)

num_sig <- nrow(filter(full_results, adj.P.Val < 0.05))
num_up <- nrow(filter((filter(full_results, adj.P.Val < 0.05)), logFC > 0))
num_down <- nrow(filter((filter(full_results, adj.P.Val < 0.05)), logFC < 0))
