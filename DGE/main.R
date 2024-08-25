library(tximport)
library(rhdf5)
library(tidyr)
library(dplyr)
library(GenomicFeatures)


dir <- list.dirs("D:/KYS/Kallisto quant on collection 766 and collection 805_ Abundances (tabular)")
n <- list.files("D:/KYS/Kallisto quant on collection 766 and collection 805_ Abundances (tabular)")
files <- file.path(dir,n)
files<-files[grepl("tabular", files)]

files<-files[14:25]
n<-n[14:25]

names(files) <- files

# txdb <- makeTxDbFromGFF(file="D:/KYS/Mus_musculus.GRCm39.110.chr.gff3", format="gff3", dataSource="ENSMBL")
# saveDb(x=txdb, file = "genomic.TxDb")
txdb <- loadDb("genomic.TxDb")
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME") 

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
# library(annotables)

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

# mart <- useEnsembl("ensembl","mmusculus_gene_ensembl", version = 110, host="https://asia.ensembl.org")

# Martgenes <- getBM(attributes=c("ensembl_gene_id", "external_gene_name"),
#            filters = 'ensembl_gene_id', values = ensembl.ids, mart = mart, bmHeader = TRUE)
# 
# setdiff(ensembl.ids, Martgenes$`Gene stable ID`)

# df <- as.data.frame(counts)
# df$`Gene stable ID` <- rownames(df)
# 
# counts_merged <- merge(df, Martgenes, by = "Gene stable ID", 
#                   all.x = TRUE)

# counts <- na.omit(counts_merged)
# counts$Symbol <- counts$`Gene name`
# counts$`Gene stable ID` <- NULL
# counts$`Gene name` <- NULL

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

count <- count[!(is.na(count$Symbol) | count$Symbol==""), ]
count$id <- Martgenes[match(count$Symbol, Martgenes$external_gene_name),]$ensembl_gene_id

saveRDS(count, "count_with_sym_and_id.rds")

# library(AnnotationHub)
# hub <- AnnotationHub()
# query(hub, c("mus musculus","ensdb"))
# ensdb <- hub[["AH113713"]]

# gns <- genes(ensdb)
# txs <- transcriptsBy(ensdb)

# txwid <- sapply(width(txs), max)
# gnwid <- setNames(width(gns), names(gns))
# id2len <- data.frame(TxLen = txwid, GeneLen = gnwid[names(txwid)])

# count$length <- id2len[count$id,]$TxLen

# colnames(count) <- c(colnames(count)[1], colnames(wow), colnames(count)[14], colnames(count)[15])

# library(DESeq2) #count[,2:13]
# df <- count[,2:13]
# rownames(df) <- count$Symbol

# col <- as.data.frame(c("c","t","c","c","c","c","c","t","t","t","t","t"))
# colnames(col) <- "cond"
# rownames(col) <- colnames(df)

# dds <- DESeqDataSetFromTximport(txi.kallisto.tsv, col, ~cond)
# fpkmres <- fpkm(dds)

# dds <- DESeq(dds)
# res <- results(dds, contrast=c("cond", "t", "c"))
# export <- as.data.frame(res)
