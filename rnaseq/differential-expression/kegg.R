## Install and load required packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

library(DESeq2)
library(ggplot2)
library(RColorBrewer)
library(tidyverse)

## Load data
source("data.R")

result_dir <- "output"
dir.create(result_dir, showWarnings = FALSE)

sample_names <- c(
  "HDMEC_Control_S31",
  "HDMEC_Poly_I_C_S32",
  "HDMEC_OCS_S35",
  "HDMEC_OCS_Poly_I_C_S36",
  "HDMEC_Pep2_S33",
  "HDMEC_Pep2_Poly_I_C_S34"
)

simple_sample_names <- c(
  "Media Control",
  "Poly I:C",
  "OCS",
  "OCS + Poly I:C",
  "Peptide 2",
  "Peptide 2 + Poly I:C"
)

conditions <- c(
  "Control",
  "PIC",
  "Control",
  "OCS_PIC",
  "Control",
  "P2_PIC"
)

dds <- load_data(file.path(data_dir, "counts.csv"), sample_names, conditions)
res <- results(dds, contrast = c("condition", "P2_PIC", "PIC"))

# Select significantly DE genes
res_sign <- subset(res, padj < 0.05) %>%
  as.data.frame %>%
  rownames_to_column(var = "Gene")

# Match gene names to Entrez ID
library(biomaRt)

mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

gene_id <- getBM(
  attributes = c("external_gene_name", "entrezgene_id"),
  filters = "external_gene_name",
  values = res_sign$Gene,
  mart = mart
)

manual_names <- (
  read.csv("unmatched.csv") %>%
    filter(Match.type != "Unmatched" | !is.na(EntrezID))
)[, c("Input", "Approved.symbol", "EntrezID")]

manual_ids <- getBM(
  attributes = c("external_gene_name", "entrezgene_id"),
  filters = "external_gene_name",
  values = (manual_names %>%
    filter(Approved.symbol != "" & is.na(EntrezID))
  )[, "Approved.symbol"],
  mart = mart
)

manual_entrez <- manual_names %>%
  mutate(EntrezID = ifelse(Approved.symbol %in% manual_ids$external_gene_name,
                           manual_ids$entrezgene_id, EntrezID)) %>%
  rename(Gene = Approved.symbol)

res_sign <- res_sign %>%
  mutate(EntrezID = gene_id[match(Gene, gene_id$external_gene_name), "entrezgene_id"]) %>% # nolint
  mutate(EntrezID = ifelse(
    Gene %in% manual_entrez$Input, manual_entrez$EntrezID, EntrezID
  )) %>%
  drop_na


# Overall changed pathways
library(clusterProfiler)

kegg_all <- enrichKEGG(
  gene = res_sign$EntrezID,
  organism = "hsa",
  pvalueCutoff = 0.01
)

dotplot(kegg_all) +
  theme(
    axis.text.x = element_text(size = 6, angle = 45, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(size = 6),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 6),
    plot.background = element_rect(fill = "white"),
    panel.background = element_blank()
  )

ggsave(file.path(result_dir, "dot.pdf"), scale = 1)