## Install and load required packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

library(DESeq2)
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

gene_names <- c(
  "IL6", "CXCL8", "IL1B", "IFNB1", "IFNA2", "STAT1", "STAT2",
  "IFIT2", "IFIT3", "CXCL1", "CXCL10", "OASL", "GBP4", "IF130", "PLEKHA4",
  "RSAD2", "MX2", "CMPK2", "HERC6", "BIRC3", "ISG20", "CSF3"
)

## Prepare heatmap
rld <- rlog(dds, blind = FALSE)

counts <- assay(rld) %>%
    as.data.frame %>%
    rownames_to_column(var = "Gene") %>%
    rename_at(vars(all_of(sample_names)), ~simple_sample_names)

selected_samples <- c(
  "Media Control",
  "OCS",
  "Peptide 2",
  "Poly I:C",
  "OCS + Poly I:C",
  "Peptide 2 + Poly I:C"
)

`%notin%` <- Negate(`%in%`)

target_genes <- counts %>%
  select(c(Gene, selected_samples)) %>%
  filter(Gene %in% gene_names) %>%
  drop_na

write.csv(target_genes, file.path(result_dir, "heatmap.csv"))
