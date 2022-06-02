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

# P2 + PIC vs PIC
res <- results(dds, contrast = c("condition", "P2_PIC", "PIC"))
write.csv(res %>% as.data.frame %>% rownames_to_column("Gene"),
          file.path(result_dir, "p2pic-pic.de.csv"))

# OCS + PIC vs PIC
res <- results(dds, contrast = c("condition", "OCS_PIC", "PIC"))
write.csv(res %>% as.data.frame %>% rownames_to_column("Gene"),
          file.path(result_dir, "ocspic-pic.de.csv"))

# PIC vs Control
res <- results(dds, contrast = c("condition", "PIC", "Control"))
write.csv(res %>% as.data.frame %>% rownames_to_column("Gene"),
          file.path(result_dir, "pic-control.de.csv"))
