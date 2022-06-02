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

dds <- load_data(file.path("data", "counts.csv"), sample_names, conditions)
vsd <- vst(dds, blind = FALSE)

# the dot between P2 + PI:C and OCS + PI:C is HDMEC_Poly_I_C_S32
pca_data <- plotPCA(vsd, intgroup = c("condition"), returnData = TRUE)
pca_data$simple <- simple_sample_names

write.csv(pca_data, file.path(result_dir, "pca.data.csv"))
