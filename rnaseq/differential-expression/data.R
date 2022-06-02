library(DESeq2)

load_data <- function(counts_file_path, sample_names, conditions) {
    counts <- as.matrix(read.csv(counts_file_path, row.names = "Gene"))
    counts <- counts[, sample_names]
    col.data <- data.frame("condition" = factor(conditions))
    row.names(col.data) <- sample_names

    ddsm <- DESeqDataSetFromMatrix(
        countData = round(counts),
        colData = col.data,
        design = ~condition
    )
    ddsm <- ddsm[rowSums(counts(ddsm)) > 0, ]

    return(DESeq(ddsm))
}
