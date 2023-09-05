# Test 6.5 before limma
# use pseudobulk samples to compute variance partition

# init ####
library(Seurat)
library(Signac)
library(rtracklayer)
library(GenomicRanges)
library(readr)
library(BiocParallel)
library(edgeR)
library(variancePartition)

# load RNAseq data
setwd('~/Data/Alexi_Duhe/test1')
load("~/Data/Alexi_Duhe/test1/integrated_labeled.RData")

# RNAseq ####
# make pseudobulk samples
DefaultAssay(integrated.labeled) <- "RNA"
lines <- unique(integrated.labeled$cell.line.ident)
Idents(integrated.labeled)
times <- unique(integrated.labeled$time.ident)
integrated.labeled <- subset(integrated.labeled, cell.type != "unidentified" & cell.type != 'immature')
cell_type <- as.vector(unique(integrated.labeled$cell.type))

# count number of cells in each line, type, and time
i <- 1
j <- 1
k <- 1

for (j in 1:length(cell_type)) {
  for (k in 1:length(times)) {
    for (i in 1:length(lines)) {
      print(paste(cell_type[j], times[k], lines[i], sep = ";"))
      temp_count <-
        rowSums(subset(x = integrated.labeled,
                       subset = (cell.line.ident == lines[i]) & cell.type == cell_type[j] &
                         (time.ident == times[k])))
      if (i == 1 & j == 1 & k == 1) {
        export_df <- data.frame(temp_count,
                                stringsAsFactors = F)
        colnames(export_df) <- paste(cell_type[j],
                                     times[k],
                                     lines[i],
                                     sep = "_")
      } else {
        temp_colnames <- colnames(export_df)
        export_df <- data.frame(cbind(export_df,
                                      temp_count),
                                stringsAsFactors = F)
        colnames(export_df) <-
          c(temp_colnames, paste(cell_type[j],
                                 times[k],
                                 lines[i],
                                 sep = "_"))
      }
    }
  }
}

for (i in 1:length(colnames(export_df))) {
  colnames(export_df)[i] <- str_replace(colnames(export_df)[i], ' ', '')
}



colnames(export_df)
rownames(export_df)
save(export_df, file = "025_raw_pseudobulk_matrix_for_varpartition.RData")



