


library(dreamlet)
library(muscat)
library(ExperimentHub)
library(zenith)
library(scater)
library(readr)
library(stringr)

library(loomR)

# load RNAseq data
setwd('~/Data/Alexi_Duhe/test1')
load("~/Data/Alexi_Duhe/test1/integrated_labeled.RData")

# chnage active assay from SCT to RNA
DefaultAssay(integrated.labeled) <- "RNA"

# Create dataframe from Seurat object for pseudobulk ----
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
  colnames(export_df)[i] <- str_replace(colnames(export_df)[i], ' ', '_')
}

#only keep singlet cells with suffucuent reads
test1 <- as.SingleCellExperiment(integrated.labeled) # create sce from seurat obj
test1 <- test1[rowSums(counts(test1) > 0) < 0, ]
test1 <- test1[ ,colData(test1)$multiplets == 'singlet']

# comput QC metrics
qc <- perCellQCMetrics(test1)

# remove cells with few or many detected genes
ol <- isOutlier(metric = qc$detected, nmads = 2)
test1 <- test1[ , !ol]

#compute normalized data
test1 <- test1[rowSums(counts(test1) > 1) >=10, ]
test1 <- computeLibraryFactors(test1)
test1 <- logNormCounts(test1)

# set variable indicating stim or ctrl
test1$StimStatus = test1$stim






# # Xandy Dreamlet Vignette
# # https://rdrr.io/github/GabrielHoffman/dreamlet/f/vignettes/mashr.Rmd
# # https://satijalab.org/seurat/archive/v3.1/conversion_vignette
#
# library(dreamlet)
# library(muscat)
# library(ExperimentHub)
# library(zenith)
# library(scater)
#
# #Preprocess ----
# # download data
# eh <- ExperimentHub()
# sce <- eh[['EH2259']]
#
#
# #only keep singlet cells with suffucuent reads
# sce <- sce[rowSums(counts(sce) > 0) < 0, ]
# sce <- sce[ ,colData(sce)$multiplets == 'singlet']
#
# # comput QC metrics
# qc <- perCellQCMetrics(sce)
#
# # remove cells with few or many detected genes
# ol <- isOutlier(metric = qc$detected, nmads = 2)
# sce <- sce[ , !ol]
#
# #compute normalized data
# sce <- sce[rowSums(counts(sce) > 1) >=10, ]
# sce <- computeLibraryFactors(sce)
# sce <- logNormCounts(sce)
#
# # set variable indicating stim or ctrl
# sce$StimStatus = sce$stim



