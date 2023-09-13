


library(dreamlet)
library(muscat)
library(ExperimentHub)
library(zenith)
library(scater)
library(readr)
library(stringr)
library(readxl)

#library(loomR)

# load RNAseq data
setwd('~/Data/Alexi_Duhe/test1')
load("~/Data/Alexi_Duhe/test1/integrated_labeled.RData")
load("~/Data/Alexi_Duhe/test1/025_covariates_full_df.RData")

#Covariate table
MGS_table <- read_excel("/nvmefs/scARC_Duan_018/MGS_iPSC lines_60samples_scRNA_ATAC-seq_bulkATAC-seq status_AK_HZ.xlsx")
MGS_table$`new iPSC line ID` <- str_replace(MGS_table$`new iPSC line ID`,
                                            "CD00000",
                                            "CD_")

# 1 = cell line id code, 3 = affect, 4 = sex, 5 = batch number
covar_table <- MGS_table[, c(1, 3, 4, 5, 9)]

# in column `sex (1=M, 2=F)` replace row values where 1 with M and 2 with F
covar_table$`sex (1=M, 2=F)` <- str_replace(str_replace(covar_table$`sex (1=M, 2=F)`,
                                                        "1", "M"), "2", "F")
# rename columns
colnames(covar_table) <- c("cell_line", "aff", "sex", "age", "batch")

# remove rows where batch = 'NA'
covar_table <- covar_table[!is.na(covar_table$batch), ]

# remove values before '&' symbol and '& ' in column batch
covar_table$batch <- str_remove(covar_table$batch, "[0-9]+ & ")

#remove values before '(need redo) &' as well as '(need redo) &' from batch col
covar_table$batch <- str_remove(covar_table$batch, "[0-9]+ \\(need redo\\)  & ")

for (i in length(integrated.labeled)) {
  for (j in length(rownames(covar_table))) {
    if (integrated.labeled$cell.line.ident[i] == covar_table$cell_line[j]){
      print(covar_table[j], integrated.labeled$cell.line.ident[i])
  }
  }
}

for (i in length(rownames(covar_table))) {
  for (j in length(integrated.labeled)) {
    if (integrated.labeled@meta.data$cell.line.ident[i] == covar_table$cell_line[j]){
      print(covar_table$cell_line[j], integrated.labeled$cell.line.ident[i])
    }
  }
}


integrated.labeled$aff <- 'NA'
for (i in 1:length(rownames(covar_table))) {
  if (covar_table$cell_line[i] %in% unique(cts2$cell.line.ident)){
    print(covar_table$cell_line[i])
    for (j in 1:length(cts2$cell.line.ident)) {
      if (covar_table$cell_line[i] == cts2$cell.line.ident[j]) {
        cts2$aff[j] <- covar_table$aff[i]
  }
}}}







for (i in 1:length(rownames(covar_table))) {
  if (covar_table$cell_line[i] %in% unique(cts2$cell.line.ident)){
    line_ <- covar_table$cell_line[i]
    print(covar_table$cell_line[i])
    aff_rep <- covar_table$aff[i]
    print(covar_table$aff[i])

    for (line_ in cts2$cell.line.ident) {
      #print(line_)
      cts2$aff <- aff_rep
      #if (covar_table$cell_line[i] == cts2$cell.line.ident[j]) {
        #cts2$aff[j] <- covar_table$aff[j]
      }}}



for (i in 1:length(covar_table)) {
  if (covar_table$cell_line[i] %in% unique(cts2$cell.line.ident)){
    aff_rep <- covar_table$aff[i]
    print(covar_table$cell_line[i])
    }}





    for (j in cts2$cell.line.ident) {
      if (covar_table$cell_line[i] == cts2$cell.line.ident[j]) {
        cts2$aff[j] <- covar_table$aff[i]
      }
    }}}

for (i in cts2$cell.line.ident){
  print(i)
}







# chnage active assay from SCT to RNA
DefaultAssay(integrated.labeled) <- "RNA"

# Create dataframe from Seurat object for pseudobulk ----
Idents(integrated.labeled)
integrated.labeled <- subset(integrated.labeled, cell.type != "unidentified" & cell.type != 'immature')

# replace space in time ident with _
integrated.labeled$time.ident <- str_replace(integrated.labeled$time.ident, ' ', '_')
# create column containing celltype and time id
integrated.labeled$typextime <- paste0(integrated.labeled$cell.type, '_' , integrated.labeled$time.ident)

# pseudobulk method that does work (Large Matrix)
cts <- AggregateExpression(integrated.labeled,
                    group.by = c('typextime', 'cell.line.ident'),
                    assays = 'RNA',
                    slot = 'counts',
                    return.seurat = F)
# extract counts data from matrix
cts <- cts$RNA


#Add informations from covariate table




cts2 <- as.SingleCellExperiment(integrated.labeled)
factor(cts2$aff)
cts2$id <- paste0(cts2$aff, cts2$cell.line.ident)

cts2 <- aggregateToPseudoBulk(cts2,
                              assay = 'counts',
                              cluster_id = 'typextime',
                              sample_id = 'id')

res.proc <- processAssays(cts2, ~ aff, min.count = 5)




































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



