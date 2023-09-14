


library(dreamlet)
library(muscat)
library(ExperimentHub)
library(zenith)
library(scater)
library(readr)
library(stringr)
library(readxl)
library(mashr)
library(zenith)

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


# chnage active assay from SCT to RNA
DefaultAssay(integrated.labeled) <- "RNA"
# Create dataframe from Seurat object for pseudobulk ----
Idents(integrated.labeled)
integrated.labeled <- subset(integrated.labeled, cell.type != "unidentified" & cell.type != 'immature')
# replace space in time ident with _
integrated.labeled$time.ident <- str_replace(integrated.labeled$time.ident, ' ', '_')

# Subset cell types to analyse differential expression separtly 
GABA_cts <- subset(integrated.labeled, cell.type == 'GABA')
nmglut_cts <- subset(integrated.labeled, cell.type == 'nmglut')
npglut_cts <- subset(integrated.labeled, cell.type == 'npglut')


apply_dream <- function(cts) {
  # transpose to sce
  cat('Transposing to single cell experement\n')
  cts <- as.SingleCellExperiment(cts)
  
  # pseudobulk
  cat('Aggregating to pseudobulk\n')
  cts <- aggregateToPseudoBulk(cts, 
                               assay = 'counts',
                               cluster_id = 'time.ident', # comparison variable 
                               sample_id = 'cell.line.ident') # cell ideintification 
  # normalize and apply voom 
  cat('Normalizing and applying voom with dreamlet weights\n')
  res.proc = processAssays(cts, ~ time.ident, min.count = 5)
  
  # differential expression analysis within each assay
  # evaluated on the voom normalized data
  cat('Evaluating the voom-normalized data\n')
  res.d1 = dreamlet(res.proc, ~ time.ident)
  
  # run mashr model to borrow information across genes and
  # cell tpes in estimating coefficients' posterior distribution
  res_mash = run_mash(res.d1, coef = '(Intercept)')
  return(res_mash)
  cat('########\n')
}


mashr_summary <- function(res_mash) {
  head(res_mash$logFC.original)
  
  # how many gene-by-celltype are sig
  cat('\nSignificant gene-by-celltype')
  print(table(get_lfsr(res_mash$model) < 0.05, useNA = 'ifany'))
  
  # how many genes are significant in at least 1 cell type
  cat('\nGenes that are significant in at least 1 cell type')
  print(table( apply(get_lfsr(res_mash$model), 1, min, na.rm=TRUE) < 0.05))
  
  # how many genes are significant in each cell type
  cat('\nSignificant genes in each time point\n')
  print(apply(get_lfsr(res_mash$model), 2, function(x) sum(x < 0.05, na.rm=TRUE)))
  
  # examine top set of genes
  # which genes are significant in at least 1 cell type
  cat('\nTop 10 genes are significant in at least 1 cell type\n')
  print(sort(names(get_significant_results(res_mash$model)))[1:10])
  
  # There is a lot of variation in the raw logFC
  cat('\nOriginal raw logFC \n')
  print(res_mash$logFC.original["ISG20",])
  
  # posterior mean after borrowing across cell type and genes
  cat('\nPosterior mean after borrowing accross cell types and genes\n')
  print(get_pm(res_mash$model)["ISG20",])
  
  # gene set analysis using mashr results
  
  # Load Gene Ontology database
  # use gene 'SYMBOL', or 'ENSEMBL' id
  # use get_MSigDB() to load MSigDB
  go.gs = get_GeneOntology("CC", to="SYMBOL")
  
  # valid values for statistic:
  # "tstatistic", "abs(tstatistic)", "logFC", "abs(logFC)"
  df_gs = zenith_gsa(res_mash, go.gs) 
  return(df_gs)
}
  

GABA_cts <- apply_dream(GABA_cts)
nmglut_cts <- apply_dream(nmglut_cts)
npglut_cts <- apply_dream(npglut_cts)

GABA_mashr <- mashr_summary(GABA_cts)

# Heatmap of results
plotZenithResults(GABA_mashr, 5, 1)

# forest plot based on mashr results
plotForest(GABA_cts, "ISG20")

# volcano plot based on mashr results
# yaxis uses local false sign rate (lfsr)
plotVolcano(GABA_cts)








































library(mashr)
head(GABA_cts$logFC.original)

head(res_mash$logFC.original)

# how many gene-by-celltype are sig
table(get_lfsr(res_mash$model) < 0.05, useNA = 'ifany')
# FALSE  TRUE  <NA>
#   5873 64747  6678

# how many genes are significant in at least 1 cell type
table( apply(get_lfsr(res_mash$model), 1, min, na.rm=TRUE) < 0.05)
# FALSE  TRUE
# 419    25347


# how many genes are significant in each cell type
apply(get_lfsr(res_mash$model), 2, function(x) sum(x < 0.05, na.rm=TRUE))
# GABA  nmglut  npglut
# 21833  20337  22577

# examine top set of genes
# which genes are significant in at least 1 cell type
sort(names(get_significant_results(res_mash$model)))[1:10]
# [1] "7SK.65"    "A1BG"      "A1BG-AS1"  "A1CF"      "A2m"       "A2M"
# [7] "A2M-AS1"   "A2ML1"     "A2ML1-AS1" "A2MP1"

# There is a lot of variation in the raw logFC
res_mash$logFC.original["ISG20",]
#       GABA     nmglut     npglut
# 2.02171075 0.62082978 0.09145955

# posterior mean after borrowing across cell type and genes
get_pm(res_mash$model)["ISG20",]
#      GABA    nmglut    npglut
# 2.0036516 0.6125717 0.1185369

# gene set analysis using mashr results

# Load Gene Ontology database
# use gene 'SYMBOL', or 'ENSEMBL' id
# use get_MSigDB() to load MSigDB
go.gs = get_GeneOntology("CC", to="SYMBOL")

# valid values for statistic:
# "tstatistic", "abs(tstatistic)", "logFC", "abs(logFC)"
df_gs = zenith_gsa(res_mash, go.gs)

# Heatmap of results
plotZenithResults(df_gs, 5, 1)

# forest plot based on mashr results
plotForest(res_mash, "ISG20")

# volcano plot based on mashr results
# yaxis uses local false sign rate (lfsr)
plotVolcano(res_mash)












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



