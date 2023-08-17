# Chuxuan Li 11/07/2022
# demultiplex, remove nonhuman, QC

# init #### XXXX
library(Seurat)
library(Signac)
library(RColorBrewer)
library(viridis)
library(ggplot2)
library(stringr)
library(readr)
library(readxl)
library(patchwork)

setwd('/data/Alexi_Duhe/test1_ATAC')

load("/data/Alexi_Duhe/test1/GRCh38_mapped_matched_human_demux.RData")
load("/data/Alexi_Duhe/test1/GRCh38_mapped_removed_rat_assigned_demux.RData")
load('/data/Alexi_Duhe/test1/GRCh38_mapped_after_QC_list.RData')
load("~/NVME/scARC_Duan_018/R_scARC_Duan_018/codes/EnsDb_UCSC_hg38_annotation.RData")

# read h5, remove rat cells and genes ####
h5.list <- list.files('/data/FASTQ/Duan_Project_025_hybrid',
                      'filtered_feature_bc_matrix.h5',
                      recursive = TRUE,
                      include.dirs = FALSE,
                      full.names = TRUE)


frag.list <- sort(list.files(path = '/data/FASTQ/Duan_Project_025_hybrid',
                             pattern = "atac_fragments.tsv.gz$",
                             recursive = T,
                             include.dirs = F,
                             full.names = T))

ATAC.list <- vector(mode = "list", length = length(h5.list))

name.id <- str_extract(h5.list, "[0-9]+-[0-6]")

df.cell.counts <- data.frame(library = name.id,
                             total = rep(0, by = length(name.id)),
                             rat = rep(0, by = length(name.id)),
                             nonrat = rep(0, by = length(name.id)),
                             human = rep(0, by = length(name.id)),
                             nonhuman = rep(0, by = length(name.id)))

time.id <- sub('.*-', '', name.id, paste0('hr'))
# create group variable: extract the group identifier from path list created in h5.list at ith file
group.id <- sub('-.*', '', name.id)
# create names vector to store names variables in vector


#----
# for (i in 1){
#   h5.file <- Read10X_h5(filename = h5.list[i])
#   cat("h5: ", str_extract(string = h5.list[i],
#                           pattern = "[0-9]+-[0-6]"))
#   frag.path <- frag.list[i]
#   
#   
#   chromAssay <- CreateChromatinAssay(counts = h5.file,
#                                      sep = c(":", "-"),
#                                      fragments = frag,
#                                      annotation = ens_use)
#   
#   chrom.obj <- CreateSeuratObject(counts = chromAssay,
#                                   assay = 'peaks',
#                                   project = name.id[i])
#   
# }
# 
# 
# chrom.obj <-CreateSeuratObject(counts = h5.file$`Gene Expression`,
#                                assay = 'RNA',
#                                project = name.id[i])
# 
# chrom.obj[['ATAC']] <- CreateChromatinAssay(counts = h5.file$Peaks,
#                                             sep = c(':', '-'),
#                                             fragments = frag.path,
#                                             annotation = ens_use)
#----


for (i in 1:length(ATAC.list)){
  h5.file <- Read10X_h5(filename = h5.list[i])
  cat("h5: ", str_extract(string = h5.list[i],
                          pattern = "[0-9]+-[0-6]"))
  frag.path <- frag.list[i]
  
  rna.dat <- CreateSeuratObject(counts = h5.file$`Gene Expression`,
                                assay = 'RNA')
  
  rna.dat[['ATAC']] <- CreateChromatinAssay(counts = h5.file$Peaks,
                                                     sep = c(':', '-'),
                                                     fragments = frag.path,
                                                     annotation = ens_use)
  ATAC.list[[i]] <- rna.dat
}



# for (i in 1:length(ATAC.list)){
#   cat("h5: ", str_extract(string = h5.list[i],
#                           pattern = "[0-9]+-[0-6]"))
#   h5.file <- Read10X_h5(filename = h5.list[i])
#   
#   frag.path <- frag.list[i]
#   
#   DefaultAssay(QCed.list[[i]]) <- 'RNA'
#   
#   
#   rna.dat[['ATAC']] <- CreateChromatinAssay(counts = h5.file$Peaks,
#                                             sep = c(':', '-'),
#                                             fragments = frag.path,
#                                             annotation = ens_use)
#   
#   QCed.list[[i]] <- rna.dat
#   
# }



# Quality control ----


for (i in 1:length(ATAC.list)) {
  DefaultAssay(ATAC.list[[i]]) <- 'ATAC'
  ATAC.list[[i]] <- NucleosomeSignal(ATAC.list[[i]])
  ATAC.list[[i]] <- TSSEnrichment(ATAC.list[[i]])
}

# filter out low quality cells 


for (i in 1:length(ATAC.list)) {
  
  ATAC.list[[i]] <- PercentageFeatureSet(ATAC.list[[i]],
                                       pattern = "^MT-",
                                       col.name = 'percent.mt')
  
  ATAC.list[[i]] <- subset(ATAC.list[[i]], 
                           subset = nFeature_ATAC > 300 & 
                             nFeature_ATAC < 7500 &
                             nCount_ATAC > 500 &
                             nCount_ATAC < 40000 &
                             percent.mt < 20 &
                             nucleosome_signal < 2 &
                             TSS.enrichment > 2)
}

# Peak Calling ----

for (i in 1:length()) {
  peaks <- CallPeaks(ATAC.list[[i]])
  
  
  
  
  
  
  
}







use.intv <-
  chrom.obj@assays$ATAC@counts@Dimnames[[1]][str_detect(chrom.obj@assays$ATAC@counts@Dimnames[[1]], "^chr")]


counts <-
  chrom.obj@assays$ATAC@counts[(which(rownames(chrom.obj@assays$ATAC@counts) %in% use.intv)), ]



# chrom.obj <- subset(chrom.obj@assays$ATAC, features = rownames(counts))


df.cell.counts$total[i] <- ncol(chrom.obj)
# assign library identity
chrom.obj$time.ident <- paste0(time.id[i], "hr")

print(unique(chrom.obj$time.ident))

# assign rat identity
humanbc <- colnames(human.list[[i]])


chrom.obj$rat.ident <- "rat"
chrom.obj$rat.ident[colnames(chrom.obj) %in% humanbc] <- "human"


df.cell.counts$rat[i] <- sum(chrom.obj$rat.ident == "rat")
df.cell.counts$nonrat[i] <- sum(chrom.obj$rat.ident == "human")


chrom.obj$cell.line.ident <- "unmatched"
lines <- unique(human.list[[i]]$cell.line.ident)


for (l in lines) {
  
  chrom.obj$cell.line.ident[colnames(chrom.obj) %in% colnames(human.list[[i]])[human.list[[i]]$cell.line.ident == l]] <- l
}

df.cell.counts$human[i] <- sum(chrom.obj$cell.line.ident != "unmatched")
df.cell.counts$nonhuman[i] <- sum(chrom.obj$cell.line.ident == "unmatched")

ATAC.list[[i]] <- chrom.obj

}


write.table(df.cell.counts,
            file = "ATAC_df.cell.counts_matched_to_RNAseq.csv",
            sep = ",",
            quote = F,
            col.names = T,
            row.names = F)


# 1. summarize number of cells and counts per library ####
plotCellCountsRectangles <- function(df) {
  p <- ggplot(df,
              aes(x = library,
                  y = rep_len(2, length(name.id)),
                  fill = cell.counts)) +
    geom_col(width = 1,
             position = position_dodge(width = 0.5)) +
    geom_text(aes(label = cell.counts),
              nudge_y = -1,
              hjust = "middle") +
    scale_fill_continuous(low = "#a2c1f2",
                          high = "#5c96f2") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5,
                                    size = 11),
          axis.title = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks = element_blank()) +
    ggtitle("# Nuclei") +
    coord_flip() +
    NoLegend()
  print(p)
}


plotUMIBars <- function(df) {
  p <- ggplot(df,
              aes(x = library,
                  y = reads.per.cell)) +
    ylab("# UMIs per nuclei") +
    geom_col(position = "dodge",
             color = "gray27",
             fill = "gray90") +
    theme_minimal() +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.x = element_text(size = 11),
          axis.text.x = element_text(size = 8)) +
    coord_flip()
  print(p)
}

plotTotalUMIBars <- function(df) {
  p <- ggplot(df,
              aes(x = library,
                  y = n.reads)) +
    ylab("# total UMIs") +
    geom_col(position = "dodge",
             color = "gray27",
             fill = "gray90") +
    theme_minimal() +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.x = element_text(size = 11),
          axis.text.x = element_text(size = 8)) +
    coord_flip()
  print(p)
}

df_raw <- data.frame(library = name.id,
                     cell.counts = rep_len(0, length(name.id)),
                     n.reads = rep_len(0, length(name.id)),
                     reads.per.cell = rep_len(0, length(name.id)))
df_human <- data.frame(library = name.id,
                       cell.counts = rep_len(0, length(name.id)),
                       n.reads = rep_len(0, length(name.id)),
                       reads.per.cell = rep_len(0, length(name.id)))


for (i in 1:length(name.id)){
  print(name.id[i])
  rawobj <- ATAC.list[[i]]
  humanonlyobj <-
    subset(rawobj, cell.line.ident != "unmatched")
  
  df_raw$cell.counts[i] <-
    length(rawobj@assays$Peaks@counts@Dimnames[[2]])
  
  df_raw$n.reads[i] <-
    sum(rawobj@assays$Peaks@counts@x)
  
  df_raw$reads.per.cell[i] <-
    sum(rawobj@assays$Peaks@counts@x)/df_raw$cell.counts[i]
  
  df_human$cell.counts[i] <-
    length(humanonlyobj@assays$Peaks@counts@Dimnames[[2]])
  
  df_human$n.reads[i] <-
    sum(humanonlyobj@assays$Peaks@counts@x)
  
  df_human$reads.per.cell[i] <-
    sum(humanonlyobj@assays$Peaks@counts@x)/df_raw$cell.counts[i]
}
write.table(df_raw,
            file = "./raw_df_cell_counts_UMI_counts_per_lib.csv",
            quote = F,
            sep = ",",
            row.names = F,
            col.names = T)

write.table(df_human,
            file = "./human_df_cell_counts_UMI_counts_per_lib.csv",
            quote = F,
            sep = ",",
            row.names = F,
            col.names = T)

df_raw$library <- factor(df_raw$library,
                         levels = rev(name.id))

df_human$library <- factor(df_human$library,
                           levels = rev(name.id))

plotCellCountsRectangles(df_human)
plotUMIBars(df_human)
plotTotalUMIBars(df_human)

save(ATAC.list,
     file = "ATAC_uncombined_obj_list_for_QC.RData")

# plot TSS and nucleosome signal for each library ####
# add annotation
for (i in 1:length(ATAC.list)){
  Signac::Annotation(ATAC.list[[i]]) <- ens_use
}
# compute QC metrics
for (i in 1:length(ATAC.list)){
  # ATAC.list[[i]] <- NucleosomeSignal(ATAC.list[[i]])
  # ATAC.list[[i]] <- TSSEnrichment(object = ATAC.list[[i]], fast = FALSE)
  #
  # pdf(file = paste0("./QC/TSS_plots/",
  #                   str_extract(h5.list[[i]],
  #                               "[0-9]+-[0|1|6]"),
  #                   "_TSS_plot.pdf"))
  
  p <- TSSPlot(ATAC.list[[i]]) +
    NoLegend() +
    theme(text = element_text(size = 14))
  print(p)
  # dev.off()
  
  # pdf(file = paste0("./QC/nuc_signal_histograms/",
  #                   str_extract(h5.list[[i]],
  #                               "[0-9]+-[0|1|6]"),
  #                   "_nucleosome_signal_hist.pdf"))
  
  p <- FragmentHistogram(ATAC.list[[i]]) +
    theme(text = element_text(size = 14))
  print(p)
  # dev.off()
}


TSS_sum <- 0
num_cells <- 0

for (i in cleanobj_lst){
  TSS_sum <- TSS_sum + sum(i$TSS.enrichment)
  num_cells <- num_cells + length(colnames(i))
}
TSS_sum / num_cells

nuc_sig_sum <- 0
num_cells <- 0
for (i in cleanobj_lst){
  print(sum(is.infinite(i$nucleosome_signal)))
  print(sum(is.nan(i$nucleosome_signal)))
  
  ns <- i$nucleosome_signal[!is.infinite(i$nucleosome_signal)]
  ns <- ns[!is.nan(ns)]
  nuc_sig_sum <- nuc_sig_sum + sum(ns)
  num_cells <- num_cells + length(colnames(i))
}
nuc_sig_sum / num_cells
