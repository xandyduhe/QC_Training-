# /nvmefs/scARC_Duan_018/Duan_project_025_RNA/Analysis_part2_GRCh38/20230507_demultiplex_QC_on_025_17_46.R
# /nvmefs/scARC_Duan_018/Duan_project_025_RNA/Analysis_part2_GRCh38/20230507_demultiplex_QC_on_025_17_46.R

library(Seurat)
library(RColorBrewer)
library(viridis)
library(ggplot2)
library(stringr)
library(readr)

# set working directories
setwd("/data/Alexi_Duhe/test1")

# load data file created in QC_test_1_separate_cells.R
load('removed_rat_mapped_hybrid.RData')

# locate data ----
# list of the paths to the files containing gene expression data
h5.list <- list.files('/data/FASTQ/Duan_Project_025_hybrid',
                      'filtered_feature_bc_matrix.h5',
                      recursive = TRUE,
                      include.dirs = FALSE,
                      full.names = TRUE)

# Initialize variables and data structures
objlist <- vector(mode = 'list', length = length(h5.list))
name.id <- str_extract(h5.list, '[0-9]+-[0-6]')#group+time
# extract the time identifier from path list created in h5.list at ith file
time.id <- sub('.*-', '', name.id)
# extract the group identifier from path list created in h5.list at ith file
group.id <- sub('-.*', '', name.id)

cell.counts <- data.frame(library = name.id,
                          total = rep(0, by = length(h5.list)),
                          rat = rep(0, by = length(h5.list)),
                          nonrat = rep(0, by = length(h5.list)),
                          human = rep(0, by = length(h5.list)),
                          nonhuman = rep(0, by = length(h5.list)))

# create empty list to store data
no.rat.list <- vector(mode = 'list', length = length(h5.list))

for (i in 1:length(objlist)) {

  h5.file <- Read10X_h5(filename = h5.list[i])
  #
  obj <- CreateSeuratObject(counts = h5.file$`Gene Expression`,
                            project = name.id[i])

  #assign rat identity
  # get column names (barcodes) of human cells from the cleaned Seurat object
  human.bc <- colnames(clean.objlist[[i]])

  obj$time.ident <- paste0(time.id[i], ' hr')
  obj$rat.ident <- 'rat'
  obj$rat.ident[colnames(obj) %in% human.bc] <- 'human' ##

  cell.counts$total[i] <- ncol(obj)
  cell.counts$rat[i] <- sum(obj$rat.ident == 'rat')
  cell.counts$nonrat[i] <- sum(obj$rat.ident == 'human')
  # cell.counts$human[i] <- sum([colnames(obj) %in% human.bc])

  objlist[[i]] <- obj

  no.rat.obj <- subset(obj, rat.ident == 'human')

  no.rat.list[[i]] <- no.rat.obj
}

path.list <- sort(list.files(path = '/nvmefs/scARC_Duan_025_GRCh38/barcodes_demuxed_by_library',
                             full.names = TRUE,
                             pattern = 'best.tsv',
                             recursive = TRUE))

barcode.list <- vector(mode = 'list', length = length(path.list))

for (i in 1:length(path.list)) {
  #
  #   barcode.list[i] <- read.delim(path.list[i],
  #                                 header = FALSE,
  #                                 row.names = NULL)
  barcode.list.read <- read.delim(path.list[i],
                                  header = FALSE,
                                  row.names = NULL)

  barcode.list[[i]] <- barcode.list.read

  colnames(barcode.list[[i]]) <- c('barcode', 'line')
}

names(barcode.list) <-
  paste0(group.id)
# names(barcode.list) <-
#   paste0(str_remove(str_extract(path.list, '[0-2].best'),
#                     '\\.best'))

names(no.rat.list) <- name.id
# names(no.rat.list) <- str_extract(h5.list,
#                                   '[0-9]+-[0-6]')

for (i in 1:length(barcode.list)) {

  lines <- unique(barcode.list[[i]]$line)

  no.rat.list[[i]]$cell.line.ident <- 'unmatched'

  for (k in 1:length(lines)) {
    line.specif.barcode <-
      barcode.list[[i]]$barcode[barcode.list[[i]]$line == lines[k]]

    no.rat.list[[i]]$cell.line.ident[no.rat.list[[i]]@assays$RNA@counts@Dimnames[[2]]
                                     %in% line.specif.barcode] <- lines[k]
  }
  cell.counts$human[i] <-
    sum(no.rat.list[[i]]$cell.line.ident != 'unmatched')

  cell.counts$nonhuman[i] <-
    sum(no.rat.list[[i]]$cell.line.ident == 'unmatched')
}

write.table(cell.counts,
            file = 'cell_counts_from_GRCh38_mapped_data.csv',
            sep = ',',
            quote = FALSE,
            col.names = TRUE,
            row.names = FALSE)



# remove cells not matched to cell line barcodes

human.list <- vector(mode = 'list',
                     length = length(objlist))

for (i in 1:length(objlist)) {
  human.list[[i]] <-
    subset(no.rat.list[[i]],
           cell.line.ident != 'unmatched')
}

# QC ----

# merge w out correcting

for (i in 1:length(human.list)) {

  human.list[[i]][['percent.mt']] <- PercentageFeatureSet(human.list[[i]],
                                                          pattern = '^MT-')
}

for (i in 1:(length(human.list) - 1)) {

  if (i == 1) {
    rough.merged <- merge(human.list[[i]],
                          human.list[[i + 1]])
  }else {
    rough.merged <- merge(rough.merged,
                          human.list[[i + 1]])
  }
}

VlnPlot(rough.merged,
        features = c('nFeature_RNA',
                     'nCount_RNA',
                     'percent.mt'),
        ncol = 3,
        fill.by = feature,
        pt.size = 0)

# TSS Enrichment
sum(rough.merged$percent.mt > 20)
sum(rough.merged$nFeature_RNA > 7500)
sum(rough.merged$nCount_RNA > 40000)
sum(rough.merged$nFeature_RNA < 300)
sum(rough.merged$nCount_RNA < 500)

QCed.list <- vector(mode = 'list',
                    length = length(human.list))

for (i in 1:length(human.list)) {

  QCed.list[[i]] <- subset(human.list[[i]],
                           subset = nFeature_RNA > 300 &
                             nFeature_RNA < 7500 &
                             nCount_RNA > 500 &
                             nCount_RNA < 40000 &
                             percent.mt < 20)
}

save(objlist, file = 'GRCH38_mapped_raw_list.RData')
save(no.rat.list,
     file = 'GRCh38_mapped_removed_rat_assigned_demux.RData')
save(human.list, file = 'GRCh38_mapped_matched_human_demux.RData')
save(QCed.list, file = 'GRCh38_mapped_after_QC_list.RData')


#---- Process Notes -----
# This R code processes scRNA-seq data using the Seurat package.
# It performs quality control (QC) and data filtering,
# demultiplexes the data based on barcodes, and applies further QC.

# 1. **Loading Libraries and Setting Working Directory:**
#    - Libraries including Seurat, RColorBrewer, viridis, ggplot2,
#      stringr, and readr are loaded.
#    - The working directory is set to a specific path.
#
# 2. **Loading Cleaned Data from Previous Analysis:**
#    - The cleaned Seurat objects from a previous analysis (from the file
#      'removed_rat_mapped_hybrid.RData') are loaded.
#
# 3. **Locating Data Files:**
#    - A list of paths to data files (in this case,
#      'filtered_feature_bc_matrix.h5' files) is created.
#
# 4. **Creating Seurat Objects and Assigning Cell Identities:**
#    - For each data file, a Seurat object is created using `Read10X_h5`.
#    - Cell identities (rat, human) are assigned based on a list of human
#      barcodes from the previous analysis.
#
# 5. **Storing Cell Counts:**
#    - A data frame called `cell.counts` is created to store various cell counts,
#      including total, rat, non-rat, human, and non-human cells.
#
# 6. **Subsetting Data for Non-Rat Cells:**
#    - A new list called `no.rat.list` is created to store Seurat objects
#      containing only non-rat cells.
#    - Looping through each Seurat object, a subset containing only human cells
#      is created and added to `no.rat.list`.
#
# 7. **Saving Processed Data:**
#    - The list of Seurat objects is saved to an RData file
#     ('GRCH38_mapped_raw_list.RData').
#
# 8. **Processing Demultiplexed Barcodes:**
#    - A list of paths to barcode files is created.
#    - Barcodes are read from each file and stored in `barcode.list`, with the
#      corresponding group identifiers as names.
#
# 9. **Matching Barcodes to Cell Lines:**
#    - Cell lines are matched to barcodes in the `no.rat.list` Seurat objects
#      based on specific barcodes for each line.
#    - Cell counts are updated based on the matched barcodes.
#
# 10. **Writing Cell Counts to a CSV File:**
#     - The cell count information is written to a CSV file
#       ('cell_counts_from_GRCh38_mapped_data.csv').
#
# 11. **Saving Processed Data with Rat Cells Removed:**
#     - The list of Seurat objects containing only matched human cells is saved
#       to an RData file ('GRCh38_mapped_removed_rat_assigned_demux.RData').
#
# 12. **Further Quality Control:**
#     - For each Seurat object in `human.list`, the percentage of MT reads is
#       calculated and added as a feature.
#     - A loop merges Seurat objects using the `merge` function, which is not
#       further explained.
#
# 13. **Generating QC Plots:**
#     - Violin plots of different features (nFeature_RNA, nCount_RNA, percent.mt)
#       are created using `VlnPlot`.
#     - Various QC thresholds are calculated using the `sum` function.
#
# 14. **Applying Additional Quality Control:**
#     - A new list called `QCed.list` is created to store Seurat objects that
#       pass additional QC criteria.
#     - The Seurat objects in `human.list` are subsetted based on specified QC
#       thresholds and added to `QCed.list`.
#
# 15. **Saving Processed Data after QC:**
#     - The list of Seurat objects passing QC is saved to an RData file
#       ('GRCh38_mapped_after_QC_list.RData').
#
# This code performs data filtering, cell demultiplexing, and quality control
# on scRNA-seq data. It involves several steps of processing and filtering,
# with the final output being a list of Seurat objects representing human
# cells that have undergone QC.







