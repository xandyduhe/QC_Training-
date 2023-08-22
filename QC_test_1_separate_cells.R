# Alexandra Duhe 8/3/2023 9:00 AM
#~/NVME/scARC_Duan_018/Duan_project_025_RNA/Analysis_part1_hybrid_genome/20230503_separate_human_from_rat_in_025_17_and_46_new_libraries-copy.R
#weds

# load libraries
library(Seurat)
library(glmGamPoi)
library(sctransform)

library(ggplot2)
library(RColorBrewer)
library(viridis)
library(graphics)
library(ggrepel)

library(dplyr)
library(readr)
library(stringr)
library(readxl)

library(future)
library(tidyverse)
library(gridExtra)
library(scCustomize)
library(umap)


plan("multisession", workers = 8)
options(expressions = 20000)
options(future.globals.maxSize = 207374182400)



setwd("/data/Alexi_Duhe/test1")
# locate data ----
h5.list <- list.files('/data/FASTQ/Duan_Project_025_hybrid',
                      'filtered_feature_bc_matrix.h5',
                      recursive = TRUE,
                      include.dirs = FALSE,
                      full.names = TRUE)

# Merge data files using for-loop ----

raw.hybrid <- vector(mode = 'list', length = length(h5.list))
name.id <- str_extract(h5.list,
                       pattern = '[0-9]+-[0-6]') #group+time
time.id <- sub('.*-', '', name.id)
group.id <- sub('-.*', '', name.id)

for (i in 1:length(h5.list)) {

  h5file <- Read10X_h5(filename = h5.list[i],
                       use.names = T,
                       unique.features = T)

  raw.hybrid[[i]] <- CreateSeuratObject(counts = h5file$`Gene Expression`,
                                        project = name.id)
  raw.hybrid[[i]]$time.ident <- time.id[[i]]
  raw.hybrid[[i]]$group.ident <- group.id[[i]]

  print(paste0(i, " number of genes: ", nrow(raw.hybrid[[i]]),
               ", number of cells: ", ncol(raw.hybrid[[i]])))
}

# 1. %MT read ----
trans.obj <-raw.hybrid

for (i in 1:length(trans.obj)) {

print(paste0('Now at ', i))

  # 1. % MT reads
  trans.obj[[i]] <- PercentageFeatureSet(trans.obj[[i]],
                                       pattern = "^MT-",
                                       col.name = 'percent.mt')

  trans.obj[[i]] <- FindVariableFeatures(trans.obj[[i]], nfeatures = 3000)

  # 3. Normalize & Identify highly variable features
  trans.obj[[i]] <- SCTransform(trans.obj[[i]],
                           vars.to.regress = 'percent.mt',
                           method = 'glmGamPoi',
                           variable.features.n = 8000,
                           return.only.var.genes = FALSE,
                           seed.use = 2022,
                           verbose = TRUE)
}

# 5. Clustering ----

ffeatures <- c('Gfap', 'S100b', 'GAD1', 'SLC17A6')
scale.obj <-trans.obj

for (i in 1:length(scale.obj)) {

  # 4. Principal Component Analysis
  scale.obj[[i]] <- RunPCA(scale.obj[[i]],
                         verbose = TRUE,
                         seed.use = 2022)

  scale.obj[[i]] <- RunUMAP(object = scale.obj[[i]],
                          reduction = 'pca',
                          dims = 1:30,
                          seed.use = 2022)

  scale.obj[[i]] <- FindNeighbors(scale.obj[[i]],
                                reduction = 'pca',
                                dims = 1:30)

  scale.obj[[i]] <- FindClusters(scale.obj[[i]],
                              resolution = c(0.5),
                              random.seed = 2022)

    # Dimension Plot
  png(paste0(name.id[[i]], '_dimplot_by_cluster.png'))
  plot.dim <- DimPlot(scale.obj[[i]],
                           group.by = 'seurat_clusters',
                           label = TRUE,
                           reduction = 'umap') +
              ggtitle(name.id[i], ' plot by cluster') +
              theme(test = element_text(size = 12))
  print(plot.dim)
  dev.off()

    # Feature Plot
  colorpal <- viridis(n = 10,
                        option = 'C',
                        direction = -1)
  png(paste0(name.id[[i]], '_specific_genes.png'))
  plot.feature <- FeaturePlot(scale.obj[[i]],
                                features = ffeatures, ncol = 2)
  print(plot.feature)
  dev.off()

    # Violin plot
  png(paste0(name.id[[i]], '_Gfap_vln.png'))
  plot.violin <- VlnPlot_scCustom(scale.obj[[i]],
                                    features = 'Gfap') +
                 scale_y_continuous(trans = 'log10') +
                 ggtitle(paste0(name.id[[i]], ' Gfap expression by cluster'))
  print(plot.violin)
  dev.off

  print(plot.feature + plot.dim + plot.violin)
}

rrc.at.group <- list(c(10,8,13),  #objlist[[1]]  9-0
                     c(6,12),     #objlist[[2]]  9-1
                     c(6),        #objlist[[3]]  9-6

                     c(9),   #objlist[[4]]  13-0
                     c(10),  #objlist[[5]]  13-1
                     c(11),  #objlist[[6]]  13-6

                     c(7),   #objlist[[7]]  21-0
                     c(9),   #objlist[[8]]  21-1
                     c(9),   #objlist[[9]]  21-6

                     c(10),  #objlist[[10]] 23-0
                     c(10),  #objlist[[11]] 23-1
                     c(9),   #objlist[[12]] 23-6

                     c(7,14,2),  #objlist[[13]] 53-0
                     c(5,8,14),  #objlist[[14]] 53-1
                     c(9,7))     #objlist[[15]] 53-6

# Check and remove rat cells manually by looking at plots----

rat.obj <-scale.obj
clean.objlist <- vector(mode = 'list', length = length(rat.obj))

for (i in 1:length(rat.obj)) {

  rrc.at.group <- remove.rat.clusters[[i]]
  cat('rat cluster number: ', rrc.at.group[[i]])

	rat.obj[[i]]$rat.ident <- 'human'
	rat.obj[[i]]$rat.ident[rat.obj[[i]]$seurat_clusters %in% rrc.at.group[[i]]] <- 'rat'

	human.barcodes <- colnames(objlist[[i]])[objlist[[i]]$rat.ident == 'human']

	write.table(human.barcodes,
				file = paste0('.',
				              name.id[i],
				              '_human_barcodes.txt'),
				sep = '\t',
				quote = FALSE,
				row.names = FALSE,
				col.names = FALSE)

	clean.objlist[[i]] <- subset(rat.obj[[i]], subset = rat.ident != 'rat')
}

rat.cell.counts <- matrix(nrow = length(rat.obj),
                          ncol = 3,
                          dimnames = list(name.id,
                                          c('number of rat astrocytes',
                                            'total number of cells',
                                            'percentage of astrocyte cells (%)')))

for (i in 1:length(rat.obj)) {

  # calculate the number fo rat cells
  rat.cell.counts[i , 1] <- (sum(rat.obj[[i]]$rat.ident == 'rat'))
  # calculate the number of cells
  rat.cell.counts[i , 2] <- (length(rat.obj[[i]]$rat.ident))
  # calculate the % of cells
  rat.cell.counts[i , 3] <- ((rat.cell.counts[i,1]/rat.cell.counts[i,2])*100)

}

write.table(rat.cell.counts,
            file = 'rat_cell_count_by_lib.csv',
            sep = ',',
            quote = FALSE)

save(raw.hybrid,
     file = 'raw_data_mapped_to_hybrid.RData')

save(trans.obj,
     file = 'transformed_data_hybrid_genome.RData')

save(scale.obj,
     file = 'scaled_cluster_list_hybrid.RData')

save(clean.objlist,
     file = 'removed_rat_mapped_hybrid.RData')


#---- Process Notes -----
# This R code performs a series of data preprocessing,
# quality control, and analysis steps on single-cell RNA sequencing
# (scRNA-seq) data using the Seurat package. The code is organized into
# different sections, each focused on a specific task.
# Detailed overview of the code:
#
# 1. **Loading Libraries and Setting Up Future Processing:**
#   - Several R libraries are loaded, including Seurat, sctransform,
#     ggplot2, and others.
#   - The processing plan is set to "multisession" mode with one worker,
#     which enables parallel processing.
#
# 2. **Setting Options:**
#   - The maximum number of expressions to evaluate in a single
#     future is adjusted.
#   - The maximum size of future globals is set.
#
# 3. **Setting Working Directory and Locating Data Files:**
#   - The working directory is set to a specific path.
#   - A list of paths to files with a specific naming pattern
#     is created using `list.files`.
#
# 4. **Merging Data Files and Creating Seurat Objects:**
#   - A list called `objlist` is created to store Seurat objects.
#   - For each data file, the gene expression data is read using
#     `Read10X_h5` and converted into a Seurat object.
#   - Time and group identifiers are extracted from the
#     file names and assigned to the Seurat object.
#   - Each Seurat object is stored in the `objlist`.
#
# 5. **Saving Raw Data:**
#   - The list of Seurat objects is saved as an RData file.
#
# 6. **Quality Control and Preprocessing:**
#   - A loop processes each Seurat object in the `objlist`.
#   - The percentage of mitochondrial (MT) reads is calculated
#     and added as a feature.
#   - Filtering is performed based on the number of expressed
#     features and the percentage of MT reads.
#   - Data is normalized and highly variable features are
#     identified using `SCTransform`.
#
# 7. **Clustering and Dimensionality Reduction:**
#   - Another loop processes each Seurat object.
#   - Principal Component Analysis (PCA) is performed.
#   - Nearest neighbors and clusters are identified using
#     `FindNeighbors` and `FindClusters`.
#   - UMAP dimensionality reduction is applied.
#
# 8. **Removing Rat Cells:**
#   - Rat cell clusters are identified and labeled for removal
#     based on predetermined indices.
#   - Rat cells are labeled as "rat" in the Seurat objects.
#   - Human cell barcodes are written to a file.
#   - Rat cell clusters are removed from the Seurat objects,
#     resulting in a new list called `clean.objlist`.
#
# 9. **Calculating and Saving Rat Cell Counts:**
#   - The number of rat cells, total cell counts, and the percentage of rat
#     astrocyte cells are calculated for each dataset.
#   - The results are written to a CSV file.
#
# 10. **Saving Cleaned Data:**
#   - The list of cleaned Seurat objects is saved as an RData file.
#
# 11. **Visualizations:**
#   - Plots are generated to visualize the data, clusters, and features
#     for each dataset.
#
# This code demonstrates a comprehensive analysis pipeline for scRNA-seq data,
# including data loading, preprocessing, quality control, clustering,
# dimensionality reduction, cell type annotation, and visualization.
# The code focuses on removing rat cells and summarizing their proportions in
# each dataset.





# ---- Code with comments ----
# plan("multisession", workers = 8)
# # plan("sequential")
# # plan()
#
# #Set max number of expressions to evaluate in a single future.
# #defines the size of each processing 'chunk'
# # high values improve performance but consume memory
# options(expressions = 20000)
#
# options(future.globals.maxSize = 207374182400)
#
#
# #set working directory
# # setwd("/data/FASTQ/Duan_Project_025_hybrid")
# setwd("/data/Alexi_Duhe/test1")
# # locate data
# #return list of the paths to the files
# h5.list <- list.files('/data/FASTQ/Duan_Project_025_hybrid',
#                       'filtered_feature_bc_matrix.h5',
#                       recursive = TRUE,
#                       include.dirs = FALSE,
#                       full.names = TRUE)
#
# # Merge data files using for-loop
#
# # create vector to store list of objects
# raw.hybrid <- vector(mode = 'list', length = length(h5.list))
#
# # extract name identifier from path list created in h5.list at ith h5 file
# name.id <- str_extract(h5.list,
#                        pattern = '[0-9]+-[0-6]') #group+time
# # extract the time identifier from path list created in h5.list at ith file
# time.id <- sub('.*-', '', name.id)
# # extract the group identifier from path list created in h5.list at ith file
# group.id <- sub('-.*', '', name.id)
# # create names vector to store names variables in vector
#
# for (i in 1:length(h5.list)) {
#
#   # create variable w the read the ith h5 file
#   h5file <- Read10X_h5(filename = h5.list[i],
#                        use.names = T,
#                        unique.features = T)
#
#   # create seurat object from gene expression. project name = name variable at i-th h5 file
#   raw.hybrid[[i]] <- CreateSeuratObject(counts = h5file$`Gene Expression`,
#                                         project = name.id)
#
#   # assign time.ident to seurat object based on time variable for  ith h5 file path
#   raw.hybrid[[i]]$time.ident <- time.id[[i]]
#
#   # assign group.ident to seurat object based on group variable for ith h5 file path
#   raw.hybrid[[i]]$group.ident <- group.id[[i]]
#
#   # # store seurat object in previously created object vector at current index i
#   # objlist[[i]] <- gobj
#
#   print(paste0(i, " number of genes: ", nrow(raw.hybrid[[i]]),
#                ", number of cells: ", ncol(raw.hybrid[[i]])))
# }
#
# # 1. %MT read
# # transformed object
# trans.obj <-raw.hybrid
#
# for (i in 1:length(trans.obj)) {
#
#   print(paste0('Now at ', i))
#
#   # 1. % MT reads
#   # gobj = miscelaneous object used as intermediary to save memory
#   trans.obj[[i]] <- PercentageFeatureSet(trans.obj[[i]],
#                                          pattern = "^MT-",
#                                          col.name = 'percent.mt')
#
#   trans.obj[[i]] <- FindVariableFeatures(trans.obj[[i]], nfeatures = 3000)
#
#   # 2. Filtering
#   # trans.obj[[i]] <- subset(objlist[[i]],
#   #                        subset = nFeature_RNA > 300 &
#   #                        nFeature_RNA < 7500 &
#   #                        percent.mt < 20)
#
#   # 3. Normalize & Identify highly variable features
#   # SCTransform() performs NormalizeData(), ScaleData(), FindVariableFeatures()
#   trans.obj[[i]] <- SCTransform(trans.obj[[i]],
#                                 vars.to.regress = 'percent.mt',
#                                 method = 'glmGamPoi',
#                                 variable.features.n = 8000,
#                                 return.only.var.genes = FALSE,
#                                 seed.use = 2022,
#                                 verbose = TRUE)
# }
#
# # 5. Clustering
#
# ffeatures <- c('Gfap', 'S100b', 'GAD1', 'SLC17A6')
#
# scale.obj <-trans.obj
#
# for (i in 1:length(scale.obj)) {
#
#   # 4. Principal Component Analysis: capturing heterogeneity
#   scale.obj[[i]] <- RunPCA(scale.obj[[i]],
#                            verbose = TRUE,
#                            seed.use = 2022)
#
#   scale.obj[[i]] <- RunUMAP(object = scale.obj[[i]],
#                             reduction = 'pca',
#                             dims = 1:30,
#                             seed.use = 2022)
#
#   scale.obj[[i]] <- FindNeighbors(scale.obj[[i]],
#                                   reduction = 'pca',
#                                   dims = 1:30) #dims = # of PC
#
#   scale.obj[[i]] <- FindClusters(scale.obj[[i]],
#                                  resolution = c(0.5),
#                                  random.seed = 2022)
#
#   # Dimension Plot
#   plot.dim <- DimPlot(scale.obj[[i]],
#                       group.by = 'seurat_clusters',
#                       label = TRUE,
#                       reduction = 'umap') +
#     ggtitle(name.id[i], ' plot by cluster') +
#     theme(test = element_text(size = 12))
#
#   # Feature Plot
#   colorpal <- viridis(n = 10,
#                       option = 'C',
#                       direction = -1) #set color palate
#
#   plot.feature <- FeaturePlot(scale.obj[[i]],
#                               features = ffeatures)
#
#
#   # Violin plot
#   plot.violin <- VlnPlot_scCustom(scale.obj[[i]],
#                                   features = 'Gfap') +
#     scale_y_continuous(trans = 'log10') +
#     ggtitle(paste0(name.id[[i]], ' Gfap expression by cluster'))
#
#   # png(paste0('./test_QC1_plots/',
#   #           name.id), '_DimFeatVln_plots_by_cluster.png')
#   print(plot.feature + plot.dim + plot.violin)
#   # dev.off()
# }
#
# # # setwd('../../..')
# # setwd('/data/Alexi_Duhe/test1')
#
# rrc.at.group <- list(c(10,8,13),  #objlist[[1]]  9-0
#                      c(6,12),     #objlist[[2]]  9-1
#                      c(6),        #objlist[[3]]  9-6
#
#                      c(9),   #objlist[[4]]  13-0
#                      c(10),  #objlist[[5]]  13-1
#                      c(11),  #objlist[[6]]  13-6
#
#                      c(7),   #objlist[[7]]  21-0
#                      c(9),   #objlist[[8]]  21-1
#                      c(9),   #objlist[[9]]  21-6
#
#                      c(10),  #objlist[[10]] 23-0
#                      c(10),  #objlist[[11]] 23-1
#                      c(9),   #objlist[[12]] 23-6
#
#                      c(7,14,2),  #objlist[[13]] 53-0
#                      c(5,8,14),  #objlist[[14]] 53-1
#                      c(9,7))     #objlist[[15]] 53-6
#
# # Check and remove rat cells manually by looking at plots
#
# rat.obj <-scale.obj
# clean.objlist <- vector(mode = 'list', length = length(rat.obj))
#
# for (i in 1:length(rat.obj)) {
#
#   rrc.at.group <- remove.rat.clusters[[i]]
#
#   cat('rat cluster number: ', rrc.at.group[[i]])
#
#   rat.obj[[i]]$rat.ident <- 'human'
#   rat.obj[[i]]$rat.ident[rat.obj[[i]]$seurat_clusters %in% rrc.at.group[[i]]] <- 'rat'
#
#   human.barcodes <- colnames(objlist[[i]])[objlist[[i]]$rat.ident == 'human']
#
#   write.table(human.barcodes,
#               file = paste0('.',
#                             name.id[i],
#                             '_human_barcodes.txt'),
#               sep = '\t',
#               quote = FALSE,
#               row.names = FALSE,
#               col.names = FALSE)
#
#   clean.objlist[[i]] <- subset(rat.obj[[i]], subset = rat.ident != 'rat')
# }
#
# rat.cell.counts <- matrix(nrow = length(rat.obj),
#                           ncol = 3,
#                           dimnames = list(name.id,
#                                           c('number of rat astrocytes',
#                                             'total number of cells',
#                                             'percentage of astrocyte cells (%)')))
#
# for (i in 1:length(rat.obj)) {
#
#   # calculate the number fo rat cells
#   rat.cell.counts[i , 1] <- (sum(rat.obj[[i]]$rat.ident == 'rat'))
#
#   # calculate the number of cells
#   rat.cell.counts[i , 2] <- (length(rat.obj[[i]]$rat.ident))
#
#   # calculate the % of cells
#   rat.cell.counts[i , 3] <- ((rat.cell.counts[i,1]/rat.cell.counts[i,2])*100)
#
# }
#
# write.table(rat.cell.counts,
#             file = 'rat_cell_count_by_lib.csv',
#             sep = ',',
#             quote = FALSE)
#
# save(raw.hybrid,
#      file = 'raw_data_mapped_to_hybrid.RData')
#
# save(trans.obj,
#      file = 'transformed_data_hybrid_genome.RData')
#
# save(scale.obj,
#      file = 'scaled_cluster_list_hybrid.RData')
#
# save(clean.objlist,
#      file = 'removed_rat_mapped_hybrid.RData')
