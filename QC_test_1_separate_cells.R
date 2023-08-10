# Alexandra Duhe 8/3/2023 9:00 AM
#~/NVME/scARC_Duan_018/Duan_project_025_RNA/Analysis_part1_hybrid_genome/20230503_separate_human_from_rat_in_025_17_and_46_new_libraries-copy.R


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
# plan("sequential")
# plan()

#Set max number of expressions to evaluate in a single future.
#defines the size of each processing 'chunk'
# high values improve performance but consume memory
options(expressions = 20000)

options(future.globals.maxSize = 207374182400)


#set working directory
# setwd("/data/FASTQ/Duan_Project_025_hybrid")
setwd("/data/Alexi_Duhe/test1")
# locate data ----
#return list of the paths to the files
h5.list <- list.files('/data/FASTQ/Duan_Project_025_hybrid',
                      'filtered_feature_bc_matrix.h5',
                      recursive = TRUE,
                      include.dirs = FALSE,
                      full.names = TRUE)

# Merge data files using for-loop ----

# create vector to store list of objects
objlist <- vector(mode = 'list', length = length(h5.list))
h5names <- c()

for (i in 1:length(h5.list)) {

  # create name variable; extract name identifier from path list created in h5.list at ith h5 file
  name.id <- str_extract(h5.list[i],
                         pattern = '[0-9]+-[0-6]') #group+time
  # create time variable; extract the time identifier from path list created in h5.list at ith file
  time.id <- sub('.*-', '', name.id)
  # create group variable: extract the group identifier from path list created in h5.list at ith file
  group.id <- sub('-.*', '', name.id)
  # create names vector to store names variables in vector
  h5names <- c(h5names, as.vector(name.id))

  # create variable w the read the ith h5 file
  h5file <- Read10X_h5(filename = h5.list[i],
                       use.names = T,
                       unique.features = T)

  # create seurat object from gene expression. project name = name variable at i-th h5 file
  gobj <- CreateSeuratObject(counts = h5file$`Gene Expression`,
                             project = name.id)

  # assign time.ident to seurat object based on time variable for  ith h5 file path
  gobj$time.ident <- time.id

  # assign group.ident to seurat object based on group variable for ith h5 file path
  gobj$group.ident <- group.id

  # store seurat object in previously created object vector at current index i
  objlist[[i]] <- gobj
}

# save raw data
save(objlist,
     file = 'raw_data_mapped_to_hybrid.RData')

# 1. %MT read ----
# transformed object
trans.obj <- vector(mode = 'list',
                    length = length(objlist))

for (i in 1:length(objlist)) {

  # 1. % MT reads
  # gobj = miscelaneous object used as intermediary to save memory
  objlist[[i]] <- PercentageFeatureSet(objlist[[i]],
                                       pattern = "^MT-",
                                       col.name = 'percent.mt')

  # 2. Filtering
  objlist[[i]] <- subset(objlist[[i]],
                         subset = nFeature_RNA > 300 &
                         nFeature_RNA < 7500 &
                         percent.mt < 20)

  # 3. Normalize & Identify highly variable features
  # SCTransform() performs NormalizeData(), ScaleData(), FindVariableFeatures()
  objlist[[i]] <- SCTransform(objlist[[i]],
                           vars.to.regress = 'percent.mt',
                           variable.features.n = 8000,
                           return.only.var.genes = FALSE,
                           verbose = TRUE)
  trans.obj[[i]] <- objlist[[i]]
}

# 5. Clustering ----
scale.obj <- vector(mode = 'list',
                    length = length(objlist))

for (i in 1:length(objlist)) {

  # 4. Principal Component Analysis: capturing heterogeneity
  objlist[[i]] <- RunPCA(objlist[[i]],
                         verbose = TRUE)

  objlist[[i]] <- FindNeighbors(objlist[[i]],
                                dims = 1:15) #dims = # of PC

  objlist[[i]] <- FindClusters(objlist[[i]],
                              resolution = c(0.6))

  objlist[[i]] <- RunUMAP(object = objlist[[i]],
                       dims = 1:15)

  scale.obj[[i]] <- objlist[[i]]
}

# Check and remove rat cells manually by looking at plots----
clean.objlist <- vector(mode = 'list', length = length(objlist))

ffeatures <- c('Gfap', 'S100b', 'GAD1', 'SLC17A6')

# plots
for (i in 1:length(objlist)) {

  name.id <- str_extract(h5.list[i],
                         pattern = '[0-9]+-[0-6]') #group+time

    # Dimension Plot
    plot.dim <- DimPlot(objlist[[i]],
                           group.by = 'SCT_snn_res.0.6',
                           label = TRUE,
                           reduction = 'umap') +
      ggtitle(h5names[i], ' plot by cluster')

    # Feature Plot
    colorpal <- viridis(n = 10,
                        option = 'C',
                        direction = -1) #set color palate

    plot.feature <- FeaturePlot(objlist[[i]],
                                features = ffeatures)


    # Violin plot
    plot.violin <- VlnPlot_scCustom(objlist[[i]],
                                    features = 'Gfap') +
      scale_y_continuous(trans = 'log10')

    # png(paste0('./test_QC1_plots/',
    #            name.id), '_DimFeatVln_plots_by_cluster.png')
    print(plot.feature + plot.dim + plot.violin)
    # dev.off()
}

# setwd('../../..')
setwd('/data/Alexi_Duhe/test1')

remove.rat.clusters <- list(c(10,8,13),  #objlist[[1]]  9-0
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

for (i in 1:length(objlist)) {

  rrc.at.group <- remove.rat.clusters[[i]]

  cat('rat cluster number: ', rrc.at.group)

	objlist[[i]]$rat.ident <- 'human'
	objlist[[i]]$rat.ident[objlist[[i]]$SCT_snn_res.0.6 %in% rrc.at.group] <- 'rat'

	human.barcodes <- colnames(objlist[[i]])[objlist[[i]]$rat.ident == 'human']

	write.table(human.barcodes,
				file = paste0('.',
				              h5names[i],
				              'ALEXI_TEST1_human_barcodes.txt'),
				sep = '\t',
				quote = FALSE,
				row.names = FALSE,
				col.names = FALSE)

	clean.objlist[[i]] <- subset(objlist[[i]], subset = rat.ident != 'rat')
}

rat.cell.counts <- matrix(nrow = length(objlist),
                          ncol = 3,
                          dimnames = list(h5names,
                                          c('number of rat astrocytes',
                                            'total number of cells',
                                            'percentage of astrocyte cells (%)')))

for (i in 1:length(objlist)) {

  # calculate the number fo rat cells
  rat.cell.counts[i , 1] <- (sum(objlist[[i]]$rat.ident == 'rat'))

  # calculate the number of cells
  rat.cell.counts[i , 2] <- (length(objlist[[i]]$rat.ident))

  # calculate the % of cells
  rat.cell.counts[i , 3] <- ((rat.cell.counts[i,1]/rat.cell.counts[i,2])*100)

}

write.table(rat.cell.counts,
            file = 'ALEXI_TEST1_rat_count_by_lib.csv',
            sep = ',',
            quote = FALSE)


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

