# 8/4/2023 4:30 PM

# /nvmefs/scARC_Duan_018/Duan_project_025_RNA/Analysis_part2_GRCh38/20230507_normalize_integrate_cluster_025_17_46.R


library(Seurat)
library(ggplot2)
library(pals)
library(stringr)
library(future)
library(glmGamPoi)
library(cowplot)

setwd("/data/Alexi_Duhe/test1/")


plan("multisession", workers = 8)
options(expressions = 20000)
options(future.globals.maxSize = 207374182400)

# load previous data file created in QC_test_2_demultiplex.R
load('GRCh38_mapped_after_QC_list.RData')


h5.list <- list.files('/data/FASTQ/Duan_Project_025_hybrid',
                      'filtered_feature_bc_matrix.h5',
                      recursive = TRUE,
                      include.dirs = FALSE,
                      full.names = TRUE)

# Initialize variables and data structures
name.id <- str_extract(h5.list, '[0-9]+-[0-6]')#group+time
# extract the time identifier from path list created in h5.list at ith file
time.id <- sub('.*-', '', name.id)
# extract the group identifier from path list created in h5.list at ith file
group.id <- sub('-.*', '', name.id)




# Normalize and transform the data ----

# create empty list to store data following transformation
trans.list <- vector(mode = 'list', length(QCed.list))

# Iterate through each samle in the preprocessed data
for (i in 1:length(QCed.list)) {

  # Print progress message indicating the current sample being processed
  print(paste0('Working on ', i, ' of ', length(QCed.list)))

  # Calculate the percentage of mt gene expression for each cell
  QCed.list[[i]] <- PercentageFeatureSet(QCed.list[[i]],
                                         pattern = c('^MT-'),
                                         col.name = 'percent.mt')

  # Identify variable features using a gamma-Poisson distribution
  QCed.list[[i]] <- FindVariableFeatures(QCed.list[[i]], nfeatures = 8000,)

  # Apply SCTransform to normalize, scale, find variable features
  QCed.list[[i]] <- SCTransform(QCed.list[[i]],
                                vars.to.regress = 'percent.mt',
                                method = 'glmGamPoi',
                                return.only.var.genes = F,
                                variable.features.n = 8000,
                                verbose = T,
                                seed.use = 42)

  trans.list[[i]] <- QCed.list[[i]]
}


# Select a common set of features across all samples for integration
ffeatures <- SelectIntegrationFeatures(object.list = trans.list,
                                       nfeatures = 3000,
                                       fvf.nfeatures = 3000)

# Prepare data for integration by aligning the selected features
trans.list <- PrepSCTIntegration(trans.list,
                                 anchor.features = ffeatures)

# Apply the scale.pca function to each sample in the transformed list
trans.list <- lapply(X = trans.list,
                     FUN = function(x) {
                       x <- ScaleData(x,
                                      features = ffeatures)
                       x <- RunPCA(x,
                                   features = ffeatures)
                     })

# Find integration anchors across samples for joint analysis
anchors <- FindIntegrationAnchors(object.list = trans.list,
                                  anchor.features = ffeatures,
                                  reference = c(1,2,3),
                                  reduction = 'rpca',
                                  normalization.method = 'SCT',
                                  scale = FALSE,
                                  dims = 1:50)

# Integrate data from different samples using the integration anchors
integrated.data <- IntegrateData(anchorset = anchors,
                                 verbose = TRUE)
save(integrated.data, file = 'integrated_object.RData')

# Calculate the percentage of mt gene expression in each cell
integrated.data[['percent.mt']] <- PercentageFeatureSet(integrated.data,
                                                        assay = 'RNA',
                                                        pattern = '^MT-')
# Set the default assay for the integrated data to 'RNA'
DefaultAssay(integrated.data) <- 'RNA'

# Set cell identities to the orig identifiers
Idents(integrated.data) <- 'orig.ident'

# Creat violin plot to visualize the key QC metrics
VlnPlot(integrated.data,
        features = c('nFeature_RNA',
                     'nCount_RNA',
                     'percent.mt'),
        ncol = 3, fill.by = 'feature', pt.size = 0)
# Set the default assay to 'integrated' for downstream analysis

unique(integrated.data$orig.ident)
unique(integrated.data$time.ident)

DefaultAssay(integrated.data) <- 'integrated'


# Scale data for downstream analysis
integrated.data <-  ScaleData(integrated.data)

# Run PCA on scaled data
integrated.data <- RunPCA(integrated.data, verbose = TRUE,
                          seed.use = 11)

# Run UMAP on PCA results
integrated.data <- RunUMAP(integrated.data, reduction = 'pca',
                           dims = 1:30,
                           seed.use = 11)

# Find nearest neighbors based on UMAP coordinates
integrated.data <- FindNeighbors(integrated.data, reduction = 'pca',
                                 dims = 1:30)

# Perform clustering using Louvain algorithm
integrated.data <- FindClusters(integrated.data, resolution = 0.6,
                                random.seed = 11)

library.colors <- DiscretePalette(n = length(unique(Idents(integrated.data))),
                                  'alphabet')
# Visualize clustering results using a dimentionality reduction plot
p1 <- DimPlot(integrated.data,
              cols = library.colors,
              label = T,
              group.by = 'seurat_clusters',
              raster = F) +
  ggtitle('Clustering RNAseq data') +
  theme(text = element_text(size = 12),
        axis.text = element_text(size = 12))
print(p1)

# remove cluster 17
integrated.data <- subset(integrated.data, seurat_clusters != '17')
Idents(integrated.data) <- 'seurat_clusters'
clust <- c('1', '2', '3', '4', '5', '6',
                    '7', '8', '9', '10', '11', '12', '13',
                    '14', '15', '16', '17', '18', '19', '20')
names(clust) <- levels(integrated.data)
integrated.data <- RenameIdents(integrated.data, clust)

DimPlot(integrated.data,
        reduction = 'umap',
        label = TRUE,
        repel = FALSE,
        pt.size = 0.3,
        raster = F,
        cols = library.colors) +
  ggtitle('Clustiering RNAseq Data') +
  theme(text = element_text(size = 12),
        axis.text = element_text(size = 12))

DefaultAssay(integrated.data) <- 'integrated'
# Create DimPlot based on orig identifiers
DimPlot(integrated.data,
        label = FALSE,
        cols = library.colors,
        group.by = 'orig.ident',
        raster = F) +
  ggtitle('By Library') +
  theme(text = element_text(size = 12),
        axis.text = element_text(size = 12))

DimPlot(integrated.data,
        label = FALSE,
        group.by = 'time.ident',
        raster = F) +
  ggtitle('By Time Point') +
  theme(text = element_text(size = 12),
        axis.text = element_text(size = 12))

lib.nine <- subset(integrated.data, orig.ident == '09-0' |
                     orig.ident == '09-1' |
                     orig.ident == '09-6')

DimPlot(lib.nine,
        label = FALSE,
        cols = library.colors,
        group.by = 'orig.ident') +
  ggtitle('Library 09') +
  theme(text = element_text(size = 12),
        axis.text = element_text(size = 12))

no.lib.nine <- subset(integrated.data, orig.ident != '09-0' &
                        orig.ident != '09-1' &
                        orig.ident != '09-6')

DimPlot(no.lib.nine,
        label = FALSE,
        cols = library.colors,
        group.by = 'orig.ident',
        raster = F) +
  ggtitle('Without Library 09') +
  theme(text = element_text(size = 12),
        axis.text = element_text(size = 12))




# for (i in 1:length(name.id) - 3) {
#
#
#
#   lib <- subset(integrated.data, orig.ident == (name.id[i]) |
#                        orig.ident == (name.id[i + 1]) |
#                        orig.ident == (name.id[i + 2]))
#
#   DimPlot(lib,
#           label = FALSE,
#           cols = library.colors,
#           group.by = 'orig.ident') +
#     ggtitle('Library ', paste0(group.id[i])) +
#     theme(text = element_text(size = 12),
#           axis.text = element_text(size = 12))
#
#   no.lib <- subset(integrated.data, orig.ident != (name.id[i]) &
#                           orig.ident != (name.id[i + 1]) &
#                           orig.ident != (name.id[i + 2]))
#
#   DimPlot(no.lib,
#           label = FALSE,
#           cols = library.colors,
#           group.by = 'orig.ident',
#           raster = F) +
#     ggtitle('Without Library ', paste0(group.id[i])) +
#     theme(text = element_text(size = 12),
#           axis.text = element_text(size = 12))
#
#
# }



# Cell type identification ----

# Set the default assay to SCT for cell type identification
DefaultAssay(integrated.data) <- 'SCT'

# Feature plot for specific marker genes
FeaturePlot(integrated.data,
            features = c('SOX2',
                         "VIM"),
            raster = F)
FeaturePlot(integrated.data,
            features = c('GAD1',
                         'GAD2',
                         "SLCC17A6"),
            ncol = 2,
            raster = F)
FeaturePlot(integrated.data,
            features = c('POU5F1',
                         'NANOG'),
            raster = F)
FeaturePlot(integrated.data,
            features = c('NEFM',
                         'MAP2'),
            raster = F)

# Define new cluster itentifiers
# List of gene markers
trimmed_markers <- c(
  "GAD1", "GAD2", "SLC17A6", "SLC17A7", "EBF1", # Striatal and various
  "SEMA3E", # Subcerebral marker
  "BCL11B",  # Cortical marker
  "SST", # Inhibitory marker
  "SATB2",  "NEFM", # Excitatory markers
  "VIM", "SOX2",  # NPC (Neural Progenitor Cell) markers
  "SERTAD4", "rna_FOXG1",  # Forebrain markers
  "POU3F2", "LHX2", # General cortex markers
  "ADCYAP1", "CUX1", "CUX2", "MAP2", "DCX" # Various neuronal markers
)

# Convert the gene markers vector to a factor with specified levels
#trimmed_markers <- factor(trimmed_markers, levels = trimmed_markers)


# Create a stacked violin plot using trimmed maker genes
VlnPlot(integrated.data,
        features = trimmed_markers,
        stack = TRUE,
        sort = F,
        flip = TRUE) +
  theme(legend.position = 'none') +
  ggtitle('Expression by Gene Identity')


# save(integrated.data, file = '08_25_2023_work.RData')
# save(anchors, file = '08_25_2023_anchors.RData')
#
# load('08_25_2023_work.RData')
# load('08_25_2023_anchors.RData')



# Set identities of the integrated data to 'seurat_clusters'
Idents(integrated.data) <- 'seurat_clusters'

# Define new cluster identifiers
# new.cluster.ids
#     NEFM+ glut
#     GABA
#     SST+ GABA
#     NEFM- glut
#     SEMA3E+ GABA
#     NEFM- glut
#     immature
#     unknown


new.cluster.ids <-
  c("NEFM+ glut", #1
    "NEFM- glut", #2
    "GABA", #3
    "NEFM- glut", #4
    "SEMA3E+ GABA", #5
    "GABA", #6
    "NEFM+ glut", #7
    "NEFM+ glut", #8
    "GABA",#9
    "GABA",#10
    "SST+ GABA", #11
    "immature", #12
    "NEFM- glut", #13
    "unknown", #14
    "NEFM- glut", #15
    "NEFM- glut", #16
    "GABA", #17
    "NEFM- glut", #18
    'GABA',#19
    'immature') #20


# length(new.cluster.ids)
#length(unique(integrated.data$seurat_clusters))

# Rename cluster levels with new cluster identifiers
names(new.cluster.ids) <- levels(integrated.data)
integrated.labeled <- RenameIdents(integrated.data,
                                   new.cluster.ids)

# Create DimPlot to visualize labeled clusters using UMAP
DimPlot(integrated.labeled,
        reduction = 'umap',
        label = TRUE,
        repel = FALSE,
        pt.size = 0.3,
        raster = F,
        cols = c("#E6D2B8", #nmglut
                 "#CCAA7A", #npglut
                 "#B33E52", #GABA
                 "#E6B8BF", #SST GABA
                 "#CC7A88", #SEMA3E GABA
                 "#0075DC", #immature neuron
                 "#4C005C",
                 "#993F00")) + #unknown
  ggtitle('Labeled by Cell Type') +
  theme(text = element_text(size = 12),
        axis.text = element_text(size = 12))

cell.ctn <- table(integrated.data@meta.data$seurat_clusters)

DimPlot(integrated.labeled,
        reduction = 'umap',
        label = TRUE,
        repel = FALSE,
        pt.size = 0.3,
        raster = F,
        cols = c("#E6D2B8", #nmglut
                 "#CCAA7A", #npglut
                 "#B33E52", #GABA
                 "#E6B8BF", #SST GABA
                 "#CC7A88", #SEMA3E GABA
                 "#0075DC", #immature neuron
                 "#4C005C",
                 "#993F00")) + #unknown
  ggtitle('Labeled by Cell Type') +
  theme(text = element_text(size = 12),
        axis.text = element_text(size = 12))

# ### Siwei test code ######
# DimPlot(integrated.labeled,
#         reduction = 'umap',
#         label = TRUE,
#         repel = FALSE,
#         pt.size = 0.3,
#         cols = c("#E6D2B8", #nmglut
#                  "#CCAA7A", #npglut
#                  "#B33E52", #GABA
#                  "#E6B8BF", #SST GABA
#                  "#CC7A88", #SEMA3E GABA
#                  "#0075DC", #immature neuron
#                  "#4C005C", #unknown neuron
#                  "#993F00", #unknown
#                  "#996F00")) + #unknown
#   ggtitle('labeled by cell type') +
#   theme(text = element_text(size = 10),
#         axis.text = element_text(size = 10))

###

# Cell type labeling and visualization ----

# Set 'fine.cell.type' to the active identifier in 'integrated.labeled'
integrated.labeled$fine.cells.type <- integrated.labeled@active.ident

# Convert 'fine.cell.type' to character and assign to 'cell.type'
integrated.labeled$cell.type <- as.character(integrated.labeled$fine.cells.type)

# Assign specific cell types to the 'cell.type' column
integrated.labeled$cell.type[integrated.labeled$fine.cells.type %in%
                               c("SEMA3E+ GABA",
                                 "SST+ GABA",
                                 "GABA")] <- "GABA"

# Assign 'unidentified' to 'cell.type' for certain categories
integrated.labeled$cell.type[integrated.labeled$cell.type %in%
                               c("unknown neuron",
                                 "immature neuron",
                                 "unknown")] <- "unidentified"

# Copy 'cell.type' to 'cell.type.forplot'
integrated.labeled$cell.type.forplot <- integrated.labeled$cell.type

# Assign specific cell types to 'cell.type' based on 'fine.cell.type'
integrated.labeled$cell.type[integrated.labeled$fine.cells.type == "NEFM- glut"] <- "nmglut"
integrated.labeled$cell.type[integrated.labeled$fine.cells.type == "NEFM+ glut"] <- "npglut"

# unique(integrated.labeled$cell.type)
# unique(integrated.labeled$cell.type.forplot)

# Create new column 'cell.type.counts' and populate with 'cell.type.forplot'
integrated.labeled$cell.type.counts <- integrated.labeled$cell.type.forplot

# Get unique cells types
types <- unique(integrated.labeled$cell.type.forplot)

# Loop through each cell type to count occurrences and update 'cell.type.counts'
for (i in 1:length(types)) {
  count <- sum(integrated.labeled$cell.type.forplot == types[i])
  print(count)
  integrated.labeled$cell.type.counts[integrated.labeled$cell.type.counts == types[i]] <-
    paste0(types[i], "\n", count)
}

# unique(integrated.labeled$cell.type.counts)

# Create a DimPlot to visualize clusters labeled by cell type counts
DimPlot(integrated.labeled, reduction = "umap", group.by = "cell.type.counts",
        raster = F,
        label = TRUE, repel = F, pt.size = 0.3,
        cols = c("#B33E52",
                 "#E6D2B8",
                 "#CCAA7A",
                 "#54990F",
                 "#0075DC")) +
  ggtitle("labeled by cell type") +
  theme(text = element_text(size = 10),
        axis.text = element_text(size = 10))

save(trans.list,
     file = 'transformed_prepared_for_integration.RData')
save(integrated.labeled,
     file = "integrated_labeled.RData")
save(integrated.data,
     file = 'integrated_clustered_obj.RData')

# integrated.data <-
load("integrated_clustered_obj.RData")



#---- Process Notes -----
#This code is part of the analysis process for single-cell RNA sequencing
# (scRNA-seq) data. It involves normalization, integration, and clustering of
# scRNA-seq datasets using the Seurat package.
#
# 1. **Libraries**: The required R libraries are loaded, including Seurat for
#      single-cell analysis, ggplot2 for visualization, pals for color palettes,
#      stringr for string manipulation, future for parallel processing,
#      and glmGamPoi for a particular normalization method.
#
# 2. **Parallel Processing Setup**: The code sets up parallel processing
#      using two worker threads. It also configures options related to expression
#      evaluation and memory usage.
#
# 3. **Load Preprocessed Data**: The preprocessed scRNA-seq data, which has
#      gone through quality control and filtering, is loaded from a previous
#      file (`GRCh38_mapped_after_QC_list.RData`).
#
# 4. **Normalization**: The code proceeds with the normalization of the data.
#      It iterates through the list of preprocessed datasets (`QCed.list`).
#      For each dataset, the following steps are performed:
#   - Calculate the percentage of mitochondrial reads for each cell.
#   - Find variable features using the `FindVariableFeatures` function,
#     including regression on the percentage of mitochondrial reads.
#   - Perform SCTransform normalization using the `SCTransform` function,
#     again with regression on the percentage of mitochondrial reads.
#
# 5. **PCA and Integration Anchors**: After normalization, the code performs
#      principal component analysis (PCA) on the integrated datasets. Integration
#      anchors are found using the `FindIntegrationAnchors` function,
#      which identifies shared features and data patterns among the datasets.
#
# 6. **Integration**: The data is integrated using the integration anchors
#      obtained in the previous step with the `IntegrateData` function.
#
# 7. **Quality Control (QC)**: The integrated data's mitochondrial content is
#      calculated using the `PercentageFeatureSet` function. The default assay is
#      set to 'RNA', and the `Vlnplot` function is used to visualize the
#      distribution of features such as the number of RNA features, RNA counts,
#      and mitochondrial percentages.
#
# 8. **Clustering**: The integrated data is clustered using the standard
#      Seurat workflow. Principal component analysis (PCA), UMAP dimension reduction,
#      finding neighbors, and clustering are performed. The resulting clusters are
#      visualized using `DimPlot`.
#
# 9. **Cell Type Identification**: Different cell types are identified based on
#      the expression of specific marker genes. `FeaturePlot` is used to visualize
#      the expression of marker genes for various cell types. The code also defines
#      a list of trimmed marker genes and creates a stacked violin plot to
#      visualize their expression.
#
# 10. **Cell Type Labeling and Visualization**: Clusters are assigned
#       human-readable cell type labels. The `DimPlot` function is used to
#       visualize the labeled clusters in the UMAP space.
#
# 11. **Cell Type Counts**: The counts of each cell type are summarized and
#       appended to the cell type names. The `DimPlot` function is used again to
#       visualize the clusters based on cell type counts.
#
# 12. **Save Data**: The labeled integrated data is saved to an RData file
#       named "integrated_labeled.RData".
#
# This code performs several important steps in scRNA-seq data analysis,
# including normalization, integration, clustering, and cell type identification.
# The Seurat package and associated functions are extensively utilized
# to accomplish these tasks.
#
#


# Notes ----
#
# # /nvmefs/scARC_Duan_018/Duan_project_025_RNA/Analysis_part2_GRCh38/20230509_025_17_46_plot_QC_stats.R
# # Load the required libraries for data processing and analysis
# library(Seurat)     # For single-cell data analysis
# library(ggplot2)    # For creating visualizations
# library(pals)       # Palette generation for plots
# library(stringr)    # String manipulation functions
# library(future)     # Parallel and asynchronous processing
# library(glmGamPoi)  # Generalized linear modeling for gene expression
#
# # Configure parallel processing for efficiency using multiple workers
# plan("multisession", workers = 2)
#
# # Increase the number of expressions to be evaluated at once for memory optimization
# options(expressions = 20000)
#
# # Set the maximum size for global variables, important for memory control
# options(future.globals.maxSize = 207374182400)  # 200 GB
#
# # Load the previously quality-controlled and preprocessed data
# load('GRCh38_mapped_after_QC_list.RData')
#
# # Normalize and transform the data using SCTransform algorithm
#
# # Create an empty list to store the transformed data for each sample
# trans.list <- vector(mode = 'list', length = length(QCed.list))
#
# # Iterate through each sample in the preprocessed data
# for (i in 1:length(QCed.list)) {
#   # Print a progress message indicating the current sample being processed
#   print(paste0('Currently working on sample: ', i))
#
#   # Calculate the percentage of mitochondrial gene expression for each cell
#   obj <- PercentageFeatureSet(QCed.list[[i]],
#                               pattern = c('^MT-'),
#                               col.name = 'percent.mt')
#
#   # Identify variable features using a generalized linear model with a gamma-Poisson distribution
#   obj <- FindVariableFeatures(obj,
#                               nfeatures = 8000,
#                               vars.to.regress = 'percent.mt',
#                               method = 'glmGamPoi',
#                               return.only.var.genes = FALSE,
#                               variable.features.n = 8000,
#                               verbose = TRUE)
#
#   # Apply the SCTransform method to normalize and transform the gene expression data
#   trans.list[[i]] <- SCTransform(obj,
#                                  vars.to.regress = 'percent.mt',
#                                  method = 'glmGamPoi',
#                                  return.only.var.genes = FALSE,
#                                  variable.features.n = 8000,
#                                  verbose = TRUE)
# }
#
# # Principal Component Analysis (PCA) for data integration
#
# # Select a common set of features across all samples for integration
# ffeatures <- SelectIntegrationFeatures(object.list = trans.list,
#                                        nfeatures = 3000,
#                                        fvf.nfeatures = 3000)
#
# # Prepare the data for integration by aligning the selected features
# trans.list <- PrepSCTIntegration(trans.list,
#                                  anchor.features = ffeatures)
#
# # Define a function to scale and perform PCA on the input data
# scale.pca <- function(x) {
#   x <- ScaleData(x,
#                  features = ffeatures)
#   x <- RunPCA(x,
#               features = ffeatures)
# }
#
# # Apply the scale.pca function to each sample in the transformed list
# trans.list <- lapply(X = trans.list,
#                      FUN = scale.pca)
#
# # Find integration anchors across samples for joint analysis
# anchors <- FindIntegrationAnchors(object.list = trans.list,
#                                   anchor.features = ffeatures,
#                                   reference = c(1,2,3),
#                                   reduction = 'rpca',
#                                   normalization.method = 'SCT',
#                                   scale = FALSE,
#                                   dims = 1:50)
#
# # Integrate data from different samples using the integration anchors
# integrated.data <- IntegrateData(anchorset = anchors,
#                                  verbose = TRUE)
#
# # Perform quality control (QC) on the integrated data
#
# # Calculate the percentage of mitochondrial gene expression for each cell
# integrated.data[['percent.mt']] <- PercentageFeatureSet(integrated.data,
#                                                         assay = 'RNA',
#                                                         pattern = '^MT-')
#
# # Set the default assay for the integrated data to 'RNA'
# DefaultAssay(integrated.data) <- 'RNA'
#
# # Set cell identities to the original identifiers
# Idents(integrated.data) <- 'original.ident'
#
# # Create a violin plot to visualize key QC metrics
# VlnPlot(integrated.data,
#         features = c('nFeature_RNA',
#                      'nCount_RNA',
#                      'percent.mt'),
#         ncol = 3, fill.by = 'feature', pt.size = 0)
#
#
# # Clustering and Visualization
#
# # Identify unique original identifiers and time identifiers
# unique_original_idents <- unique(integrated.data$original.ident)
# unique_time_idents <- unique(integrated.data$time.ident)
#
# # Set the default assay to 'integrated' for downstream analysis
# DefaultAssay(integrated.data) <- 'integrated'
#
# # Standard workflow for visualization and clustering using Seurat
#
# # Scale the data for downstream analysis
# scaled_data <- integrated.data %>%
#   ScaleData()
#
# # Run Principal Component Analysis (PCA) on the scaled data
# pca_result <- scaled_data %>%
#   RunPCA(verbose = TRUE,
#          seed.use = 2023)
#
# # Run Uniform Manifold Approximation and Projection (UMAP) using PCA results
# umap_result <- pca_result %>%
#   RunUMAP(reduction = 'pca',
#           dims = 1:30,
#           seed.use = 2023)
#
# # Find nearest neighbors based on the UMAP coordinates
# neighbors_result <- umap_result %>%
#   FindNeighbors(reduction = 'pca',
#                 dims = 1:30)
#
# # Perform clustering using the Louvain algorithm
# cluster_result <- neighbors_result %>%
#   FindClusters(resolution = 0.6,
#                random.seed = 2023)
#
# # Visualize clustering results using a dimensionality reduction plot
# DimPlot(integrated.data,
#         label = TRUE,
#         group.by = 'seurat_clusters') +
#   NoLegend() +
#   ggtitle('Clustering RNAseq data') +
#   theme(text = element_text(size = 10),
#         axis.text = element_text(size = 10))
#
# # Define a color palette for library-based grouping
# library.colors <- DiscretePalette(n = length(unique(Idents(integrated.data))),
#                                   'alphabet')
#
# # Create a dimensionality reduction plot with colors based on original identifiers
# DimPlot(integrated.data,
#         label = FALSE,
#         cols = library.colors,
#         group_by = 'original.ident') +
#   ggtitle('Grouped by Library') +
#   theme(text = element_text(size = 10),
#         axis.text = element_text(size = 10))
#
# # Cell Type Identification
#
# # Set the default assay to 'SCT' for cell type identification
# DefaultAssay(integrated.data) <- 'SCT'
#
# # Cell Type Visualization and Labeling
#
# # Create a feature plot for specific marker genes
# # Features include SOX2 and VIM
# FeaturePlot(integrated.data, features = c('SOX2', 'VIM'))
#
# # Features include GAD1, GAD2, and SLCC17A6 with 2 columns
# FeaturePlot(integrated.data, features = c('GAD1', 'GAD2', 'SLCC17A6'), ncol = 2)
#
# # Features include POU5F1 and NANOG
# FeaturePlot(integrated.data, features = c('POU5F1', 'NANOG'))
#
# # Features include NEFM and MAP2
# FeaturePlot(integrated.data, features = c('NEFM', 'MAP2'))
#
# # Define a list of trimmed marker genes for StackedVlnPlot
# trimmed_markers <- c("GAD1", "GAD2", "SLC17A6", "SLC17A7",
#                      "EBF1", "SEMA3E", "BCL11B", "SST", "SATB2", "NEFM",
#                      "VIM", "SOX2", "SLC17A7", "SERTAD4", "FOXG1",
#                      "POU3F2", "LHX2", "ADCYAP1", "CUX1", "CUX2",
#                      "MAP2", "DCX")
#
# # Create a StackedVlnPlot using trimmed marker genes
# StackedVlnPlot(obj = integrated.data, features = trimmed_markers) +
#   coord_flip()
#
# # Set the identities of the integrated data to 'seurat_clusters'
# Idents(integrated.data) <- 'seurat_clusters'
#
# # Define new cluster identifiers
# new.cluster.ids <- c("NEFM- glut", "NEFM+ glut", "NEFM- glut", "NEFM+ glut", "GABA",
#                      "GABA", "NEFM- glut", "SST+ GABA", "NEFM- glut", "NEFM- glut",
#                      "SEMA3E+ GABA", "NEFM- glut", "NEFM- glut", "GABA", "GABA",
#                      "unknown neuron", "immature neuron", "unknown")
#
# # Print the length of new cluster identifiers
# print(length(new.cluster.ids))
#
# # Print the number of unique clusters in integrated data
# print(length(unique(integrated.data$seurat_clusters)))
#
# # Rename cluster levels with new cluster identifiers
# names(new.cluster.ids) <- levels(integrated.data)
# integrated.labeled <- RenameIdents(integrated.data, new.cluster.ids)
#
# # Create a DimPlot to visualize labeled clusters using UMAP
# DimPlot(integrated.labeled,
#         reduction = 'umap',
#         label = TRUE,
#         repel = FALSE,
#         pt.size = 0.3,
#         cols = c("#E6D2B8", "#CCAA7A", "#B33E52", "#E6B8BF", "#CC7A88",
#                  "#0075DC", "#4C005C", "#993F00")) + # Color codes for different clusters
#   ggtitle('Labeled by Cell Type') +
#   theme(text = element_text(size = 10),
#         axis.text = element_text(size = 10))
#
# # Cell Type Labeling and Visualization
#
# # Set 'fine.cells.type' to the active identifier in 'integrated.labeled'
# integrated.labeled$fine.cells.type <- integrated.labeled@active.ident
#
# # Convert 'fine.cells.type' to character and assign to 'cell.type'
# integrated.labeled$cell.type <- as.character(integrated.labeled$fine.cell.type)
#
# # Assign specific cell types to the 'cell.type' column
# integrated.labeled$cell.type[integrated.labeled$fine.cell.type %in%
#                                c("SEMA3E+ GABA",
#                                  "SST+ GABA",
#                                  "GABA")] <- "GABA"
#
# # Assign 'unidentified' to 'cell.type' for certain categories
# integrated.labeled$cell.type[integrated.labeled$cell.type %in%
#                                c("unknown neuron",
#                                  "immature neuron",
#                                  "unknown")] <- "unidentified"
#
# # Copy 'cell.type' to 'cell.type.forplot'
# integrated.labeled$cell.type.forplot <- integrated.labeled$cell.type
#
# # Assign specific cell types to 'cell.type' based on 'fine.cell.type'
# integrated.labeled$cell.type[integrated.labeled$fine.cell.type == "NEFM- glut"] <- "nmglut"
# integrated.labeled$cell.type[integrated.labeled$fine.cell.type == "NEFM+ glut"] <- "npglut"
#
# # Print unique cell types
# print(unique(integrated.labeled$cell.type))
#
# # Print unique cell types for plotting
# print(unique(integrated.labeled$cell.type.forplot))
#
# # Create a new column 'cell.type.counts' and populate with 'cell.type.forplot'
# integrated.labeled$cell.type.counts <- integrated.labeled$cell.type.forplot
#
# # Get unique cell types
# types <- unique(integrated.labeled$cell.type.forplot)
#
# # Loop through each cell type to count occurrences and update 'cell.type.counts'
# for (i in 1:length(types)) {
#   count <- sum(integrated.labeled$cell.type.forplot == types[i])
#   print(count)
#   integrated.labeled$cell.type.counts[integrated.labeled$cell.type.counts == types[i]] <-
#     paste0(types[i], "\n", count)
# }
#
# # Print unique cell types with counts
# print(unique(integrated.labeled$cell.type.counts))
#
# # Create a DimPlot to visualize clusters labeled by cell type counts
# DimPlot(integrated.labeled, reduction = "umap", group.by = "cell.type.counts",
#         label = TRUE, repel = FALSE, pt.size = 0.3,
#         cols = c("#B33E52", "#E6D2B8", "#CCAA7A", "#54990F")) +
#   ggtitle("Labeled by Cell Type") +
#   theme(text = element_text(size = 10),
#         axis.text = element_text(size = 10))
#
# # Save the labeled integrated data to a file
# save(integrated_labeled, file = "integrated_labeled.RData")
# #
