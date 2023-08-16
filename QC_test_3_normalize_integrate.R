# 8/4/2023 4:30 PM

# /nvmefs/scARC_Duan_018/Duan_project_025_RNA/Analysis_part2_GRCh38/20230507_normalize_integrate_cluster_025_17_46.R


library(Seurat)
library(ggplot2)
library(pals)
library(stringr)
library(future)
library(glmGamPoi)
library(cowplot)

library(RColorBrewer)

setwd("/data/Alexi_Duhe/test1/")


plan("multisession", workers = 8)
options(expressions = 20000)
options(future.globals.maxSize = 207374182400)

# load previous data file created in QC_test_2_demultiplex.R
load('GRCh38_mapped_after_QC_list.RData')



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
                                seed.use = 2023)

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

# # Scale data for downstream analysis
# integrated.data <-  ScaleData(integrated.data)
#
# # Run PCA on scaled data
# integrated.data <- RunPCA(integrated.data, verbose = TRUE,
#                           seed.use = 2023)
#
# # Run UMAP on PCA results
# integrated.data <- RunUMAP(integrated.data, reduction = 'pca',
#                            dims = 1:30,
#                            seed.use = 2023)
#
# # Find nearest neighbors based on UMAP coordinates
# integrated.data <- FindNeighbors(integrated.data, reduction = 'pca',
#                                  dims = 1:30)
#
# # Perform clustering using Louvain algorithm
# integrated.data <- FindClusters(integrated.data, resolution = 0.5,
#                                 random.seed = 2023)

# PCA for data integration ----


# Save integrated.data before transformati
# load(integrated_object.RData)
# Perform Quality control on the integrated data ----

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

# Clustering and visualization ----

# Identify unique orig identifiers and time identifiers
unique(integrated.data$orig.ident)
unique(integrated.data$time.ident)

# Set the default assay to 'integrated' for downstream analysis
DefaultAssay(integrated.data) <- 'integrated'
#
# Scale data for downstream analysis
integrated.data <-  ScaleData(integrated.data)

# Run PCA on scaled data
integrated.data <- RunPCA(integrated.data, verbose = TRUE,
                          seed.use = 2023)

# Run UMAP on PCA results
integrated.data <- RunUMAP(integrated.data, reduction = 'pca',
                           dims = 1:30,
                           seed.use = 2023)

# Find nearest neighbors based on UMAP coordinates
integrated.data <- FindNeighbors(integrated.data, reduction = 'pca',
                                 dims = 1:30)

# Perform clustering using Louvain algorithm
integrated.data <- FindClusters(integrated.data, resolution = 0.5,
                                random.seed = 2023)

# Define a color pallet for library-based grouping
library.colors <- DiscretePalette(n = length(unique(Idents(integrated.data))),
                                  'alphabet')
# Visualize clustering results using a dimentionality reduction plot
p1 <- DimPlot(integrated.data,
              cols = library.colors,
              label = T,
              group.by = 'integrated_snn_res.0.5',
              raster = F) +
  NoLegend() +
  ggtitle('Clustering RNAseq data') +
  theme(text = element_text(size = 10),
        axis.text = element_text(size = 10))
print(p1)

# Create DimPlot based on orig identifiers
DimPlot(integrated.data,
        label = FALSE,
        cols = library.colors,
        group.by = 'orig.ident',
        raster = F) +
  ggtitle('by library') +
  theme(text = element_text(size = 10),
        axis.text = element_text(size = 10))

DimPlot(integrated.data,
        label = FALSE,
        group.by = 'time.ident',
        raster = F) +
  ggtitle('By time point') +
  theme(text = element_text(size = 10),
        axis.text = element_text(size = 10))

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
            raster = F,
            ncol = 2)
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
  "SLC17A7", "SERTAD4", "rna_FOXG1",  # Forebrain markers
  "POU3F2", "LHX2", # General cortex markers
  "ADCYAP1", "CUX1", "CUX2", "MAP2", "DCX" # Various neuronal markers
)

# Convert the gene markers vector to a factor with specified levels
trimmed_markers <- factor(trimmed_markers, levels = trimmed_markers)


# Create a stacked violin plot using trimmed maker genes
VlnPlot(integrated.data,
        features = trimmed_markers,
        stack = TRUE,
        sort = F,
        flip = TRUE) +
  theme(legend.position = 'none') +
  ggtitle('Expression by Gene Identity')



# Set identities of the integrated data to 'seurat_clusters'
Idents(integrated.data) <- 'seurat_clusters'

# Define new cluster identifiers
new.cluster.ids <-
  c("NEFM- glut", "NEFM+ glut", "NEFM- glut", "NEFM+ glut", "GABA",
    "GABA", "NEFM- glut", "SST+ GABA", "NEFM- glut", "NEFM- glut",
    "SEMA3E+ GABA", "NEFM- glut", "NEFM- glut", "GABA", "GABA",
    "unknown neuron", "immature neuron", "unknown", 'unknown', 'unknown')

# length(new.cluster.ids)
#length(unique(integrated.data$seurat_clusters))

### Siwei test code #####
DimPlot(integrated.data,
        cols = library.colors,
        label = T,
        group.by = 'integrated_snn_res.0.5',
        raster = F)

FeaturePlot(integrated.data,
            features = c('NEFM'),
            raster = F)

new.cluster.ids[integrated.data$seurat_clusters]
integrated.data$new_cluster_ids <-
  "unidentified"
integrated.data$new_cluster_ids[integrated.data$seurat_clusters %in% c("0", "6", "7")] <-
  "NEFM+ glut"
integrated.data$new_cluster_ids[integrated.data$seurat_clusters %in% c("2", "3")] <-
  "NEFM- glut"
integrated.data$new_cluster_ids[integrated.data$seurat_clusters %in% c("1", "4", "5",
                                                                       "9", "10", "11", "12")] <-
  "GABA"
integrated.data$new_cluster_ids <-
  factor(integrated.data$new_cluster_ids,
         levels = c("NEFM+ glut",
                    "NEFM- glut",
                    "GABA",
                    "unidentified"))

Idents(integrated.data) <- "new_cluster_ids"

# Create DimPlot to visualize labeled clusters using UMAP
DimPlot(integrated.data,
        reduction = 'umap',
        label = TRUE,
        repel = FALSE,
        raster = F,
        pt.size = 0.3,
        cols = c("#B33E52",
                 "#E6D2B8",
                 "#CCAA7A",
                 "#54990F",
                 "#0075DC")) +
  ggtitle('labeled by cell type') +
  theme(text = element_text(size = 10),
        axis.text = element_text(size = 10))
###

# Rename cluster levels with new cluster identifiers
names(new.cluster.ids) <- levels(integrated.data)
integrated.labeled <- RenameIdents(integrated.data,
                                   new.cluster.ids)

# Create DimPlot to visualize labeled clusters using UMAP
DimPlot(integrated.labeled,
        reduction = 'umap',
        label = TRUE,
        repel = FALSE,
        raster = F,
        pt.size = 0.3,
        cols = c("#E6D2B8", #nmglut
                 "#CCAA7A", #npglut
                 "#B33E52", #GABA
                 "#E6B8BF", #SST GABA
                 "#CC7A88", #SEMA3E GABA
                 "#0075DC", #immature neuron
                 "#4C005C", #unknown neuron
                 "#993F00", #unknown
                 "#996F00")) + #unknown
  ggtitle('labeled by cell type') +
  theme(text = element_text(size = 10),
        axis.text = element_text(size = 10))

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
        label = TRUE, repel = F, pt.size = 0.3, raster = F,
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
