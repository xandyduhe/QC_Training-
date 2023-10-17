# 10 Oct 2023
# script to ntegrate human_only_lst

# init ####
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(patchwork)
library(stringr)
# library(ggplot2)
# library(readr)
# library(RColorBrewer)
library(dplyr)
# library(viridis)
# library(graphics)
# library(readxl)
# library(readr)

library(future)

plan("multisession", workers = 2)
options(expressions = 20000)
options(future.globals.maxSize = 1000 * (1024 ^ 3))

## 10 Oct 2023 #####
setwd("/data/Alexi_Duhe/test1")
print("Start loading data")
load("Duan_all_lines_correct_raw_transformed_unmatched_cells_removed_29Sept2023.RData")
rm(transformed_lst_mapped_to_demux_bc)

print("Finished loading data")
# Will run a QC and remove some low qual cells to reduce the total library size
print("Running QC and remove low-qual cells")
QCed.list <- vector(mode = 'list',
                    length = length(human_only_lst))

for (i in 1:length(human_only_lst)) {
	print(paste("i = ",
		    i, ",",
		    names(human_only_lst)[i]))

	  QCed.list[[i]] <- subset(human_only_lst[[i]],
                                   subset = nFeature_RNA > 500 & # get some more stringent QC here
                                   nFeature_RNA < 7500 &
                                   nCount_RNA > 750 & # more stringent QC here
                                   nCount_RNA < 40000 &
                                   percent.mt < 20)
	print(paste(ncol(human_only_lst[[i]]) - ncol(QCed.list[[i]]),
		    "cells removed."))
}

print("Start SelectIntegrationFeatures")
# prepare reciprocal PCA integration ####
# find anchors
# ! change nfeatures = 2000, 3000 is too large (12 Oct 2023
features <- SelectIntegrationFeatures(object.list = QCed.list,
                                      nfeatures = 3000,
                                      fvf.nfeatures = 3000)
print("SelectIntegrationFeatures done, start PrepSCTIntegration")

human_only_lst <- PrepSCTIntegration(QCed.list,
                                     anchor.features = features)
print("PrepSCTIntegration Done, start ScaleData + RunPCA")

human_only_lst <-
  lapply(X = human_only_lst, FUN = function(x) {
  x <- ScaleData(x, features = features)
  x <- RunPCA(x, features = features)
})

print("Saving RDs...")

saveRDS(human_only_lst,
	file = "human_only_lst_scaled_PCAed_4_FindIntegrationAnchors.RDs")
## check if all members of human_only_lst has unique barcodes
# for (i in 1:length(human_only_lst)) {
#   print(paste(i,
#               names(human_only_lst)[i],
#               length(colnames(human_only_lst[[i]])),
#               length(unique(colnames(human_only_lst[[i]]))),
#               (length(colnames(human_only_lst[[i]])) == length(unique(colnames(human_only_lst[[i]])))),
#               sep = ","))
#
# }

print("ScaleData + RunPCA Done, start FindIntegrationAnchors")

anchors <- FindIntegrationAnchors(object.list = human_only_lst,
                                  anchor.features = features,
                                  reference = c(1, 2, 3),
                                  reduction = "rpca",
                                  normalization.method = "SCT",
                                  scale = F,
                                  dims = 1:50)
print("FindIntegrationAnchors done, saving data")
#Found 9036 anchors
#all.genes <- transformed_lst[[1]]@assays$RNA@counts@Dimnames[[1]]
# save(anchors, file = "anchors_after_cleaned_lst_QC.RData")
saveRDS(anchors,
        file = "anchors_after_cleaned_lst_QC_10Oct2023.RDs")
print("Data saved")

# anchors <-
#   readRDS(file = "anchors_after_cleaned_lst_QC_10Oct2023.RDs")
# # rm(cleaned_lst)
# # integrate ####
# integrated <- IntegrateData(anchorset = anchors,
#                             dims = 1:50,
#                             normalization.method = "SCT",
#                             verbose = T)
# saveRDS(integrated,
#         file = "integrated_obj_after_cleaned_lst_10Oct2023.RDs")

# remove library 22----
human_only_4_integration <- vector("list", (length(human_only_lst_scaled_PCAed_4_FindIntegrationAnchors) - 3))
k = 1
for (i in 1:length(human_only_lst_scaled_PCAed_4_FindIntegrationAnchors)){
  #print(i)
  #print(unique(human_only_lst_scaled_PCAed_4_FindIntegrationAnchors[[i]]@meta.data$lib.ident))
  if (i == 7 ){
    next}
  if (i == 8) {
    next }
  if (i == 9) {
    next
    # print('no')
    # obj <- human_only_lst_scaled_PCAed_4_FindIntegrationAnchors[[i]]
    # new_test_obj[[i]] <- obj
  }else {
    print(i)
    obj <- human_only_lst_scaled_PCAed_4_FindIntegrationAnchors[[i]]
    human_only_4_integration[[k]] <- obj
    k=k+1
  }
}

for (i in 1:length(human_only_4_integration)){
  print(unique(human_only_4_integration[[i]]@meta.data$lib.ident))
}

save(human_only_4_integration, file = 'human_only_4_integration_16Oct23.RData')
