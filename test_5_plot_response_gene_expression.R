#Xandy
# ~/NVME/scARC_Duan_018/Duan_project_025_RNA/Analysis_part2_GRCh38/20221114_plot_response_gene_expression.R
# plot response gene expression

library(ggplot2)
library(Seurat)
library(cowplot)

setwd('data/Alexi_Duhe/test1')
load("~/Data/Alexi_Duhe/test1/integrated_labeled.RData")

RNAobj <- subset(integrated.labeled, cell.type != "unidentified")
RNAobj$typextime <- ""
types <- unique(RNAobj$cell.type)
times <- unique(RNAobj$time.ident)
for (i in 1:length(types)) {
  for (j in 1:length(times)) {
    RNAobj$typextime[RNAobj$cell.type == types[i] &
                       RNAobj$time.ident == times[j]] <- paste(types[i], times[j], sep = " ")
  }
}

GABA <- subset(integrated.labeled, cell.type == "GABA")
npglut <- subset(integrated.labeled, cell.type == "npglut")
nmglut <- subset(integrated.labeled, cell.type == "nmglut")
unique(RNAobj$typextime)
resp_genes <- c("FOS", "NPAS4", "BDNF", "VGF")
names(resp_genes) <- c("darkred", "darkred", "deepskyblue3", "deepskyblue3")

# feature plot ####
for (i in 1:length(resp_genes)) {
  g <- resp_genes[i]
  c <- names(resp_genes)[i]
  cat(g, c)
  p1 <- FeaturePlot(GABA, features = g, cols = c("grey", c), split.by = "time.ident") +
    theme(legend.position = "right")
  p2 <- FeaturePlot(nmglut, features = g, cols = c("grey", c), split.by = "time.ident") +
    theme(legend.position = "right")
  p3 <- FeaturePlot(npglut, features = g, cols = c("grey", c), split.by = "time.ident") +
    theme(legend.position = "right")
  p <- plot_grid(p1, p2, p3, align = "v", nrow = 3, label_size = 10, labels = c("GABA", "NEFM- glut", "NEFM+ glut"))
  print(p)
}
