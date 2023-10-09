# QC_Training-
Training code for scRNA-seq Quality Control

Part 1: Analysis of Individual sequencing batches 
I.	Raw scRNA-seq data processing 
  1.	For each library: raw single-cell multi-omic sequencing data is organized into FASTQ files
    i.	Then input into 10x Cell Ranger ARC pipeline for analysis using Cell Ranger ARC command “count”
    ii.	Output of Cell Ranger ARC pipeline includes a filtered_feature_bc_matrix.h5 file
  2.	Separate cells by species: Data is mapped to hybrid genome of Hg38 and mm10/Rn6 
    i.	(rat vs human)
  3.	The filtered count matrix (filtered_feature_bc_matrix.h5 ) mapped to the hybrid genome was imported into R
  4.	Preliminary analyses: extract a list of barcodes for putative human neurons
    i.	k-nearest neighbor (KNN) clustering performed on data 
  5.	This list then input into demuxlet: demultiplex the cell lines co-cultured in the library, returns list of barcodes for human neurons -> used in the next step.
II.	Quality Control 
  1.	scRNA-seq library QC metrics examined for raw data and human neurons after separation from mouse/rat astrocytes 
    i.	# of cells in the library
    ii.	# of read counts (unique molecular identifier, UMI) in the library
    iii.	# of read counts (UMI) per cell
    iv.	distribution of # of UMI in cells (nCount)
    v.	distribution of # of genes detected in cells (nFeature)
    vi.	distribution of % of reads mapped to mitochondrial genes in cells (percent.mt).
  2.	scATAC-seq library QC metrics examined for raw data
    i.	% of fragment in peaks (fRiP) 
    ii.	# of reads in the library
    iii.	nucleosome signal 
  1.	(ratio of mononucleosome fragments to nucleosome-free fragments – cal by Signac) 
    iv.	TSS enrichment
  3.	Cells were filtered by
    i.	# of genes detected in each cell (nFeature)
    ii.	# of read counts in each cell (nCount)
    iii.	% of reads mapped to mitochondrial genes (percent.mt)
    iv.	nucleosome signal
    v.	TSS enrichment. 
      1.	Cutoff: 
      2.	300 < nFeature < 7500
      3.	500 < nCount < 40000
      4.	percent.mt < 20
      5.	nucleosome_signal < 2
      6.	TSS enrichment > 2
III.	Cell type assignment using scRNA-seq data 
  1.	One integrated object was generated after “Quality Control”
    i.	For each sequencing batch
  2.	Integrated object is scaled and clustered using Seurat. 
  3.	The results were illustrated on the UMAP space
    i.	a list of marker genes was used to assign cell types to each cluster (GABA, nmglut, npglut)
  4.	Each cluster assigned a cell-type identity 
    i.	based on marker gene expression,
    ii.	main identities: inhibitory neurons and excitatory neurons. 
    iii.	Each cell within cluster labeled by the identity of cluster. 
  5.	IF current list of marker genes was insufficient to identify clusters
    i.	marked as "unidentified" -> removed from downstream analysis
IV.	Integration of scATAC-seq data
  1.	Import scATAC-seq data for pre-comparison-manipulations  
    i.	scATAC-seq data separately integrated using the 10x -aggr command in Linux directly on 10x Cell Ranger output files. 
    ii.	Resulting integrated count matrix imported into R as a ChromatinAssay object
    iii.	ChromatinAssay Object made into a Seurat object
      1.	Integrated scRNA-seq object was added as the RNA assay to the scATAC-seq object resulting in a multi-omic object
    iv.	The scATAC-seq assay was normalized using TF-IDF and SVD in Signac 
  2.	Compare scATAC-seq clusters to scRNA-seq clusters 
    i.	Clusters found for the scATAC-seq data
    ii.	Cell type labels from scRNA-seq analysis transferred to scATAC-seq data by matching barcodes
    iii.	Clustering result of scATAC-seq data examined against cell type labels to detect any discrepancies between scATAC-seq and scRNA-seq clustering results. 
V.	Calling cell-type and time point specific peaks (scATAC-seq)
  1.	After obtaining labels from scRNA-seq data analysis, scATAC-seq count matrix was separated into groups that corresponded to each main cell type 
    i.	two excitatory neuron subtypes, one inhibitory neuron cell type, and other cell types
    ii.	each time point (0hr, 1hr, and 6hrs of stimulation). 
  2.	For each of these groups, MACS2 was used to call peaks on the subsetted count matrices, 
  3.	the resulting peak sets were then combined by setting the argument "combine.peaks" to true in the Signac function CallPeaks. 
  4.	The resulting union peak set was applied to the multi-omic object from “Integration of scATAC-seq data”
  5.	New ChromatinAssay object was created and added to the multi-omic object by counting fragments using the new union peak set.


Part 2: Integration of Sequencing batches and analysis of integrated data 
  I.	Raw RNA-seq data processing 
    a.	Similar to the analysis of individual sequencing batches described in part I, 
    b.	the integration of all libraries from all 6 sequencing batches starts from the normalization of individual scRNA-seq libraries. 
    c.	From each library, the filtered_feature_bc_matrix.h5 from mapping the raw sequencing data to Hg38 genome was imported into R and normalized by SCTransform. 
    d.	Using demuxlet-generated barcodes corresponding to each human cell line, the raw data from each library was matched against the barcodes, and 
    e.	any unmatched cells were discarded. 
      i  .	* note on sequencing batch 022: out of 20 cell lines in this batch, 15 cell lines (CD_22, CD_23, CD_31, CD_36, CD_37, CD_38, CD_39, CD_40, CD_44, CD_45, CD_46, CD_47, CD_48, CD_58, CD_61) were repeated in later sequencing batches, so the data from these   cell lines were removed from the count matrices from batch 022. Due to the removal of these lines, one library in batch 022 (library 22) was entirely removed because all cell lines in this library were repeated.
  II.	Integration of scRNA-seq and scATAC-seq libraries 
    a.	All scRNA-seq libraries were prepared for integration using the Seurat package; 
      i.	integration anchors were found for the list of libraries using reciprocal PCA, which 
        1.	project datasets into each other’s PCA space and obtain cells that are mutual nearest neighbors. 
      ii.	Finally, all libraries were integrated using the IntegrateData() function in Seurat. 
    b.	Similar to the method described in part I, scATAC-seq data was separately integrated using the 10x -aggr command in Linux on 10x Cell Ranger output files from all libraries.
      i.	The resulting integrated count matrix was then imported into R 
      ii.	the integrated RNA assay was added to the scATAC-seq object, resulting in a multi-omic object. 
    c.	Integrated scRNA-seq data was clustered. Clusters were labeled as described in part 1 
      i.	The cell type labels from scRNA-seq analysis were transferred to scATAC-seq data by matching barcodes. 
III.	Quality Control 
  a.	The same QC steps as part I were applied to the integrated multi-omic object. 
  b.	The cutoff 
    i.	300 < nFeature < 7500
    ii.	500 < nCount < 40000
    iii.	percent.mt < 15
    iv.	nucleosome_signal < 2
    v.	TSS enrichment > 2

    
Part 3: Differential gene expression and chromatin accessibility in stimulated vs. not stimulated groups 
  I.	The goal of this analysis: to identify genes and peaks differentially expressed/accessible in stimulated (1hr or 6hrs of stimulation) versus unstimulated (0hr of stimulation) cells.
  II.	Pseudobulk count matrices were created 
    a.	summing up the counts across all cells in one cell line for each gene/peak. 
  III.	Combat-seq was then used to correct the effect of co-cultured batch. 
  IV.	Genes/peaks were then filtered by the following criterion: 
    a.	cpm > 1 in half of all samples in any of the 3 time points. 
   V.	Normalization of pseudobulk counts using Limma-voom 
    a.	for each cell type and each pair of time points (1vs0hr, 6vs0hr). 
  VI.	A linear model was fit for gene expression
    a.	covariates 
    b.	age
    c.	sex
    d.	disease status
    e.	co-cultured group
    f.	cellular composition
   VII.	Differential gene expression/chromatin accessibility was tested using the Limma function contrasts.fit
    a.	then smoothed using eBayes. 
  III.	The p-values from the test results were adjusted using the BH method
    a.	the results were filtered by FDR < 0.05.

    
Part 4: Differential gene expression in case vs control cell lines 
  I.	Model-based Analysis of Single Cell Transcriptomics (MAST). 
    a.	Differential gene expression analysis:
      i.	comparing case and control cell lines
  II.	
  III.	
  IV.	Starting from the integrated Seurat object that contains all libraries
    a.	covariates assigned to every cell: age, sex, disease status, co-cultured group, and cellular composition 
  V.	The data separated by cell type and time point into subsets
    a.	For each subset, a SingleCellAssay object was created using the FromMatrix function
  VI.	Genes were filtered by the criterion:
    a.	greater than 5% of cells contain nonzero values for the given gene
  VII.	Zero-inflated linear regression model from MAST was fitted to the expression value using the model
    a.	y = b1*affected status + b2*age + b3*sex + b4*co-cultured batch + b5*cell type fraction. 
  VIII.	Likelihood ratio test was used to test the difference between SZ case and control
    a.	where control was used as the baseline. 
    b.	The p-values from the test results were adjusted using the BH method, and the results were filtered by FDR < 0.05.

