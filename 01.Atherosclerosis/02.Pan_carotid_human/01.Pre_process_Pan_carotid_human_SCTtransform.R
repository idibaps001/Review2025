#pre-process Pan_carotid_human scRNAseq data | SCTtransform 
library(Seurat)
library(tidyverse)
library(scDblFinder)
library(celda)
library(sctransform)
library(data.table)
library(cluster)
library(scRNAutils)
# Set seed for doublet removal reproducibility
set.seed(1)
# PROCESS HUMAN CAROTID MATRICES USING SCTransform
setwd("/Users/lucia/Documents/34.immunometabolism/data/Pan_carotid_human")

#######################################################################
# Read Sample 1: GSM4705589 RPE004 (Caucasian, Male, 83, Symptomatic) #
#######################################################################
# Load Matrix
rpe004_matrix = fread("GSE155512_RAW/GSM4705589_RPE004_matrix.txt.gz", sep="\t", header = TRUE)
rownames(rpe004_matrix) = rpe004_matrix$gene
head(rownames(rpe004_matrix))
dim(rpe004_matrix)#[1] 15796  2615

# Make sparse matrix and remove "gene" column
rpe004_sparse_mtx = as.sparse(rpe004_matrix)
rpe004_sparse_mtx = rpe004_sparse_mtx[, -1]

# 15796 genes x 2614 cells
dim(rpe004_sparse_mtx)

##################################
# STEP 1: Pre-processing (Create Seurat obj containing doublets)
# Dims of seurat object: 15179 genes x 2614 cells
rpe004_seurat_sct = CreateSeuratObject(counts = rpe004_sparse_mtx, 
                                       project = "pan_rpe004_carotid", 
                                       min.cells = 10, 
                                       min.features = 200)

###################################
# STEP 2 before/after removing doublets

Seurat_SCT_process = function(seurat_obj, seurat_filter=FALSE, 
                              sample_id, study_name, artery, disease_status){
  seurat_obj$sample = sample_id
  seurat_obj$study = study_name
  # Check for mt, hb percentage and other quality metrics 
  seurat_obj[["percent.mt"]] = PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  
  # NOTE: Skip this step during the first processing round to be used for doublet detection
  if (seurat_filter) {
      seurat_obj =  subset(seurat_obj, subset= nFeature_RNA > 200 & nCount_RNA > 300)
  }
  
  # Normalize data, find variable genes, scale data and regress out percent.mt variance
  # SCT enables extraction of meaningful insights from more PCs so we'll set dims=1:30
  # SCTransform Arg internally uses glmGamPoi
  seurat_obj = SCTransform(seurat_obj, method = "glmGamPoi", 
                           vars.to.regress = "percent.mt") %>%
    RunPCA( verbose = FALSE) %>% 
    RunUMAP(dims = 1:30, verbose = FALSE) %>%
    FindNeighbors(reduction = "pca", dims = 1:30, k.param = 20, verbose = FALSE) %>%
    FindClusters( verbose = FALSE)
    return(seurat_obj)
}

# SCT normalize
rpe004_seurat_sct = Seurat_SCT_process(rpe004_seurat_sct, seurat_filter = FALSE,
                                       sample_id = "pan_rpe004", 
                                       study_name = "pan_et_al")
                                       
# Clean counts matrix (without doublets) from ambient RNA using the decontX wrapper
   decontX_remove = function(seurat_obj) {
                                         decontX_results = celda::decontX(seurat_obj@assays$RNA$counts)
                                         decontaminated_matrix = decontX_results$decontXcounts
                                         seurat_obj@assays$RNA$counts = decontaminated_matrix
                                         return(seurat_obj)
                                       }
                                       
    rpe004_seurat_sct =  decontX_remove(rpe004_seurat_sct)
    head(rpe004_seurat_sct@assays$RNA$counts)
                                                                            
# Finding the right resolution is quite important for doublet detection and removal
rpe004_seurat_sct =  FindClusters(rpe004_seurat_sct, resolution = 0.5)   

# Visualize clusters
p1_before_QC = DimPlot(rpe004_seurat_sct, reduction = "umap", label = TRUE) +
  ggtitle("Before removing doublets")

##################################
# STEP 3
# Detect doublets with scDblFinder
# Called 27 doublets

sce <- as.SingleCellExperiment(rpe004_seurat_sct)
sce <- scDblFinder(sce,clusters = "SCT_snn_res.0.5")
table(sce$scDblFinder.class) #2587 singlet, 27 doublet
metadata <-  sce@colData %>% as.data.frame()
rpe004_seurat_sct@meta.data <- metadata

# First filtering step: Use the vector of consensus doublet IDs output below
rpe004_seurat_sct <- subset(rpe004_seurat_sct, subset = scDblFinder.class == "singlet")
ncol(rpe004_seurat_sct)#2587

# New dims after removing doublets: 15174 genes x 2587 cells

###################################
# STEP 4 after removing doublets

# SCT normalize
rpe004_seurat_sct = Seurat_SCT_process(rpe004_seurat_sct, seurat_filter = TRUE,
                                       sample_id = "pan_rpe004", 
                                       study_name = "pan_et_al") #15179 genes * 2587 cells
ncol(rpe004_seurat_sct)

p1_after_QC = DimPlot(rpe004_seurat_sct, reduction = "umap", label = TRUE) +
  ggtitle("After removing doublets + decontX")

FeaturePlot(rpe004_seurat_sct, features = c("MYH11", "CNN1", "TNFRSF11B", "ACTA2", 
                                            "CRTAC1", "LUM"), order = TRUE,  pt.size = 0.1)
FeaturePlot(rpe004_seurat_sct, features = c("CD68", "CD14", "CD163", "APOE", 
                                            "S100A8", "IL1B"), order = TRUE,  pt.size = 0.1)
FeaturePlot(rpe004_seurat_sct, features = c("MYH11", "CNN1", "FHL5", "FOXC1"), order = TRUE,  pt.size = 0.1)
# Save RDS object
# NUmber of cells after processing: 2587
saveRDS(rpe004_seurat_sct, "pan_rpe004_processed_SCTtransformed.rds")

#########################################################################
# Read Sample 2: GSM4705590 RPE005 (Caucasian, Male, 67, Asymptomatic)  #
#########################################################################

# Load Matrix
rpe005_matrix = fread("GSE155512_RAW/GSM4705590_RPE005_matrix.txt.gz", sep="\t", header = TRUE)
rownames(rpe005_matrix) = rpe005_matrix$gene
head(rownames(rpe005_matrix))
dim(rpe005_matrix)#[1] 17397  3487
# Make sparse matrix and remove "gene" column
rpe005_sparse_mtx = as.sparse(rpe005_matrix)
rpe005_sparse_mtx = rpe005_sparse_mtx[, -1]

# 17397 genes x 3486 cells
dim(rpe005_sparse_mtx)

##################################
# STEP 1: Pre-processing (Create Seurat obj containing doublets)
# Dims of seurat object: 16670 genes x 3486 cells
rpe005_seurat_sct = CreateSeuratObject(counts = rpe005_sparse_mtx, 
                                       project = "pan_rpe005", 
                                       min.cells = 10, 
                                       min.features = 200)

###################################
# STEP 2 before/after removing doublets
# SCT normalize
rpe005_seurat_sct = Seurat_SCT_process(rpe005_seurat_sct, seurat_filter = FALSE,
                                       sample_id = "pan_rpe005", 
                                       study_name = "pan_et_al") #16670 genes * 3486 cells


# Clean counts matrix (without doublets) from ambient RNA using the decontX wrapper
rpe005_seurat_sct =  decontX_remove(rpe005_seurat_sct)
head(rpe005_seurat_sct@assays$RNA$counts)

# Finding the right resolution is quite important for doublet detection and removal
rpe005_seurat_sct =  FindClusters(rpe005_seurat_sct, resolution = 0.5)   

# Visualize clusters
p1_before_QC = DimPlot(rpe005_seurat_sct, reduction = "umap", label = TRUE) +
  ggtitle("Before removing doublets")

##################################
# STEP 3
# Detect doublets with scDblFinder
# Called 17 doublets

sce <- as.SingleCellExperiment(rpe005_seurat_sct)
sce <- scDblFinder(sce,clusters = "SCT_snn_res.0.5")
table(sce$scDblFinder.class) #3469 singlet, 17 doublet
metadata <-  sce@colData %>% as.data.frame()
rpe005_seurat_sct@meta.data <- metadata

# First filtering step: Use the vector of consensus doublet IDs output below
rpe005_seurat_sct <- subset(rpe005_seurat_sct, subset = scDblFinder.class == "singlet")
ncol(rpe005_seurat_sct)#3469

# New dims after removing doublets: 16670 genes x 3469 cells

###################################
# STEP 4 after removing doublets

# SCT normalize
rpe005_seurat_sct = Seurat_SCT_process(rpe005_seurat_sct, seurat_filter = TRUE,
                                       sample_id = "pan_rpe005", 
                                       study_name = "pan_et_al") # 16670 genes *   3469 cells

p1_after_QC = DimPlot(rpe005_seurat_sct, reduction = "umap", label = TRUE) +
  ggtitle("After removing doublets + decontX")

FeaturePlot(rpe005_seurat_sct, features = c("MYH11", "CNN1", "TNFRSF11B", "ACTA2", 
                                            "CRTAC1", "LUM"), order = TRUE,  pt.size = 0.1)
FeaturePlot(rpe005_seurat_sct, features = c("CD68", "CD14", "CD163", "APOE", 
                                            "S100A8", "IL1B"), order = TRUE,  pt.size = 0.1)
FeaturePlot(rpe005_seurat_sct, features = c("MYH11", "CNN1", "FHL5", "FOXC1"), order = TRUE,  pt.size = 0.1)
# Save RDS object
# NUmber of cells after processing: 3469
saveRDS(rpe005_seurat_sct, "pan_rpe005_processed_SCTtransformed.rds")

##########################################################################
# Read Sample 3: GSM4705591 RPE006 (Caucasian, Female, 76, Asymptomatic) #
##########################################################################
# Load Matrix
rpe006_matrix = fread("GSE155512_RAW/GSM4705591_RPE006_matrix.txt.gz", sep="\t", header = TRUE)

rownames(rpe006_matrix) = rpe006_matrix$gene
head(rownames(rpe006_matrix))
dim(rpe006_matrix)

# Make sparse matrix and remove "gene" column
rpe006_sparse_mtx = as.sparse(rpe006_matrix)
rpe006_sparse_mtx = rpe006_sparse_mtx[, -1]

# 15687 genes x 2767 cells
dim(rpe006_sparse_mtx)

##################################
# STEP 1: Pre-processing (Create Seurat obj containing doublets)
# Dims of seurat object: 15445 genes x 2767 cells
rpe006_seurat_sct = CreateSeuratObject(counts = rpe006_sparse_mtx, 
                                       project = "pan_rpe006", 
                                       min.cells = 10, 
                                       min.features = 200)

###################################
# STEP 2 before/after removing doublets
# SCT normalize
rpe006_seurat_sct = Seurat_SCT_process(rpe006_seurat_sct, seurat_filter = FALSE,
                                       sample_id = "pan_rpe006", 
                                       study_name = "pan_et_al") #15445 genes * 2767 cells


# Clean counts matrix (without doublets) from ambient RNA using the decontX wrapper
rpe006_seurat_sct =  decontX_remove(rpe006_seurat_sct)
head(rpe006_seurat_sct@assays$RNA$counts)

# Finding the right resolution is quite important for doublet detection and removal
rpe006_seurat_sct =  FindClusters(rpe006_seurat_sct, resolution = 0.2)   

# Visualize clusters
p1_before_QC = DimPlot(rpe006_seurat_sct, reduction = "umap", label = TRUE) +
  ggtitle("Before removing doublets")

##################################
# STEP 3
# Detect doublets with scDblFinder
# Called 40 doublets

sce <- as.SingleCellExperiment(rpe006_seurat_sct)
sce <- scDblFinder(sce,clusters = "SCT_snn_res.0.2")
table(sce$scDblFinder.class) #2727 singlet, 40 doublet
metadata <-  sce@colData %>% as.data.frame()
rpe006_seurat_sct@meta.data <- metadata

# First filtering step: Use the vector of consensus doublet IDs output below
rpe006_seurat_sct <- subset(rpe006_seurat_sct, subset = scDblFinder.class == "singlet")
dim(rpe006_seurat_sct)#2727

# New dims after removing doublets: 15445 genes x 2727 cells

###################################
# STEP 4 after removing doublets

# SCT normalize
rpe006_seurat_sct = Seurat_SCT_process(rpe006_seurat_sct, seurat_filter = TRUE,
                                       sample_id = "pan_rpe006", 
                                       study_name = "pan_et_al") # 15445 genes *   2663 cells

ncol(rpe006_seurat_sct)

p1_after_QC = DimPlot(rpe006_seurat_sct, reduction = "umap", label = TRUE) +
  ggtitle("After removing doublets + decontX")

FeaturePlot(rpe006_seurat_sct, features = c("MYH11", "CNN1", "TNFRSF11B", "ACTA2", 
                                            "CRTAC1", "LUM"), order = TRUE,  pt.size = 0.1)
FeaturePlot(rpe006_seurat_sct, features = c("CD68", "CD14", "CD163", "APOE", 
                                            "S100A8", "IL1B"), order = TRUE,  pt.size = 0.1)
FeaturePlot(rpe006_seurat_sct, features = c("MYH11", "CNN1", "FHL5", "FOXC1"), order = TRUE,  pt.size = 0.1)
# Save RDS object
# NUmber of cells after processing: 2727
saveRDS(rpe006_seurat_sct, "pan_rpe006_processed_SCTtransformed.rds")

######################################################################################################
#                                                                                                    #
#       Add metadata &  Prefix_barcode & SingleR annotation & Subset monocytes and macrophages       #
#                                        rpe_004                                                     #
######################################################################################################
library(Seurat)
library(tidyverse)
setwd("/Users/lucia/Documents/34.immunometabolism/data/Pan_carotid_human")
rpe004_seurat_sct <- readRDS("pan_rpe004_processed_SCTtransformed.rds")

rpe004_seurat_sct$Project <- "Pan_rpe004_carotid"
rpe004_seurat_sct$GSE_ID <- "GSE155512"
rpe004_seurat_sct$GSM_ID <- "GSM4705589"
rpe004_seurat_sct$Tissue <- "carotid artery atherosclerotic core (AC)"
rpe004_seurat_sct$Study_name <- "Pan_et_al"
rpe004_seurat_sct$Species <- "human"
rpe004_seurat_sct$Age_group <- "old"
rpe004_seurat_sct$orig.ident <- "Pan_rpe004_carotid"
View(rpe004_seurat_sct@meta.data)

rpe004_seurat_sct <- RenameCells(rpe004_seurat_sct,add.cell.id = "Pan_rpe004")

library(SingleR)
library(celldex)
surveyReferences()
ref.data <- HumanPrimaryCellAtlasData()
#ensembl = TRUE
#ref.data
#table(ref.data$label.main)
#table(ref.data$label.fine)
#table(ref.data$label.ont)

sce <- LayerData(rpe004_seurat_sct, assay = "RNA", layer = "counts")
predictions <- SingleR(test=sce, assay.type.test=1, 
                       ref=ref.data, labels=ref.data$label.main)


table(predictions$labels)
plotScoreHeatmap(predictions)#picture size: 840*600
rpe004_seurat_sct$SingleR <- predictions$labels
saveRDS(rpe004_seurat_sct, "pan_rpe004_processed_SCTtransformed.rds")

Idents(rpe004_seurat_sct) = "SingleR"
rpe004_seurat_sct.monomac <- subset(rpe004_seurat_sct, idents = c("Monocyte","Macrophage"))
Idents(rpe004_seurat_sct.monomac) = "SingleR"
rpe004_seurat_sct.monomac <- subset(rpe004_seurat_sct.monomac, idents = c("Monocyte","Macrophage"))
table(rpe004_seurat_sct.monomac$SingleR)
# Macrophage   Monocyte 
#   162        51 
dim(rpe004_seurat_sct.monomac)#[1] 15179   213

# SCT normalize
rpe004_seurat_sct.monomac = Seurat_SCT_process(rpe004_seurat_sct.monomac, seurat_filter = TRUE,
                                       sample_id = "pan_rpe004", 
                                       study_name = "pan_et_al")

# SCT normalize
rpe004_seurat_sct.monomac = Seurat_SCT_process(rpe004_seurat_sct.monomac, seurat_filter = TRUE,
                                               sample_id = "pan_rpe004", 
                                               study_name = "pan_et_al")
dim(rpe004_seurat_sct.monomac)#[1] 11118   213

colnames(rpe004_seurat_sct.monomac@meta.data)
metadata <- rpe004_seurat_sct.monomac@meta.data
metadata$barcode <- rownames(metadata)
#View(metadata)
metadata2 <- metadata[,c("orig.ident","Project","GSE_ID", "GSM_ID", "Tissue",  "Study_name","Species", "Age_group", "SingleR","nCount_RNA" , "nFeature_RNA","percent.mt" , "nCount_SCT", "nFeature_SCT", "scDblFinder.class" )]
#View(metadata2)
rpe004_seurat_sct.monomac@meta.data <- metadata2
Idents(rpe004_seurat_sct.monomac) = "SingleR"
DimPlot(rpe004_seurat_sct.monomac)

#11118 genes * 213 cells
saveRDS(rpe004_seurat_sct.monomac, "pan_rpe004_processed_SCTtransformed_monomac.rds")

rpe004_seurat_sct.monomac <- readRDS("pan_rpe004_processed_SCTtransformed_monomac.rds")

######################################################################################################
#                                                                                                    #
#       Add metadata &  Prefix_barcode & SingleR annotation & Subset monocytes and macrophages       #
#                                        rpe_005                                                     #
######################################################################################################
library(Seurat)
library(tidyverse)
setwd("/Users/lucia/Documents/34.immunometabolism/data/Pan_carotid_human")
rpe005_seurat_sct <- readRDS("pan_rpe005_processed_SCTtransformed.rds")

rpe005_seurat_sct$Project <- "Pan_rpe005_carotid"
rpe005_seurat_sct$GSE_ID <- "GSE155512"
rpe005_seurat_sct$GSM_ID <- "GSM4705590"
rpe005_seurat_sct$Tissue <- "carotid artery atherosclerotic core (AC)"
rpe005_seurat_sct$Study_name <- "Pan_et_al"
rpe005_seurat_sct$Species <- "human"
rpe005_seurat_sct$Age_group <- "old"
rpe005_seurat_sct$orig.ident <- "Pan_rpe005_carotid"
View(rpe005_seurat_sct@meta.data)

rpe005_seurat_sct <- RenameCells(rpe005_seurat_sct,add.cell.id = "Pan_rpe005")

library(SingleR)
library(celldex)
surveyReferences()
ref.data <- HumanPrimaryCellAtlasData()
#ensembl = TRUE
#ref.data
#table(ref.data$label.main)
#table(ref.data$label.fine)
#table(ref.data$label.ont)

sce <- LayerData(rpe005_seurat_sct, assay = "RNA", layer = "counts")
predictions <- SingleR(test=sce, assay.type.test=1, 
                       ref=ref.data, labels=ref.data$label.main)

table(predictions$labels)
plotScoreHeatmap(predictions)#picture size: 840*600
rpe005_seurat_sct$SingleR <- predictions$labels
saveRDS(rpe005_seurat_sct, "pan_rpe005_processed_SCTtransformed.rds")

Idents(rpe005_seurat_sct) = "SingleR"
rpe005_seurat_sct.monomac <- subset(rpe005_seurat_sct, idents = c("Monocyte","Macrophage"))
Idents(rpe005_seurat_sct.monomac) = "SingleR"
rpe005_seurat_sct.monomac <- subset(rpe005_seurat_sct.monomac, idents = c("Monocyte","Macrophage"))
table(rpe005_seurat_sct.monomac$SingleR)
# Macrophage   Monocyte 
#    605        495
dim(rpe005_seurat_sct.monomac)#[1] 16670  1100

# SCT normalize
rpe005_seurat_sct.monomac = Seurat_SCT_process(rpe005_seurat_sct.monomac, seurat_filter = TRUE,
                                               sample_id = "pan_rpe005", 
                                               study_name = "pan_et_al")

# SCT normalize
rpe005_seurat_sct.monomac = Seurat_SCT_process(rpe005_seurat_sct.monomac, seurat_filter = TRUE,
                                               sample_id = "pan_rpe005", 
                                               study_name = "pan_et_al")
dim(rpe005_seurat_sct.monomac)#[1]  14160  1100

colnames(rpe005_seurat_sct.monomac@meta.data)
metadata <- rpe005_seurat_sct.monomac@meta.data
metadata$barcode <- rownames(metadata)
#View(metadata)
metadata2 <- metadata[,c("orig.ident","Project","GSE_ID", "GSM_ID", "Tissue",  "Study_name","Species", "Age_group", "SingleR","nCount_RNA" , "nFeature_RNA","percent.mt" , "nCount_SCT", "nFeature_SCT", "scDblFinder.class" )]
#View(metadata2)
rpe005_seurat_sct.monomac@meta.data <- metadata2
Idents(rpe005_seurat_sct.monomac) = "SingleR"
DimPlot(rpe005_seurat_sct.monomac)

#14160 genes * 1100 cells
saveRDS(rpe005_seurat_sct.monomac, "pan_rpe005_processed_SCTtransformed_monomac.rds")

######################################################################################################
#                                                                                                    #
#       Add metadata &  Prefix_barcode & SingleR annotation & Subset monocytes and macrophages       #
#                                        rpe_006                                                     #
######################################################################################################
library(Seurat)
library(tidyverse)
setwd("/Users/lucia/Documents/34.immunometabolism/data/Pan_carotid_human")
rpe006_seurat_sct <- readRDS("pan_rpe006_processed_SCTtransformed.rds")

rpe006_seurat_sct$Project <- "Pan_rpe006_carotid"
rpe006_seurat_sct$GSE_ID <- "GSE155512"
rpe006_seurat_sct$GSM_ID <- "GSM4705591"
rpe006_seurat_sct$Tissue <- "carotid artery atherosclerotic core (AC)"
rpe006_seurat_sct$Study_name <- "Pan_et_al"
rpe006_seurat_sct$Species <- "human"
rpe006_seurat_sct$Age_group <- "old"
rpe006_seurat_sct$orig.ident <- "Pan_rpe006_carotid"
View(rpe006_seurat_sct@meta.data)

rpe006_seurat_sct <- RenameCells(rpe006_seurat_sct,add.cell.id = "Pan_rpe006")

library(SingleR)
library(celldex)
surveyReferences()
ref.data <- HumanPrimaryCellAtlasData()
#ensembl = TRUE
#ref.data
#table(ref.data$label.main)
#table(ref.data$label.fine)
#table(ref.data$label.ont)

sce <- LayerData(rpe006_seurat_sct, assay = "RNA", layer = "counts")
predictions <- SingleR(test=sce, assay.type.test=1, 
                       ref=ref.data, labels=ref.data$label.main)

table(predictions$labels)
plotScoreHeatmap(predictions)#picture size: 840*600
rpe006_seurat_sct$SingleR <- predictions$labels
saveRDS(rpe006_seurat_sct, "pan_rpe006_processed_SCTtransformed.rds")

Idents(rpe006_seurat_sct) = "SingleR"
rpe006_seurat_sct.monomac <- subset(rpe006_seurat_sct, idents = c("Monocyte","Macrophage"))
Idents(rpe006_seurat_sct.monomac) = "SingleR"
rpe006_seurat_sct.monomac <- subset(rpe006_seurat_sct.monomac, idents = c("Monocyte","Macrophage"))
table(rpe006_seurat_sct.monomac$SingleR)
# Macrophage   Monocyte 
#   391         335 
dim(rpe006_seurat_sct.monomac)#[1] 15445   726

# SCT normalize
rpe006_seurat_sct.monomac = Seurat_SCT_process(rpe006_seurat_sct.monomac, seurat_filter = TRUE,
                                               sample_id = "pan_rpe006", 
                                               study_name = "pan_et_al")

# SCT normalize
rpe006_seurat_sct.monomac = Seurat_SCT_process(rpe006_seurat_sct.monomac, seurat_filter = TRUE,
                                               sample_id = "pan_rpe006", 
                                               study_name = "pan_et_al")
dim(rpe006_seurat_sct.monomac)#[1] 12727   726

colnames(rpe006_seurat_sct.monomac@meta.data)
metadata <- rpe006_seurat_sct.monomac@meta.data
metadata$barcode <- rownames(metadata)
#View(metadata)
metadata2 <- metadata[,c("orig.ident","Project","GSE_ID", "GSM_ID", "Tissue",  "Study_name","Species", "Age_group", "SingleR","nCount_RNA" , "nFeature_RNA","percent.mt" , "nCount_SCT", "nFeature_SCT", "scDblFinder.class" )]
#View(metadata2)
rpe006_seurat_sct.monomac@meta.data <- metadata2
Idents(rpe006_seurat_sct.monomac) = "SingleR"
DimPlot(rpe006_seurat_sct.monomac)

#12727 genes * 726 cells
saveRDS(rpe006_seurat_sct.monomac, "pan_rpe006_processed_SCTtransformed_monomac.rds")
