#pre-process Normal_human Ascending aortic wall tissue | different age | scRNAseq data | SCTtransform 
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
setwd("/Users/lucia/Documents/34.immunometabolism/data/Pan_normal_human/")

##############################################################################################
# Read sample 6: GSM6696847	Normal human, replicate 6, scRNAseq Ascending aortic wall tissue #
##############################################################################################
ascending_aortic_n6 = Seurat::Read10X("GSE216860_RAW/GSM6696847/outs/filtered_feature_bc_matrix/")
# Raw, unfiltered matrix is 58825 genes x 10484 cells 
dim(ascending_aortic_n6)

###############################################################
# STEP 1: Pre-processing (Create Seurat obj containing doublets)
# Dims of seurat object: 22137 genes x 10422 cells
ascending_aortic_n6 = CreateSeuratObject(counts = ascending_aortic_n6, 
                                      project = "Pan_ascending_aortic_n6", 
                                      min.cells = 10, 
                                      min.features = 200)
dim(ascending_aortic_n6)

###################################
# STEP 2 before/after removing doublets

Seurat_SCT_process = function(seurat_obj, seurat_filter=FALSE, 
                              sample_id, study_name, artery, disease_status){
  seurat_obj$sample = sample_id
  seurat_obj$study = study_name
  # Check for mt percentage 
  seurat_obj[["percent.mt"]] = PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  
  # Filter low quality cells (start with parameters used in the paper)
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
options(future.globals.maxSize= 8912896000)
ascending_aortic_n6 = Seurat_SCT_process(ascending_aortic_n6, seurat_filter = FALSE,
                                       sample_id = "Pan_ascending_aortic_n6", 
                                       study_name = "Pan_et_al_normal")

# Clean counts matrix (without doublets) from ambient RNA using the decontX wrapper
decontX_remove = function(seurat_obj) {
  decontX_results = celda::decontX(seurat_obj@assays$RNA$counts)
  decontaminated_matrix = decontX_results$decontXcounts
  seurat_obj@assays$RNA$counts = decontaminated_matrix
  return(seurat_obj)
}

ascending_aortic_n6 =  decontX_remove(ascending_aortic_n6)
head(ascending_aortic_n6@assays$RNA$counts)

# Finding the right resolution is quite important for doublet detection and removal
ascending_aortic_n6 =  FindClusters(ascending_aortic_n6, resolution = 0.2)   

# Visualize clusters
n6_before_QC = DimPlot(ascending_aortic_n6, reduction = "umap", label = TRUE) +
  ggtitle("Before removing doublets")

##################################
# STEP 3
# Detect doublets with scDblFinder
# Called 363 doublets

sce <- as.SingleCellExperiment(ascending_aortic_n6)
sce <- scDblFinder(sce,clusters = "SCT_snn_res.0.2")
table(sce$scDblFinder.class) #10059 singlet, 363 doublet
metadata <-  sce@colData %>% as.data.frame()
ascending_aortic_n6@meta.data <- metadata

# First filtering step: Use the vector of consensus doublet IDs output below
ascending_aortic_n6 <- subset(ascending_aortic_n6, subset = scDblFinder.class == "singlet")
dim(ascending_aortic_n6)#10059

# New dims after removing doublets: 22137 genes x 10059 cells
###################################
# STEP 4 after removing doublets
# SCT normalize
ascending_aortic_n6 = Seurat_SCT_process(ascending_aortic_n6, seurat_filter = TRUE,
                                         sample_id = "Pan_ascending_aortic_n6", 
                                         study_name = "Pan_et_al_normal") #22136 genes * 10059 cells
dim(ascending_aortic_n6)

n6_after_QC = DimPlot(ascending_aortic_n6, reduction = "umap", label = TRUE) +
  ggtitle("After removing doublets + decontX")

FeaturePlot(ascending_aortic_n6, features = c("MYH11", "CNN1", "TNFRSF11B", "ACTA2", 
                                            "CRTAC1", "LUM"), order = TRUE,  pt.size = 0.1)
FeaturePlot(ascending_aortic_n6, features = c("CD68", "CD14", "CD163", "APOE", 
                                            "S100A8", "IL1B"), order = TRUE,  pt.size = 0.1)
FeaturePlot(ascending_aortic_n6, features = c("MYH11", "CNN1", "FHL5", "FOXC1"), order = TRUE,  pt.size = 0.1)
# Save RDS object
# NUmber of cells after processing: 10512
saveRDS(ascending_aortic_n6, "Pan_ascending_aortic_n6_processed_SCTtransformed.rds")

######################################################################################################
#                                                                                                    #
#       Add metadata &  Prefix_barcode & SingleR annotation & Subset monocytes and macrophages       #
#                                                                                                    #
######################################################################################################
library(Seurat)
library(tidyverse)
setwd("/Users/lucia/Documents/34.immunometabolism/data/Pan_normal_human/")

ascending_aortic_n6 <- readRDS("Pan_ascending_aortic_n6_processed_SCTtransformed.rds")
ascending_aortic_n6$Project <- "Zhang_ascending_aortic_n6"
ascending_aortic_n6$GSE_ID <- "GSE216860"
ascending_aortic_n6$GSM_ID <- "GSM6696847"
ascending_aortic_n6$Tissue <- "ascending aortic wall tissue"
ascending_aortic_n6$Study_name <- "Zhang_et_al_normal"
ascending_aortic_n6$Species <- "human"
ascending_aortic_n6$Age_group <- "old"
ascending_aortic_n6$orig.ident <- "Zhang_ascending_aortic_n6"
View(ascending_aortic_n6@meta.data)

ascending_aortic_n6 <- RenameCells(ascending_aortic_n6,add.cell.id = "GSM6696847")

library(SingleR)
library(celldex)
surveyReferences()
ref.data <- HumanPrimaryCellAtlasData()
#ensembl = TRUE
#ref.data
#table(ref.data$label.main)
#table(ref.data$label.fine)
#table(ref.data$label.ont)

sce <- LayerData(ascending_aortic_n6, assay = "RNA", layer = "counts")
predictions <- SingleR(test=sce, assay.type.test=1, 
                       ref=ref.data, labels=ref.data$label.main)
table(predictions$labels)
plotScoreHeatmap(predictions)#picture size: 840*600
ascending_aortic_n6$SingleR <- predictions$labels
table(ascending_aortic_n6$SingleR)
saveRDS(ascending_aortic_n6, "Pan_ascending_aortic_n6_processed_SCTtransformed.rds")

Idents(ascending_aortic_n6) = "SingleR"
ascending_aortic_n6.monomac <- subset(ascending_aortic_n6, idents = c("Monocyte","Macrophage"))
Idents(ascending_aortic_n6.monomac) = "SingleR"
ascending_aortic_n6.monomac <- subset(ascending_aortic_n6.monomac, idents = c("Monocyte","Macrophage"))
table(ascending_aortic_n6.monomac$SingleR)
# Macrophage   Monocyte 
#   2373        693 

# SCT normalize
ascending_aortic_n6.monomac = Seurat_SCT_process(ascending_aortic_n6.monomac, seurat_filter = TRUE,
                                         sample_id = "Pan_ascending_aortic_n6", 
                                         study_name = "Pan_et_al_normal")
# SCT normalize
ascending_aortic_n6.monomac = Seurat_SCT_process(ascending_aortic_n6.monomac, seurat_filter = TRUE,
                                                 sample_id = "Pan_ascending_aortic_n6", 
                                                 study_name = "Pan_et_al_normal")
dim(ascending_aortic_n6.monomac)#[1] 19480  3066

colnames(ascending_aortic_n6.monomac@meta.data)
metadata <- ascending_aortic_n6.monomac@meta.data
metadata$barcode <- rownames(metadata)
#View(metadata)
metadata2 <- metadata[,c("orig.ident","Project","GSE_ID", "GSM_ID", "Tissue",  "Study_name","Species", "Age_group", "SingleR","nCount_RNA" , "nFeature_RNA","percent.mt" , "nCount_SCT", "nFeature_SCT", "scDblFinder.class" )]
#View(metadata2)
ascending_aortic_n6.monomac@meta.data <- metadata2
DimPlot(ascending_aortic_n6.monomac,group.by = "SingleR")

#11716  genes * 440 cells
saveRDS(ascending_aortic_n6.monomac, "Pan_ascending_aortic_n6_processed_SCTtransformed_monomac.rds")
