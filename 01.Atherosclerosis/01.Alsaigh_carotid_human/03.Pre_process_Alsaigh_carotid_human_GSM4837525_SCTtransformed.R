#pre-process Alsaigh_carotid_human scRNAseq data | SCTtransform 
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
setwd("/Users/lucia/Documents/34.immunometabolism/data/Alsaigh_carotid_human/")

###############################################################################
# Read sample 3: GSM4837525 Patient2 Carotid Artery Atherosclerotic Core (AC) #
###############################################################################
carotid_ac_p2 = Seurat::Read10X("GSE159677_RAW/GSM4837525/outs/filtered_feature_bc_matrix/")
# Raw, unfiltered matrix is 33538 genes x 15960 cells
dim(carotid_ac_p2)

###############################################################
# STEP 1: Pre-processing (Create Seurat obj containing doublets)
# Dims of seurat object: 17901 genes x 15923 cells
ac_p2_seurat_sct = CreateSeuratObject(counts = carotid_ac_p2, 
                                      project = "alsaigh_p2_carotid_AC", 
                                      min.cells = 10, 
                                      min.features = 200)
dim(ac_p2_seurat_sct)

###################################
# STEP 2 before/after removing doublets

Seurat_SCT_process = function(seurat_obj, seurat_filter=FALSE, 
                              sample_id, study_name, artery, disease_status){
  seurat_obj$sample = sample_id
  seurat_obj$study = study_name
  # Check for mt percentage
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
options(future.globals.maxSize= 891289600)
ac_p2_seurat_sct = Seurat_SCT_process(ac_p2_seurat_sct, seurat_filter = FALSE,
                                       sample_id = "alsaigh_ac_p2", 
                                       study_name = "alsaigh_et_al")

# Clean counts matrix (without doublets) from ambient RNA using the decontX wrapper
decontX_remove = function(seurat_obj) {
  decontX_results = celda::decontX(seurat_obj@assays$RNA$counts)
  decontaminated_matrix = decontX_results$decontXcounts
  seurat_obj@assays$RNA$counts = decontaminated_matrix
  return(seurat_obj)
}

ac_p2_seurat_sct =  decontX_remove(ac_p2_seurat_sct)
head(ac_p2_seurat_sct@assays$RNA$counts)

# Finding the right resolution is quite important for doublet detection and removal
ac_p2_seurat_sct =  FindClusters(ac_p2_seurat_sct, resolution = 0.1)   

# Visualize clusters
ac_p2_before_QC = DimPlot(ac_p2_seurat_sct, reduction = "umap", label = TRUE) +
  ggtitle("Before removing doublets")

##################################
# STEP 3
# Detect doublets with scDblFinder
# Called 909 doublets

sce <- as.SingleCellExperiment(ac_p2_seurat_sct)
sce <- scDblFinder(sce,clusters = "SCT_snn_res.0.1")
table(sce$scDblFinder.class) #15014 singlet, 909 doublet
metadata <-  sce@colData %>% as.data.frame()
ac_p2_seurat_sct@meta.data <- metadata

# First filtering step: Use the vector of consensus doublet IDs output below
ac_p2_seurat_sct <- subset(ac_p2_seurat_sct, subset = scDblFinder.class == "singlet")
dim(ac_p2_seurat_sct)#15014

# New dims after removing doublets: 17901 genes x 15014 cells
###################################
# STEP 4 after removing doublets

# SCT normalize
ac_p2_seurat_sct = Seurat_SCT_process(ac_p2_seurat_sct, seurat_filter = TRUE,
                                       sample_id = "alsaigh_ac_p2", 
                                       study_name = "alsaigh_et_al") #17901 genes x 14829 cells
dim(ac_p2_seurat_sct)

ac_p2_after_QC = DimPlot(ac_p2_seurat_sct, reduction = "umap", label = TRUE) +
  ggtitle("After removing doublets + decontX")



FeaturePlot(ac_p2_seurat_sct, features = c("MYH11", "CNN1", "TNFRSF11B", "ACTA2", 
                                            "CRTAC1", "LUM"), order = TRUE,  pt.size = 0.1)
FeaturePlot(ac_p2_seurat_sct, features = c("CD68", "CD14", "CD163", "APOE", 
                                            "S100A8", "IL1B"), order = TRUE,  pt.size = 0.1)
FeaturePlot(ac_p2_seurat_sct, features = c("MYH11", "CNN1", "FHL5", "FOXC1"), order = TRUE,  pt.size = 0.1)
# Save RDS object
# NUmber of cells after processing: 14829
saveRDS(ac_p2_seurat_sct, "alsaigh_ac_p2_processed_SCTtransformed.rds")

######################################################################################################
#                                                                                                    #
#       Add metadata &  Prefix_barcode & SingleR annotation & Subset monocytes and macrophages       #
#                                                                                                    #
######################################################################################################
library(Seurat)
library(tidyverse)
setwd("/Users/lucia/Documents/34.immunometabolism/data/Alsaigh_carotid_human/")

ac_p2_seurat_sct <- readRDS("alsaigh_ac_p2_processed_SCTtransformed.rds")
ac_p2_seurat_sct$Project <- "Alsaigh_p2_carotid_AC"
ac_p2_seurat_sct$GSE_ID <- "GSE159677"
ac_p2_seurat_sct$GSM_ID <- "GSM4837525"
ac_p2_seurat_sct$Tissue <- "carotid artery atherosclerotic core (AC)"
ac_p2_seurat_sct$Study_name <- "Alsaigh_et_al"
ac_p2_seurat_sct$Species <- "human"
ac_p2_seurat_sct$Age_group <- "old"
ac_p2_seurat_sct$orig.ident <- "Alsaigh_p2_carotid_AC"
View(ac_p2_seurat_sct@meta.data)

ac_p2_seurat_sct <- RenameCells(ac_p2_seurat_sct,add.cell.id = "Alsaigh_acp2")

library(SingleR)
library(celldex)
surveyReferences()
ref.data <- HumanPrimaryCellAtlasData()
#ensembl = TRUE
#ref.data
#table(ref.data$label.main)
#table(ref.data$label.fine)
#table(ref.data$label.ont)

sce <- LayerData(ac_p2_seurat_sct, assay = "RNA", layer = "counts")
predictions <- SingleR(test=sce, assay.type.test=1, 
                       ref=ref.data, labels=ref.data$label.main)

table(predictions$labels)
plotScoreHeatmap(predictions)#picture size: 840*600
ac_p2_seurat_sct$SingleR <- predictions$labels
table(ac_p2_seurat_sct$SingleR)
saveRDS(ac_p2_seurat_sct, "alsaigh_ac_p2_processed_SCTtransformed.rds")

Idents(ac_p2_seurat_sct) = "SingleR"
ac_p2_seurat_sct.monomac <- subset(ac_p2_seurat_sct, idents = c("Monocyte","Macrophage"))
Idents(ac_p2_seurat_sct.monomac) = "SingleR"
ac_p2_seurat_sct.monomac <- subset(ac_p2_seurat_sct.monomac, idents = c("Monocyte","Macrophage"))
table(ac_p2_seurat_sct.monomac$SingleR)
# Macrophage   Monocyte 
#    817       1027  


# SCT normalize
ac_p2_seurat_sct.monomac = Seurat_SCT_process(ac_p2_seurat_sct.monomac, seurat_filter = TRUE,
                                              sample_id = "alsaigh_ac_p2", 
                                              study_name = "alsaigh_et_al") 

# SCT normalize
ac_p2_seurat_sct.monomac = Seurat_SCT_process(ac_p2_seurat_sct.monomac, seurat_filter = TRUE,
                                              sample_id = "alsaigh_ac_p2", 
                                              study_name = "alsaigh_et_al") 

dim(ac_p2_seurat_sct.monomac)#[1] 14419  1844

colnames(ac_p2_seurat_sct.monomac@meta.data)
metadata <- ac_p2_seurat_sct.monomac@meta.data
metadata$barcode <- rownames(metadata)
#View(metadata)
metadata2 <- metadata[,c("orig.ident","Project","GSE_ID", "GSM_ID", "Tissue",  "Study_name","Species", "Age_group", "SingleR","nCount_RNA" , "nFeature_RNA","percent.mt" , "nCount_SCT", "nFeature_SCT", "scDblFinder.class" )]
#View(metadata2)
ac_p2_seurat_sct.monomac@meta.data <- metadata2
DimPlot(ac_p2_seurat_sct.monomac,group.by = "SingleR")

#14419 genes * 1844 cells
saveRDS(ac_p2_seurat_sct.monomac, "alsaigh_ac_p2_processed_SCTtransformed_monomac.rds")
