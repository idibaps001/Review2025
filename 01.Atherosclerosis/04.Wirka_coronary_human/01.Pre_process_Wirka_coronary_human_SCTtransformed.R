#pre-process Wirka_coronary_human scRNAseq data | SCTtransform 
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
# PROCESS HUMAN CORONARY MATRICES USING SCTransform
setwd("/Users/lucia/Documents/34.immunometabolism/data/Wirka_coronary_human/")
coronary <- read.table(file="GSE131778_human_coronary_scRNAseq_wirka_et_al_GEO.txt.gz",sep="\t", header = T, row.names=1)
coronary = as.sparse(coronary)

##################################
# STEP 1: Pre-processing (Create Seurat obj containing doublets)
# Dims of seurat object: 18194 genes x 11756 cells
coronary  <- CreateSeuratObject(counts=coronary , project = "wirka.human", min.cells = 10, min.features = 200)
dim(coronary)#[1] 18194 11756

###################################
# STEP 2 before/after removing doublets
Seurat_SCT_process = function(seurat_obj, seurat_filter=FALSE, 
                              sample_id, study_name, artery, disease_status){
  seurat_obj$sample = sample_id
  seurat_obj$study = study_name
  # Check for mt percentage and other quality metrics 
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
coronary = Seurat_SCT_process(coronary, seurat_filter = FALSE,
                                       sample_id = "coronary", 
                                       study_name = "Wirka_et_al")

# Clean counts matrix (without doublets) from ambient RNA using the decontX wrapper
decontX_remove = function(seurat_obj) {
  decontX_results = celda::decontX(seurat_obj@assays$RNA$counts)
  decontaminated_matrix = decontX_results$decontXcounts
  seurat_obj@assays$RNA$counts = decontaminated_matrix
  return(seurat_obj)
}

coronary =  decontX_remove(coronary)
head(coronary@assays$RNA$counts)

# Finding the right resolution is quite important for doublet detection and removal
coronary =  FindClusters(coronary, resolution = 0.1)   

# Visualize clusters
p1_before_QC = DimPlot(coronary, reduction = "umap", label = TRUE) +
  ggtitle("Before removing doublets")

##################################
# STEP 3
# Detect doublets with scDblFinder
# Called 469  doublets

sce <- as.SingleCellExperiment(coronary)
sce <- scDblFinder(sce,clusters = "SCT_snn_res.0.1")
table(sce$scDblFinder.class) #11287 singlet, 469  doublet
metadata <-  sce@colData %>% as.data.frame()
coronary@meta.data <- metadata

# First filtering step: Use the vector of consensus doublet IDs output below
coronary <- subset(coronary, subset = scDblFinder.class == "singlet")
dim(coronary)#11287

# New dims after removing doublets: 18194 genes x 11287 cells

###################################
# STEP 4 after removing doublets

# SCT normalize
coronary = Seurat_SCT_process(coronary, seurat_filter = TRUE,
                              sample_id = "coronary", 
                              study_name = "Wirka_et_al") #18194 genes * 11287 cells

ncol(coronary)
nrow(coronary)
p1_after_QC = DimPlot(coronary, reduction = "umap", label = TRUE) +
  ggtitle("After removing doublets + decontX")

FeaturePlot(coronary, features = c("MYH11", "CNN1", "TNFRSF11B", "ACTA2", 
                                            "CRTAC1", "LUM"), order = TRUE,  pt.size = 0.1)
FeaturePlot(coronary, features = c("CD68", "CD14", "CD163", "APOE", 
                                            "S100A8", "IL1B"), order = TRUE,  pt.size = 0.1)
FeaturePlot(coronary, features = c("MYH11", "CNN1", "FHL5", "FOXC1"), order = TRUE,  pt.size = 0.1)
# Save RDS object
# NUmber of cells after processing: 11287
saveRDS(coronary, "Wirka_processed_SCTtransformed.rds")

######################################################################################################
#                                                                                                    #
#       Add metadata &  Prefix_barcode & SingleR annotation & Subset monocytes and macrophages       #
#                                                                                                    #
######################################################################################################
library(Seurat)
library(tidyverse)
setwd("/Users/lucia/Documents/34.immunometabolism/data/Wirka_coronary_human/")

coronary <- readRDS("Wirka_processed_SCTtransformed.rds")
coronary$Project <- "Wirka_coronary"
coronary$GSE_ID <- "GSE131778"
coronary$GSM_ID <- "GSM3819856_63"
coronary$Tissue <- "coronary atherosclerotic core (AC)"
coronary$Study_name <- "Wirka_et_al"
coronary$Species <- "human"
coronary$Age_group <- "old"
coronary$orig.ident <- "Wirka_coronary"
View(coronary@meta.data)

coronary <- RenameCells(coronary,add.cell.id = "Wirka")

library(SingleR)
library(celldex)
surveyReferences()
ref.data <- HumanPrimaryCellAtlasData()
#ensembl = TRUE
#ref.data
#table(ref.data$label.main)
#table(ref.data$label.fine)
#table(ref.data$label.ont)

sce <- LayerData(coronary, assay = "RNA", layer = "counts")
predictions <- SingleR(test=sce, assay.type.test=1, 
                       ref=ref.data, labels=ref.data$label.main)

table(predictions$labels)
plotScoreHeatmap(predictions)#picture size: 840*600
coronary$SingleR <- predictions$labels
saveRDS(coronary, "Wirka_processed_SCTtransformed.rds")

Idents(coronary) = "SingleR"
coronary.monomac <- subset(coronary, idents = c("Monocyte","Macrophage"))
Idents(coronary.monomac) = "SingleR"
coronary.monomac <- subset(coronary.monomac, idents = c("Monocyte","Macrophage"))
table(coronary.monomac$SingleR)
# Macrophage   Monocyte 
#   1086        734 

# SCT normalize
coronary.monomac = Seurat_SCT_process(coronary.monomac, seurat_filter = TRUE,
                              sample_id = "coronary", 
                              study_name = "Wirka_et_al") #14398 genes * 1820 cells



# SCT normalize
coronary.monomac = Seurat_SCT_process(coronary.monomac, seurat_filter = TRUE,
                                      sample_id = "coronary", 
                                      study_name = "Wirka_et_al") #14398 genes * 1820 cells

dim(coronary.monomac)
colnames(coronary.monomac@meta.data)
metadata <- coronary.monomac@meta.data
metadata$barcode <- rownames(metadata)
#View(metadata)
metadata2 <- metadata[,c("orig.ident","Project","GSE_ID", "GSM_ID", "Tissue",  "Study_name","Species", "Age_group", "SingleR","nCount_RNA" , "nFeature_RNA","percent.mt" , "nCount_SCT", "nFeature_SCT", "scDblFinder.class" )]
#View(metadata2)
coronary.monomac@meta.data <- metadata2
Idents(coronary.monomac) = "SingleR"
DimPlot(coronary.monomac)

#14398 genes * 1820 cells
saveRDS(coronary.monomac, "Wirka_processed_SCTtransformed_monomac.rds")
