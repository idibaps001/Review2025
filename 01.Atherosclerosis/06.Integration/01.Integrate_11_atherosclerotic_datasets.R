#Integrate 11 datasets about atherosclerotic
library(Seurat)
library(tidyverse)
library(celda)
library(sctransform)
library(data.table)
library(cluster)
#remotes::install_github('satijalab/seurat-wrappers')
library(SeuratWrappers)
library(patchwork)
options(future.globals.maxSize = 3e10)
setwd("/Users/lucia/Documents/34.immunometabolism/data")

#Load 11 datasets
coronary.monomac <- readRDS("Wirka_coronary_human/Wirka_processed_SCTtransformed_monomac.rds")
rpe004_seurat_sct.monomac <- readRDS("Pan_carotid_human/pan_rpe004_processed_SCTtransformed_monomac.rds")
rpe005_seurat_sct.monomac <- readRDS("Pan_carotid_human/pan_rpe005_processed_SCTtransformed_monomac.rds")
rpe006_seurat_sct.monomac <- readRDS("Pan_carotid_human/pan_rpe006_processed_SCTtransformed_monomac.rds")
ac_p1_seurat_sct.monomac <- readRDS("Alsaigh_carotid_human/alsaigh_ac_p1_processed_SCTtransformed_monomac.rds")
pa_p1_seurat_sct.monomac <- readRDS("Alsaigh_carotid_human/alsaigh_pa_p1_processed_SCTtransformed_monomac.rds")
ac_p2_seurat_sct.monomac <- readRDS("Alsaigh_carotid_human/alsaigh_ac_p2_processed_SCTtransformed_monomac.rds")
pa_p2_seurat_sct.monomac <- readRDS("Alsaigh_carotid_human/alsaigh_pa_p2_processed_SCTtransformed_monomac.rds")
ac_p3_seurat_sct.monomac <- readRDS("Alsaigh_carotid_human/alsaigh_ac_p3_processed_SCTtransformed_monomac.rds")
pa_p3_seurat_sct.monomac <- readRDS("Alsaigh_carotid_human/alsaigh_pa_p3_processed_SCTtransformed_monomac.rds")
ascending_aortic_n6.monomac <- readRDS("Pan_normal_human/Pan_ascending_aortic_n6_processed_SCTtransformed_monomac.rds")

#Merge datasets
merged_obj <- merge(x = coronary.monomac, y = list(rpe004_seurat_sct.monomac,rpe005_seurat_sct.monomac,
                                                   rpe006_seurat_sct.monomac,ac_p1_seurat_sct.monomac,
                                                   pa_p1_seurat_sct.monomac, ac_p2_seurat_sct.monomac,
                                                   pa_p2_seurat_sct.monomac, ac_p3_seurat_sct.monomac, 
                                                   pa_p3_seurat_sct.monomac, ascending_aortic_n6.monomac      ))

merged_obj <- SCTransform(merged_obj, method = "glmGamPoi", 
                          vars.to.regress = "percent.mt")
merged_obj <- RunPCA(merged_obj, npcs = 30, verbose = F)

#Do integration for SCTtransformed datasets
Integration_atherosclerotic <- IntegrateLayers(
  object = merged_obj,
  method = RPCAIntegration,
  normalization.method = "SCT",
  verbose = F
)

Integration_atherosclerotic <- FindNeighbors(Integration_atherosclerotic, dims = 1:30, reduction = "integrated.dr")
Integration_atherosclerotic <- FindClusters(Integration_atherosclerotic, resolution = 2)
Integration_atherosclerotic <- RunUMAP(Integration_atherosclerotic, dims = 1:30, reduction = "integrated.dr")
# Visualization
DimPlot(Integration_atherosclerotic, reduction = "umap", group.by = "Tissue")

#Integration_atherosclerotic <- JoinLayers(Integration_atherosclerotic) 运行失败

dim(Integration_atherosclerotic)
#[1] 22024 16668

saveRDS(Integration_atherosclerotic,"Integration_atherosclerotic.rds"  )

#################################### Do Violin Plots ######################################################
Integration_atherosclerotic <- read_rds("Integration_atherosclerotic.rds")

Integration_atherosclerotic@meta.data$celltype2 <- "Others"
length(which(Integration_atherosclerotic@meta.data$SingleR_map_to_TAM == "16_ECMHomeoMac"    ))#2310
table( Integration_atherosclerotic@meta.data$SingleR_map_to_TAM   )#2310
Integration_atherosclerotic@meta.data$celltype2[which(   Integration_atherosclerotic@meta.data$SingleR_map_to_TAM == "16_ECMHomeoMac"  )] <- "16_ECMHomeoMac" 
table(Integration_atherosclerotic@meta.data$celltype2)
#16_ECMHomeoMac         Others 
# 2310          14358  
table(Integration_atherosclerotic@meta.data$celltype2, Integration_atherosclerotic@meta.data$Tissue)
#               ascending aortic wall tissue  carotid artery atherosclerotic core (AC) carotid artery proximal adjacent (PA) coronary atherosclerotic core (AC)
#16_ECMHomeoMac                            4                                     2302                                     0                                  4
#Others                                 3062                                     8491                                   989                               1816

Idents(Integration_atherosclerotic ) <- "Tissue"
p=VlnPlot(Integration_atherosclerotic ,idents = "carotid artery atherosclerotic core (AC)" ,   features = c("FABP5","CD36","CD68","PLIN2","LIPA","APOE","LPL","LAMP1","TREM2","ABCA1","MARCO" ),group.by = "celltype2", cols = c("#AD867E","#a47053","#7b5965","#95594c","#533068","#66597b","#50609f","#4e6980","#678171","#7bab77","#979c70"),raster=FALSE, stack = TRUE, sort = FALSE, flip = TRUE)+ 
  theme(legend.position = "none") + 
  ggtitle("Atherosclerotic core") + 
  theme(plot.margin = unit(c(0.1,0.1,0.1,1), "inches"))+
  scale_x_discrete(labels = c("16_ECMHomeoMac" = "Cluster16", "Others" = "Others"))

setwd("/Users/lucia/Documents/34.immunometabolism")
ggsave(p,filename = "plots/Lipid_markers_of_16_ECMHomeoMac_in_Carotid_Artery_AC.jpeg",height = 6,width = 4) 


p=VlnPlot(Integration_atherosclerotic ,idents = "carotid artery atherosclerotic core (AC)" ,   features = c("LGALS1", "LGALS3", "SIRPA", "MMP9","MRC1", "HLA-DPA1", "HLA-DRB1" ), cols = c("#95594c","#533068","#66597b","#50609f","#4e6980","#678171","#7bab77"), group.by = "celltype2",raster=FALSE, stack = TRUE, sort = FALSE, flip = TRUE)+ 
  theme(legend.position = "none") + 
  ggtitle("Atherosclerotic core") + 
  theme(plot.margin = unit(c(0.1,0.1,0.1,1), "inches"))+
  scale_x_discrete(labels = c("16_ECMHomeoMac" = "Cluster16", "Others" = "Others"))

setwd("/Users/lucia/Documents/34.immunometabolism")
ggsave(p,filename = "plots/Immuno_markers_of_16_ECMHomeoMac_in_Carotid_Artery_AC.jpeg",height = 4,width = 4) 
