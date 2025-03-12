#Map atherosclerotic dataset to the pan-cancer dataset
library(Seurat)
library(tidyverse)
library(celda)
library(sctransform)
library(data.table)
library(cluster)
library(SeuratWrappers)
library(patchwork)
library(SingleR)
options(future.globals.maxSize = 3e10)

#data.ref
Load the data follow "Immunometabolism.R"
Idents(mac.atlas.zenodo) <- "short.label"
#mac.atlas.zenodo <- subset(mac.atlas.zenodo, idents = "23_NA", invert = TRUE)
dim(mac.atlas.zenodo) #[1]  64082 363315

#DefaultAssay(mac.atlas.zenodo) <- "RNA"    #[1]   2000 363315
#colnames(GetAssayData(mac.atlas.zenodo, assay = "RNA", slot = "counts"))
#rownames(mac.atlas.zenodo@meta.data)

#The ref dataset is too large, I have to do a sampling. I keep 20000 cells.
set.seed(1111)
# Downsample the number of cells per identity class
Idents(mac.atlas.zenodo) <- "short.label"
#ref.sub <- mac.atlas.zenodo[, sample(colnames(mac.atlas.zenodo), size =20000, replace=F)]
#ref.sub <- subset(ref.sub, idents = "23_NA", invert = TRUE)

###########################################################################################################
#  No.01                                                                                                  #
#               Map Wirka_processed_SCTtransformed_monomac.rds to mac.atlas.zenodo                        #
#                                                                                                         #
###########################################################################################################

#data.query
setwd("/Users/lucia/Documents/34.immunometabolism/data")
coronary.monomac <- readRDS("Wirka_coronary_human/Wirka_processed_SCTtransformed_monomac.rds")

#setwd("/Users/lucia/Documents/34.immunometabolism/data")
#Integration_atherosclerotic <- readRDS("Integration_atherosclerotic.RDS")#这数据集太大了，1.5小时都没注释上

#********map
#https://github.com/SingleR-inc/SingleR/issues/98  Do not use SCTtransformed data
#https://bioconductor.org/books/release/SingleRBook/classic-mode.html#choices-of-assay-data

#sce.test <- LayerData(Integration_atherosclerotic, assay = "RNA", layer = "counts")
sce.test <- LayerData(coronary.monomac, assay = "RNA", layer = "counts")
sce.ref <- LayerData(mac.atlas.zenodo, assay = "RNA", layer = "data")
label.ref <- mac.atlas.zenodo$short.label
pred <- SingleR(test = sce.test, ref = sce.ref, labels = label.ref, de.method = "wilcox")

#******************test dataset 1820 cells, ref 20000 cells,need 15 min calculation
table(pred$labels)
#0_AlvMac     1_MetM2Mac   10_InflamMac  11_MetalloMac      12_MBMMac  13_CalciumMac    14_ProliMac      15_LYZMac 
#    1            632            169              3             39              4             33              6 
#16_ECMHomeoMac     17_IFNMac3      18_ECMMac   19_ClassMono        2_C3Mac     21_HemeMac     22_IFNMac4      3_ICIMac1 
#        4              8              1            182            425             13              4             49 
#4_ICIMac2  6_SPP1AREGMac       7_IFNMac     9_AngioMac 
#      101              1             71             74 


#******************test dataset 1820 cells, ref 363315 cells,need 2 hours calculation
table(pred$labels)

#0_AlvMac     1_MetM2Mac   10_InflamMac  11_MetalloMac      12_MBMMac  13_CalciumMac    14_ProliMac      15_LYZMac 
#    3            613            147              4              8              7             40             20 
#16_ECMHomeoMac     17_IFNMac3      18_ECMMac   19_ClassMono        2_C3Mac       20_TDoub     21_HemeMac     22_IFNMac4 
#    4              9                  2            181             404              2             56             18 
#23_NA      3_ICIMac1      4_ICIMac2  6_SPP1AREGMac       7_IFNMac     9_AngioMac 
#   9             61             94              2             66             70 
 
coronary.monomac$SingleR_map_to_TAM <- pred$labels
setwd("/Users/lucia/Documents/34.immunometabolism/data/Wirka_coronary_human/")
saveRDS(coronary.monomac, "Wirka_processed_SCTtransformed_monomac.rds")
    
###########################################################################################################
#  No.02                                                                                                  #
#               Map pan_rpe004_processed_SCTtransformed_monomac.rds to mac.atlas.zenodo                   #
#                                                                                                         #
###########################################################################################################

#data.query
setwd("/Users/lucia/Documents/34.immunometabolism/data/Pan_carotid_human/")
rpe004_seurat_sct.monomac <- readRDS("pan_rpe004_processed_SCTtransformed_monomac.rds")


#********map
sce.test <- LayerData(rpe004_seurat_sct.monomac, assay = "RNA", layer = "counts")
sce.ref <- LayerData(mac.atlas.zenodo, assay = "RNA", layer = "data")
label.ref <- mac.atlas.zenodo$short.label


library(BiocParallel)
pred <- SingleR(test = sce.test, ref = sce.ref, labels = label.ref, de.method = "wilcox", BPPARAM=MulticoreParam(8))

#******************test dataset 213 cells, ref 363315 cells,need 12:20- before 14:50  calculation
table(pred$labels)
# 1_MetM2Mac    14_ProliMac      15_LYZMac 16_ECMHomeoMac     17_IFNMac3   19_ClassMono        2_C3Mac       20_TDoub 
#     173              3              3              1              1             12              7              1 
# 23_NA      3_ICIMac1      4_ICIMac2       7_IFNMac     9_AngioMac 
#  1              3              3              3              2 
rpe004_seurat_sct.monomac$SingleR_map_to_TAM <- pred$labels
saveRDS(rpe004_seurat_sct.monomac, "pan_rpe004_processed_SCTtransformed_monomac.rds")

###########################################################################################################
#   No.03                                                                                                 #
#               Map pan_rpe005_processed_SCTtransformed_monomac.rds to mac.atlas.zenodo                   #
#                                                                                                         #
###########################################################################################################

#data.query
setwd("/Users/lucia/Documents/34.immunometabolism/data/Pan_carotid_human/")
rpe005_seurat_sct.monomac <- readRDS("pan_rpe005_processed_SCTtransformed_monomac.rds")

#********map
sce.test <- LayerData(rpe005_seurat_sct.monomac, assay = "RNA", layer = "counts")
sce.ref <- LayerData(mac.atlas.zenodo, assay = "RNA", layer = "data")
label.ref <- mac.atlas.zenodo$short.label

pred <- SingleR(test = sce.test, ref = sce.ref, labels = label.ref, de.method = "wilcox")

#******************test dataset 1100  cells, ref 363315 cells,need 102min calculation
table(pred$labels)
#1_MetM2Mac   10_InflamMac      12_MBMMac  13_CalciumMac    14_ProliMac      15_LYZMac 16_ECMHomeoMac     17_IFNMac3 
#   644             20              2              4              5             32              2              8 
#19_ClassMono        2_C3Mac       20_TDoub     22_IFNMac4          23_NA      3_ICIMac1      4_ICIMac2       7_IFNMac 
#        62            174              4             15              3              7             34              7 
#8_IFNGMac     9_AngioMac 
#    1             76 

rpe005_seurat_sct.monomac$SingleR_map_to_TAM <- pred$labels
saveRDS(rpe005_seurat_sct.monomac, "pan_rpe005_processed_SCTtransformed_monomac.rds")

###########################################################################################################
#   No.04                                                                                                 #
#               Map pan_rpe006_processed_SCTtransformed_monomac.rds to mac.atlas.zenodo                   #
#                                                                                                         #
###########################################################################################################

#data.query
setwd("/Users/lucia/Documents/34.immunometabolism/data/Pan_carotid_human/")
rpe006_seurat_sct.monomac <- readRDS("pan_rpe006_processed_SCTtransformed_monomac.rds")


#********map
sce.test <- LayerData(rpe006_seurat_sct.monomac, assay = "RNA", layer = "counts")
sce.ref <- LayerData(mac.atlas.zenodo, assay = "RNA", layer = "data")
label.ref <- mac.atlas.zenodo$short.label

pred <- SingleR(test = sce.test, ref = sce.ref, labels = label.ref, de.method = "wilcox")

#******************test dataset 726  cells, ref 363315 cells,need 18:22 -  calculation
table(pred$labels)
# 1_MetM2Mac   10_InflamMac      12_MBMMac  13_CalciumMac    14_ProliMac      15_LYZMac 16_ECMHomeoMac     17_IFNMac3 
#    345            124              4              1              7              4              4              2 
# 18_ECMMac   19_ClassMono        2_C3Mac       20_TDoub     21_HemeMac     22_IFNMac4      3_ICIMac1      4_ICIMac2 
#      1             26            148              2              7              8              3              7 
# 7_IFNMac     9_AngioMac 
#   27              6 
rpe006_seurat_sct.monomac$SingleR_map_to_TAM <- pred$labels
saveRDS(rpe006_seurat_sct.monomac, "pan_rpe006_processed_SCTtransformed_monomac.rds")


###########################################################################################################
#   No.05                                                                                                 #
#               Map alsaigh_ac_p1_processed_SCTtransformed_monomac.rds to mac.atlas.zenodo                #
#                                                                                                         #
###########################################################################################################

#data.query
setwd("/Users/lucia/Documents/34.immunometabolism/data/Alsaigh_carotid_human/")
ac_p1_seurat_sct.monomac <- readRDS("alsaigh_ac_p1_processed_SCTtransformed_monomac.rds")


#********map
sce.test <- LayerData(ac_p1_seurat_sct.monomac, assay = "RNA", layer = "counts")
sce.ref <- LayerData(mac.atlas.zenodo, assay = "RNA", layer = "data")
label.ref <- mac.atlas.zenodo$short.label

pred <- SingleR(test = sce.test, ref = sce.ref, labels = label.ref, de.method = "wilcox")

#******************test dataset 4803  cells, ref 363315 cells,need  -  calculation
table(pred$labels)
# 0_AlvMac     1_MetM2Mac   10_InflamMac  11_MetalloMac      12_MBMMac  13_CalciumMac    14_ProliMac      15_LYZMac 
#    1            341            110            104             18              6             41             51 
# 16_ECMHomeoMac     17_IFNMac3      18_ECMMac   19_ClassMono        2_C3Mac       20_TDoub     21_HemeMac     22_IFNMac4 
#       2157             24              7            541            218              7             35             70 
# 23_NA      3_ICIMac1      4_ICIMac2    5_StressMac  6_SPP1AREGMac       7_IFNMac      8_IFNGMac     9_AngioMac 
#   61             73            206             23            153            187              6            363 

ac_p1_seurat_sct.monomac$SingleR_map_to_TAM <- pred$labels
saveRDS(ac_p1_seurat_sct.monomac, "alsaigh_ac_p1_processed_SCTtransformed_monomac.rds")


###########################################################################################################
#   No.06                                                                                                 #
#               Map alsaigh_ac_p1_processed_SCTtransformed_monomac.rds to mac.atlas.zenodo                #
#                                                                                                         #
###########################################################################################################

#data.query
setwd("/Users/lucia/Documents/34.immunometabolism/data/Alsaigh_carotid_human/")
pa_p1_seurat_sct.monomac <- readRDS("alsaigh_pa_p1_processed_SCTtransformed_monomac.rds")


#********map
sce.test <- LayerData(pa_p1_seurat_sct.monomac, assay = "RNA", layer = "counts")
sce.ref <- LayerData(mac.atlas.zenodo, assay = "RNA", layer = "data")
label.ref <- mac.atlas.zenodo$short.label

pred <- SingleR(test = sce.test, ref = sce.ref, labels = label.ref, de.method = "wilcox")

#******************test dataset 266  cells, ref 363315 cells,need  1 hour  calculation
table(pred$labels)
#1_MetM2Mac 10_InflamMac    12_MBMMac  14_ProliMac    15_LYZMac   17_IFNMac3    18_ECMMac 19_ClassMono      2_C3Mac 
#    68            8            1            2            3            2            1           94           30 
#21_HemeMac   22_IFNMac4    3_ICIMac1    4_ICIMac2  5_StressMac     7_IFNMac   9_AngioMac 
#     2           13            2            1            3           16           20 

plotScoreHeatmap(pred)
pa_p1_seurat_sct.monomac$SingleR_map_to_TAM <- pred$labels
saveRDS(pa_p1_seurat_sct.monomac, "alsaigh_pa_p1_processed_SCTtransformed_monomac.rds")

###########################################################################################################
#   No.07                                                                                                 #
#               Map alsaigh_ac_p2_processed_SCTtransformed_monomac.rds to mac.atlas.zenodo                #
#                                                                                                         #
###########################################################################################################

#data.query
setwd("/Users/lucia/Documents/34.immunometabolism/data/Alsaigh_carotid_human/")
ac_p2_seurat_sct.monomac <- readRDS("alsaigh_ac_p2_processed_SCTtransformed_monomac.rds")


#********map
sce.test <- LayerData(ac_p2_seurat_sct.monomac, assay = "RNA", layer = "counts")
sce.ref <- LayerData(mac.atlas.zenodo, assay = "RNA", layer = "data")
label.ref <- mac.atlas.zenodo$short.label

pred <- SingleR(test = sce.test, ref = sce.ref, labels = label.ref, de.method = "wilcox")

#******************test dataset 1844  cells, ref 363315 cells,need  3 hours  calculation
table(pred$labels)
#0_AlvMac     1_MetM2Mac   10_InflamMac  11_MetalloMac      12_MBMMac  13_CalciumMac    14_ProliMac      15_LYZMac 
#   1            499            109             18              9             36             14             30 
#16_ECMHomeoMac     17_IFNMac3      18_ECMMac   19_ClassMono        2_C3Mac       20_TDoub     21_HemeMac     22_IFNMac4 
#   120                 16              4            241            166             10             18             55 
#23_NA      3_ICIMac1      4_ICIMac2    5_StressMac  6_SPP1AREGMac       7_IFNMac      8_IFNGMac     9_AngioMac 
#   44           28             53             42              9             74             18            230 
plotScoreHeatmap(pred)
ac_p2_seurat_sct.monomac$SingleR_map_to_TAM <- pred$labels
saveRDS(ac_p2_seurat_sct.monomac, "alsaigh_ac_p2_processed_SCTtransformed_monomac.rds")


###########################################################################################################
#   No.08                                                                                                 #
#               Map alsaigh_pa_p2_processed_SCTtransformed_monomac.rds to mac.atlas.zenodo                #
#                                                                                                         #
###########################################################################################################

#data.query
setwd("/Users/lucia/Documents/34.immunometabolism/data/Alsaigh_carotid_human/")
pa_p2_seurat_sct.monomac <- readRDS("alsaigh_pa_p2_processed_SCTtransformed_monomac.rds")


#********map
sce.test <- LayerData(pa_p2_seurat_sct.monomac, assay = "RNA", layer = "counts")
sce.ref <- LayerData(mac.atlas.zenodo, assay = "RNA", layer = "data")
label.ref <- mac.atlas.zenodo$short.label

pred <- SingleR(test = sce.test, ref = sce.ref, labels = label.ref, de.method = "wilcox")

#******************test dataset 283   cells, ref 363315 cells,need  1 hour  calculation
table(pred$labels)
#1_MetM2Mac 10_InflamMac    12_MBMMac  14_ProliMac    15_LYZMac   17_IFNMac3    18_ECMMac 19_ClassMono      2_C3Mac 
#45           26            1            1            2            2            1          122           29 
#20_TDoub   22_IFNMac4        23_NA    3_ICIMac1  5_StressMac     7_IFNMac    8_IFNGMac   9_AngioMac 
#1            8            3            3            3           14            1           21 

plotScoreHeatmap(pred)
pa_p2_seurat_sct.monomac$SingleR_map_to_TAM <- pred$labels
saveRDS(pa_p2_seurat_sct.monomac, "alsaigh_pa_p2_processed_SCTtransformed_monomac.rds")

###########################################################################################################
#   No.09                                                                                                 #
#               Map alsaigh_ac_p3_processed_SCTtransformed_monomac.rds to mac.atlas.zenodo                #
#                                                                                                         #
###########################################################################################################

#data.query
setwd("/Users/lucia/Documents/34.immunometabolism/data/Alsaigh_carotid_human/")
ac_p3_seurat_sct.monomac <- readRDS("alsaigh_ac_p3_processed_SCTtransformed_monomac.rds")


#********map
sce.test <- LayerData(ac_p3_seurat_sct.monomac, assay = "RNA", layer = "counts")
sce.ref <- LayerData(mac.atlas.zenodo, assay = "RNA", layer = "data")
label.ref <- mac.atlas.zenodo$short.label

pred <- SingleR(test = sce.test, ref = sce.ref, labels = label.ref, de.method = "wilcox")

#******************test dataset 2107   cells, ref 363315 cells,need 2.5  hours  calculation
table(pred$labels)
#0_AlvMac     1_MetM2Mac   10_InflamMac  11_MetalloMac      12_MBMMac  13_CalciumMac    14_ProliMac      15_LYZMac 
#   2            649            124              1              1              1             32             40 
#16_ECMHomeoMac     17_IFNMac3    18_ECMMac   19_ClassMono      2_C3Mac       20_TDoub     21_HemeMac     22_IFNMac4 
#      18              7              8            374            122              7              5             82 
#23_NA      3_ICIMac1      4_ICIMac2    5_StressMac  6_SPP1AREGMac       7_IFNMac      8_IFNGMac     9_AngioMac 
# 54             11            135             64             15             23              8            324 


plotScoreHeatmap(pred)
ac_p3_seurat_sct.monomac$SingleR_map_to_TAM <- pred$labels
saveRDS(ac_p3_seurat_sct.monomac, "alsaigh_ac_p3_processed_SCTtransformed_monomac.rds")

###########################################################################################################
#   No.10                                                                                                 #
#               Map alsaigh_pa_p3_processed_SCTtransformed_monomac.rds to mac.atlas.zenodo                #
#                                                                                                         #
###########################################################################################################

#data.query
setwd("/Users/lucia/Documents/34.immunometabolism/data/Alsaigh_carotid_human/")
pa_p3_seurat_sct.monomac <- readRDS("alsaigh_pa_p3_processed_SCTtransformed_monomac.rds")

#********map
sce.test <- LayerData(pa_p3_seurat_sct.monomac, assay = "RNA", layer = "counts")
sce.ref <- LayerData(mac.atlas.zenodo, assay = "RNA", layer = "data")
label.ref <- mac.atlas.zenodo$short.label

pred <- SingleR(test = sce.test, ref = sce.ref, labels = label.ref, de.method = "wilcox")

#******************test dataset 440   cells, ref 363315 cells,need 2  hours  calculation
table(pred$labels)
#1_MetM2Mac  10_InflamMac 13_CalciumMac   14_ProliMac     15_LYZMac     18_ECMMac  19_ClassMono       2_C3Mac 
#    119             6             2             4             2             2           214             9 
#20_TDoub    21_HemeMac    22_IFNMac4     4_ICIMac2    5_StressMac 6_SPP1AREGMac      7_IFNMac    9_AngioMac 
#     2             7            38             1            11             3             1            19 

plotScoreHeatmap(pred)
pa_p3_seurat_sct.monomac$SingleR_map_to_TAM <- pred$labels
saveRDS(pa_p3_seurat_sct.monomac, "alsaigh_pa_p3_processed_SCTtransformed_monomac.rds")

###########################################################################################################
#   No.11                                                                                                 #
#               Map Pan_ascending_aortic_n6_processed_SCTtransformed_monomac.rds to mac.atlas.zenodo                #
#                                                                                                         #
###########################################################################################################

#data.query
setwd("/Users/lucia/Documents/34.immunometabolism/data/Pan_normal_human/")
ascending_aortic_n6.monomac <- readRDS("Pan_ascending_aortic_n6_processed_SCTtransformed_monomac.rds")

#********map
sce.test <- LayerData(ascending_aortic_n6.monomac, assay = "RNA", layer = "counts")
sce.ref <- LayerData(mac.atlas.zenodo, assay = "RNA", layer = "data")
label.ref <- mac.atlas.zenodo$short.label

pred <- SingleR(test = sce.test, ref = sce.ref, labels = label.ref, de.method = "wilcox")

#******************test dataset 3066   cells, ref 363315 cells,need   hours  calculation
table(pred$labels)
# 1_MetM2Mac   10_InflamMac  11_MetalloMac      12_MBMMac  13_CalciumMac    14_ProliMac 16_ECMHomeoMac     17_IFNMac3 
# 1518            144             32              1             17             37              4              2 
# 18_ECMMac   19_ClassMono        2_C3Mac       20_TDoub     21_HemeMac     22_IFNMac4          23_NA      3_ICIMac1 
# 55             60             36              6            276             57             82             41 
# 4_ICIMac2    5_StressMac  6_SPP1AREGMac       7_IFNMac      8_IFNGMac     9_AngioMac 
# 130            550              1              1              1             15 

plotScoreHeatmap(pred)
ascending_aortic_n6.monomac$SingleR_map_to_TAM <- pred$labels
saveRDS(ascending_aortic_n6.monomac, "Pan_ascending_aortic_n6_processed_SCTtransformed_monomac.rds")
