#Review on in immunometabolism in macrophages - cancer/atherosclerosis
#https://zenodo.org/records/11222158
#Note this atlas does not include raw count data as these can be download from the original sources, 
#but does include log-normalized count data, SCTransformed data and integrated data via RPCA.

library(Seurat)
library(SeuratObject)
library(tidyverse)
setwd("/Users/lucia/Documents/34.immunometabolism")
mac.atlas.zenodo <- readRDS("mac.atlas.zenodo.200524.rds")

mac.atlas.zenodo.metadata <- mac.atlas.zenodo@meta.data
View(mac.atlas.zenodo.metadata)
rownames(mac.atlas.zenodo.metadata) <- mac.atlas.zenodo.metadata$cellid
mac.atlas.zenodo.metadata$orig.ident <- mac.atlas.zenodo.metadata$cellid
mac.atlas.zenodo@meta.data <- mac.atlas.zenodo.metadata
# Assuming 'seurat_obj' is your Seurat object and 'umap_1' and 'umap_2' are in the metadata
UMAP_1 = mac.atlas.zenodo@meta.data$UMAP_1
UMAP_2 = mac.atlas.zenodo@meta.data$UMAP_2
mac.atlas.zenodo.umap <- data.frame(UMAP_1 = UMAP_1, UMAP_2 = UMAP_2)
rownames(mac.atlas.zenodo.umap) <- rownames(mac.atlas.zenodo.metadata)
mac.atlas.zenodo.umap <- as.matrix(mac.atlas.zenodo.umap)
# Create a new 'DimReduc' object (UMAP)
umap_reduction <- CreateDimReducObject(
  embeddings = mac.atlas.zenodo.umap,
  key = "UMAP_",
  assay = "RNA"  # Or specify the assay you are using, if it's different
)
# Add the UMAP reduction to the Seurat object
mac.atlas.zenodo[["umap"]] <- umap_reduction
UMAPPlot(mac.atlas.zenodo,group.by = "short.label",label = TRUE)

#Extract the colors of the UMAP plot from the pan-cancer analysis review article
library(png)
library(cluster)
library(ggplot2)

img <- readPNG("plot1.png")
# Obtain the width, height, and number of color channels of the image (typically 3 channels: R, G, B)
dim(img)
# Convert image data into a DataFrame, where each row represents the RGB values of a pixel
# img is a three-dimensional array with dimensions (width, height, 3)
img_data <- data.frame(
  R = as.vector(img[,,1]),
  G = as.vector(img[,,2]),
  B = as.vector(img[,,3])
)

head(img_data)

set.seed(123)  
kmeans_result <- kmeans(img_data, centers = 50)  # Extract 30 colors

# Obtain the colors (RGB values) of the clustering results
cluster_colors <- kmeans_result$centers
# Display the extracted colors (RGB values)
print("main color (RGB values):")
print(cluster_colors)
# Obtain the color category for each pixel
img_data$cluster <- factor(kmeans_result$cluster)

# Calculate the frequency of color occurrences in each cluster
color_counts <- img_data %>%
  group_by(cluster) %>%
  summarise(count = n())

# Obtain the colors (RGB) of the clusters
cluster_colors_rgb <- rgb(cluster_colors[,1], cluster_colors[,2], cluster_colors[,3], maxColorValue = 1)

color_counts$color <- cluster_colors_rgb

# Visualize the frequency of occurrence for each color
ggplot(color_counts[1:25,], aes(x = cluster, y = 1, fill = color)) +
  geom_bar(stat = "identity") +
  scale_fill_identity() +  # Fill the bars using RGB color values
  theme_void() +
  ggtitle("Top 5 Colors in Image")  +  
  geom_text(aes(label = color), 
            position = position_stack(vjust = 0.5),  
            color = "black",  
            size = 5,
            angle = 90
            )  


ggplot(color_counts[26:50,], aes(x = cluster, y = 1, fill = color)) +
  geom_bar(stat = "identity") +
  scale_fill_identity() +  # Fill the bars using RGB color values
  theme_void() +
  ggtitle("Top 5 Colors in Image")  +  
  geom_text(aes(label = color), 
            position = position_stack(vjust = 0.5),  
            color = "black",  
            size = 5,
            angle = 90
  )  

#check the color of each cell type
table( mac.atlas.zenodo$short.label )
#  0_AlvMac         1_MetM2Mac        2_C3Mac      3_ICIMac1            4_ICIMac2 
#  #50609F          #BF493E           #769765      #E7E291              #4E6980
#  5_StressMac      6_SPP1AREGMac     7_IFNMac     8_IFNGMac            9_AngioMac 
#  #B16743          #64A5D0           #7B334C      #7BAB77              #C1909A
#  10_InflamMac     11_MetalloMac     12_MBMMac    13_CalciumMac        14_ProliMac 
#  #854E30          #7B7788           #AD4C45      #C18D66              #6C6497
#  15_LYZMac        16_ECMHomeoMac    17_IFNMac3   18_ECMMac            19_ClassMono 
#  #CEA876          #35224D           #C1CFB1      #533068               #992F5F
#  20_TDoub         21_HemeMac        22_IFNMac4          23_NA 
#  #D7C182          #672456           #BA963F             #FFFFFF

custom_colors <- c(
  "#50609F", "#BF493E", "#769765", "#E7E291", "#4E6980", 
  "#B16743","#64A5D0", "#7B334C","#7BAB77", "#C1909A",
  "#854E30", "#7B7788","#AD4C45","#C18D66", "#6C6497",
  "#CEA876", "#35224D", "#C1CFB1","#533068", "#992F5F",
  "#D7C182", "#672456", "#BA963F", "#FFFFFF"
  )
#"#50609F""#6A7BA6""#5C6B8E""#7A8AAB""#4E5A79""#5B6D8A"
UMAPPlot(mac.atlas.zenodo,group.by = "short.label",label = TRUE,cols = custom_colors, raster = FALSE)


library(ape)
library(ggnewscale)
library(ggtree)
library(ggrepel)
############################
#FUNCTIONS 
############################
adjust.label = function(cluster, x, y){
  label.coords$x.mean[label.coords$short.label == cluster] =
    label.coords$x.mean[label.coords$short.label == cluster] + x
  label.coords$y.mean[label.coords$short.label == cluster] = 
    label.coords$y.mean[label.coords$short.label == cluster] + y
  label.coords
}
############################
#CLUSTER UMAP 
############################
set.seed(1)
mac.atlas.zenodo.metadata = mac.atlas.zenodo.metadata[order(mac.atlas.zenodo.metadata$short.label), ]
#randomize order for better cluster visualization
mac.atlas.zenodo.metadata = mac.atlas.zenodo.metadata[sample(1:nrow(mac.atlas.zenodo.metadata), nrow(mac.atlas.zenodo.metadata)), ]

mac.atlas.zenodo.metadata$cluster = as.factor(mac.atlas.zenodo.metadata$cluster)

mac.atlas.zenodo.metadata = mac.atlas.zenodo.metadata %>%
  filter(cluster != 23)


############################
#GENERATE LABEL COORDS 
############################

label.coords = mac.atlas.zenodo.metadata %>%
  group_by(short.label) %>%
  summarize(
    x.mean = mean(UMAP_1),
    y.mean = mean(UMAP_2)
  )

label.coords$lab.text.color = "dark"
label.coords$lab.text.color[
  label.coords$short.label %in% c(
    '16_ECMHomeoMac',
    '18_ECMMac',
    '4_ICIMac2',
    '0_AlvMac',
    '7_IFNMac',
    '21_HemeMac'
  )
] = 'light'

#adjustments to coords
label.coords = adjust.label('0_AlvMac', -1, 0)
label.coords = adjust.label('1_MetM2Mac', 0.3, -2)
label.coords = adjust.label('4_ICIMac2', -2, 0)
label.coords = adjust.label('18_ECMMac', 0, 0.5)
label.coords = adjust.label('17_IFNMac3', 0, 1)
label.coords = adjust.label('16_ECMHomeoMac', 2, 1)
label.coords = adjust.label('3_ICIMac1', 1, 0)
label.coords = adjust.label('5_StressMac', 1, -1)
label.coords = adjust.label('12_MBMMac', -0.5, 1)
label.coords = adjust.label('13_CalciumMac', -1, 0)
label.coords = adjust.label('11_MetalloMac', 0, 1)
label.coords = adjust.label('6_SPP1AREGMac', 1.3, -0.8)
label.coords = adjust.label('2_C3Mac', 0, -0.5)
label.coords = adjust.label('10_InflamMac', 1, -1)
label.coords = adjust.label('20_TDoub', 0, 0.6)
label.coords = adjust.label('15_LYZMac', 0, 0.3)
label.coords = adjust.label('7_IFNMac', 0, 2)
label.coords = adjust.label('21_HemeMac', 3.5, -1)


library(ggsci)

clus.umap.raw = ggplot(mac.atlas.zenodo.metadata, aes(x = UMAP_1, y = UMAP_2, color = short.label)) +
  geom_point() +
  #scale_color_igv() +
  scale_color_manual( values = custom_colors    ) +
  new_scale_color() +
  geom_label(
    inherit.aes = F,
    data = label.coords,
    aes(
      label = short.label,
      x = x.mean,
      y = y.mean,
      fill = short.label,
      color = lab.text.color
    ),
    size = 8,
    label.size = 1.2
  ) +
  scale_color_manual(values = c('#000000', '#ffffff')) +
  #scale_fill_igv() +
  scale_fill_manual( values = custom_colors    ) +
  theme_classic() +
  xlab('UMAP1') +
  ylab('UMAP2') +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = 'none',
    axis.title = element_text(size = 18)
  ) 
#save as 1400*1400 picture

###################################################################
#Process each dataset individually, perform SCTtransform on each dataset, and keep the filtering strategy consistent with the pan-cancer review.
DefaultAssay(mac.atlas.zenodo) <- "SCT"
VlnPlot(mac.atlas.zenodo,features = c("ZEB1","ZEB2","CXCR2","CX3CR1","CD68","CD14","ADGRE1","HLA-DRA","GAPDH","SNAI1","SNAI2","TWIST1","TWIST2","FOXC2"),group.by = "short.label",raster=FALSE, stack = TRUE, sort = FALSE, flip = TRUE)+ theme(legend.position = "none")+ggtitle("") + theme(plot.margin = unit(c(0.1,0.1,0.1,1), "inches"))
DotPlot(mac.atlas.zenodo, idents = levels(mac.atlas.zenodo$short.label)[1:23]   ,  features = c("ZEB1","ZEB2","CXCR2","CX3CR1","CD68","CD14","ADGRE1","HLA-DRA","GAPDH","SNAI1","SNAI2","TWIST1","TWIST2","FOXC2"),group.by = "short.label" )

#############Classify the cells into the “16_ECMHomeoMac” group and “Others”.  
mac.atlas.zenodo@meta.data$celltype2 <- "Others"

length(which(mac.atlas.zenodo@meta.data$short.label == "16_ECMHomeoMac"    ))#10011
table( mac.atlas.zenodo@meta.data$short.label   )#10011

mac.atlas.zenodo@meta.data$celltype2[which(   mac.atlas.zenodo@meta.data$short.label == "16_ECMHomeoMac"  )] <- "16_ECMHomeoMac" 
table(mac.atlas.zenodo@meta.data$celltype2)
#16_ECMHomeoMac         Others 
#10011         353304 
table(mac.atlas.zenodo@meta.data$celltype2, mac.atlas.zenodo@meta.data$tissue)
#                Blood Effusion LymphNode Normal  Tumor
#16_ECMHomeoMac      0       51       184    429   9347
#Others            252     2736      6006  74553 269757


#Include only Tumor cells and compare '16_ECMHomeoMac' with 'Others'.
Idents(mac.atlas.zenodo) <- "tissue"
#VlnPlot(mac.atlas.zenodo,features = c("FABP5","CD36","CD63","CD68","PLIN2","LIPA","APOE","IL1B","HSPA5","SIRPA","LPL"),group.by = "short.label",raster=FALSE, stack = TRUE, sort = FALSE, flip = TRUE)+ theme(legend.position = "none")+ggtitle("") + theme(plot.margin = unit(c(0.1,0.1,0.1,1), "inches"))
p = VlnPlot(mac.atlas.zenodo, idents = "Tumor" , features = c("FABP5","CD36","CD68","PLIN2","LIPA","APOE","LPL","LAMP1","TREM2","ABCA1","MARCO" ),group.by = "celltype2", cols = c("#AD867E","#a47053","#7b5965","#95594c","#533068","#66597b","#50609f","#4e6980","#678171","#7bab77","#979c70") ,raster=FALSE, stack = TRUE, sort = FALSE, flip = TRUE)+ 
  theme(legend.position = "none") + 
  ggtitle("Pan-cancer") + 
  theme(plot.margin = unit(c(0.1,0.1,0.1,1), "inches")) +
  scale_x_discrete(labels = c("16_ECMHomeoMac" = "Cluster16", "Others" = "Others"))
  
ggsave(p,filename = "plots/Lipid_markers_of_16_ECMHomeoMac_in_Tumor.jpeg",height = 6,width = 4) 


p = VlnPlot(mac.atlas.zenodo, idents = "Tumor" ,features = c("LGALS1", "LGALS3", "SIRPA", "MMP9","MRC1", "HLA-DPA1", "HLA-DRB1"), cols = c("#95594c","#533068","#66597b","#50609f","#4e6980","#678171","#7bab77"),group.by = "celltype2",raster=FALSE, stack = TRUE, sort = FALSE, flip = TRUE)+ 
  theme(legend.position = "none") + 
  ggtitle("Pan-cancer") + 
  theme(plot.margin = unit(c(0.1,0.1,0.1,1), "inches"))+
  scale_x_discrete(labels = c("16_ECMHomeoMac" = "Cluster16", "Others" = "Others"))

ggsave(p,filename = "plots/Immumo_markers_of_16_ECMHomeoMac_in_Tumor.jpeg",height = 4,width = 4) 
