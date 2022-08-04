###SPatial Genomics####
##set working directory##
getwd()

##install/load needed packages##
library(Seurat)
library(dplyr) # data manipulation
library(ggplot2)
library(patchwork)
install.packages('BiocManager')
BiocManager::install('limma')
library(limma) # optional
install.packages("hdf5r")
library(hdf5r)
getwd()
#load the data
# specify your working directory
root_dir <- "/Volumes/GoogleDrive/My Drive/Spatial Genomics" 
setwd(root_dir)
# url prefix to download files
url_prefix <- "https://drive.google.com/uc?export=download&id="
# url IDs are extracted from Google Drive
spatial_folder_id <- "13B6ulZwz9XmWBD1RsnTHTMn3ySgIpFCY"
h5_mat_id <- "1nnch2ctksJ4rXYG5FzlPh3bIgnylTicV"
# file names
spatial_folder_name <- "./V1_Mouse_Brain_Sagittal_Posterior_Section_2_spatial.tar.gz"
h5_mat_name <- "./V1_Mouse_Brain_Sagittal_Posterior_Section_2_filtered_feature_bc_matrix.h5"
# download files from the resource URL to a destination path where the file is saved 
download.file(url = paste0(url_prefix, spatial_folder_id), 
              destfile = file.path(spatial_folder_name))
download.file(url = paste0(url_prefix, h5_mat_id), 
              destfile = file.path(h5_mat_name))
# extract (or list) contents from the tar archives
untar(spatial_folder_name)
untar(spatial_folder_name, list = TRUE)  # list contents
getwd()
# Load a 10x Genomics Visium Spatial Experiment into a Seurat object
brain_data <- Seurat::Load10X_Spatial(
  # The directory contains the read count matrix H5 file and the image data in a subdirectory called `spatial`. 
  data.dir = root_dir, 
  filename = h5_mat_name,
  assay = "Spatial", # specify name of the initial assay
  slice = "slice1", # specify name of the stored image
  filter.matrix = TRUE, 
  to.upper = FALSE
)
#Extract dimension of the active assay using dim,nrow or ncol.
dim(x =brain_data)
brain_data
dim(brain_data)
names(x=brain_data)
nrow(x=brain_data)
ncol(x=brain_data)
#feature and sample names
#Extract the feature and sample names using rownames or colnames.
head(x = rownames(brain_data), n = 10)
tail(rownames(brain_data), n = 5)
rownames(brain_data)
#Sample-level metadata
class(brain_data[[]])
colnames(brain_data[[]]) # automatically calculated while creating the Seurat object
rownames(brain_data[[]])
head(brain_data@meta.data)
brain_data$nCount_Spatial[1:3]
# nFeature_Spatial: the number of unique genes in each sample
sum(brain_data$nFeature_Spatial ==  colSums(brain_data@assays$Spatial@counts > 0))
# nCount_Spatial: the total number of detected molecules in each sample
sum(brain_data$nCount_Spatial ==  colSums(brain_data@assays$Spatial@counts))
##objects(e.g. Assay) together with feature-level metadata
# A vector of names of associated objects can be had with the names function
# These can be passed to the double [[ extract operator to pull them from the Seurat object
names(x = brain_data)
brain_data[['Spatial']]
brain_data[['slice1']]
brain_data@assays$Spatial@counts[5:10, 1:3]
#Feature-level metadata is associated with each individual assay. Feature-level metadata can be accessed through the double bracket [[ extract operator on the Assay objects, or the meta.features slot.
brain_data[['Spatial']]@meta.features
head(brain_data[['Spatial']][[]])                                                                                                                                     
##Analysis pipeline in Seurat
#The general steps to perform preprocessing, dimensiona reduction and clustering for spatial transcriptomics data are quite similar to the Seurat workflow analyzing single-cell RNA sequencing data
#Quality control
brain_data$percent.mt<- PercentageFeatureSet(brain_data, pattern = "^mt-")
VlnPlot(
  brain_data, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), 
  pt.size = 0.1, ncol = 3) & 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
# Jointly (rather than separately) consider the QC metrics when filtering
plot1 <- FeatureScatter(
  brain_data, feature1 = "nCount_Spatial", feature2 = "percent.mt") + NoLegend()
plot2 <- FeatureScatter(
  brain_data, feature1 = "nCount_Spatial", feature2 = "nFeature_Spatial") +
  NoLegend()
plot1 + plot2
brain_subset <- subset(
  brain_data, 
  subset = nFeature_Spatial < 8000 & nFeature_Spatial > 1000 & 
    nCount_Spatial < 50000 & percent.mt < 30)
# filterout 
print(paste("Filter out", ncol(brain_data) - ncol(brain_subset), 
            "samples because of the outlier QC metrics, with", ncol(brain_subset),
            "samples left."))
#normalization
SpatialFeaturePlot(
  brain_subset, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt")) &
  theme(legend.position = "bottom")  
brain_norm <- SCTransform(brain_subset, assay = "Spatial", verbose = FALSE)

names(brain_norm)  # new names
dim(brain_norm@assays$SCT@counts)  # check dimension after normalization
dim(brain_norm@assays$SCT@data) 
dim(brain_norm@assays$SCT@scale.data)
#Downstream tasks
brain_obj <- RunPCA(brain_norm, assay = "SCT", verbose = FALSE)
# compute K nearest neighbors (KNN)
brain_obj <- FindNeighbors(brain_obj, reduction = "pca", dims = 1:30)
# Leiden algorithm for community detection
brain_obj <- FindClusters(brain_obj, verbose = FALSE)
# PCA result is the default UMAP input, use dimensions 1:30 as input features
brain_obj <- RunUMAP(brain_obj, reduction = "pca", dims = 1:30)

plot3 <- DimPlot(brain_obj, reduction = "umap", label = TRUE) + NoLegend()
plot4 <- SpatialDimPlot(brain_obj, label = TRUE, label.size = 3) + NoLegend()
plot3 + plot4
brain_obj@reductions # show all the results in the seurat
# identity class of each sample
table(brain_obj@active.ident)
# find all markers of cluster 10
cluster10_markers <- FindMarkers(brain_obj, ident.1 = 10, min.pct = 0.25)
head(cluster10_markers, n = 5)
VlnPlot(brain_obj, features = c("Esr1", "Trh", "Ccn3"))
SpatialFeaturePlot(object = brain_obj, 
                   features = rownames(cluster10_markers)[1:3], 
                   alpha = c(0.1, 1), ncol = 3)
# find markers for every cluster compared to all remaining cells, 
# report only the positive ones
# this code chunk is not evaluated for now because of time constraints
brain_obj_markers <- FindAllMarkers(brain_obj, only.pos = TRUE, min.pct = 0.25, 
                                    logfc.threshold = 0.25)
brain_obj_markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
top10 <- brain_obj_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)
DoHeatmap(brain_obj, features = top10$gene) + NoLegend()
brain_moransi <- FindSpatiallyVariableFeatures(
  brain_obj, assay = "SCT", 
  features = VariableFeatures(brain_obj)[1:10],
  selection.method = "moransi") 
install.packages('Rfast2')
moransi_output_df <- brain_moransi@assays$SCT@meta.features %>%
  na.exclude
head(moransi_output_df[order(moransi_output_df$MoransI_observed, decreasing = T), ])
top_features_moransi <- head(
  SpatiallyVariableFeatures(brain_moransi, 
                            selection.method = "moransi"), 3)
SpatialFeaturePlot(brain_moransi, 
                   features = top_features_moransi, ncol = 3, alpha = c(0.1, 1)) + 
  plot_annotation(
    title = "Top 3 genes with the largest Moran's I",
    subtitle = "among 10 top variable genes for illustration purposes")
bottom_features_moransi <- tail(
  SpatiallyVariableFeatures(brain_moransi, 
                            selection.method = "moransi"), 3)
SpatialFeaturePlot(brain_moransi, 
                   features = bottom_features_moransi, ncol = 3, alpha = c(0.1, 1)) + 
  plot_annotation(
    title = "Bottom 3 genes with the smallest Moran's I",
    subtitle = "among 10 top variable genes for illustration purposes")
brain_variogram <- FindSpatiallyVariableFeatures(
  brain_obj, assay = "SCT", 
  features = VariableFeatures(brain_obj)[1:10],
  selection.method = "markvariogram")  
variogram_output_df <- brain_variogram@assays$SCT@meta.features %>%
  na.exclude # there are NA rows b/c we only calculated the variogram for 10 genes
head(variogram_output_df[order(variogram_output_df$r.metric.5), ])
top_features_variogram <- head(
  SpatiallyVariableFeatures(brain_variogram, 
                            selection.method = "markvariogram"), 3)
SpatialFeaturePlot(brain_variogram, 
                   features = top_features_variogram, ncol = 3, alpha = c(0.1, 1)) + 
  plot_annotation(
    title = "3 genes with the top spatially variable rank (by mark-variogram)",
    subtitle = "among 10 top variable genes for illustration purposes")
bottom_features_variogram <- tail(
  SpatiallyVariableFeatures(brain_variogram, 
                            selection.method = "markvariogram"), 3)
SpatialFeaturePlot(brain_variogram, 
                   features = bottom_features_variogram, ncol = 3, alpha = c(0.1, 1)) + 
  plot_annotation(
    title = "3 genes with the bottom spatially variale rank (by mark-variogram)",
    subtitle = "among 10 top variable genes for illustration purposes")
sessionInfo()
