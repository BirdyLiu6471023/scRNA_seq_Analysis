library(Matrix)
library(MatrixModels)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library (stringr)
library(tidyverse)
library(RCurl)
library(cowplot)
library(cluster)
library (factoextra)
library (metap)
library(ensembldb)
library(biomaRt)
library (kableExtra)
library(leiden)
library(igraph)
library(doParallel)

# ==================================== Note ====================================

# Generally, every time you rerun this code, you need to modify or check 2 place:
# 1) First part: Where you need to modify
# 2) Last part:  Assign Cell Type to Clusters
# (you could identify the part with the title beside the orange # between script and console parts)

# If other changes are made, would be:
# 1. 
# 2. 

# ========================== Where you need to modify ==========================
# give an id for this run, so that we could distinguish the results: 
id <- "sl06202023_res05_20000"

# resolution value: 
resolution_value <- 0.5

# where you store everything 
based_directory <- "/Users/macbook/Desktop/SGcell_Evolution/R_work_env" 

# where your data store: 
data_path <- "/Users/macbook/Desktop/Bio-Research\ /PS050_cellranger_count_outs/filtered_feature_bc_matrix/"

# if find all of biomarker and store them: 
store_biomarkers <- TRUE

# nUMI upper bound:
nUMI_upper <- 100000

# =======================Set up environment (don't look) ======================

# create a place to store the result from this run: 
setwd(glue::glue("{based_directory}")) 
dir.create(glue::glue("{id}_resolution_{resolution_value}"))

## setting the place you defined as working directory: 
setwd(glue::glue("{based_directory}/{id}_resolution_{resolution_value}")) 

## create directories to store biomarkers.csv, plots and so on: 
biomarkers_csv_directory <- glue::glue("biomarkers")
if (!dir.exists(biomarkers_csv_directory)) {
  print("creating biomarkers csv directory in your working environment...")
  dir.create(biomarkers_csv_directory)
  print("created!")
} else{
  print(":) biomarkers csv directory already exists!")
}

cell_type_biomarker_plot_directory <- "plot_celltype_with_biomarkers"
if (!dir.exists(cell_type_biomarker_plot_directory)) {
  print("creating cell type biomarker plot directory in your working environment...")
  dir.create(cell_type_biomarker_plot_directory)
  print("created!")
} else{
  print(":) cell type biomarker plot directory already exists!")
}

filtered_Rdata <- "filtered_Rdata"
if (!dir.exists(filtered_Rdata)) {
  print("creating filtered Rdata directory in your working environment...")
  dir.create(filtered_Rdata)
  print("created!")
} else{
  print(":) filtered Rdata directory already exists!")
}

# ====================== Read data & Calculate Metrics==========================
EpCAM.data <- Read10X(data.dir = data_path)

# Create Seurat object
EpCAM <- CreateSeuratObject(counts = EpCAM.data, project = "mSG_EpCAM", save.SNN = TRUE)

# Add number of genes per UMI for each cell to metadata
EpCAM$log10GenesPerUMI <- log10(EpCAM$nFeature_RNA) / log10(EpCAM$nCount_RNA)

# Compute percent mito ratio
EpCAM$mitoRatio <- PercentageFeatureSet(object = EpCAM, pattern = "^mt-")
EpCAM$mitoRatio <- EpCAM@meta.data$mitoRatio / 100

# Create metadata dataframe
metadata <- EpCAM@meta.data

# Add cell IDs to metadata
metadata$cells <- rownames(metadata)

# Rename columns
metadata <- metadata %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)
# Create sample column
metadata$sample <- "mSG"

# Add metadata back to Seurat object
EpCAM@meta.data <- metadata


## ===================== Cell-Level Filtering ==================================
# Cell-level FILTERING
set.seed (12)
filtered_seurat <- subset(x = EpCAM, 
                          subset= (nUMI >= 1) & 
                            (nUMI <= nUMI_upper) &
                            (nGene >= 100) & 
                            (log10GenesPerUMI > 0.5) &
                            (mitoRatio < 0.20))

# Output a logical vector for every gene on whether the more than zero counts per cell
# Extract counts
counts <- GetAssayData(object = filtered_seurat, slot = "counts")
nonzero <- counts > 0
# Sums all TRUE values and returns TRUE if more than 1 TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= 1
# Only keeping those genes expressed in more than 1 cell (include rare populations e.g. myoepithelial)
filtered_counts <- counts[keep_genes, ]
# Reassign to filtered Seurat object
filtered_seurat <- CreateSeuratObject(filtered_counts, meta.data = filtered_seurat@meta.data)
# Save filtered subset to new metadata
#metadata_clean <- filtered_seurat@meta.data
filtered_seurat[["percent.mt"]] <- PercentageFeatureSet(filtered_seurat, pattern = "^mt-")

# ===========================Normalize Count====================================

# Normalize the counts, regress out mitoRatio 
set.seed(12)
filtered_seurat_norm <- SCTransform(filtered_seurat, vars.to.regress = "mitoRatio")

# Create .RData object to load at any time
#save(filtered_seurat_norm, file="E:/Columbia/Piera lab/R/Data/RNA velocity/PS050/analysis/filtered_seurat_norm.RData")

# Identification of highly variable features (feature selection)
set.seed(12)
filtered_seurat_norm <- FindVariableFeatures(filtered_seurat_norm, selection.method = "vst", nfeatures = 2000)

#=============================Run PCA===========================================
# Run PCA
set.seed (12)
filtered_seurat_norm_PCA <- RunPCA(object = filtered_seurat_norm)
# Plot PCA
PCAPlot(filtered_seurat_norm_PCA, pt.size = 3)+theme(text = element_text(size=30))
ggsave("PCAPlot.png", h = 5000, w = 7000, units = "px")

#=============================Run UMAP==========================================
# Run UMAP
set.seed(12)
filtered_seurat_norm_UMAP <- RunUMAP(filtered_seurat_norm_PCA, 
                                     dims = 1:15,
                                     reduction = "pca")
# Plot UMAP                             
DimPlot(filtered_seurat_norm_UMAP, pt.size = 3)+theme(text = element_text(size=30))
ggsave("UMAPPlot.png", h = 5000, w = 7000, units = "px")

# ==========================Clustering =========================================
set.seed (12)
filtered_seurat_norm_UMAP <- FindNeighbors(object = filtered_seurat_norm_UMAP, reduction = "pca", dims = 1:15, verbose = TRUE, graph.name = "mSG")

set.seed (12)
filtered_seurat_norm_UMAP <- FindClusters(object = filtered_seurat_norm_UMAP, verbose = TRUE, algorithm = 4, resolution = resolution_value, graph.name = "mSG")


set.seed(12)
embed_matrix <- Embeddings(filtered_seurat_norm_UMAP[['umap']])
embed_matrix <- as.matrix(embed_matrix)
distance_matrix <- dist(embed_matrix)

#Isolate cluster name
clusters <- filtered_seurat_norm_UMAP@active.ident

#Compute silhouette score for each cluster
silhouette <- silhouette(as.numeric(clusters), dist = distance_matrix)
filtered_seurat_norm_UMAP@meta.data$silhouette_score <- silhouette[,3]
fviz_silhouette(silhouette, label = FALSE, print.summary = TRUE)
filename<- glue::glue("silhouette_score_{resolution_value}.png")
ggsave(filename, h = 2000, w = 4000, units = "px")


# =================================Find BioMarker===============================
# how many cluster:

num_cluster <- max(as.numeric(clusters))

#num_cluster <- (filtered_seurat_norm_UMAP@meta.data %>%
#  dplyr::distinct(clusters) %>%
#  dplyr::count())[1,1]

# if find and store biomarkers: 
if (store_biomarkers){
for (i in 1:num_cluster){
  set.seed (12)
  markers <- FindMarkers(filtered_seurat_norm_UMAP, ident.1 = i, logfc.threshold = 0.25)
  filename<-glue::glue("biomarkers/cluster{i}_markers_{resolution_value}.csv")
  write.csv2 (markers, file = filename)
  print("file saved!")
}}

# ==============================Loading Data====================================





# ===============================Visualization==================================
# Normalize RNA data for visualization purposes
set.seed(12)
filtered_seurat_norm_UMAP_RNA <- NormalizeData(filtered_seurat_norm_UMAP, verbose = FALSE)

#--------------------------Acinar1: Car 6, Lpo, Bhlha15, Pip--------------------
FeaturePlot(filtered_seurat_norm_UMAP_RNA, 
            reduction = "umap", 
            features = c("Car6", "Lpo", "Bhlha15", "Pip"), 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE,
            repel = TRUE, pt.size = 3, label.size = 8)
ggsave("plot_celltype_with_biomarkers/scRNAseq_mSG_UMAP_Acinar1.png", h = 5000, w = 8000, units = "px")

# --------- Acinar2: Pip, Smgc, Prlr, Bhlha15 ----------------------------------
FeaturePlot(filtered_seurat_norm_UMAP_RNA, 
            reduction = "umap", 
            features = c("Pip", "Smgc", "Prlr", "Bhlha15"), 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE,
            repel = TRUE, pt.size = 3, label.size = 8)
ggsave("plot_celltype_with_biomarkers/scRNAseq_mSG_UMAP_Acinar2.png", h = 5000, w = 8000, units = "px")

#----------Ductal1: Klk1, Ngf, Clcnkb, Bsnd -----------------------------------
FeaturePlot(filtered_seurat_norm_UMAP_RNA, 
            reduction = "umap", 
            features = c ("Klk1", "Ngf", "Clcnkb", "Bsnd"), 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE,
            repel = TRUE, pt.size = 3, label.size = 8)
ggsave("plot_celltype_with_biomarkers/scRNAseq_mSG_UMAP_Ductal1.png", h = 5000, w = 8000, units = "px")

#----------Ductal2:Foxi1, Ascl3,Clcnkb, Bsnd, Kit, Elf5-------------------------
FeaturePlot(filtered_seurat_norm_UMAP_RNA, 
            reduction = "umap", 
            features = c ("Foxi1", "Ascl3", "Clcnkb", "Bsnd", "Kit", "Elf5"), 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE,
            repel = TRUE, pt.size = 3, label.size = 8)
ggsave("plot_celltype_with_biomarkers/scRNAseq_mSG_UMAP_Ductal2.png", h = 5000, w = 8000, units = "px")

#----------Ductal3:Slc9a3, Clic6,Kit, Elf5--------------------------------------
FeaturePlot(filtered_seurat_norm_UMAP_RNA, 
            reduction = "umap", 
            features = c ("Foxi1", "Ascl3", "Clcnkb", "Bsnd"), 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE,
            repel = TRUE, pt.size = 3, label.size = 8)
ggsave("plot_celltype_with_biomarkers/scRNAseq_mSG_UMAP_Ductal3.png", h = 5000, w = 8000, units = "px")


#----------Basal: Krt14, Krt15--------------------------------------------------
FeaturePlot(filtered_seurat_norm_UMAP_RNA, 
                  reduction = "umap", 
                  features = c ("Krt14", "Krt15"), 
                  order = TRUE,
                  min.cutoff = 'q10', 
                  label = TRUE,
                  repel = TRUE, pt.size = 3, label.size = 8)
ggsave("plot_celltype_with_biomarkers/scRNAseq_mSG_UMAP_Basal.png", h = 5000, w = 8000, units = "px")


#======================= Assign Cell Type to Clusters============================
# Create a lookup table for cluster renaming
filtered_seurat_norm_UMAP <- RenameIdents(filtered_seurat_norm_UMAP, 
                                          `1` = "Ductal2", 
                                          `2` = "Ductal1", 
                                          `3` = "?", 
                                          `4` = "Ductal3_1", 
                                          `5` = "Acinar2", 
                                          `6` = "Ductal3_2", 
                                          `7` = "Acinar1", 
                                          `8` = "Basal", 
                                          `9` = "Endothelial")


filtered_data <- glue::glue("filtered_data/filtered_seurat_norm_{resolution_value}.RData")
save(filtered_seurat_norm_UMAP, file="filtered_Rdata/filtered_seurat_norm.RData")


