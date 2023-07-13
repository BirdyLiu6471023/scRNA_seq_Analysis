# ==================== Load Packages and information =============

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
library(ensembldb)
library(biomaRt)
library (kableExtra)
library(leiden)
library(igraph)
library(doParallel)

convertString <- function(s) {
  firstChar <- substr(s, 1, 1)
  restChars <- tolower(substr(s, 2, nchar(s)))
  converted <- paste(firstChar, restChars, sep = "")
  return(converted)
}

biomarkers_dictionary_dataset1 <- list(
  "Acinar1" = c("Car6", "Lpo", "Bhlha15", "Pip"),
  "Acinar2" = c("Serpinb11", "Pip", "Smgc", "Prlr", "Bhlha15"), 
  "Ductal1" = c("Klk1", "Ngf", "Clcnkb", "Bsnd"),
  "Ductal2" = c("Foxi1", "Ascl3", "Clcnkb", "Bsnd", "Kit", "Elf5", "Hepacam2"),
  "Ductal3" = c("Slc9a3","Clic6","Kit","Elf5","Esp4","Gjb2"),
  "Basal" = c("Ucma","Krt5", "Krt14", "Krt15", "Ptn", "Tacstd2"))

biomarkers_dictionary_dataset2 <- list(
  "Basal" = c("Ucma","Krt5", "Krt14", "Krt15", "Ptn", "Tacstd2"), 
  "Cycling" = c('Mki67'),
  "Myoepithelial" = c("Acta2", "Cnn1"),
  "Endothelial" = c("Cdh5", "Pecam1", "Epcam"),
  "Intermediate" = c("Orm1", "Orm3", "Ifi202b", "Defb1","Cat", "Tacstd2"),
  "Krt19high" = c("Krt19", "Duox2","Hsd11b1"),
  "Ductal3" = c("Slc9a3","Clic6","Kit","Elf5","Esp4","Gjb2"),
  "Stromal-like" = c("Pdgfrb", "Twist1")
)

plot_genes <- function(object, feature_list, cell_name, labeled = FALSE){
  FeaturePlot(object, 
              reduction = "umap", 
              features = feature_list, 
              order = TRUE,
              min.cutoff = 'q10', 
              label = TRUE,
              repel = TRUE, pt.size = 3, label.size = 8)
  
  if (labeled){
    file_name <- glue::glue("plot_celltype_with_biomarkers/labeled_combined_UMAP_{cell_name}.png")
  }else{
    file_name <- glue::glue("plot_celltype_with_biomarkers/combined_UMAP_{cell_name}.png")
  }
  ggsave(file_name, h = 5000, w = 8000, units = "px")
  
  VlnPlot(object = object,features = feature_list)
  if (labeled){
    file_name <- glue::glue("plot_celltype_with_biomarkers/labeled_combined_Vln_{cell_name}.png")
  }else{
    file_name <- glue::glue("plot_celltype_with_biomarkers/combined_Vln_{cell_name}.png")
  }
  ggsave(file_name, h = 5000, w = 8000, units = "px")
  
}

# ========================== Where you need to modify ==========================
# give an id for this run, so that we could distinguish the results: 
id <- "Whole_sl0713_UMI20000_log0.5_remove_ribosome"

# resolution value: 
resolution_value <- 0.5

# where you store everything 
based_directory <- "/Users/macbook/Desktop/SGcell_Evolution/R_work_env" 

# if find all of biomarker and store them: 
store_biomarkers <- TRUE

dataset_1 = TRUE


# =======================Set up environment (don't look) ======================
if (dataset_1) {
  biomarkers_dic <- biomarkers_dictionary_dataset1
  data_path <-"/Users/macbook/Desktop/Bio-Research\ /PS050_cellranger_count_outs/filtered_feature_bc_matrix"
}else{
  biomarkers_dic <- biomarkers_dictionary_dataset2
  data_path <- "/Users/macbook/Desktop/Bio-Research\ /filtered_feature_bc_matrix"
}

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


##====Read Data====

library(Seurat)

Obj.data <- Read10X(data.dir = data_path)
Obj <- CreateSeuratObject(counts = Obj.data, project = "mSG")

## ==== Calculate Basic Info and Plot ==== 

# Add number of genes per UMI for each cell to metadata
Obj$log10GenesPerUMI <- log10(Obj$nFeature_RNA) / log10(Obj$nCount_RNA)

# Compute percent mito ratio
Obj$mitoRatio <- PercentageFeatureSet(object = Obj, pattern = "^mt-")
Obj$mitoRatio <- Obj@meta.data$mitoRatio / 100

# Compute percent Rps ratio
Obj$rpsRatio <- PercentageFeatureSet(object =Obj, pattern = "^Rps")
Obj$rpsRatio <- Obj@meta.data$rpsRatio / 100

# Compute percent Rpl ratio
Obj$rplRatio <- PercentageFeatureSet(object = Obj, pattern = "^Rpl")
Obj$rplRatio <- Obj@meta.data$rplRatio / 100

Obj@meta.data <- Obj@meta.data %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)

VlnPlot(object = Obj,features = c("nGene", "nUMI", "mitoRatio", "log10GenesPerUMI")) + theme(text = element_text(size = 12))
ggsave("dataset_info.png", h = 5000, w = 8000, units = "px")

VlnPlot(object = Obj,features = c("rpsRatio", "rplRatio", "mitoRatio")) + theme(text = element_text(size = 12))
ggsave("dataset_info_ribosome_mitoRatio.png", h = 5000, w = 8000, units = "px")

## ====== Filter cells =====
set.seed (12)
filtered_seurat <- subset(x =Obj, 
                          subset= (nUMI >= 1) & 
                            (nUMI <=20000) &
                            (nGene >= 100) & 
                            (log10GenesPerUMI > 0.5) &
                            (mitoRatio < 0.20))

rp_genes <- grep(pattern = "^Rps|^Rpl", x = rownames(x = Obj@assays$RNA@counts), value = TRUE)
Obj <- subset(Obj, features = rownames(Obj)[!rownames(Obj) %in% rp_genes])

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

# Normalize the counts, regress out mitoRatio, "rpsRatio", "rplRatio"
set.seed(12)
filtered_seurat_norm <- SCTransform(filtered_seurat, vars.to.regress = c("mitoRatio"))


# Identification of highly variable features (feature selection)
set.seed(12)
filtered_seurat_norm <- FindVariableFeatures(filtered_seurat_norm, selection.method = "vst", nfeatures = 2000)

## ====================== Check cell cycle=========================
#s.genes <- sapply(cc.genes$s.genes, convertString)
#g2m.genes <- sapply(cc.genes$g2m.genes, convertString)
#filtered_seurat_norm <- CellCycleScoring(filtered_seurat_norm, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

#Cell_Cycle_Assignment <- RunPCA(filtered_seurat_norm, features = c(s.genes, g2m.genes))
#DimPlot(Cell_Cycle_Assignment)
#ggsave("Cell_Cycle_Plot.png", h = 5000, w = 7000, units = "px")


#filtered_seurat_norm$CC.Difference <- filtered_seurat_norm$S.Score - filtered_seurat_norm$G2M.Score
#filtered_seurat_norm_regress_cell_cycle <- ScaleData(filtered_seurat_norm, vars.to.regress = "CC.Difference", features = rownames(filtered_seurat_norm))

#Cell_Cycle_Assignment <- RunPCA(filtered_seurat_norm_regress_cell_cycle, features = c(s.genes, g2m.genes))
#DimPlot(Cell_Cycle_Assignment)
#ggsave("Cell_Cycle_Plot_Regress_Out.png", h = 5000, w = 7000, units = "px")


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

filtered_data <- glue::glue("filtered_Rdata/filtered_seurat_norm_{resolution_value}.RData")
save(filtered_seurat_norm_UMAP, file=filtered_data)

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

# ===============================Visualization==================================
# Normalize RNA data for visualization purposes
set.seed(12)
filtered_seurat_norm_UMAP_RNA <- NormalizeData(filtered_seurat_norm_UMAP, verbose = FALSE)

# check basic info 
metrics <-  c("nUMI", "nGene", "mitoRatio", "rpsRatio", "rplRatio")
FeaturePlot(filtered_seurat_norm_UMAP_RNA, 
            reduction = "umap", 
            features = metrics,
            pt.size = 0.4, 
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE)
ggsave("scRNAseq_mSG_infomation.png", h = 5000, w = 7000, units = "px")


for (name in names(biomarkers_dic)){
  print(name)
  plot_genes(filtered_seurat_norm_UMAP_RNA, feature_list = biomarkers_dic[[name]], cell_name = name)
}

###################################################################################

# ========= Rename Indets =============

filtered_seurat_norm_UMAP_rename <- RenameIdents(filtered_seurat_norm_UMAP, 
                                          `1` = "basal", 
                                          `2` = "basal_w/ribogene", 
                                          `3` = "myoepithelial", 
                                          `4` = "Intermediate/Krt19high", 
                                          `5` = "ductal3", 
                                          `6` = "endothelial", 
                                          `7` = "immune_cell"
)

DimPlot(filtered_seurat_norm_UMAP_rename, pt.size = 3)+theme(text = element_text(size=30))
ggsave("UMAPPlot.png", h = 5000, w = 7000, units = "px")

filtered_data <- glue::glue("filtered_Rdata/filtered_seurat_norm_renamed.RData")
save(filtered_seurat_norm_UMAP_rename, file=filtered_data)

# ======== Replot with Name ======== 
set.seed(12)
filtered_seurat_norm_UMAP_rename_plot <- NormalizeData(filtered_seurat_norm_UMAP_rename, verbose = FALSE)

for (name in names(biomarkers_dic)){
  print(name)
  plot_genes(filtered_seurat_norm_UMAP_rename_plot, feature_list = biomarkers_dic[[name]], cell_name = name, labeled = TRUE)
}

# ========= Getting Stored Data ========
data_name<-glue::glue("{based_directory}/{id}_resolution_{resolution_value}/filtered_seurat_norm_renamed.RData")
filtered_seurat_norm_UMAP <- get(load(data_name))

