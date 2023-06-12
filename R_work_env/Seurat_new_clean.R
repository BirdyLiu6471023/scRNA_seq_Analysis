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

setwd("E:/Columbia/Piera lab/R/Data/RNA velocity/PS050/analysis")

EpCAM.data <- Read10X(data.dir = "E:/Columbia/Piera lab/R/Data/RNA velocity/PS050/analysis/")

# Create Seurat object
EpCAM <- CreateSeuratObject(counts = EpCAM.data, project = "mSG_EpCAM", save.SNN = TRUE)

EpCAM

# Explore merged metadata
View(EpCAM@meta.data)

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

# Create .RData object to load at any time
save(EpCAM, file="E:/Columbia/Piera lab/R/Data/RNA velocity/PS050/analysis/EpCAM_merged.RData")

# Visualize the number of cell counts per sample
p7<- metadata %>% 
  ggplot(aes(x=sample)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")

# Visualize the number UMIs/transcripts per cell. Usually more than 10,000 UMIs are cell doublets
p8 <- metadata %>% 
  ggplot(aes(x=nUMI)) + 
  geom_density(alpha = 0.2, color = "deepskyblue", fill = "grey") + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 10, color = "tomato", lwd = 1.5) +
  geom_vline(xintercept = 20000, color = "tomato", lwd = 1.5)

# Visualize the distribution of genes detected per cell via histogram
p9<- metadata %>% 
  ggplot(aes(x=nGene)) + 
  geom_density(alpha = 0.2, color = "deepskyblue", fill = "grey") + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 100, color = "tomato", lwd = 1.5)
# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
p10<- metadata %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 100) +
  facet_wrap(~sample)

# Visualize the distribution of mitochondrial gene expression detected per cell. We define poor quality samples for mitochondrial counts as cells which surpass the 0.2 mitochondrial ratio mark
p11<- metadata %>% 
  ggplot(aes(color=sample, x=mitoRatio)) + 
  geom_density(alpha = 0.2, color = "deepskyblue", fill = "grey") + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.2, color = "tomato", lwd = 1.5)

# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI 
p12<- metadata %>%
  ggplot(aes(x=log10GenesPerUMI)) +
  geom_density(alpha = 0.2, color = "deepskyblue", fill = "grey") +
  theme_classic() +
  geom_vline(xintercept = 0.8, color = "tomato", lwd = 1.5)
p7 + p8 + p9 + p10 + p11 + p12 +
  plot_layout(ncol = 3)
ggsave("scRNAseq_mSG_not_filtered.tiff", h = 5000, w = 7000, units = "px")

# Cell-level FILTERING
set.seed (12)
filtered_seurat <- subset(x = EpCAM, 
                          subset= (nUMI >= 1) & 
                            (nUMI <=20000) &
                            (nGene >= 100) & 
                            (log10GenesPerUMI > 0.5) &
                            (log10GenesPerUMI < 0.9) &
                            (mitoRatio < 0.20))

View(filtered_seurat@meta.data)

# Output a logical vector for every gene on whether the more than zero counts per cell
# Extract counts
counts <- GetAssayData(object = filtered_seurat, slot = "counts")

# Output a logical vector for every gene on whether the more than zero counts per cell
nonzero <- counts > 0

# Sums all TRUE values and returns TRUE if more than 1 TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= 1

# Only keeping those genes expressed in more than 1 cell (include rare populations e.g. myoepithelial)
filtered_counts <- counts[keep_genes, ]

# Reassign to filtered Seurat object
filtered_seurat <- CreateSeuratObject(filtered_counts, meta.data = filtered_seurat@meta.data)

# Save filtered subset to new metadata
metadata_clean <- filtered_seurat@meta.data

# Create .RData object to load at any time
save(filtered_seurat, file="E:/Columbia/Piera lab/R/Data/RNA velocity/PS050/analysis/filtered_seurat.RData")


# Visualize the number of cell counts per sample
p1 <- 
  metadata_clean %>% 
  ggplot(aes(x=sample)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")

# Visualize the number UMIs/transcripts per cell. Usually more than 10,000 UMIs are cell doublets
p2 <- 
  metadata_clean %>% 
  ggplot(aes(x=nUMI)) + 
  geom_density(alpha = 0.2, color = "deepskyblue", fill = "grey") + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 1, color = "tomato", lwd = 1.5) +
  geom_vline(xintercept = 20000, color = "tomato", lwd = 1.5)

# Visualize the distribution of genes detected per cell via histogram
p3 <-
  metadata_clean %>% 
  ggplot(aes(x=nGene)) + 
  geom_density(alpha = 0.2, color = "deepskyblue", fill = "grey") + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 100, color = "tomato", lwd = 1.5)

# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
p4 <-
  metadata_clean %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 1000) +
  geom_hline(yintercept = 300) +
  facet_wrap(~sample)

# Visualize the distribution of mitochondrial gene expression detected per cell. We define poor quality samples for mitochondrial counts as cells which surpass the 0.2 mitochondrial ratio mark
p5 <-
  metadata_clean %>% 
  ggplot(aes(color=sample, x=mitoRatio)) + 
  geom_density(alpha = 0.2, color = "deepskyblue", fill = "grey") + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.2, color = "tomato", lwd = 1.5)

# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI 
p6 <-
  metadata_clean %>%
  ggplot(aes(x=log10GenesPerUMI)) +
  geom_density(alpha = 0.2, color = "deepskyblue", fill = "grey") +
  theme_classic() +
  geom_vline(xintercept = 0.5, color = "tomato", lwd = 1.5) +
  geom_vline(xintercept = 0.9, color = "tomato", lwd = 1.5)

#Patchwork assemble
p1 + p2 + p3 + p4 + p5 + p6 +
  plot_layout(ncol = 3)
ggsave("scRNAseq_mSG_filtered.tiff", h = 5000, w = 7000, units = "px")

filtered_seurat[["percent.mt"]] <- PercentageFeatureSet(filtered_seurat, pattern = "^mt-")
plot1 <- FeatureScatter(filtered_seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(filtered_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
ggsave("scRNAseq_mSG_filtered_counts_features_mt.tiff", h = 2000, w = 4000, units = "px")

# Normalize the counts
set.seed(12)
filtered_seurat_norm <- SCTransform(filtered_seurat, vars.to.regress = "mitoRatio")

# Create .RData object to load at any time
save(filtered_seurat_norm, file="E:/Columbia/Piera lab/R/Data/RNA velocity/PS050/analysis/filtered_seurat_norm.RData")

# Check which assays are stored in objects
filtered_seurat_norm@assays

# Identification of highly variable features (feature selection)
filtered_seurat_norm <- FindVariableFeatures(filtered_seurat_norm, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(filtered_seurat_norm), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(filtered_seurat_norm)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

ggsave("scRNAseq_mSG_top10_genes.tiff", h = 5000, w = 7000, units = "px")

# Run PCA
set.seed (12)
filtered_seurat_norm_PCA <- RunPCA(object = filtered_seurat_norm)
# Plot PCA
PCAPlot(filtered_seurat_norm_PCA)  

# Run UMAP
set.seed(12)
filtered_seurat_norm_UMAP <- RunUMAP(filtered_seurat_norm_PCA, 
                                dims = 1:15,
                                reduction = "pca")
# Plot UMAP                             
DimPlot(filtered_seurat_norm_UMAP)

# Run t-SNE
set.seed(12)
filtered_seurat_norm_tSNE <- RunTSNE(
  filtered_seurat_norm_PCA, reduction = "pca",
  assay = NULL,
  seed.use = 12,
  tsne.method = "Rtsne",
  dim.embed = 2,
  dims = 1:15,
  reduction.key = "tSNE_")

# Plot t-SNE
TSNEPlot(filtered_seurat_norm_tSNE)

# CLUSTERING CELLS BASED ON TOP PCs (METAGENES)
# Identify significant PCs
DimHeatmap(filtered_seurat_norm_UMAP, 
           dims = 1:12, 
           cells = 500, 
           balanced = TRUE)

# Printing out the most variable genes driving PCs
print(x = filtered_seurat_norm_UMAP[["pca"]], 
      dims = 1:12, 
      nfeatures = 15)
VizDimLoadings(filtered_seurat_norm_UMAP, dims = 1:3, reduction = "pca")


# Elbow plot: threshold for identifying the majority of the variation. However, this method can be quite subjective.
ElbowPlot(object = filtered_seurat_norm_UMAP, 
          ndims = 40)


# CLUSTER THE CELLS
# Determine the K-nearest neighbor graph 
set.seed (12)
filtered_seurat_norm_UMAP <- FindNeighbors(object = filtered_seurat_norm_UMAP, reduction = "pca", dims = 1:15, verbose = TRUE, graph.name = "mSG")

# Determine the clusters with Leiden algorithm                          
set.seed (12)
filtered_seurat_norm_UMAP_0.5 <- FindClusters(object = filtered_seurat_norm_UMAP, verbose = TRUE, algorithm = 4, resolution = 0.5, graph.name = "mSG")
set.seed (12)
filtered_seurat_norm_UMAP_0.45 <- FindClusters(object = filtered_seurat_norm_UMAP, verbose = TRUE, algorithm = 4, resolution = 0.45, graph.name = "mSG")

# Explore resolutions
filtered_seurat_norm_UMAP_0.5@meta.data %>% 
  View()

# Subclustering (if needed)
basal_subcluster <- FindSubCluster(filtered_seurat_norm_UMAP_0.5, cluster = "5", graph.name ="mSG", subcluster.name = "Acinar_new",  resolution = 0.2, algorithm = 4)

DimPlot(acinar_subcluster,
        reduction = "umap", group.by = "Acinar_new", label = TRUE, label.size = 6)

# Clustering quality using SILHOUETTE plot
#Compute distance matrix to UMAP coordinates
set.seed(12)
distance_matrix <- dist(Embeddings(filtered_seurat_norm_UMAP_0.5[['umap']])[, 1:2])

#Isolate cluster name
clusters <- filtered_seurat_norm_UMAP_0.5@active.ident

#Compute silhouette score for each cluster
silhouette <- silhouette(as.numeric(clusters), dist = distance_matrix)
filtered_seurat_norm_UMAP_0.5@meta.data$silhouette_score <- silhouette[,3]

#Plot
fviz_silhouette(silhouette, label = FALSE, print.summary = TRUE)

ggsave("scRNAseq_mSG_filtered_silhouette_0.5_dims_15.tiff", h = 2000, w = 4000, units = "px")


# Determine metrics to plot present in seurat_integrated@meta.data
metrics <-  c("nUMI", "nGene", "mitoRatio")

FeaturePlot(filtered_seurat_norm_UMAP_0.5, 
            reduction = "umap", 
            features = metrics,
            pt.size = 0.4, 
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE)
            
ggsave("scRNAseq_mSG_feature_plot_0.5.tiff", h = 5000, w = 7000, units = "px")

# Defining the information in the seurat object of interest
columns <- c(paste0("PC_", 1:16),
             "ident",
             "UMAP_1", "UMAP_2")

# Extracting this data from the seurat object
pc_data <- FetchData(filtered_seurat_norm_UMAP_0.5, 
                     vars = columns)

# Extract the UMAP coordinates for the first 10 cells
filtered_seurat_norm_UMAP_0.5@reductions$umap@cell.embeddings[1:10, 1:2]


# Adding cluster label to center of cluster on UMAP
umap_label <- FetchData(filtered_seurat_norm_UMAP_0.5, 
                        vars = c("ident", "UMAP_1", "UMAP_2"))  %>%
  group_by(ident) %>%
  summarise(x=mean(UMAP_1), y=mean(UMAP_2))

# Plotting a UMAP plot for each of the PCs
map(paste0("PC_", 1:16), function(pc){
  ggplot(pc_data, 
         aes(UMAP_1, UMAP_2)) +
    geom_point(aes_string(color=pc), 
               alpha = 0.7) +
    scale_color_gradient(guide = FALSE, 
                         low = "grey90", 
                         high = "blue")  +
    geom_text(data=umap_label, 
              aes(label=ident, x, y)) +
    ggtitle(pc)
}) %>% 
  plot_grid(plotlist = .)

ggsave ("scRNAseq_mSG_UMAP_plot_PC_0.5.tiff", h = 5000, w = 6500, units = "px")

# Examine PCA results 
print(filtered_seurat_norm_UMAP_0.5["pca"], dims = 1:5, nfeatures = 10)

# Select the RNA counts slot to be the default assay
DefaultAssay(filtered_seurat_norm_UMAP_0.5) <- "RNA"

# Normalize RNA data for visualization purposes
filtered_seurat_norm_UMAP_0.5_RNA <- NormalizeData(filtered_seurat_norm_UMAP_0.5, verbose = FALSE)

# Visualize BASAL+ACINAR cluster-specific markers expression with UMAP
FeaturePlot(filtered_seurat_norm_UMAP_0.5_RNA, 
            reduction = "umap", 
            features = c("Actb", "Krt15", "Bhlha15", "Aqp5", "Pip", "Smgc", "Scgb2b27", "Prlr", "Mki67", "Lpo", "Acta2"), 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE,
            repel = TRUE)
ggsave("scRNAseq_mSG_UMAP_acinar_PC_0.5.tiff", h = 5000, w = 8000, units = "px")

# Visualize BASAL+ACINAR cluster-specific markers expression with violin plots
VlnPlot(filtered_seurat_norm_UMAP_0.5_RNA, features = c("Krt15", "Bhlha15", "Aqp5", "Pip", "Smgc", "Scgb2b27", "Prlr", "Mki67", "Lpo", "Acta2"))

ggsave("scRNAseq_mSG_vln_acinar_PC_0.5.tiff", h = 5000, w = 8000, units = "px")

# Visualize DUCTAL cluster-specific markers expression with UMAP
FeaturePlot(filtered_seurat_norm_UMAP_0.5_RNA, 
            reduction = "umap", 
            features = c("Krt7", "Klk1", "Ngf", "Bsnd", "Ascl3", "Kit", "Elf5", "Clic6"), 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)
ggsave("scRNAseq_mSG_UMAP_ductal_PC_0.5.tiff", h = 5000, w = 7000, units = "px")

# Visualize DUCTAL cluster-specific markers expression with violin plots
VlnPlot(filtered_seurat_norm_UMAP_0.5_RNA, features = c("Krt7", "Klk1", "Ngf", "Bsnd", "Ascl3", "Kit", "Elf5", "Clic6"))

ggsave("scRNAseq_mSG_vln_ductal_PC_0.5.tiff", h = 5000, w = 7000, units = "px")


#identification of each cluster-specific marker

set.seed (12)
cluster1.markers <- FindMarkers(filtered_seurat_norm_UMAP_0.5, ident.1 = 1, logfc.threshold = 0.25)

write.csv2 (cluster1.markers, file = "cluster1_markers_0.5.csv")

set.seed (12)
cluster2.markers <- FindMarkers(filtered_seurat_norm_UMAP_0.5, ident.1 = 2, logfc.threshold = 0.25)

write.csv2 (cluster2.markers, file = "cluster2_markers_0.5.csv")

set.seed (12)
cluster3.markers <- FindMarkers(filtered_seurat_norm_UMAP_0.5, ident.1 = 3, logfc.threshold = 0.25)

write.csv2 (cluster3.markers, file = "cluster3_markers_0.5.csv")

set.seed (12)
cluster4.markers <- FindMarkers(filtered_seurat_norm_UMAP_0.5, ident.1 = 4, logfc.threshold = 0.25)

write.csv2 (cluster4.markers, file = "cluster4_markers_0.5.csv")

set.seed (12)
cluster5.markers <- FindMarkers(filtered_seurat_norm_UMAP_0.5, ident.1 = 5, logfc.threshold = 0.25)

write.csv2 (cluster5.markers, file = "cluster5_markers_0.5.csv")

set.seed (12)
cluster6.markers <- FindMarkers(filtered_seurat_norm_UMAP_0.5, ident.1 = 6, logfc.threshold = 0.25)

write.csv2 (cluster6.markers, file = "cluster6_markers_0.5.csv")

set.seed (12)
cluster7.markers <- FindMarkers(filtered_seurat_norm_UMAP_0.5, ident.1 = 7, logfc.threshold = 0.25)

write.csv2 (cluster7.markers, file = "cluster7_markers_0.5.csv")

set.seed (12)
cluster8.markers <- FindMarkers(filtered_seurat_norm_UMAP_0.5, ident.1 = 8, logfc.threshold = 0.25)

write.csv2 (cluster8.markers, file = "cluster8_markers_0.5.csv")

set.seed (12)
cluster9.markers <- FindMarkers(filtered_seurat_norm_UMAP_0.5, ident.1 = 9, logfc.threshold = 0.1)
write.csv2 (cluster9.markers, file = "cluster9_markers_0.5.csv")

set.seed (12)
cluster10.markers <- FindMarkers(filtered_seurat_norm_UMAP_0.5, ident.1 = 10, logfc.threshold = 0.1)
write.csv2 (cluster10.markers, file = "cluster10_markers_0.5.csv")


# Assign name to clusters
mSG.combined <- RenameIdents(filtered_seurat_norm_UMAP_0.5, `1` = "Ductal 2", `2` = "Ductal 1_1", `3` = "Acinar 2", `4` = "Ductal 3_1", `5` = "Junk", `6` = "Ductal 3_2", `7` = "Acinar 1", `8` = "Basal", `9` = "Ductal 1_2", `10` = "Mesenchymal")
save(mSG.combined, file="E:/Columbia/Piera lab/R/Data/RNA velocity/PS050/analysis/filtered_seurat_norm.RData")

# Find differentially expressed markers among two specific clusters
FindMarkers(filtered_seurat_norm_UMAP_0.5, ident.1 = "4", ident.2 = "6")

# Find differentially expressed markers among one cluster and all the others
FindMarkers(filtered_seurat_norm_UMAP_0.5, ident.1 = "10", ident.2 = NULL, only.pos = TRUE)

# Plot UMAP with cluster names
DimPlot(mSG.combined, label = TRUE, label.size = 8)
ggsave("scRNAseq_mSG_UMAP_0.5_annotated.tiff", h = 5000, w = 6000, units = "px")

markers.to.plot <- c("Tspan33", "Gpc6", "Slc13a2", "Slc9a3", "Esp4", "Gjb2", "Car2", "Ceacam1", "Elf5", "Kit", "Trp53i11", "Slc4a11", "Stap1", "Hepacam2", "Ascl3", "Nrip3", "Foxi1", "Rcan2", "Aldh1a3", "Clcnkb", "Igfbp5", "Esrrg", "Ccdc141", "Ngf", "Scgb1c1", "Adamts2", "Slc2a4", "Krt14", "Krt5", "Ptn", "Col4a1", "Robo2", "Il18r1", "Pdpn", "Clca3a1", "Wnt10b", "Tg", "St3gal1", "Pigr", "Cd44", "Frmpd1os", "Kcnn4", "Folr1", "Agt", "Oit1", "Lpo", "Snhg18", "Esp8", "Etv1", "Gjb1", "Creb3l1", "Slc12a8", "Bhlha15", "Snca", "Lman1l", "Prlr", "Gdpd1", "Serpinb11", "Gm44410", "Tll1", "Tmem59l", "Dcdc2a", "Gm26917", "Lars2", "Kctd12" , "Mafb", "Cracr2a")

DotPlot(mSG.combined, features = rev(markers.to.plot), cols = c("deepskyblue", "coral1"), dot.scale = 10) + FontSize(x.text = 15, y.text = 15, x.title =25, y.title = 25) + RotatedAxis()

ggsave("scRNAseq_mSG_dotplot_genes_0.5.tiff", h = 5000, w = 8000, units = "px")

