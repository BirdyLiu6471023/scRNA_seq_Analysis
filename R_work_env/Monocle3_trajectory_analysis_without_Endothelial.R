#devtools::install_github('cole-trapnell-lab/monocle3')
library(monocle3)
#remotes::install_github('satijalab/seurat-wrappers')
library(SeuratWrappers)
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


#MONOCLE3 WORKFLOW --------------------

# loading the Rdata: filtered_seurat_norm.RData
setwd("/Users/macbook/Desktop/SGcell_Evolution/R_work_env/sl06202023_20000_resolution_0.5")
mSG.combined <- get(load("/Users/macbook/Desktop/SGcell_Evolution/R_work_env/sl06202023_20000_resolution_0.5/filtered_Rdata/filtered_seurat_norm_without_Endothelial.RData"))
resolution_value = 0.5

previous_cluster_name <- mSG.combined@active.ident
Idents(mSG.combined) <- "new_identity"

# =======================Remove Previous PCA and UMAP results ==================
# Clear PCA results
mSG.combined[["pca"]] <- NULL

# Clear UMAP results
mSG.combined[["umap"]] <- NULL

# Clear clustering results
mSG.combined[["seurat_clusters"]] <- NULL

#=============================Run PCA===========================================
# Run PCA
set.seed (12)
mSG.combined <- RunPCA(object = mSG.combined)
# Plot PCA
PCAPlot(mSG.combined, pt.size = 3)+theme(text = element_text(size=30))
ggsave("PCAPlot_without_Endothelial2.png", h = 5000, w = 7000, units = "px")

#=============================Run UMAP==========================================
# Run UMAP
set.seed(12)
mSG.combined <- RunUMAP(mSG.combined, 
                                     dims = 1:15,
                                     reduction = "pca")
# Plot UMAP                             
DimPlot(mSG.combined, pt.size = 3)+theme(text = element_text(size=30))
ggsave("UMAPPlot_without_Endothelial2.png", h = 5000, w = 7000, units = "px")

# ==========================Clustering =========================================
set.seed (12)
mSG.combined <- FindNeighbors(object = mSG.combined, reduction = "pca", dims = 1:15, verbose = TRUE, graph.name = "mSG")

set.seed (12)
mSG.combined <- FindClusters(object = mSG.combined, verbose = TRUE, algorithm = 4, resolution = resolution_value, graph.name = "mSG")


set.seed(12)
embed_matrix <- Embeddings(mSG.combined[['umap']])
embed_matrix <- as.matrix(embed_matrix)
distance_matrix <- dist(embed_matrix)

#Isolate cluster name
clusters <- mSG.combined@active.ident

#Compute silhouette score for each cluster
silhouette <- silhouette(as.numeric(clusters), dist = distance_matrix)
mSG.combined@meta.data$silhouette_score <- silhouette[,3]
fviz_silhouette(silhouette, label = FALSE, print.summary = TRUE)
filename<- glue::glue("silhouette_score_{resolution_value}_without_Endothelial.png")
ggsave(filename, h = 2000, w = 4000, units = "px")


# =======================Replot UMAP ===========================================
DimPlot(mSG.combined, pt.size = 3)+theme(text = element_text(size=30))
ggsave("UMAPPlot_without_Endothelial2_reclustering.png", h = 5000, w = 7000, units = "px")

mSG.combined@active.ident<-previous_cluster_name
DimPlot(mSG.combined, pt.size = 3)+theme(text = element_text(size=30))
ggsave("UMAPPlot_without_Endothelial2_previous_clusters.png", h = 5000, w = 7000, units = "px")

# ======================pseudotime analysis=====================================

# Converting Seurat obj into Monocle3 obj
cds <- as.cell_data_set(mSG.combined)        
### Warning: Monocle 3 trajectories require cluster partitions, which Seurat does not calculate. Please run 'cluster_cells' on your cell_data_set object
#cds <- cluster_cells(cds)
#plot_cells(cds, color_cells_by = "partition")

# Assign gene name
fData(cds)$gene_short_name <- rownames(fData(cds))

# Retrieve clustering information from Surat object
## 1. Assign partitions
recreate.partitions <- c(rep(1, length(cds@colData@rownames)))
names(recreate.partitions) <- cds@colData@rownames
recreate.partitions <- as.factor(recreate.partitions)
cds@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partitions

## 2. Assign cluster information
list.cluster <- mSG.combined@active.ident
cds@clusters@listData[["UMAP"]][["clusters"]] <- list.cluster

 ##3. Assign UMAP coordinates
cds@int_colData@listData[["reducedDims"]]@listData[["UMAP"]] <- mSG.combined@reductions$umap@cell.embeddings

#Plot
set.seed (12)
cds <- learn_graph(cds)
cluster.before.traj <- plot_cells(cds, color_cells_by = "ident", label_groups_by_cluster = F, label_leaves=F, group_label_size = 10, cell_size = 1.5)
cluster.before.traj
ggsave("scRNAseq_mSG.combined_trajectory_filtered.png", h = 5000, w = 6000, units = "px")

#============================Order cells in pseudotime==========================

set.seed (12)
cds <- order_cells(cds, reduction_method = "UMAP", root_cells = colnames(cds[, cds@clusters@listData[["UMAP"]][["clusters"]] == 'Basal']))
plot_cells(cds, color_cells_by = "pseudotime", group_cells_by = "cluster", label_groups_by_cluster = F, label_branch_points = F, label_roots = F, label_leaves = T, graph_label_size=5, cell_size = 2, trajectory_graph_color = "grey")
ggsave("scRNAseq_mSG.combined_pseudotime_filtered.png", h = 5000, w = 6000, units = "px")


#===========plot the aggregate module scores within each group of cell==========

# Monocle’s graph_test() function detects genes that vary over a trajectory. This may run very slowly. Adjust the number of cores as needed
set.seed (12)
cds_graph_test_results <- graph_test(cds,
                                     neighbor_graph = "principal_graph",
                                     cores = 4)

#rowData(cds)$gene_short_name <- row.names(rowData(cds))
# head(cds_graph_test_results, error=FALSE, message=FALSE, warning=FALSE)

deg_ids <- rownames(subset(cds_graph_test_results[order(cds_graph_test_results$morans_I, decreasing = TRUE),], q_value < 0.05))
plot_cells(cds,
           genes=head(deg_ids),
           show_trajectory_graph = FALSE,
           label_cell_groups = FALSE,
           label_leaves = FALSE,
           cell_size = 1.5) + theme(text = element_text (size = 40)) + theme(legend.key.size = unit(2, 'cm'))
ggsave("scRNAseq_mSG.combined_genes_trajectory_filtered.png", h = 5000, w = 6000, units = "px")

# colData(cds)$assigned_cell_type <- as.character(partitions(cds))

# before dimension reduction using UMAP or T-SNE or PCA, preprocess_cds 
set.seed (12)
cds <- preprocess_cds(cds, num_dim = 100)

set.seed(12)
gene_module_df <- find_gene_modules(cds[deg_ids,], reduction_method = "UMAP", resolution=c(10^seq(-6,-1)))
write.csv2 (gene_module_df, file = "gene_module_df.csv")

cell_group_df <- tibble::tibble(cell=row.names(colData(cds)), 
                                cell_group=list.cluster)

set.seed(12)
agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)

row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))

heatmap <- pheatmap::pheatmap(agg_mat, scale="column", clustering_method="ward.D2", fontsize = 18)
ggsave("gene_module_filtered.png", plot = heatmap, h = 6000, w = 4000, units = "px")

#===================KEGG pathway over-representation analysis ==================

# https://yulab-smu.top/biomedical-knowledge-mining-book/clusterprofiler-kegg.html
library(clusterProfiler)
library(org.Mm.eg.db)

# create a function to find the pathway
# mod: the number of module you are interested in 
# exclude: the gene name you find that not in database, the error would tell you which gene id to exclude
pathway_analysis <- function(mod, exclude = NULL){
  ## check if database exists: 
  if (!exists("database")){
    database <<- as.list(org.Mm.egALIAS2EG) 
    print("loading database...")
  }
  significant_module <- subset(gene_module_df, module == mod)
  significant_genes <-significant_module$id
  if (!is.null(exclude)){
    significant_genes <- significant_genes[-which(significant_genes %in% exclude)]
  }
  significant_ids <- unlist(unname(database[significant_genes]))
  kk <- enrichKEGG(gene         = significant_ids,
                   organism     = 'mmu',
                   pvalueCutoff = 0.05)
  return(kk)
}

num_module <- max(as.numeric(gene_module_df$module))
for (i in 1:num_module) {
  analysis <- pathway_analysis(i)
  filtered_data <- glue::glue('gene_module_{i}.RData')
  save(analysis, file=filtered_data)
}


