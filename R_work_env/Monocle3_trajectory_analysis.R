devtools::install_github('cole-trapnell-lab/monocle3')
library(monocle3)
remotes::install_github('satijalab/seurat-wrappers')
library(SeuratWrappers)

install.packages ("R.utils")
library(R.utils)

#MONOCLE3 WORKFLOW --------------------
# 1.Convert Seurat object into an object of cell dataset (cds) class 
cds<- as.cell_data_set(S.obj)
cds
#to get cell metadata
colData(cds)
#to get gene data
fData(cds)
rownames(fData(cds)) [1:10]
fData(cds)$gene_short_name <- rownames(fData(cds)
#to get counts
counts(cds)

#2.Cluster cells using UMAP info from Seurat
#assign partitions
recreate.partition <-rep(1, lenght (cds@colData@rownames))
names(recreate.partition) <- cds@colData@rownames
recreate.partition <- as.factor(recreate.partition)

cds@clusters$UMAP$partitions <- recreate.partition

#Assign the cluster info
list_cluster <- S.obj@active.ident
cds@clusters$UMAP$clusters <- list_cluster

#Assign UMAP coordinates - cell embeddings

cds

#plot
 
setwd("E:/Columbia/Piera lab/R/Data/RNA velocity/PS050/filtered_feature_bc_matrix_PS050")
cds <- load_mm_data(mat_path = "E:/Columbia/Piera lab/R/Data/RNA velocity/PS050/filtered_feature_bc_matrix_PS050/matrix.mtx", 
                    feature_anno_path = "E:/Columbia/Piera lab/R/Data/RNA velocity/PS050/filtered_feature_bc_matrix_PS050/features.tsv", 
                    cell_anno_path = "E:/Columbia/Piera lab/R/Data/RNA velocity/PS050/filtered_feature_bc_matrix_PS050/barcodes.tsv")
# Converting Seurat obj into Monocle3 obj
cds <- as.cell_data_set(mSG.combined)        
### Warning: Monocle 3 trajectories require cluster partitions, which Seurat does not calculate. Please run 'cluster_cells' on your cell_data_set object
cluster_cells(cds)

# Get cell metadata
head(colData(cds))
# Get feature/gene metadata
fData(cds)
rownames(fData(cds))[1:10]
fData(cds)$gene_short_name <- rownames(fData(cds))
head(fData(cds))

# Get counts 
head(counts(cds))

# Retrieve clustering information from Surat object
## 1. Assign partitions
recreate.partitions <- c(rep(1, length(cds@colData@rownames)))
names(recreate.partitions) <- cds@colData@rownames
recreate.partitions <- as.factor(recreate.partitions)
recreate.partitions

cds@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partitions

## 2. Assign cluster information
list.cluster <- mSG.combined@active.ident
cds@clusters@listData[["UMAP"]][["clusters"]] <- list.cluster

 ##3. Assign UMAP coordinates
cds@int_colData@listData[["reducedDims"]]@listData[["UMAP"]] <- mSG.combined@reductions$umap@cell.embeddings

#Plot
cds <- learn_graph(cds)
cluster.before.traj <- plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster = F, group_label_size = 5)

cluster.before.traj

ggsave("scRNAseq_mSG_UMAP_0.5_trajectory.tiff", h = 5000, w = 6000, units = "px")

# Order cells in pseudotime
cds <- order_cells(cds, reduction_method = "UMAP", root_cells = colnames(cds[, clusters(cds) == 'Basal']))
plot_cells(cds, color_cells_by = "pseudotime", label_groups_by_cluster = T,
           label_branch_points = T, label_roots = F, label_leaves = F, graph_label_size=1.5)

head(pseudotime(cds), 10)

