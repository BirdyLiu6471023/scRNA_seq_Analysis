devtools::install_github('cole-trapnell-lab/monocle3')
library(monocle3)
remotes::install_github('satijalab/seurat-wrappers')
library(SeuratWrappers)
)

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
cluster.before.traj <- plot_cells(cds, color_cells_by = "ident", label_groups_by_cluster = F, label_leaves=F, group_label_size = 5)

cluster.before.traj

ggsave("scRNAseq_mSG.combined_trajectory.tiff", h = 5000, w = 6000, units = "px")

# Order cells in pseudotime
cds <- order_cells(cds, reduction_method = "UMAP", root_cells = colnames(cds[, clusters(cds) == 'Basal']))

plot_cells(cds, color_cells_by = "pseudotime", group_cells_by = "cluster", label_groups_by_cluster = T, label_branch_points = T, label_roots = F, label_leaves = T, graph_label_size=1.5, trajectory_graph_color = "grey", trajectory_graph_segment_size = 3, cell_size = 1.5) + theme(axis.text.x=element_text(size = 40)) + theme(axis.text.y=element_text(size = 40)) + theme(axis.title.x= element_text (size=40)) + theme(axis.title.y= element_text (size=40)) + theme(legend.title = element_text(size = 30)) + theme(legend.text = element_text(size = 30)) + theme(legend.key.size = unit(2, 'cm'))

ggsave("scRNAseq_mSG.combined_pseudotime.tiff", h = 5000, w = 6000, units = "px")

head(pseudotime(cds), 10)

# Monocle’s graph_test() function detects genes that vary over a trajectory. This may run very slowly. Adjust the number of cores as needed

cds_graph_test_results <- graph_test(cds,
                                     neighbor_graph = "principal_graph",
                                     cores = 4)

rowData(cds)$gene_short_name <- row.names(rowData(cds))

head(cds_graph_test_results, error=FALSE, message=FALSE, warning=FALSE)

deg_ids <- rownames(subset(cds_graph_test_results[order(cds_graph_test_results$morans_I, decreasing = TRUE),], q_value < 0.05))

plot_cells(cds,
           genes=head(deg_ids),
           show_trajectory_graph = FALSE,
           label_cell_groups = FALSE,
           label_leaves = FALSE,
           cell_size = 1.5) + theme(text = element_text (size = 40)) + theme(legend.key.size = unit(2, 'cm'))


ggsave("scRNAseq_mSG.combined_genes_trajectory.tiff", h = 5000, w = 6000, units = "px")

colData(cds)$assigned_cell_type <- as.character(partitions(cds))
colData(cds)$assigned_cell_type <- dplyr::recode(colData(cds)$assigned_cell_type,
                                                '1' = "Ductal 2",
                                                '2' = "Ductal 1_1",
                                                '3' = "Acinar 2",
                                                '4' = "Ductal 3_1",
                                                '5' = "Junk",
                                                '6' = "Ductal 3_2",
                                                '7' = "Acinar 1",
                                                '8' = "Basal",
                                                '9' = "Ductal 1_2",
                                                '10' = "Mesenchymal")

gene_module_df <- find_gene_modules(cds[deg_ids,], resolution=c(10^seq(-6,-1)))

cell_group_df <- tibble::tibble(cell=row.names(colData(cds)), 
                                cell_group=colData(cds)$cell.type)

agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)

row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))

pheatmap::pheatmap(agg_mat,
                   scale="column", clustering_method="ward.D2")
