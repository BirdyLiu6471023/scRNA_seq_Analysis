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


cluster_name <- "ductal2"

cluster_dataname <- glue::glue("/Users/macbook/Desktop/SGcell_Evolution/R_work_env/sl06202023_20000_resolution_0.5/filtered_Rdata/{cluster_name}_filtered_seurat_norm.RData")


# ============================ cluster inspecting ============================

setwd("/Users/macbook/Desktop/SGcell_Evolution/R_work_env/sl06202023_20000_resolution_0.5")
filtered_seurat_norm_UMAP_ <- get(load(cluster_dataname))

metadata_cluster = filtered_seurat_norm_UMAP_@meta.data
# Visualize the number of cell counts per sample
p1 <- 
  metadata_cluster %>% 
  ggplot(aes(x=sample)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")

# Visualize the number UMIs/transcripts per cell. Usually more than 10,000 UMIs are cell doublets
p2 <- 
  metadata_cluster %>% 
  ggplot(aes(x=nUMI)) + 
  geom_density(alpha = 0.2, color = "deepskyblue", fill = "grey") + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 1, color = "tomato", lwd = 1.5) +
  geom_vline(xintercept = 20000, color = "tomato", lwd = 1.5)

# Visualize the distribution of genes detected per cell via histogram
p3 <-
  metadata_cluster %>% 
  ggplot(aes(x=nGene)) + 
  geom_density(alpha = 0.2, color = "deepskyblue", fill = "grey") + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 100, color = "tomato", lwd = 1.5)

# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
p4 <-
  metadata_cluster %>% 
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
  metadata_cluster %>% 
  ggplot(aes(color=sample, x=mitoRatio)) + 
  geom_density(alpha = 0.2, color = "deepskyblue", fill = "grey") + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.2, color = "tomato", lwd = 1.5)

# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI 
p6 <-
  metadata_cluster %>%
  ggplot(aes(x=log10GenesPerUMI)) +
  geom_density(alpha = 0.2, color = "deepskyblue", fill = "grey") +
  theme_classic() +
  geom_vline(xintercept = 0.5, color = "tomato", lwd = 1.5)

#Patchwork assemble
p1 + p2 + p3 + p4 + p5 + p6 +
  plot_layout(ncol = 3)

plot_name <- glue::glue("scRNAseq_mSG_filtered_{cluster_name}.png")
ggsave(plot_name, h = 5000, w = 7000, units = "px")





