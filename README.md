# scRNA_seq_Analysis

resolution 0.5, upper bound of nUMI = 20000
Quality Control->PCA->UMAP->Clustering->Pseudotime Analysis 

## Quality Plots for Each Cluster:

**Cluster ?** \
For Cluster ?, log10GenesPerUMI is in the normal range, 0.8+, so are other filter conditions, which are all in the normal range. 

<img src="R_work_env/sl06202023_20000_resolution_0.5/cluster_quality_plot/scRNAseq_mSG_filtered_quest.png" width="700">

**Cluster Basal** \
<img src="R_work_env/sl06202023_20000_resolution_0.5/cluster_quality_plot/scRNAseq_mSG_filtered_basal.png" width="700">


**Cluster Acinar1** \
As we discussed before, the log10GenesPerUMI for cluster Acinar1 is less than 0.8 (the threshold for other clusters) , and there two peaks shown in plot.  
<img src="R_work_env/sl06202023_20000_resolution_0.5/cluster_quality_plot/scRNAseq_mSG_filtered_acinar1.png" width="700">

**Cluster Acinar2** \
<img src="R_work_env/sl06202023_20000_resolution_0.5/cluster_quality_plot/scRNAseq_mSG_filtered_acinar2.png" width="700">

**Cluster Ductal1** \
<img src="R_work_env/sl06202023_20000_resolution_0.5/cluster_quality_plot/scRNAseq_mSG_filtered_ductal1.png" width="700">

**Cluster Ductal2** \
<img src="R_work_env/sl06202023_20000_resolution_0.5/cluster_quality_plot/scRNAseq_mSG_filtered_ductal2.png" width="700">

**Cluster Ductal3_1** \
<img src="R_work_env/sl06202023_20000_resolution_0.5/cluster_quality_plot/scRNAseq_mSG_filtered_ductal3_1.png" width="700">

**Cluster Ductal3_2** \
<img src="R_work_env/sl06202023_20000_resolution_0.5/cluster_quality_plot/scRNAseq_mSG_filtered_ductal3_2.png" width="700">

**Cluster Endothelial** \
<img src="R_work_env/sl06202023_20000_resolution_0.5/cluster_quality_plot/scRNAseq_mSG_filtered_endothelial.png" width="700">

## Pseudotime after filtering out Endothelia cells 
<p float="left">
  <img src="R_work_env/sl06202023_20000_resolution_0.5/scRNAseq_mSG.combined_trajectory_filtered.png" width="700">
  <img src="R_work_env/sl06202023_20000_resolution_0.5/scRNAseq_mSG.combined_pseudotime_filtered.png" width="700" /> 
</p>




