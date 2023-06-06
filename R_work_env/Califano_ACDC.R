  #Automated Community Detection of Cell Populations (ACDC) beta version from Califano lab

install.packages("devtools")
devtools::install_github("califano-lab/acdc")
library("acdc")

# 1.using a Simulated Annealing-based optimization (stochastic search). Retrieve solutions in a user-defined amount of time (2 min) and automatically update S.obj by storing the optimal clustering solution in the seurat_cluster slot of the metadata

settings <- list(max.time=180)
S.obj <- SAClustering(S.obj=filtered_seurat_norm,
                      res.range=c(0.1,1),
                      NN.range=c(3,15),
                      reduction=TRUE,
                      control=settings)

# 2.using a Grid Search approach (deterministic search) span a 10x10 grid of NN and resolution values and store the optimal clustering solution in S.obj

S.obj <- GridSearch_pcs_fast(object = S.obj,
                             assay.name="RNA",
                             .resolutions = seq(0.1, 1.9, by = 0.2),
                             .bootstraps = 1,
                             .knns = seq(11, 101, by = 10))

# 3.Store the clustering solution obtained with a given set of parameters, e.g. NN and res, and include it in S.obj

S.obj <- getFinal(S.obj=pbmc3k.final,
                  res=1,
                  NN=30,
                  reduction=TRUE)