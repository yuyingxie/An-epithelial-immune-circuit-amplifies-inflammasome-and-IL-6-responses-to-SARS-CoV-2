library(Seurat)
library(scales)
library(scater)
library(muscat)
library(scran) 
library(sctransform)

# Check the inflammasome in C1, C9, C11, C15, C16, C20, C28, C30
# Inflammasome gene list: GSDMD, NLRP3, NLRP1, PYCARD, CASP1, CASP4, CASP5, IL1B, NLRC4, AIM2, IL18, 
inflam_g = c('GSDMD', 'NLRP3', 'NLRP1', 'PYCARD', 'CASP1', 
             'CASP4', 'CASP5', 'IL1B',   'NLRC4', 'AIM2', 
             'IL18', 'NAIP')

cluster_id = c(0, 8, 10, 14, 15, 19, 27, 29)
T.integrated = readRDS('T_SCT_integrated_11_26_2020.rds')

T.integrated$Sample_f = factor(T.integrated$Sample, 
                    levels = c('L_H1', 'L_H2', 'L_H3', 'L_H4', 'A_H1',
                              'A_H2', 'A_H3', 'A_H4', 'A_H5', 'W_H1',
                              'W_H2', 'W_H3', 'W_H4', 'W_H5', 'W_H6',
                              'L_M2', 'L_M4', 'L_M5', 'L_M7_2', 'L_M8',
                              'A_M2', 'A_M3', 'A_M7', 'W_M7', 'L_S1',
                              'L_S3_1', 'L_S3_2', 'L_S6_1', 'L_S6_2', 'L_S7_1',
                              'A_S1', 'A_S4', 'A_S5', 'A_S6', 'W_S1_1',
                              'W_S1_2', 'W_S2', 'W_S3', 'W_S4', 'W_S6'))

#DefaultAssay(T.integrated) = 'SCT'


###############################################################
## Subset cells to cluster 9
################################################################
cluster_id = 8
id = T.integrated$seurat_clusters %in% cluster_id
T_sub = subset(T.integrated, cells = colnames(T.integrated)[id])
Dis =  data.frame(Covid = T_sub$Disease, row.names = colnames(T_sub))
Dis[Dis != 'Healthy'] = 'COVID'
T_sub <- AddMetaData(object = T_sub, metadata = Dis)

VlnPlot(T_sub, features = inflam_g[i], assay = 'RNA', group.by = 'Disease',  pt.size = 0.02)

B <- as.SingleCellExperiment(T_sub, assay = 'RNA')
B$cluster = 'PBMC'

qc <- perCellQCMetrics(B)
# remove cells with few or many detected genes
ol <- isOutlier(metric = qc$detected, nmads = 2, log = TRUE)
#B <- B[, !ol]
dim(B)

n = dim(B)[2]
gene_id =  match(inflam_g, rownames(B))
tmp_id = rowSums(counts(B) > .1) > (n * .1)  # filter genes with high expression
tmp_id[gene_id] = TRUE
B <- B[tmp_id, ]
dim(B)

id = B$Disease == 'Healthy'
B = B[, !id]


(B <- prepSCE(B, 
    kid = "cluster", # subpopulation assignments
    gid = "Disease",  # group IDs (ctrl/stim)
    sid = "Sample",   # sample IDs (ctrl/stim.1234)
    drop = TRUE))  # drop all other colData columns



clusters <- rep(1, dim(B)[2])
#B <- computeSumFactors(B, clusters=clusters)
B <- computeLibraryFactors(B)
#B <- logNormCounts(B)
summary(sizeFactors(B))
#assays(B)$vstresiduals <- vst(counts(B), show_progress = FALSE)$y

gene_id =  match(inflam_g, rownames(B))
B1 = B[gene_id, ]

mm1 <- mmDS(B1, method = "dream2",  n_cells = 10, n_samples = 2,  min_count = 1, min_cells = 20)
mm <- mmDS(B1, method = "nbinom")
#mm <- mmDS(B1, method = "nbinom", n_cells = 10, n_samples = 2,  min_count = 1, min_cells = 20)

# Healthy vs COVID
Dis =  B1$group_id

##################################################
##  CLuster 11
##############################################
cluster_id = 10
id = T.integrated$seurat_clusters %in% cluster_id
T_sub = subset(T.integrated, cells = colnames(T.integrated)[id])
Dis =  data.frame(Covid = T_sub$Disease, row.names = colnames(T_sub))
Dis[Dis != 'Healthy'] = 'COVID'
T_sub <- AddMetaData(object = T_sub, metadata = Dis)

VlnPlot(T_sub, features = inflam_g[i], assay = 'RNA', group.by = 'Disease',  pt.size = 0.02)

B <- as.SingleCellExperiment(T_sub, assay = 'RNA')
B$cluster = 'PBMC'

qc <- perCellQCMetrics(B)
# remove cells with few or many detected genes
ol <- isOutlier(metric = qc$detected, nmads = 2, log = TRUE)
#B <- B[, !ol]
dim(B)

n = dim(B)[2]
gene_id =  match(inflam_g, rownames(B))
tmp_id = rowSums(counts(B) > .1) > (n * .1)  # filter genes with high expression
tmp_id[gene_id] = TRUE
B <- B[tmp_id, ]
dim(B)

id = B$Disease == 'Healthy'
B = B[, !id]


(B <- prepSCE(B, 
    kid = "cluster", # subpopulation assignments
    gid = "Disease",  # group IDs (ctrl/stim)
    sid = "Sample",   # sample IDs (ctrl/stim.1234)
    drop = TRUE))  # drop all other colData columns

(B2 <- prepSCE(B2, 
    kid = "cluster", # subpopulation assignments
    gid = "Disease",  # group IDs (ctrl/stim)
    sid = "Sample",   # sample IDs (ctrl/stim.1234)
    drop = TRUE))  # drop all other colData columns


(B <- prepSCE(B, 
    kid = "cluster", # subpopulation assignments
    gid = "Covid",  # group IDs (ctrl/stim)
    sid = "Sample",   # sample IDs (ctrl/stim.1234)
    drop = TRUE))  # drop all other colData columns


clusters <- rep(1, dim(B)[2])
#B <- computeSumFactors(B, clusters=clusters)
B <- computeLibraryFactors(B)
#B <- logNormCounts(B)
summary(sizeFactors(B))
#assays(B)$vstresiduals <- vst(counts(B), show_progress = FALSE)$y

gene_id =  match(inflam_g, rownames(B))
B1 = B[gene_id, ]

mm1 <- mmDS(B1, method = "dream2",  n_cells = 10, n_samples = 2,  min_count = 1, min_cells = 20)
mm_all <- mmDS(B1, method = "nbinom")
mm_HC <- mmDS(B1, method = "nbinom")
mm_MS <- mmDS(B1, method = "nbinom")

#mm <- mmDS(B1, method = "nbinom", n_cells = 10, n_samples = 2,  min_count = 1, min_cells = 20)

A = cbind(mm_all$PBMC[, c(1, 6)], mm_HC$PBMC[, 6], mm_MS$PBMC[, 6])
colnames(A) = c('gene', 'All', 'Health_vs_COVID', 'Moderate_vs_Severe')



##################################################
##  CLuster 9-11-15
##############################################
cluster_id = c(8, 10, 14)
id = T.integrated$seurat_clusters %in% cluster_id
T_sub = subset(T.integrated, cells = colnames(T.integrated)[id])
Dis =  data.frame(Covid = T_sub$Disease, row.names = colnames(T_sub))
Dis[Dis != 'Healthy'] = 'COVID'
T_sub <- AddMetaData(object = T_sub, metadata = Dis)


pdf('T_infl_cluster_9_11_15_Samples_Vlnplot_01_04_2021.pdf',  width = 60, height = 20)
  print(VlnPlot(T_sub, features = inflam_g, assay = 'RNA', group.by = 'Sample_f',  pt.size = 0.02))
dev.off()

pdf('T_infl_cluster_9_11_15_Disease_Vlnplot_01_04_2021.pdf', width = 40, height = 20)
  VlnPlot(T_sub, features = inflam_g, assay = 'RNA', group.by = 'Disease',  pt.size = 0.02)
dev.off()

B <- as.SingleCellExperiment(T_sub, assay = 'RNA')
B$cluster = 'PBMC'

qc <- perCellQCMetrics(B)
# remove cells with few or many detected genes
ol <- isOutlier(metric = qc$detected, nmads = 2, log = TRUE)
#B <- B[, !ol]
dim(B)

n = dim(B)[2]
gene_id =  match(inflam_g, rownames(B))
tmp_id = rowSums(counts(B) > .1) > (n * .1)  # filter genes with high expression
tmp_id[gene_id] = TRUE
B <- B[tmp_id, ]
dim(B)

id = B$Disease == 'Healthy'
B = B[, !id]


(BB <- prepSCE(B, 
    kid = "cluster", # subpopulation assignments
    gid = "Disease",  # group IDs (ctrl/stim)
    sid = "Sample",   # sample IDs (ctrl/stim.1234)
    drop = TRUE))  # drop all other colData columns

(BB1 <- prepSCE(B, 
    kid = "cluster", # subpopulation assignments
    gid = "Disease",  # group IDs (ctrl/stim)
    sid = "Sample",   # sample IDs (ctrl/stim.1234)
    drop = TRUE))  # drop all other colData columns


(BB2 <- prepSCE(B, 
    kid = "cluster", # subpopulation assignments
    gid = "Covid",  # group IDs (ctrl/stim)
    sid = "Sample",   # sample IDs (ctrl/stim.1234)
    drop = TRUE))  # drop all other colData columns

clusters <- rep(1, dim(B)[2])
#B <- computeSumFactors(B, clusters=clusters)--
BB <- computeLibraryFactors(BB)
BB1 <- computeLibraryFactors(BB1)
BB2 <- computeLibraryFactors(BB2)
#B <- logNormCounts(B)
summary(sizeFactors(B))
#assays(B)$vstresiduals <- vst(counts(B), show_progress = FALSE)$y

gene_id =  match(inflam_g, rownames(BB))
B1 = B[gene_id, ]

mm1 <- mmDS(B1, method = "dream2",  n_cells = 10, n_samples = 2,  min_count = 1, min_cells = 20)
mm_all <- mmDS(BB[gene_id, ], method = "nbinom")
mm_HC <- mmDS(BB2[gene_id, ], method = "nbinom")
mm_MS <- mmDS(BB1[gene_id, ], method = "nbinom")
A = cbind(mm_all$PBMC[, c(1, 6)], mm_HC$PBMC[, 6], mm_MS$PBMC[, 6])
colnames(A) = c('gene', 'All', 'Health_vs_COVID', 'Moderate_vs_Severe')


##################################################
##  CLuster 30
##############################################
cluster_id = 29
id = T.integrated$seurat_clusters %in% cluster_id
T_sub = subset(T.integrated, cells = colnames(T.integrated)[id])
Dis =  data.frame(Covid = T_sub$Disease, row.names = colnames(T_sub))
Dis[Dis != 'Healthy'] = 'COVID'
T_sub <- AddMetaData(object = T_sub, metadata = Dis)

VlnPlot(T_sub, features = inflam_g[i], assay = 'RNA', group.by = 'Disease',  pt.size = 0.02)

pdf('T_infl_cluster_9_11_15_Samples_Vlnplot_01_04_2021.pdf',  width = 60, height = 20)
  print(VlnPlot(T_sub, features = inflam_g, assay = 'RNA', group.by = 'Sample_f',  pt.size = 0.02))
dev.off()

pdf('T_infl_cluster_9_11_15_Disease_Vlnplot_01_04_2021.pdf', width = 40, height = 20)
  VlnPlot(T_sub, features = inflam_g, assay = 'RNA', group.by = 'Disease',  pt.size = 0.02)
dev.off()

B <- as.SingleCellExperiment(T_sub, assay = 'RNA')
B$cluster = 'PBMC'

qc <- perCellQCMetrics(B)
# remove cells with few or many detected genes
ol <- isOutlier(metric = qc$detected, nmads = 2, log = TRUE)
#B <- B[, !ol]
dim(B)

n = dim(B)[2]
gene_id =  match(inflam_g, rownames(B))
tmp_id = rowSums(counts(B) > .1) > (n * .1)  # filter genes with high expression
tmp_id[gene_id] = TRUE
B <- B[tmp_id, ]
dim(B)

id = B$Disease == 'Healthy'
B = B[, !id]


(BB <- prepSCE(B, 
    kid = "cluster", # subpopulation assignments
    gid = "Disease",  # group IDs (ctrl/stim)
    sid = "Sample",   # sample IDs (ctrl/stim.1234)
    drop = TRUE))  # drop all other colData columns

(BB1 <- prepSCE(B, 
    kid = "cluster", # subpopulation assignments
    gid = "Disease",  # group IDs (ctrl/stim)
    sid = "Sample",   # sample IDs (ctrl/stim.1234)
    drop = TRUE))  # drop all other colData columns


(BB2 <- prepSCE(B, 
    kid = "cluster", # subpopulation assignments
    gid = "Covid",  # group IDs (ctrl/stim)
    sid = "Sample",   # sample IDs (ctrl/stim.1234)
    drop = TRUE))  # drop all other colData columns

clusters <- rep(1, dim(B)[2])
#B <- computeSumFactors(B, clusters=clusters)--
BB <- computeLibraryFactors(BB)
BB1 <- computeLibraryFactors(BB1)
BB2 <- computeLibraryFactors(BB2)
#B <- logNormCounts(B)
summary(sizeFactors(B))
#assays(B)$vstresiduals <- vst(counts(B), show_progress = FALSE)$y

gene_id =  match(inflam_g, rownames(BB))
B1 = B[gene_id, ]

mm1 <- mmDS(B1, method = "dream2",  n_cells = 10, n_samples = 2,  min_count = 1, min_cells = 20)
mm_all <- mmDS(BB[gene_id, ], method = "nbinom")
mm_HC <- mmDS(BB2[gene_id, ], method = "nbinom")
mm_MS <- mmDS(BB1[gene_id, ], method = "nbinom")
A = cbind(mm_all$PBMC[, c(1, 6)], mm_HC$PBMC[, 6], mm_MS$PBMC[, 6])
colnames(A) = c('gene', 'All', 'Health_vs_COVID', 'Moderate_vs_Severe')


##############################################################
## Subset cells to clusters (0, 8, 10, 14, 15, 19, 27, 29)
##############################################################
cluster_id = c(0, 8, 10, 14, 15, 19, 27, 29)
id = T.integrated$seurat_clusters %in% cluster_id

#saveRDS(T_sub, file ='T_SCT_int_sub_11_28_2020.rds')
T_sub = readRDS('T_SCT_int_sub_11_28_2020.rds')

pdf(paste('T_genes_SCT_clusters_samples_03_07_2021.pdf'),  width = 60, height = 20)
    print(VlnPlot(T_sub, features = inflam_g, assay = 'SCT', group.by = 'Sample_f',  pt.size = 0.02))
dev.off()

pdf(paste('T_genes_8clusters_Vlnplot_11_29_2020.pdf'),  width = 40, height = 20)
    print(VlnPlot(T_sub, features = inflam_g, assay = 'SCT', group.by = 'Disease',  pt.size = 0.02))
dev.off()

B <- as.SingleCellExperiment(T_sub, assay = 'RNA')
#A = GetAssayData(T_sub, assay = 'SCT', slot = 'scale.data')
#id = match(rownames(A), rownames(B))
#B = B[id, ]
B$cluster = 'PBMC'

qc <- perCellQCMetrics(B)
# remove cells with few or many detected genes
ol <- isOutlier(metric = qc$detected, nmads = 2, log = TRUE)
B <- B[, !ol]
dim(B)

n = dim(B)[2]
gene_id =  match(inflam_g, rownames(B))
tmp_id = rowSums(counts(B) > .1) > (n * .1)  # filter genes with high expression
tmp_id[gene_id] = TRUE
B <- B[tmp_id, ]
dim(B)

(B <- prepSCE(B, 
    kid = "cluster", # subpopulation assignments
    gid = "Disease",  # group IDs (ctrl/stim)
    sid = "Sample",   # sample IDs (ctrl/stim.1234)
    drop = TRUE))  # drop all other colData columns

clusters <- rep(1, dim(B)[2])
B <- computeSumFactors(B, clusters=clusters)
B <- logNormCounts(B)
summary(sizeFactors(B))

gene_id =  match(inflam_g, rownames(B))
B1 = B[gene_id, ]

mm <- mmDS(B1, method = "dream",  n_cells = 10, n_samples = 2,  min_count = 1, min_cells = 20)
#mm <- mmDS(B, method = "nbinom")
#mm <- mmDS(B1, method = "nbinom")
write.csv(mm, file = '3clusters_mm.csv')





###############################################################
## Subset cells to cluster  1, 9, 11, 15, 16, 20, 28, and 30
################################################################
cluster_id = c( 0, 8, 10, 14, 15, 19, 27, 29)
id = T.integrated$seurat_clusters %in% cluster_id
T_sub = subset(T.integrated, cells = colnames(T.integrated)[id])
Dis =  data.frame(Covid = T_sub$Disease, row.names = colnames(T_sub))
Dis[Dis != 'Healthy'] = 'COVID'
T_sub <- AddMetaData(object = T_sub, metadata = Dis)

VlnPlot(T_sub, features = inflam_g[i], assay = 'RNA', group.by = 'Disease',  pt.size = 0.02)

pdf('T_infl_cluster_1_9_11_15_16_20_28_30_Samples_Vlnplot_03_03_2021.pdf',  width = 60, height = 20)
  print(VlnPlot(T_sub, features = inflam_g, assay = 'RNA', group.by = 'Sample_f',  pt.size = 0.02))
dev.off()

pdf('T_infl_cluster_1_9_11_15_16_20_28_30_Disease_Vlnplot_03_03_2021.pdf', width = 40, height = 20)
  VlnPlot(T_sub, features = inflam_g, assay = 'RNA', group.by = 'Disease',  pt.size = 0.02)
dev.off()

B <- as.SingleCellExperiment(T_sub, assay = 'RNA')
B$cluster = 'PBMC'

qc <- perCellQCMetrics(B)
# remove cells with few or many detected genes
ol <- isOutlier(metric = qc$detected, nmads = 2, log = TRUE)
#B <- B[, !ol]
dim(B)

n = dim(B)[2]
gene_id =  match(inflam_g, rownames(B))
tmp_id = rowSums(counts(B) > .1) > (n * .1)  # filter genes with high expression
tmp_id[gene_id] = TRUE
B <- B[tmp_id, ]
dim(B)

id = B$Disease == 'Healthy'
B1 = B[, !id]   # only Moderate and Severe

# BB has three groups: H, M, and S
(BB <- prepSCE(B, 
    kid = "cluster", # subpopulation assignments
    gid = "Disease",  # group IDs (ctrl/stim)
    sid = "Sample",   # sample IDs (ctrl/stim.1234)
    drop = TRUE))  # drop all other colData columns

# BB1 has two groups:  M, and S
(BB1 <- prepSCE(B1, 
    kid = "cluster", # subpopulation assignments
    gid = "Disease",  # group IDs (ctrl/stim)
    sid = "Sample",   # sample IDs (ctrl/stim.1234)
    drop = TRUE))  # drop all other colData columns

# BB2 has two groups:  Healthy vs COVID
(BB2 <- prepSCE(B, 
    kid = "cluster", # subpopulation assignments
    gid = "Covid",  # group IDs (ctrl/stim)
    sid = "Sample",   # sample IDs (ctrl/stim.1234)
    drop = TRUE))  # drop all other colData columns

#clusters <- rep(1, dim(B)[2])
#B <- computeSumFactors(B, clusters=clusters)
BB <- computeLibraryFactors(BB)
BB1 <- computeLibraryFactors(BB1)
BB2 <- computeLibraryFactors(BB2)
#B <- logNormCounts(B)
#summary(sizeFactors(B))
#assays(B)$vstresiduals <- vst(counts(B), show_progress = FALSE)$y

gene_id =  match(inflam_g, rownames(BB))
B1 = B[gene_id, ]

# Can use this to check level colnames(model.matrix(~group_id, data=as.data.frame(colData(BB))))
#mm1 <- mmDS(B1, method = "dream2",  n_cells = 10, n_samples = 2,  min_count = 1, min_cells = 20)
mm_MvsH <- mmDS(BB[gene_id, ], method = "nbinom", coef = 'group_idModerate')
mm_SvsH <- mmDS(BB[gene_id, ], method = "nbinom", coef = 'group_idSevere')
mm_HC <- mmDS(BB2[gene_id, ], method = "nbinom")
mm_MS <- mmDS(BB1[gene_id, ], method = "nbinom")
A = cbind(mm_MvsH$PBMC[, c(1, 6)], mm_SvsH$PBMC[,  6], mm_HC$PBMC[, 6], mm_MS$PBMC[, 6])
colnames(A) = c('gene', 'Moderate_vs_Healthey', 'Severe_vs_Healthy','Health_vs_COVID', 'Moderate_vs_Severe')

# want individal C1, C15, C16, C20, C28,
###############################################################
## Subset cells to cluster  1
################################################################
cluster_id = 0
id = T.integrated$seurat_clusters %in% cluster_id
T_sub = subset(T.integrated, cells = colnames(T.integrated)[id])
Dis =  data.frame(Covid = T_sub$Disease, row.names = colnames(T_sub))
Dis[Dis != 'Healthy'] = 'COVID'
T_sub <- AddMetaData(object = T_sub, metadata = Dis)

VlnPlot(T_sub, features = inflam_g[i], assay = 'RNA', group.by = 'Disease',  pt.size = 0.02)

pdf('T_infl_cluster_1_Samples_Vlnplot_03_04_2021.pdf',  width = 60, height = 20)
  print(VlnPlot(T_sub, features = inflam_g, assay = 'RNA', group.by = 'Sample_f',  pt.size = 0.02))
dev.off()

pdf('T_infl_cluster_1_Disease_Vlnplot_03_04_2021.pdf', width = 40, height = 20)
  VlnPlot(T_sub, features = inflam_g, assay = 'RNA', group.by = 'Disease',  pt.size = 0.02)
dev.off()

B <- as.SingleCellExperiment(T_sub, assay = 'RNA')
B$cluster = 'PBMC'

qc <- perCellQCMetrics(B)
# remove cells with few or many detected genes
ol <- isOutlier(metric = qc$detected, nmads = 2, log = TRUE)
#B <- B[, !ol]
dim(B)

n = dim(B)[2]
gene_id =  match(inflam_g, rownames(B))
tmp_id = rowSums(counts(B) > .1) > (n * .1)  # filter genes with high expression
tmp_id[gene_id] = TRUE
B <- B[tmp_id, ]
dim(B)

id = B$Disease == 'Healthy'
B1 = B[, !id]


(BB <- prepSCE(B, 
    kid = "cluster", # subpopulation assignments
    gid = "Disease",  # group IDs (ctrl/stim)
    sid = "Sample",   # sample IDs (ctrl/stim.1234)
    drop = TRUE))  # drop all other colData columns

(BB1 <- prepSCE(B1, 
    kid = "cluster", # subpopulation assignments
    gid = "Disease",  # group IDs (ctrl/stim)
    sid = "Sample",   # sample IDs (ctrl/stim.1234)
    drop = TRUE))  # drop all other colData columns


(BB2 <- prepSCE(B, 
    kid = "cluster", # subpopulation assignments
    gid = "Covid",  # group IDs (ctrl/stim)
    sid = "Sample",   # sample IDs (ctrl/stim.1234)
    drop = TRUE))  # drop all other colData columns

clusters <- rep(1, dim(B)[2])
#B <- computeSumFactors(B, clusters=clusters)--
BB <- computeLibraryFactors(BB)
BB1 <- computeLibraryFactors(BB1)
BB2 <- computeLibraryFactors(BB2)
#B <- logNormCounts(B)
summary(sizeFactors(B))
#assays(B)$vstresiduals <- vst(counts(B), show_progress = FALSE)$y

gene_id =  match(inflam_g, rownames(BB))
B1 = B[gene_id, ]

# Can use this to check level colnames(model.matrix(~group_id, data=as.data.frame(colData(BB))))
mm1 <- mmDS(B1, method = "dream2",  n_cells = 10, n_samples = 2,  min_count = 1, min_cells = 20)
mm_MvsH <- mmDS(BB[gene_id, ], method = "nbinom", coef = 'group_idModerate')
mm_SvsH <- mmDS(BB[gene_id, ], method = "nbinom", coef = 'group_idSevere')
mm_HC <- mmDS(BB2[gene_id, ], method = "nbinom")
mm_MS <- mmDS(BB1[gene_id, ], method = "nbinom")
A = cbind(mm_MvsH$PBMC[, c(1, 6)], mm_SvsH$PBMC[,  6], mm_HC$PBMC[, 6], mm_MS$PBMC[, 6])
colnames(A) = c('gene', 'Moderate_vs_Healthey', 'Severe_vs_Healthy','Health_vs_COVID', 'Moderate_vs_Severe')


# want individal C1, C15, C16, C20, C28,
###############################################################
## Subset cells to cluster  9
################################################################
cluster_id = 8
id = T.integrated$seurat_clusters %in% cluster_id
T_sub = subset(T.integrated, cells = colnames(T.integrated)[id])
Dis =  data.frame(Covid = T_sub$Disease, row.names = colnames(T_sub))
Dis[Dis != 'Healthy'] = 'COVID'
T_sub <- AddMetaData(object = T_sub, metadata = Dis)

VlnPlot(T_sub, features = inflam_g[i], assay = 'RNA', group.by = 'Disease',  pt.size = 0.02)

pdf('T_infl_cluster_9_Samples_Vlnplot_03_04_2021.pdf',  width = 60, height = 20)
  print(VlnPlot(T_sub, features = inflam_g, assay = 'RNA', group.by = 'Sample_f',  pt.size = 0.02))
dev.off()

pdf('T_infl_cluster_9_Disease_Vlnplot_03_04_2021.pdf', width = 40, height = 20)
  VlnPlot(T_sub, features = inflam_g, assay = 'RNA', group.by = 'Disease',  pt.size = 0.02)
dev.off()

B <- as.SingleCellExperiment(T_sub, assay = 'RNA')
B$cluster = 'PBMC'

qc <- perCellQCMetrics(B)
# remove cells with few or many detected genes
ol <- isOutlier(metric = qc$detected, nmads = 2, log = TRUE)
#B <- B[, !ol]
dim(B)

n = dim(B)[2]
gene_id =  match(inflam_g, rownames(B))
tmp_id = rowSums(counts(B) > .1) > (n * .1)  # filter genes with high expression
tmp_id[gene_id] = TRUE
B <- B[tmp_id, ]
dim(B)

id = B$Disease == 'Healthy'
B1 = B[, !id]


(BB <- prepSCE(B, 
    kid = "cluster", # subpopulation assignments
    gid = "Disease",  # group IDs (ctrl/stim)
    sid = "Sample",   # sample IDs (ctrl/stim.1234)
    drop = TRUE))  # drop all other colData columns

(BB1 <- prepSCE(B1, 
    kid = "cluster", # subpopulation assignments
    gid = "Disease",  # group IDs (ctrl/stim)
    sid = "Sample",   # sample IDs (ctrl/stim.1234)
    drop = TRUE))  # drop all other colData columns


(BB2 <- prepSCE(B, 
    kid = "cluster", # subpopulation assignments
    gid = "Covid",  # group IDs (ctrl/stim)
    sid = "Sample",   # sample IDs (ctrl/stim.1234)
    drop = TRUE))  # drop all other colData columns

clusters <- rep(1, dim(B)[2])
#B <- computeSumFactors(B, clusters=clusters)--
BB <- computeLibraryFactors(BB)
BB1 <- computeLibraryFactors(BB1)
BB2 <- computeLibraryFactors(BB2)
#B <- logNormCounts(B)
summary(sizeFactors(B))
#assays(B)$vstresiduals <- vst(counts(B), show_progress = FALSE)$y

gene_id =  match(inflam_g, rownames(BB))
B1 = B[gene_id, ]

# Can use this to check level colnames(model.matrix(~group_id, data=as.data.frame(colData(BB))))
mm1 <- mmDS(B1, method = "dream2",  n_cells = 10, n_samples = 2,  min_count = 1, min_cells = 20)
mm_MvsH <- mmDS(BB[gene_id, ], method = "nbinom", coef = 'group_idModerate')
mm_SvsH <- mmDS(BB[gene_id, ], method = "nbinom", coef = 'group_idSevere')
mm_HC <- mmDS(BB2[gene_id, ], method = "nbinom")
mm_MS <- mmDS(BB1[gene_id, ], method = "nbinom")
A = cbind(mm_MvsH$PBMC[, c(1, 6)], mm_SvsH$PBMC[,  6], mm_HC$PBMC[, 6], mm_MS$PBMC[, 6])
colnames(A) = c('gene', 'Moderate_vs_Healthey', 'Severe_vs_Healthy','Health_vs_COVID', 'Moderate_vs_Severe')



# want individal C1, C15, C16, C20, C28,
###############################################################
## Subset cells to cluster  9, 11, 15
################################################################
cluster_id = c(8, 10, 14)
id = T.integrated$seurat_clusters %in% cluster_id
T_sub = subset(T.integrated, cells = colnames(T.integrated)[id])
Dis =  data.frame(Covid = T_sub$Disease, row.names = colnames(T_sub))
Dis[Dis != 'Healthy'] = 'COVID'
T_sub <- AddMetaData(object = T_sub, metadata = Dis)

VlnPlot(T_sub, features = inflam_g[i], assay = 'RNA', group.by = 'Disease',  pt.size = 0.02)

pdf('T_infl_clusters_9_11_15_Samples_Vlnplot_03_04_2021.pdf',  width = 60, height = 20)
  print(VlnPlot(T_sub, features = inflam_g, assay = 'RNA', group.by = 'Sample_f',  pt.size = 0.02))
dev.off()

pdf('T_infl_cluster_9_11_15_Disease_Vlnplot_03_04_2021.pdf', width = 40, height = 20)
  VlnPlot(T_sub, features = inflam_g, assay = 'RNA', group.by = 'Disease',  pt.size = 0.02)
dev.off()

id = T_sub$Disease == 'Healthy'
dat = GetAssayData(T_sub['NAIP', id], assay = 'RNA', slot ='data')
mean(dat)
id = T_sub$Disease == 'Moderate'
dat = GetAssayData(T_sub['NAIP', id], assay = 'RNA', slot ='data')
mean(dat)


id = T_sub$Disease == 'Severe'
dat = GetAssayData(T_sub['NAIP', id], assay = 'RNA', slot ='data')
mean(dat)



B <- as.SingleCellExperiment(T_sub, assay = 'RNA')
B$cluster = 'PBMC'

qc <- perCellQCMetrics(B)
# remove cells with few or many detected genes
ol <- isOutlier(metric = qc$detected, nmads = 2, log = TRUE)
#B <- B[, !ol]
dim(B)

n = dim(B)[2]
gene_id =  match(inflam_g, rownames(B))
tmp_id = rowSums(counts(B) > .1) > (n * .1)  # filter genes with high expression
tmp_id[gene_id] = TRUE
B <- B[tmp_id, ]
dim(B)

id = B$Disease == 'Healthy'
B1 = B[, !id]


(BB <- prepSCE(B, 
    kid = "cluster", # subpopulation assignments
    gid = "Disease",  # group IDs (ctrl/stim)
    sid = "Sample",   # sample IDs (ctrl/stim.1234)
    drop = TRUE))  # drop all other colData columns

(BB1 <- prepSCE(B1, 
    kid = "cluster", # subpopulation assignments
    gid = "Disease",  # group IDs (ctrl/stim)
    sid = "Sample",   # sample IDs (ctrl/stim.1234)
    drop = TRUE))  # drop all other colData columns


(BB2 <- prepSCE(B, 
    kid = "cluster", # subpopulation assignments
    gid = "Covid",  # group IDs (ctrl/stim)
    sid = "Sample",   # sample IDs (ctrl/stim.1234)
    drop = TRUE))  # drop all other colData columns

clusters <- rep(1, dim(B)[2])
#B <- computeSumFactors(B, clusters=clusters)--
BB <- computeLibraryFactors(BB)
BB1 <- computeLibraryFactors(BB1)
BB2 <- computeLibraryFactors(BB2)
#B <- logNormCounts(B)
summary(sizeFactors(B))
#assays(B)$vstresiduals <- vst(counts(B), show_progress = FALSE)$y

gene_id =  match(inflam_g, rownames(BB))
B1 = B[gene_id, ]

# Can use this to check level colnames(model.matrix(~group_id, data=as.data.frame(colData(BB))))
mm1 <- mmDS(B1, method = "dream2",  n_cells = 10, n_samples = 2,  min_count = 1, min_cells = 20)
mm_MvsH <- mmDS(BB[gene_id, ], method = "nbinom", coef = 'group_idModerate')
mm_SvsH <- mmDS(BB[gene_id, ], method = "nbinom", coef = 'group_idSevere')
mm_HC <- mmDS(BB2[gene_id, ], method = "nbinom")
mm_MS <- mmDS(BB1[gene_id, ], method = "nbinom")
A = cbind(mm_MvsH$PBMC[, c(1, 6)], mm_SvsH$PBMC[,  6], mm_HC$PBMC[, 6], mm_MS$PBMC[, 6])
colnames(A) = c('gene', 'Moderate_vs_Healthey', 'Severe_vs_Healthy','Health_vs_COVID', 'Moderate_vs_Severe')



# want individal C1, C15, C16, C20, C28,
###############################################################
## Subset cells to cluster  9, 15
################################################################
cluster_id = c(8,  14)
id = T.integrated$seurat_clusters %in% cluster_id
T_sub = subset(T.integrated, cells = colnames(T.integrated)[id])
Dis =  data.frame(Covid = T_sub$Disease, row.names = colnames(T_sub))
Dis[Dis != 'Healthy'] = 'COVID'
T_sub <- AddMetaData(object = T_sub, metadata = Dis)

VlnPlot(T_sub, features = inflam_g[i], assay = 'RNA', group.by = 'Disease',  pt.size = 0.02)

pdf('T_infl_clusters_9_15_Samples_Vlnplot_03_07_2021.pdf',  width = 60, height = 20)
  print(VlnPlot(T_sub, features = inflam_g, assay = 'RNA', group.by = 'Sample_f',  pt.size = 0.02))
dev.off()

pdf('T_infl_cluster_9_15_Disease_Vlnplot_03_07_2021.pdf', width = 40, height = 20)
  VlnPlot(T_sub, features = inflam_g, assay = 'RNA', group.by = 'Disease',  pt.size = 0.02)
dev.off()

id = T_sub$Disease == 'Healthy'
dat = GetAssayData(T_sub['NAIP', id], assay = 'RNA', slot ='data')
mean(dat)
id = T_sub$Disease == 'Moderate'
dat = GetAssayData(T_sub['NAIP', id], assay = 'RNA', slot ='data')
mean(dat)


id = T_sub$Disease == 'Severe'
dat = GetAssayData(T_sub['NAIP', id], assay = 'RNA', slot ='data')
mean(dat)



B <- as.SingleCellExperiment(T_sub, assay = 'RNA')
B$cluster = 'PBMC'

qc <- perCellQCMetrics(B)
# remove cells with few or many detected genes
ol <- isOutlier(metric = qc$detected, nmads = 2, log = TRUE)
#B <- B[, !ol]
dim(B)

n = dim(B)[2]
gene_id =  match(inflam_g, rownames(B))
tmp_id = rowSums(counts(B) > .1) > (n * .1)  # filter genes with high expression
tmp_id[gene_id] = TRUE
B <- B[tmp_id, ]
dim(B)

id = B$Disease == 'Healthy'
B1 = B[, !id]


(BB <- prepSCE(B, 
    kid = "cluster", # subpopulation assignments
    gid = "Disease",  # group IDs (ctrl/stim)
    sid = "Sample",   # sample IDs (ctrl/stim.1234)
    drop = TRUE))  # drop all other colData columns

(BB1 <- prepSCE(B1, 
    kid = "cluster", # subpopulation assignments
    gid = "Disease",  # group IDs (ctrl/stim)
    sid = "Sample",   # sample IDs (ctrl/stim.1234)
    drop = TRUE))  # drop all other colData columns


(BB2 <- prepSCE(B, 
    kid = "cluster", # subpopulation assignments
    gid = "Covid",  # group IDs (ctrl/stim)
    sid = "Sample",   # sample IDs (ctrl/stim.1234)
    drop = TRUE))  # drop all other colData columns

clusters <- rep(1, dim(B)[2])
#B <- computeSumFactors(B, clusters=clusters)--
BB <- computeLibraryFactors(BB)
BB1 <- computeLibraryFactors(BB1)
BB2 <- computeLibraryFactors(BB2)
#B <- logNormCounts(B)
summary(sizeFactors(B))
#assays(B)$vstresiduals <- vst(counts(B), show_progress = FALSE)$y

gene_id =  match(inflam_g, rownames(BB))
B1 = B[gene_id, ]

# Can use this to check level colnames(model.matrix(~group_id, data=as.data.frame(colData(BB))))
mm1 <- mmDS(B1, method = "dream2",  n_cells = 10, n_samples = 2,  min_count = 1, min_cells = 20)
mm_MvsH <- mmDS(BB[gene_id, ], method = "nbinom", coef = 'group_idModerate')
mm_SvsH <- mmDS(BB[gene_id, ], method = "nbinom", coef = 'group_idSevere')
mm_HC <- mmDS(BB2[gene_id, ], method = "nbinom")
mm_MS <- mmDS(BB1[gene_id, ], method = "nbinom")
A = cbind(mm_MvsH$PBMC[, c(1, 6)], mm_SvsH$PBMC[,  6], mm_HC$PBMC[, 6], mm_MS$PBMC[, 6])
colnames(A) = c('gene', 'Moderate_vs_Healthey', 'Severe_vs_Healthy','Health_vs_COVID', 'Moderate_vs_Severe')



# want individal C1, C15, C16, C20, C28,
###############################################################
## Subset cells to cluster  11
################################################################
cluster_id = 10
id = T.integrated$seurat_clusters %in% cluster_id
T_sub = subset(T.integrated, cells = colnames(T.integrated)[id])
Dis =  data.frame(Covid = T_sub$Disease, row.names = colnames(T_sub))
Dis[Dis != 'Healthy'] = 'COVID'
T_sub <- AddMetaData(object = T_sub, metadata = Dis)

VlnPlot(T_sub, features = inflam_g[i], assay = 'RNA', group.by = 'Disease',  pt.size = 0.02)

pdf('T_infl_cluster_11_Samples_Vlnplot_03_04_2021.pdf',  width = 60, height = 20)
  print(VlnPlot(T_sub, features = inflam_g, assay = 'RNA', group.by = 'Sample_f',  pt.size = 0.02))
dev.off()

pdf('T_infl_cluster_11_Disease_Vlnplot_03_04_2021.pdf', width = 40, height = 20)
  VlnPlot(T_sub, features = inflam_g, assay = 'RNA', group.by = 'Disease',  pt.size = 0.02)
dev.off()

B <- as.SingleCellExperiment(T_sub, assay = 'RNA')
B$cluster = 'PBMC'

qc <- perCellQCMetrics(B)
# remove cells with few or many detected genes
ol <- isOutlier(metric = qc$detected, nmads = 2, log = TRUE)
#B <- B[, !ol]
dim(B)

n = dim(B)[2]
gene_id =  match(inflam_g, rownames(B))
tmp_id = rowSums(counts(B) > .1) > (n * .1)  # filter genes with high expression
tmp_id[gene_id] = TRUE
B <- B[tmp_id, ]
dim(B)

id = B$Disease == 'Healthy'
B1 = B[, !id]


(BB <- prepSCE(B, 
    kid = "cluster", # subpopulation assignments
    gid = "Disease",  # group IDs (ctrl/stim)
    sid = "Sample",   # sample IDs (ctrl/stim.1234)
    drop = TRUE))  # drop all other colData columns

(BB1 <- prepSCE(B1, 
    kid = "cluster", # subpopulation assignments
    gid = "Disease",  # group IDs (ctrl/stim)
    sid = "Sample",   # sample IDs (ctrl/stim.1234)
    drop = TRUE))  # drop all other colData columns


(BB2 <- prepSCE(B, 
    kid = "cluster", # subpopulation assignments
    gid = "Covid",  # group IDs (ctrl/stim)
    sid = "Sample",   # sample IDs (ctrl/stim.1234)
    drop = TRUE))  # drop all other colData columns

clusters <- rep(1, dim(B)[2])
#B <- computeSumFactors(B, clusters=clusters)--
BB <- computeLibraryFactors(BB)
BB1 <- computeLibraryFactors(BB1)
BB2 <- computeLibraryFactors(BB2)
#B <- logNormCounts(B)
summary(sizeFactors(B))
#assays(B)$vstresiduals <- vst(counts(B), show_progress = FALSE)$y

gene_id =  match(inflam_g, rownames(BB))
B1 = B[gene_id, ]

# Can use this to check level colnames(model.matrix(~group_id, data=as.data.frame(colData(BB))))
mm1 <- mmDS(B1, method = "dream2",  n_cells = 10, n_samples = 2,  min_count = 1, min_cells = 20)
mm_MvsH <- mmDS(BB[gene_id, ], method = "nbinom", coef = 'group_idModerate')
mm_SvsH <- mmDS(BB[gene_id, ], method = "nbinom", coef = 'group_idSevere')
mm_HC <- mmDS(BB2[gene_id, ], method = "nbinom")
mm_MS <- mmDS(BB1[gene_id, ], method = "nbinom")
A = cbind(mm_MvsH$PBMC[, c(1, 6)], mm_SvsH$PBMC[,  6], mm_HC$PBMC[, 6], mm_MS$PBMC[, 6])
colnames(A) = c('gene', 'Moderate_vs_Healthey', 'Severe_vs_Healthy','Health_vs_COVID', 'Moderate_vs_Severe')








# want individal C15, C16, C20, C28,
###############################################################
## Subset cells to cluster  15
################################################################
cluster_id = 14
id = T.integrated$seurat_clusters %in% cluster_id
T_sub = subset(T.integrated, cells = colnames(T.integrated)[id])
Dis =  data.frame(Covid = T_sub$Disease, row.names = colnames(T_sub))
Dis[Dis != 'Healthy'] = 'COVID'
T_sub <- AddMetaData(object = T_sub, metadata = Dis)


pdf('T_infl_cluster_15_Samples_Vlnplot_03_04_2021.pdf',  width = 60, height = 20)
  print(VlnPlot(T_sub, features = inflam_g, assay = 'RNA', group.by = 'Sample_f',  pt.size = 0.02))
dev.off()

pdf('T_infl_cluster_15_Disease_Vlnplot_03_04_2021.pdf', width = 40, height = 20)
  VlnPlot(T_sub, features = inflam_g, assay = 'RNA', group.by = 'Disease',  pt.size = 0.02)
dev.off()

B <- as.SingleCellExperiment(T_sub, assay = 'RNA')
B$cluster = 'PBMC'

qc <- perCellQCMetrics(B)
# remove cells with few or many detected genes
ol <- isOutlier(metric = qc$detected, nmads = 2, log = TRUE)
#B <- B[, !ol]
dim(B)

n = dim(B)[2]
gene_id =  match(inflam_g, rownames(B))
tmp_id = rowSums(counts(B) > .1) > (n * .1)  # filter genes with high expression
tmp_id[gene_id] = TRUE
B <- B[tmp_id, ]
dim(B)

id = B$Disease == 'Healthy'
B1 = B[, !id]


(BB <- prepSCE(B, 
    kid = "cluster", # subpopulation assignments
    gid = "Disease",  # group IDs (ctrl/stim)
    sid = "Sample",   # sample IDs (ctrl/stim.1234)
    drop = TRUE))  # drop all other colData columns

(BB1 <- prepSCE(B1, 
    kid = "cluster", # subpopulation assignments
    gid = "Disease",  # group IDs (ctrl/stim)
    sid = "Sample",   # sample IDs (ctrl/stim.1234)
    drop = TRUE))  # drop all other colData columns


(BB2 <- prepSCE(B, 
    kid = "cluster", # subpopulation assignments
    gid = "Covid",  # group IDs (ctrl/stim)
    sid = "Sample",   # sample IDs (ctrl/stim.1234)
    drop = TRUE))  # drop all other colData columns

clusters <- rep(1, dim(B)[2])
#B <- computeSumFactors(B, clusters=clusters)--
BB <- computeLibraryFactors(BB)
BB1 <- computeLibraryFactors(BB1)
BB2 <- computeLibraryFactors(BB2)
#B <- logNormCounts(B)
summary(sizeFactors(B))
#assays(B)$vstresiduals <- vst(counts(B), show_progress = FALSE)$y

gene_id =  match(inflam_g, rownames(BB))
B1 = B[gene_id, ]

# Can use this to check level colnames(model.matrix(~group_id, data=as.data.frame(colData(BB))))
mm1 <- mmDS(B1, method = "dream2",  n_cells = 10, n_samples = 2,  min_count = 1, min_cells = 20)
mm_MvsH <- mmDS(BB[gene_id, ], method = "nbinom", coef = 'group_idModerate')
mm_SvsH <- mmDS(BB[gene_id, ], method = "nbinom", coef = 'group_idSevere')
mm_HC <- mmDS(BB2[gene_id, ], method = "nbinom")
mm_MS <- mmDS(BB1[gene_id, ], method = "nbinom")
A = cbind(mm_MvsH$PBMC[, c(1, 6)], mm_SvsH$PBMC[,  6], mm_HC$PBMC[, 6], mm_MS$PBMC[, 6])
colnames(A) = c('gene', 'Moderate_vs_Healthey', 'Severe_vs_Healthy','Health_vs_COVID', 'Moderate_vs_Severe')


# want individal  C16, C20, C28,
###############################################################
## Subset cells to cluster  16
################################################################
cluster_id = 15
id = T.integrated$seurat_clusters %in% cluster_id
T_sub = subset(T.integrated, cells = colnames(T.integrated)[id])
Dis =  data.frame(Covid = T_sub$Disease, row.names = colnames(T_sub))
Dis[Dis != 'Healthy'] = 'COVID'
T_sub <- AddMetaData(object = T_sub, metadata = Dis)

VlnPlot(T_sub, features = inflam_g[i], assay = 'RNA', group.by = 'Disease',  pt.size = 0.02)

pdf('T_infl_cluster_16_Samples_Vlnplot_03_04_2021.pdf',  width = 60, height = 20)
  print(VlnPlot(T_sub, features = inflam_g, assay = 'RNA', group.by = 'Sample_f',  pt.size = 0.02))
dev.off()

pdf('T_infl_cluster_16_Disease_Vlnplot_03_04_2021.pdf', width = 40, height = 20)
  VlnPlot(T_sub, features = inflam_g, assay = 'RNA', group.by = 'Disease',  pt.size = 0.02)
dev.off()

B <- as.SingleCellExperiment(T_sub, assay = 'RNA')
B$cluster = 'PBMC'

qc <- perCellQCMetrics(B)
# remove cells with few or many detected genes
ol <- isOutlier(metric = qc$detected, nmads = 2, log = TRUE)
#B <- B[, !ol]
dim(B)

n = dim(B)[2]
gene_id =  match(inflam_g, rownames(B))
tmp_id = rowSums(counts(B) > .1) > (n * .1)  # filter genes with high expression
tmp_id[gene_id] = TRUE
B <- B[tmp_id, ]
dim(B)

id = B$Disease == 'Healthy'
B1 = B[, !id]


(BB <- prepSCE(B, 
    kid = "cluster", # subpopulation assignments
    gid = "Disease",  # group IDs (ctrl/stim)
    sid = "Sample",   # sample IDs (ctrl/stim.1234)
    drop = TRUE))  # drop all other colData columns

(BB1 <- prepSCE(B1, 
    kid = "cluster", # subpopulation assignments
    gid = "Disease",  # group IDs (ctrl/stim)
    sid = "Sample",   # sample IDs (ctrl/stim.1234)
    drop = TRUE))  # drop all other colData columns


(BB2 <- prepSCE(B, 
    kid = "cluster", # subpopulation assignments
    gid = "Covid",  # group IDs (ctrl/stim)
    sid = "Sample",   # sample IDs (ctrl/stim.1234)
    drop = TRUE))  # drop all other colData columns

clusters <- rep(1, dim(B)[2])
#B <- computeSumFactors(B, clusters=clusters)--
BB <- computeLibraryFactors(BB)
BB1 <- computeLibraryFactors(BB1)
BB2 <- computeLibraryFactors(BB2)
#B <- logNormCounts(B)
summary(sizeFactors(B))
#assays(B)$vstresiduals <- vst(counts(B), show_progress = FALSE)$y

gene_id =  match(inflam_g, rownames(BB))
B1 = B[gene_id, ]

# Can use this to check level colnames(model.matrix(~group_id, data=as.data.frame(colData(BB))))
mm1 <- mmDS(B1, method = "dream2",  n_cells = 10, n_samples = 2,  min_count = 1, min_cells = 20)
mm_MvsH <- mmDS(BB[gene_id, ], method = "nbinom", coef = 'group_idModerate')
mm_SvsH <- mmDS(BB[gene_id, ], method = "nbinom", coef = 'group_idSevere')
mm_HC <- mmDS(BB2[gene_id, ], method = "nbinom")
mm_MS <- mmDS(BB1[gene_id, ], method = "nbinom")
A = cbind(mm_MvsH$PBMC[, c(1, 6)], mm_SvsH$PBMC[,  6], mm_HC$PBMC[, 6], mm_MS$PBMC[, 6])
colnames(A) = c('gene', 'Moderate_vs_Healthey', 'Severe_vs_Healthy','Health_vs_COVID', 'Moderate_vs_Severe')



# want individal  C20, C28,
###############################################################
## Subset cells to cluster  20
################################################################
cluster_id = 19
id = T.integrated$seurat_clusters %in% cluster_id
T_sub = subset(T.integrated, cells = colnames(T.integrated)[id])
Dis =  data.frame(Covid = T_sub$Disease, row.names = colnames(T_sub))
Dis[Dis != 'Healthy'] = 'COVID'
T_sub <- AddMetaData(object = T_sub, metadata = Dis)


pdf('T_infl_cluster_20_Samples_Vlnplot_03_04_2021.pdf',  width = 60, height = 20)
  print(VlnPlot(T_sub, features = inflam_g, assay = 'RNA', group.by = 'Sample_f',  pt.size = 0.02))
dev.off()

pdf('T_infl_cluster_20_Disease_Vlnplot_03_04_2021.pdf', width = 40, height = 20)
  VlnPlot(T_sub, features = inflam_g, assay = 'RNA', group.by = 'Disease',  pt.size = 0.02)
dev.off()

B <- as.SingleCellExperiment(T_sub, assay = 'RNA')
B$cluster = 'PBMC'

qc <- perCellQCMetrics(B)
# remove cells with few or many detected genes
ol <- isOutlier(metric = qc$detected, nmads = 2, log = TRUE)
#B <- B[, !ol]
dim(B)

n = dim(B)[2]
gene_id =  match(inflam_g, rownames(B))
tmp_id = rowSums(counts(B) > .1) > (n * .1)  # filter genes with high expression
tmp_id[gene_id] = TRUE
B <- B[tmp_id, ]
dim(B)

id = B$Disease == 'Healthy'
B1 = B[, !id]


(BB <- prepSCE(B, 
    kid = "cluster", # subpopulation assignments
    gid = "Disease",  # group IDs (ctrl/stim)
    sid = "Sample",   # sample IDs (ctrl/stim.1234)
    drop = TRUE))  # drop all other colData columns

(BB1 <- prepSCE(B1, 
    kid = "cluster", # subpopulation assignments
    gid = "Disease",  # group IDs (ctrl/stim)
    sid = "Sample",   # sample IDs (ctrl/stim.1234)
    drop = TRUE))  # drop all other colData columns


(BB2 <- prepSCE(B, 
    kid = "cluster", # subpopulation assignments
    gid = "Covid",  # group IDs (ctrl/stim)
    sid = "Sample",   # sample IDs (ctrl/stim.1234)
    drop = TRUE))  # drop all other colData columns

clusters <- rep(1, dim(B)[2])
#B <- computeSumFactors(B, clusters=clusters)--
BB <- computeLibraryFactors(BB)
BB1 <- computeLibraryFactors(BB1)
BB2 <- computeLibraryFactors(BB2)
#B <- logNormCounts(B)
summary(sizeFactors(B))
#assays(B)$vstresiduals <- vst(counts(B), show_progress = FALSE)$y

gene_id =  match(inflam_g, rownames(BB))
B1 = B[gene_id, ]

# Can use this to check level colnames(model.matrix(~group_id, data=as.data.frame(colData(BB))))
mm1 <- mmDS(B1, method = "dream2",  n_cells = 10, n_samples = 2,  min_count = 1, min_cells = 20)
mm_MvsH <- mmDS(BB[gene_id, ], method = "nbinom", coef = 'group_idModerate')
mm_SvsH <- mmDS(BB[gene_id, ], method = "nbinom", coef = 'group_idSevere')
mm_HC <- mmDS(BB2[gene_id, ], method = "nbinom")
mm_MS <- mmDS(BB1[gene_id, ], method = "nbinom")
A = cbind(mm_MvsH$PBMC[, c(1, 6)], mm_SvsH$PBMC[,  6], mm_HC$PBMC[, 6], mm_MS$PBMC[, 6])
colnames(A) = c('gene', 'Moderate_vs_Healthey', 'Severe_vs_Healthy','Health_vs_COVID', 'Moderate_vs_Severe')



# want individal  C28,
###############################################################
## Subset cells to cluster  28
################################################################
cluster_id = 27
id = T.integrated$seurat_clusters %in% cluster_id
T_sub = subset(T.integrated, cells = colnames(T.integrated)[id])
Dis =  data.frame(Covid = T_sub$Disease, row.names = colnames(T_sub))
Dis[Dis != 'Healthy'] = 'COVID'
T_sub <- AddMetaData(object = T_sub, metadata = Dis)


pdf('T_infl_cluster_28_Samples_Vlnplot_03_04_2021.pdf',  width = 60, height = 20)
  print(VlnPlot(T_sub, features = inflam_g, assay = 'RNA', group.by = 'Sample_f',  pt.size = 0.02))
dev.off()

pdf('T_infl_cluster_28_Disease_Vlnplot_03_04_2021.pdf', width = 40, height = 20)
  VlnPlot(T_sub, features = inflam_g, assay = 'RNA', group.by = 'Disease',  pt.size = 0.02)
dev.off()

B <- as.SingleCellExperiment(T_sub, assay = 'RNA')
B$cluster = 'PBMC'

qc <- perCellQCMetrics(B)
# remove cells with few or many detected genes
ol <- isOutlier(metric = qc$detected, nmads = 2, log = TRUE)
#B <- B[, !ol]
dim(B)

n = dim(B)[2]
gene_id =  match(inflam_g, rownames(B))
tmp_id = rowSums(counts(B) > .1) > (n * .1)  # filter genes with high expression
tmp_id[gene_id] = TRUE
B <- B[tmp_id, ]
dim(B)

id = B$Disease == 'Healthy'
B1 = B[, !id]


(BB <- prepSCE(B, 
    kid = "cluster", # subpopulation assignments
    gid = "Disease",  # group IDs (ctrl/stim)
    sid = "Sample",   # sample IDs (ctrl/stim.1234)
    drop = TRUE))  # drop all other colData columns

(BB1 <- prepSCE(B1, 
    kid = "cluster", # subpopulation assignments
    gid = "Disease",  # group IDs (ctrl/stim)
    sid = "Sample",   # sample IDs (ctrl/stim.1234)
    drop = TRUE))  # drop all other colData columns


(BB2 <- prepSCE(B, 
    kid = "cluster", # subpopulation assignments
    gid = "Covid",  # group IDs (ctrl/stim)
    sid = "Sample",   # sample IDs (ctrl/stim.1234)
    drop = TRUE))  # drop all other colData columns

clusters <- rep(1, dim(B)[2])
#B <- computeSumFactors(B, clusters=clusters)--
BB <- computeLibraryFactors(BB)
BB1 <- computeLibraryFactors(BB1)
BB2 <- computeLibraryFactors(BB2)
#B <- logNormCounts(B)
summary(sizeFactors(B))
#assays(B)$vstresiduals <- vst(counts(B), show_progress = FALSE)$y

gene_id =  match(inflam_g, rownames(BB))
B1 = B[gene_id, ]

# Can use this to check level colnames(model.matrix(~group_id, data=as.data.frame(colData(BB))))
mm1 <- mmDS(B1, method = "dream2",  n_cells = 10, n_samples = 2,  min_count = 1, min_cells = 20)
mm_MvsH <- mmDS(BB[gene_id, ], method = "nbinom", coef = 'group_idModerate')
mm_SvsH <- mmDS(BB[gene_id, ], method = "nbinom", coef = 'group_idSevere')
mm_HC <- mmDS(BB2[gene_id, ], method = "nbinom")
mm_MS <- mmDS(BB1[gene_id, ], method = "nbinom")
A = cbind(mm_MvsH$PBMC[, c(1, 6)], mm_SvsH$PBMC[,  6], mm_HC$PBMC[, 6], mm_MS$PBMC[, 6])
colnames(A) = c('gene', 'Moderate_vs_Healthey', 'Severe_vs_Healthy','Health_vs_COVID', 'Moderate_vs_Severe')




# want individal  C28,
###############################################################
## Subset cells to cluster  30
################################################################
cluster_id = 29
id = T.integrated$seurat_clusters %in% cluster_id
T_sub = subset(T.integrated, cells = colnames(T.integrated)[id])
Dis =  data.frame(Covid = T_sub$Disease, row.names = colnames(T_sub))
Dis[Dis != 'Healthy'] = 'COVID'
T_sub <- AddMetaData(object = T_sub, metadata = Dis)


pdf('T_infl_cluster_30_Samples_Vlnplot_03_14_2021.pdf',  width = 60, height = 20)
  print(VlnPlot(T_sub, features = inflam_g, assay = 'RNA', group.by = 'Sample_f',  pt.size = 0.02))
dev.off()

pdf('T_infl_cluster_30_Disease_Vlnplot_03_04_2021.pdf', width = 40, height = 20)
  VlnPlot(T_sub, features = inflam_g, assay = 'RNA', group.by = 'Disease',  pt.size = 0.02)
dev.off()

B <- as.SingleCellExperiment(T_sub, assay = 'RNA')
B$cluster = 'PBMC'

qc <- perCellQCMetrics(B)
# remove cells with few or many detected genes
ol <- isOutlier(metric = qc$detected, nmads = 2, log = TRUE)
#B <- B[, !ol]
dim(B)

n = dim(B)[2]
gene_id =  match(inflam_g, rownames(B))
tmp_id = rowSums(counts(B) > .1) > (n * .1)  # filter genes with high expression
tmp_id[gene_id] = TRUE
B <- B[tmp_id, ]
dim(B)

id = B$Disease == 'Healthy'
B1 = B[, !id]


(BB <- prepSCE(B, 
    kid = "cluster", # subpopulation assignments
    gid = "Disease",  # group IDs (ctrl/stim)
    sid = "Sample",   # sample IDs (ctrl/stim.1234)
    drop = TRUE))  # drop all other colData columns

(BB1 <- prepSCE(B1, 
    kid = "cluster", # subpopulation assignments
    gid = "Disease",  # group IDs (ctrl/stim)
    sid = "Sample",   # sample IDs (ctrl/stim.1234)
    drop = TRUE))  # drop all other colData columns


(BB2 <- prepSCE(B, 
    kid = "cluster", # subpopulation assignments
    gid = "Covid",  # group IDs (ctrl/stim)
    sid = "Sample",   # sample IDs (ctrl/stim.1234)
    drop = TRUE))  # drop all other colData columns

clusters <- rep(1, dim(B)[2])
#B <- computeSumFactors(B, clusters=clusters)--
BB <- computeLibraryFactors(BB)
BB1 <- computeLibraryFactors(BB1)
BB2 <- computeLibraryFactors(BB2)
#B <- logNormCounts(B)
summary(sizeFactors(B))
#assays(B)$vstresiduals <- vst(counts(B), show_progress = FALSE)$y

gene_id =  match(inflam_g, rownames(BB))
B1 = B[gene_id, ]

# Can use this to check level colnames(model.matrix(~group_id, data=as.data.frame(colData(BB))))
mm1 <- mmDS(B1, method = "dream2",  n_cells = 10, n_samples = 2,  min_count = 1, min_cells = 20)
mm_MvsH <- mmDS(BB[gene_id, ], method = "nbinom", coef = 'group_idModerate')
mm_SvsH <- mmDS(BB[gene_id, ], method = "nbinom", coef = 'group_idSevere')
mm_HC <- mmDS(BB2[gene_id, ], method = "nbinom")
mm_MS <- mmDS(BB1[gene_id, ], method = "nbinom")
A = cbind(mm_MvsH$PBMC[, c(1, 6)], mm_SvsH$PBMC[,  6], mm_HC$PBMC[, 6], mm_MS$PBMC[, 6])
colnames(A) = c('gene', 'Moderate_vs_Healthey', 'Severe_vs_Healthy','Health_vs_COVID', 'Moderate_vs_Severe')
