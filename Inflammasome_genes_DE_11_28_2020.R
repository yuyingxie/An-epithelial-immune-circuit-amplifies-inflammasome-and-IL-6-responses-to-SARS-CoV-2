library(Seurat)
library(scales)
library(scater)
library(muscat)
library(scran)

# Check the inflammasome in C1, C9, C11, C15, C16, C20, C28, C30
# Inflammasome gene list: GSDMD, NLRP3, NLRP1, PYCARD, CASP1, CASP4, CASP5, IL1B, NLRC4, AIM2, IL18, 
inflam_g = c('GSDMD', 'NLRP3', 'NLRP1', 'PYCARD', 'CASP1', 
             'CASP4', 'CASP5', 'IL1B',   'NLRC4', 'AIM2', 
             'IL18')

cluster_id = c(0, 8, 10, 14, 15, 19, 27, 29)
T.integrated = readRDS('T_SCT_integrated_11_26_2020.rds')

#saveRDS(T.integrated, 'T_SCT_integrated_11_26_2020.rds')

DefaultAssay(T.integrated) <- "RNA"

for(i in 1:11){
     pdf(paste('T_', inflam_g[i], 'feature_plot_12_03_2020.pdf'),
     width = 8, height = 8)
    print(FeaturePlot(T.integrated, features = inflam_g[i]))
    dev.off()
    cat('i \n')
}


T.integrated$Sample_f = factor(T.integrated$Sample, 
                    levels = c('L_H1', 'L_H2', 'L_H3', 'L_H4', 'A_H1',
                              'A_H2', 'A_H3', 'A_H4', 'A_H5', 'W_H1',
                              'W_H2', 'W_H3', 'W_H4', 'W_H5', 'W_H6',
                              'L_M2', 'L_M4', 'L_M5', 'L_M7_2', 'L_M8',
                              'A_M2', 'A_M3', 'A_M7', 'W_M7', 'L_S1',
                              'L_S3_1', 'L_S3_2', 'L_S6_1', 'L_S6_2', 'L_S7_1',
                              'A_S1', 'A_S4', 'A_S5', 'A_S6', 'W_S1_1',
                              'W_S1_2', 'W_S2', 'W_S3', 'W_S4', 'W_S6'))

pdf(paste('T_genes_all_clusters_samples_12_01_2020.pdf'),  width = 60, height = 20)
print(VlnPlot(T.integrated, features = inflam_g, assay = 'RNA', group.by = 'Sample_f',  pt.size = 0.02))
dev.off()

pdf(paste('T_genes_all_clusters_disease_12_01_2020.pdf'),  width = 60, height = 20)
print(VlnPlot(T.integrated, features = inflam_g, assay = 'RNA', group.by = 'Disease',  pt.size = 0.02))
dev.off()


FeaturePlot(T.integrated, feature = inflam_g[1] )

for(i in 1:11){
    pdf(paste('T_', inflam_g[i], '_counts_Gene_11_28_2020.pdf'),
     width = 6, height = 6)
    print(FeaturePlot(T.integrated, features = inflam_g[i]))
    dev.off()
    cat(i, '\n')
}

for(i in 1:11){
    pdf(paste('T_', inflam_g[i], '_all_clusters_Vlnplot_11_28_2020.pdf'),
     width = 8, height = 8)
    print(VlnPlot(T.integrated, features = inflam_g[i], assay = 'RNA', group.by = 'Disease',  pt.size = 0.02))
    dev.off()
    cat(i, '\n')
}

##############################################################
## Subset cells to clusters (0, 8, 10, 14, 15, 19, 27, 29)
##############################################################
cluster_id = c(0, 8, 10, 14, 15, 19, 27, 29)
id = T.integrated$seurat_clusters %in% cluster_id
inflam_g = c('GSDMD', 'NLRP3', 'NLRP1', 'PYCARD', 'CASP1', 
             'CASP4', 'CASP5', 'IL1B',   'NLRC4', 'AIM2', 
             'IL18')

#T_sub = subset(T.integrated, cells = colnames(T.integrated)[id])
T_sub$Sample_f = factor(T_sub$Sample, 
                    levels = c('L_H1', 'L_H2', 'L_H3', 'L_H4', 'A_H1',
                              'A_H2', 'A_H3', 'A_H4', 'A_H5', 'W_H1',
                              'W_H2', 'W_H3', 'W_H4', 'W_H5', 'W_H6',
                              'L_M2', 'L_M4', 'L_M5', 'L_M7_2', 'L_M8',
                              'A_M2', 'A_M3', 'A_M7', 'W_M7', 'L_S1',
                              'L_S3_1', 'L_S3_2', 'L_S6_1', 'L_S6_2', 'L_S7_1',
                              'A_S1', 'A_S4', 'A_S5', 'A_S6', 'W_S1_1',
                              'W_S1_2', 'W_S2', 'W_S3', 'W_S4', 'W_S6'))


#saveRDS(T_sub, file ='T_SCT_int_sub_11_28_2020.rds')
T_sub = readRDS('T_SCT_int_sub_11_28_2020.rds')

pdf(paste('T_genes_clusters_samples_11_29_2020.pdf'),  width = 60, height = 20)
print(VlnPlot(T_sub, features = inflam_g, assay = 'RNA', group.by = 'Sample_f',  pt.size = 0.02))
dev.off()

pdf(paste('T_genes_8clusters_Vlnplot_11_29_2020.pdf'),  width = 40, height = 20)
    print(VlnPlot(T_sub, features = inflam_g, assay = 'RNA', group.by = 'Disease',  pt.size = 0.02))
dev.off()

########################################################################################################
#########  Sub 7
#######################################################################################################
# only for  Subset cells to clusters (0, 8, 10, 14, 15, 19, 27)
cluster_id = c(0, 8, 10, 14, 15, 19, 27)
id = T_sub$seurat_clusters %in% cluster_id
inflam_g = c('GSDMD', 'NLRP3', 'NLRP1', 'PYCARD', 'CASP1', 
             'CASP4', 'CASP5', 'IL1B',   'NLRC4', 'AIM2', 
             'IL18')

#T_sub7 = subset(T_sub, cells = colnames(T_sub)[id])
#saveRDS(T_sub7, file ='T_SCT_int_sub_7clusters_11_28_2020.rds')
T_sub7 = readRDS('T_SCT_int_sub_7clusters_11_28_2020.rds')

pdf(paste('T_genes_7clusters_samples_11_29_2020.pdf'),  width = 60, height = 20)
print(VlnPlot(T_sub7, features = inflam_g, assay = 'RNA', group.by = 'Sample_f',  pt.size = 0.02))
dev.off()


pdf(paste('T_genes_7clusters_Vlnplot_11_29_2020.pdf'),  width = 40, height = 20)
print(VlnPlot(T_sub7, features = inflam_g, assay = 'RNA', group.by = 'Disease',  pt.size = 0.02))
dev.off()

for(i in 1:11){
    pdf(paste('T_', inflam_g[i], '_7clusters_Vlnplot_11_29_2020.pdf'),
     width = 8, height = 8)
    print(VlnPlot(T_sub7, features = inflam_g[i], assay = 'RNA', group.by = 'Disease',  pt.size = 0.02))
    dev.off()
    cat(i, '\n')
}

###################################################################
### Only for cluster_id = c(8, 10, 14)
##################################################################

cluster_id = c(8, 10, 14)
id = T_sub7$seurat_clusters %in% cluster_id
T_sub3 = subset(T_sub7, cells = colnames(T_sub7)[id])
saveRDS(T_sub3, file ='T_SCT_int_sub_c8_10_14_11_28_2020.rds')
T_sub3 = readRDS('T_SCT_int_sub_c8_10_14_11_28_2020.rds')

pdf(paste('T_genes_3clusters_samples_11_29_2020.pdf'),  width = 60, height = 20)
print(VlnPlot(T_sub3, features = inflam_g, assay = 'RNA', group.by = 'Sample_f',  pt.size = 0.02))
dev.off()

pdf(paste('T_genes_3clusters_C9_C11_C15_Vlnplot_11_29_2020.pdf'),  width = 40, height = 20)
print(VlnPlot(T_sub7, features = inflam_g, assay = 'RNA', group.by = 'Disease',  pt.size = 0.02))
dev.off()


for(i in 1:11){
    pdf(paste('T_', inflam_g[i], '_3clusters_Vlnplot_11_29_2020.pdf'),
     width = 8, height = 8)
    print(VlnPlot(T_sub3, features = inflam_g[i], assay = 'RNA', group.by = 'Disease',  pt.size = 0.02))
    dev.off()
    cat(i, '\n')
}



###################################################################
### Only for cluster_id = 5
##################################################################

cluster_id = 5
id = T.integrated$seurat_clusters %in% cluster_id
T_5 = subset(T.integrated, cells = colnames(T.integrated)[id])

pdf(paste('T_genes_cluster_5_samples_12_01_2020.pdf'),  width = 60, height = 20)
print(VlnPlot(T_5, features = inflam_g, assay = 'RNA', group.by = 'Sample_f',  pt.size = 0.02))
dev.off()

pdf(paste('T_genes_cluster_5_Disease_12_01_2020.pdf'),  width = 40, height = 20)
print(VlnPlot(T_5, features = inflam_g, assay = 'RNA', group.by = 'Disease',  pt.size = 0.02))
dev.off()


###################################################################
### Only for cluster_id = 12
##################################################################

cluster_id = 12
id = T.integrated$seurat_clusters %in% cluster_id
T_12 = subset(T.integrated, cells = colnames(T.integrated)[id])

pdf(paste('T_genes_cluster_12_samples_12_01_2020.pdf'),  width = 60, height = 20)
print(VlnPlot(T_12, features = inflam_g, assay = 'RNA', group.by = 'Sample_f',  pt.size = 0.02))
dev.off()

pdf(paste('T_genes_cluster_12_Disease_12_01_2020.pdf'),  width = 40, height = 20)
print(VlnPlot(T_12, features = inflam_g, assay = 'RNA', group.by = 'Disease',  pt.size = 0.02))
dev.off()


####################################################
# for individual cluster
########################################################

for(i in 0:29){
    cluster_id = i
    id = T.integrated$seurat_clusters == cluster_id
    T_tmp = subset(T.integrated, cells = colnames(T.integrated)[id])

/wn    cat(i, '\n')
}3n1/=4
=3t
1dr
id = T.integrated$Experiment == 'Lee'
T_tmp = subset(T.integrated, cells = colnames(T.integrated)[id])

E = unique(T.integrated$Experiment)
for(i in 1:3){
    E_id = E[i]
    id = T.integrated$Experiment == E[i]
    T_tmp = subset(T.integrated, cells = colnames(T.integrated)[id])
    for(j in 0:29){
    cluster_id = j
    id = T_tmp$seurat_clusters == cluster_id
    T_tmp1 = subset(T_tmp, cells = colnames(T_tmp)[id])
    #pdf(paste('T_infl_E', E[i],'_cluster', j + 1,'_Samples_Vlnplot_12_01_2020.pdf'),
    # width = 60, height = 20)
    png(paste0('T_infl_E', E[i],'_cluster', j + 1,'_Samples_Vlnplot_12_01_2020.png'),
     width = 30*300, height = 10* 300)    
    print(VlnPlot(T_tmp1, features = inflam_g, assay = 'RNA', group.by = 'Sample_f',  pt.size = 0.02))
    dev.off()
#   pdf(paste('T_infl_E', E[i] ,'_cluster', j + 1,'_Disease_Vlnplot_12_01_2020.pdf'),
#     width = 4 * 300, height = 2 * 300)
    png(paste0('T_infl_E', E[i],'_cluster', j + 1,'_Disease_Vlnplot_12_01_2020.png'),
     width = 20*300, height = 10* 300)

    print(VlnPlot(T_tmp1, features = inflam_g, assay = 'RNA', group.by = 'Disease',  pt.size = 0.02))
    dev.off()
    cat(i, '\n')
    }
}


T_sub = readRDS('T_SCT_int_sub_11_28_2020.rds')
cluster_id = 29

# for individual cluster
T_sub = readRDS('T_SCT_int_sub_11_28_2020.rds')
cluster_id = c(0, 8, 10, 14, 15, 19, 27, 29)

for(i in 1:8){
    id = T_sub$seurat_clusters == cluster_id[i]
    T_tmp = subset(T_sub, cells = colnames(T_sub)[id])
    pdf(paste0('T_genes_cluster', cluster_id[i] + 1,'_Disease_Vlnplot_11_29_2020.pdf'),
     width = 40, height = 20)
    print(VlnPlot(T_tmp, features = inflam_g, assay = 'RNA', group.by = 'Disease',  pt.size = 0.02))
    dev.off()
    pdf(paste0('T_genes_cluster', cluster_id[i] + 1,'_Sample_Vlnplot_11_29_2020.pdf'),
     width = 60, height = 20)
    print(VlnPlot(T_tmp, features = inflam_g, assay = 'RNA', group.by = 'Sample_f',  pt.size = 0.02))
    dev.off()
   
    cat('(', i, ')', '\n')
}


print(VlnPlot(T_sub, features = inflam_g[j], assay = 'RNA', group.by = 'Sample_f',  pt.size = 0.02))

i = 4
print(VlnPlot(T_sub, features = inflam_g[i], assay = 'RNA', group.by = 'Sample_f',  pt.size = 0.02))

FeatureScatter(T_sub, feature1 = inflam_g[1], feature2 = inflam_g[i], group.by = 'Sample_f')

id = T_sub$seurat_clusters == cluster_id[i]


for(i in 1:8){
    id = T_sub$seurat_clusters == cluster_id[i]
    T_tmp = subset(T_sub, cells = colnames(T_sub)[id])
    for(j in 1:11){
    pdf(paste0('T_Sample', inflam_g[j], '_cluster', cluster_id[i] + 1,'_Vlnplot_11_29_2020.pdf'),
     width = 16, height = 8)
    print(VlnPlot(T_tmp, features = inflam_g[j], assay = 'RNA', group.by = 'Sample_f',  pt.size = 0.02))
    dev.off()
    cat('(', i, j, ')', '\n')
}
}

Sample_id = T_sub$Sample

for(i in 1:1){
    id = T_sub$seurat_clusters == cluster_id[i]
    T_tmp = subset(T_sub, cells = colnames(T_sub)[id])
    for(j in 2:3){
        for(g in 1:40){ 
          id1 = T_tmp$Sample == Sample_id[g]
          T_tmp1 = subset(T_tmp, cells = colnames(T_tmp)[id1])
        pdf(paste0('T_Sample_', inflam_g[j], '_cluster_', cluster_id[i] + 1,'Sample', Sample_id[g],'_Vlnplot_11_29_2020.pdf'), width = 16, height = 8)
    print(FeatureScatter(T_tmp1, feature1 = inflam_g[1],  feature2 = inflam_g[j],))
    dev.off()
    cat('(', i, j, g, ')', '\n')
        }
}
}


##################################################
## W only
######################################################
T_sub7 = readRDS('T_SCT_int_sub_7clusters_11_28_2020.rds')
id = T_sub7$Experiment == 'W'
T_W = subset(T_sub7, cells = colnames(T_sub7)[id])
T_tmp = subset(T_sub7, cells = colnames(T_sub7)[!id])

print(VlnPlot(T_W, features = inflam_g[j], assay = 'RNA', group.by = 'Disease',  pt.size = 0.02))
dat_W = GetAssayData(T_W, slot = 'data')
dat_W = dat_W[inflam_g, ]
dat_W_H = as.matrix(dat_W[, T_W$Disease == 'Healthy'  ])
dat_W_M = as.matrix(dat_W[, T_W$Disease == 'Moderate'  ])
dat_W_S = as.matrix(dat_W[, T_W$Disease == 'Severe'  ])

for(i in 1:11){
    dat_tmp = data.frame(gene = c(dat_W_H[i, ], dat_W_M[i, ],  dat_W_S[i, ]), 
    Disease = c(rep('H', dim(dat_W_H)[2]), rep('M', dim(dat_W_M)[2]), rep('S', dim(dat_W_S)[2])  ))
    boxplot(gene ~ Disease, data = dat_tmp)

}

boxplot(gene ~ Disease, data = dat_tmp)

id = T_sub7$Experiment == 'A'
T_A = subset(T_sub7, cells = colnames(T_sub7)[id])
print(VlnPlot(T_A, features = inflam_g[j], assay = 'RNA', group.by = 'Disease',  pt.size = 0.02))



dat1 = GetAssayData(T_tmp, slot = 'data')
dat1 = as.matrix(dat1[inflam_g, ])
dat1_H = as.matrix(dat1[, T_tmp$Disease == 'Healthy'  ])
dat1_M = as.matrix(dat1[, T_tmp$Disease == 'Moderate'  ])
dat1_S = as.matrix(dat1[, T_tmp$Disease == 'Severe'  ])

for(i in 1:11){
    dat_tmp = data.frame(gene = c(dat1_H[i, ], dat1_M[i, ],  dat1_S[i, ]), 
    Disease = c(rep('H', dim(dat_W_H)[2]), rep('M', dim(dat_W_M)[2]), rep('S', dim(dat_W_S)[2])  ))
    boxplot(gene ~ Disease, data = dat_tmp)
}


dat_tmp = data.frame(gene = c(apply(dat1_H, 2, sum), apply(dat1_M, 2, sum),  apply(dat1_S, 2, sum)), 
    Disease = c(rep('H', dim(dat_W_H)[2]), rep('M', dim(dat_W_M)[2]), rep('S', dim(dat_W_S)[2])  ))
  boxplot(gene ~ Disease, data = dat_tmp)


dat_tmp = data.frame(gene = c(apply(dat1_H > 0.001, 2, sum), apply(dat1_M > 0.001, 2, sum),  apply(dat1_S > 0.001, 2, sum)), 
    Disease = c(rep('H', dim(dat_W_H)[2]), rep('M', dim(dat_W_M)[2]), rep('S', dim(dat_W_S)[2])  ))
  boxplot(gene ~ Disease, data = dat_tmp)


# PCA 

pca_res <- prcomp(t(dat1), scale. = FALSE, center = T)


A = t(as.matrix(dat_W))


#B <- as.SingleCellExperiment(T_sub, assay = 'RNA')
B <- as.SingleCellExperiment(T_sub3, assay = 'RNA')
B <- as.SingleCellExperiment(T_sub7, assay = 'RNA')
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
gene_list = rownames(B) 
save(gene_list, file = 'Epi_gene_list_10_24_2020.rda')

(B <- prepSCE(B, 
    kid = "cluster", # subpopulation assignments
    gid = "Disease",  # group IDs (ctrl/stim)
    sid = "Sample",   # sample IDs (ctrl/stim.1234)
    drop = TRUE))  # drop all other colData columns

saveRDS(B, file = 'T_Sub_a_8_clusters_11_28_2020.rds')
B = readRDS('T_Sub_a_8_clusters_11_28_2020.rds')


clusters <- rep(1, dim(B)[2])
B <- computeSumFactors(B, clusters=clusters)
B <- logNormCounts(B)
summary(sizeFactors(B))

gene_id =  match(inflam_g, rownames(B))
B1 = B[gene_id, ]

#mm <- mmDS(B, method = "dream",  n_cells = 10, n_samples = 2,  min_count = 1, min_cells = 20)
mm <- mmDS(B, method = "nbinom")

mm <- mmDS(B1, method = "nbinom")


############################
# test for 3clusters
##################################
T_sub3 = readRDS('T_SCT_int_sub_c8_10_14_11_28_2020.rds')
T_sub3$Sample_f = factor(T_sub3$Sample, 
                    levels = c('L_N1', 'L_N2', 'L_N3', 'L_N4', 'A_N1',
                              'A_N2', 'A_N3', 'A_N4', 'A_N5', 'W_N1',
                              'W_N2', 'W_N3', 'W_N4', 'W_N5', 'W_N6',
                              'L_C2', 'L_C4', 'L_C5', 'L_C7_2', 'L_C8',
                              'A_C2', 'A_C3', 'A_C7', 'W_C7', 'L_C1',
                              'L_C3_1', 'L_C3_2', 'L_C6_1', 'L_C6_2', 'L_C7_1',
                              'A_C1', 'A_C4', 'A_C5', 'A_C6', 'W_C1_1',
                              'W_C1_2', 'W_C2', 'W_C3', 'W_C4', 'W_C6' ) )


B <- as.SingleCellExperiment(T_sub3, assay = 'RNA')
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

#mm <- mmDS(B, method = "dream",  n_cells = 10, n_samples = 2,  min_count = 1, min_cells = 20)
#mm <- mmDS(B, method = "nbinom")

mm <- mmDS(B1, method = "nbinom")
write.csv(mm, file = '3clusters_mm.csv')


# 3 clusters No W
id = T_sub3$Experiment == 'W'
T_sub3_noW = subset(T_sub3, cells = colnames(T_sub3)[!id])

B <- as.SingleCellExperiment(T_sub3_noW, assay = 'RNA')
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

mm <- mmDS(B1, method = "nbinom")
write.csv(mm, file = '3clusters_no_W_mm.csv')


# Only A
id = T_sub3$Experiment == 'A'
T_sub3_A = subset(T_sub3, cells = colnames(T_sub3)[id])

B <- as.SingleCellExperiment(T_sub3_A, assay = 'RNA')
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

mm <- mmDS(B1, method = "nbinom")
write.csv(mm, file = '3clusters_A_mm.csv')

dat = GetAssayData(T_sub3_A, slot = 'data')
dat = dat[inflam_g, ]

res = matrix(0, 11, 1)
rownames(res) = inflam_g
for(i in 1:11){
    tmp = data.frame(clus = dat[i, ], Dis = T_sub3_A$Disease)
    A = kruskal.test(clus ~ Dis, data = tmp)
    res[i] = A$p.value
}
write.csv(res, file = '3clusters_A_KW_test.csv')


# Only L
id = T_sub3$Experiment == 'Lee'
T_sub3_L= subset(T_sub3, cells = colnames(T_sub3)[id])

B <- as.SingleCellExperiment(T_sub3_L, assay = 'RNA')
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

mm <- mmDS(B1, method = "nbinom")
write.csv(mm, file = '3clusters_A_mm.csv')

dat = GetAssayData(T_sub3_A, slot = 'data')
dat = dat[inflam_g, ]

res = matrix(0, 11, 1)
rownames(res) = inflam_g
for(i in 1:11){
    tmp = data.frame(clus = dat[i, ], Dis = T_sub3_A$Disease)
    A = kruskal.test(clus ~ Dis, data = tmp)
    res[i] = A$p.value
}
write.csv(res, file = '3clusters_A_KW_test.csv')




###################################
## Test for cluster 9
############################
T_sub3 = readRDS('T_SCT_int_sub_c8_10_14_11_28_2020.rds')
id = T_sub3$seurat_clusters == 8
T_9 = subset(T_sub3, cells = colnames(T_sub3)[id])

B <- as.SingleCellExperiment(T_9, assay = 'RNA')
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

mm <- mmDS(B1, method = "nbinom")
write.csv(mm, file = '3clusters_mm.csv')

###############################
# Test for cluster 30

T_sub = readRDS('T_SCT_int_sub_11_28_2020.rds')

id = T_sub$seurat_clusters == 29

dat = GetAssayData(T_30, slot = 'data')
dat = dat[inflam_g, ]
res = rep(0, 30)
for(i in 1:30){
    tmp = data.frame(clus = cluster_per[, i], Dis = D)
    A = kruskal.test(clus ~ Dis, data = tmp)
    res[i] = A$p.value
}

res = matrix(0, 11, 1)
rownames(res) = inflam_g
for(i in 1:11){
    tmp = data.frame(clus = dat[i, ], Dis = T_30$Disease)
    A = kruskal.test(clus ~ Dis, data = tmp)
    res[i] = A$p.value
}


