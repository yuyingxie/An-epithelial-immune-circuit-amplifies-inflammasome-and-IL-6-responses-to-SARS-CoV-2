library(Seurat)
library(dplyr)
library(Matrix)
library(data.table)
library(gridExtra)
library(ggplot2)
library(cowplot)
library(ggrastr)
library(matrixStats)
library(parallel)
library(reshape2)
library(ComplexHeatmap)
library(slingshot)
library(qlcMatrix)
set.seed(1)
source( "/mnt/ufs18/home-133/xyy/COVID19_B/seuratWrapper.R")
source('../COVID_Function.R')

#The suffix of each barcode sequence represents patient information.
#Sample  Age  Gender Disease
# 1  cov01	 75	   F	Severe
# 2  cov02  53	   F	Moderate
# 3  cov03  75	   F	Moderate
# 4  cov04  59	   M	Severe
# 5  cov07  84	   F	Healthy

# 6  cov08  68	   F	Healthy
# 7  cov09  38	   M	Healthy
# 8  cov10  60	   F	Severe
# 9  cov11  48	   M	Severe
# 10 cov12  47	   F	Moderate

# 11 cov17  90	   M	Healthy
# 12 cov18  70	   F	Healthy

p_id = c('cov01', 'cov02', 'cov03', 'cov04', 'cov07', 'cov08', 'cov09', 'cov10', 'cov11', 'cov12', 'cov17', 'cov18')
disease_type = c('Severe', 'Moderate', 'Moderate', 'Severe', 'Healthy',
              'Healthy', 'Healthy', 'Severe', 'Severe', 'Moderate',
              'Healthy', 'Healthy')  

Patient_id = c('A_C1', 'A_C2', 'A_C3', 'A_C4', 'A_N1',
               'A_N2', 'A_N3', 'A_C5', 'A_C6', 'A_C7',
                'A_N4', 'A_N5')

Sample_id = c('A_C1', 'A_C2', 'A_C3', 'A_C4', 'A_N1',
               'A_N2', 'A_N3', 'A_C5', 'A_C6', 'A_C7',
                'A_N4', 'A_N5')

Gender_tmp = c('F', 'F', 'F', 'M', 'F', 
               'F', 'M', 'F', 'M', 'F',
               'M', 'F')

Age_tmp = c(75, 53, 75, 59, 84,
            68, 38, 60, 48, 47,
            90, 70)


result.matrix = matrix(0, 12,3)
i = 0
dpi = 300
A.list = list()
nFeature.low = 200 
nFeature.hi = 10000
mt.hi = 20


for(sample_tmp in p_id){
    i = i + 1
    barcode.path <- paste0('Data/GSE155673_', sample_tmp, '_barcodes.tsv.gz')
    features.path <- paste0("Data/GSE155673_features.tsv.gz")
    matrix.path <- paste0('Data/GSE155673_', sample_tmp, '_matrix.mtx.gz')
    mat <- readMM(file = matrix.path)
    feature.names = read.delim(features.path, header = FALSE, stringsAsFactors = FALSE)
    barcode.names = read.delim(barcode.path, header = FALSE, stringsAsFactors = FALSE)
    colnames(mat) = barcode.names$V1
    rownames(mat) = feature.names$V2
    n = dim(mat)[2]
    print(sample_tmp)
    sample.tmp.seurat <- CreateSeuratObject(counts = mat, min.cells = 10)
    sample.tmp.seurat[['percent.mt']] <- PercentageFeatureSet(sample.tmp.seurat, pattern = "^MT-")

    meta = data.frame(Disease = rep(disease_type[i], n), Age = rep(Age_tmp[i], n), Gender = rep(Gender_tmp[i], n),
                     Sample = rep(Sample_id[i], n), Patient = rep(Patient_id[i], n), row.names = colnames(sample.tmp.seurat),
                     Experiment = rep('A', n))
    sample.tmp.seurat = AddMetaData(sample.tmp.seurat, metadata = meta)
   
   png(file = paste('QC/',sample_tmp,"_qc.png",sep=''), width = dpi*16, height = dpi*8, units = "px",res = dpi,type='cairo')
   print(VlnPlot(object = sample.tmp.seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
   dev.off()
  
   png(file=paste('QC/',sample_tmp,"_umi-mito.png",sep=''), width = dpi*6, height = dpi*5, units = "px",res = dpi,type='cairo')
   print(FeatureScatter(object = sample.tmp.seurat, feature1 = "nCount_RNA", feature2 = "percent.mt"))
   dev.off()
  
   png(file=paste('QC/',sample_tmp,"_umi-gene.png",sep=''), width = dpi*6, height = dpi*5, units = "px",res = dpi,type='cairo')
   print(FeatureScatter(object = sample.tmp.seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA"))
   dev.off()
  
  result.matrix[i,1] = sample_tmp
  result.matrix[i,2] = dim(sample.tmp.seurat)[2]
  sample.tmp.seurat <- subset(x = sample.tmp.seurat, subset =  nCount_RNA > 500  &  nFeature_RNA > nFeature.low & 
                        nFeature_RNA < nFeature.hi & percent.mt < mt.hi) #& percent.mt < 0.2
  result.matrix[i,3] = dim(sample.tmp.seurat)[2]
  sample.tmp.seurat <- NormalizeData(sample.tmp.seurat, verbose = FALSE)
  sample.tmp.seurat <- FindVariableFeatures(sample.tmp.seurat, selection.method = "vst", nfeatures = 4000,verbose = FALSE)
  A.list[Sample_id[i]] = sample.tmp.seurat
   png(file = paste('QC/filtered_',sample_tmp,"_qc.png",sep=''), width = dpi*16, height = dpi*8, units = "px",res = dpi,type='cairo')
   print(VlnPlot(object = sample.tmp.seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
   dev.off()
  
   png(file=paste('QC/filtered_',sample_tmp,"_umi-mito.png",sep=''), width = dpi*6, height = dpi*5, units = "px",res = dpi,type='cairo')
   print(FeatureScatter(object = sample.tmp.seurat, feature1 = "nCount_RNA", feature2 = "percent.mt"))
   dev.off()
  
   png(file=paste('QC/filtered_',sample_tmp,"_umi-gene.png",sep=''), width = dpi*6, height = dpi*5, units = "px",res = dpi,type='cairo')
   print(FeatureScatter(object = sample.tmp.seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA"))
   dev.off()
}

i = 1
plot(A.list[[i]]$nCount_RNA, (A.list[[i]]$nFeature_RNA))
abline(v = 50000, col = 'red')
abline(v = 125000, col = 'blue')
abline(h = 5500, col = 'red')

# These data is very noisy
# remove cell with nCount_RNA >125000
max_nCount = 125000
for(i in 1:12){
    A.list[[i]] = subset(x = A.list[[i]], subset =  nCount_RNA < max_nCount)
    print(i)
}
i = 1
plot(A.list[[i]]$nCount_RNA, (A.list[[i]]$nFeature_RNA))
abline(v = 50000, col = 'red')
abline(h = 5500, col = 'red')

max_nCount = 125000

for(i in 1:12){
    id = A.list[[i]]$nCount_RNA > 50000 & A.list[[i]]$nFeature_RNA < 5500
    A.list[[i]] = subset(x = A.list[[i]], cells = colnames(A.list[[i]])[!id]  )
    print(i)
}

i = 1
plot(A.list[[i]]$nCount_RNA, (A.list[[i]]$nFeature_RNA))
abline(v = 20000, col = 'red')
abline(h = 3000, col = 'red')

for(i in 1:12){
    id = A.list[[i]]$nCount_RNA > 20000 & A.list[[i]]$nFeature_RNA < 3000
    A.list[[i]] = subset(x = A.list[[i]], cells = colnames(A.list[[i]])[!id]  )
    print(i)
}

i = 1
plot(A.list[[i]]$nCount_RNA, (A.list[[i]]$nFeature_RNA))
abline(v = 7000, col = 'red')
abline(h = 1000, col = 'red')

for(i in 1:12){
    sample_tmp = p_id[i]
    id = A.list[[i]]$nCount_RNA > 7000 & A.list[[i]]$nFeature_RNA < 1000
    A.list[[i]] = subset(x = A.list[[i]], cells = colnames(A.list[[i]])[!id]  )
    print(i)
    sample.tmp.seurat = A.list[[i]]
    result.matrix[i,1] = Sample_id[i]
    result.matrix[i,2] = dim(sample.tmp.seurat)[2]
    result.matrix[i,3] = dim(sample.tmp.seurat)[2]
    sample.tmp.seurat <- NormalizeData(sample.tmp.seurat, verbose = FALSE)
    sample.tmp.seurat <- FindVariableFeatures(sample.tmp.seurat, selection.method = "vst", nfeatures = 4000,verbose = FALSE)
    A.list[Sample_id[i]] = sample.tmp.seurat
    png(file = paste('QC/filtered_',sample_tmp,"_qc.png",sep=''), width = dpi*16, height = dpi*8, units = "px",res = dpi,type='cairo')
    print(VlnPlot(object = sample.tmp.seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
    dev.off()
  
    png(file=paste('QC/filtered_',sample_tmp,"_umi-mito.png",sep=''), width = dpi*6, height = dpi*5, units = "px",res = dpi,type='cairo')
    print(FeatureScatter(object = sample.tmp.seurat, feature1 = "nCount_RNA", feature2 = "percent.mt"))
    dev.off()
  
    png(file=paste('QC/filtered_',sample_tmp,"_umi-gene.png",sep=''), width = dpi*6, height = dpi*5, units = "px",res = dpi,type='cairo')
    print(FeatureScatter(object = sample.tmp.seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA"))
    dev.off()
}

i = 1
plot(A.list[[i]]$nCount_RNA, (A.list[[i]]$nFeature_RNA))
abline(v = 7000, col = 'red')
abline(h = 1000, col = 'red')

result.dataframe = as.data.frame(result.matrix)
colnames(result.dataframe) = c('sample','before','filter')
write.table(result.dataframe, file='statistics_filter.txt',row.names = FALSE,quote = FALSE,sep='\t')


gene2[duplicated(gene2)]

gene2[is.na(match(gene2, rownames(A.list$A_C7)))]

P.features <- SelectIntegrationFeatures(object.list = A.list, nfeatures = 6000)
saveRDS(P.features, 'A_gene6000_11_25_2020.rds')

saveRDS(A.list, file = 'A_list.rds')


srt <- readRDS("unfiltered_srt_object.RDS")

# the unfiltered srt object has all the output of the below code included
# but should anyone decide to reprocess this is the script
srt = NormalizeData(srt, display.progress = FALSE)
srt = FindVariableFeatures(srt, do.plot = F, display.progress = FALSE)
srt = ScaleData(srt, display.progress = FALSE)
srt <- RunPCA(object = srt, verbose = FALSE)
peb <- ElbowPlot(srt, ndims = 50)
srt <- FindNeighbors(object = srt, dims = 1:25)
srt <- FindClusters(object = srt, resolution = 0.4)
srt <- RunUMAP(object = srt, dims = 1:25)

# plot the umap rep
UMAPPlot(srt, label = T) 

# look at some important qc metrics
p_mito <- ggplot(srt@meta.data, aes(x=seurat_clusters, y = mito)) + 
  geom_jitter(height = 0, width = 0.2) +
  geom_boxplot() +
  theme_cowplot()
p_ncount <- ggplot(srt@meta.data, aes(x=seurat_clusters, y = ncount)) + 
  geom_jitter(height = 0, width = 0.2) +
  geom_boxplot() +
  theme_cowplot() + ggtitle("Count RNA Per Cell")
p_nfeat <- ggplot(srt@meta.data, aes(x=seurat_clusters, y = nfeat)) + 
  geom_jitter(height = 0, width = 0.2) +
  geom_boxplot() +
  theme_cowplot() + ggtitle("Unique Feat Per Cell")


# given the expression and awful qc values of cluster 8
# we are going to get rid of it entirely
# the marker genes for the cluster are also suggestive of dead cells 
srt_clean <- subset(srt, idents = "8", invert = T)
# and also get rid of cells that have more than 25% mito RNA
srt_clean <- subset(srt_clean, subset = mito < 25)

# new plot
UMAPPlot(srt_clean, label = T) 

# remove the old object
rm(srt)

#################################
# Correlations

# sub cluster by monocyte related-clusters
srt_sub_clust <- c("3", "4", "5", "10", "11", "16",  "22")
Idents(srt_clean) <- "seurat_clusters"
sub_srt <- subset(srt_clean, idents = srt_sub_clust)

# pull out the gene data
gene_dat <- sub_srt@assays$RNA@data
dim(gene_dat)
# its awfully big
# so we better be careful correlating

# pull out two genes of interest
# S100A12 we found to be important in the Olink analysis
# IFI27 is a very important viral response gene 
sub_genedat <- gene_dat[which(rownames(gene_dat) %in% c("S100A12", "IFI27")),]

# transform it for easy correlations
gene_dat_t <- t(gene_dat)

# it would have been better to swap zeros for NAs here
# but this matrix is very big, and corSparse (from qlcMatrix) is very fast
corres <- corSparse(gene_dat_t, t(sub_genedat))

# configure for later use
corres <- as.data.frame(corres)
corres$gene <- rownames(gene_dat)

# dont want IFI27 correlations anymore
corres$V2 <- NULL

# keep only correlations that are above .3
# just given the number of points these are likely to be sig
# we will confirm that later
corres_sig <- corres[which(abs(corres$V1) > .3),]
colnames(corres_sig)[1] <- "corval"

# pick out the top contenders that correlate with S100A12 (EN-RAGE)
# this is manual but they are the top 5 most and least correlated genes
IFN_genes_cor <- c("HLA-DPA1", "HLA-DPB1", "RPS19", "HLA-DRA", "CD74", 
                    "S100A8", "S100A9", "VCAN", "CD14", "PLBD1")
# subset the cor list
corres_sig_best <- corres_sig[which(corres_sig$gene %in% IFN_genes_cor),]
corres_sig_best <- corres_sig_best[order(corres_sig_best$corval, decreasing = F),]

# factor for plotting
corres_sig_best$gene <- factor(corres_sig_best$gene, levels = corres_sig_best$gene)

# plot the results
p_cor <- ggplot(data=corres_sig_best, aes(x=gene, y=corval)) +
  geom_bar(stat="identity",  fill = "#478A5E") + coord_flip() + theme_cowplot() +
  ggtitle("Correlation EN-RAGE") + ylab("Pearson Correlation") + xlab("")
p_cor
