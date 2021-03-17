library(readr)
library(viridis)
library(gridExtra)
source( "seuratWrapper.R")
source('COVID_Function.R')
library(phateR) 
library(dplyr)   # define funciton %>%
library(hdf5r)
###pbmc integration and reduction
##integrate data with seurat v3
library(Seurat)
library(Matrix)
library(ggplot2)
library(MAST)
library(SingleR)
setwd('/mnt/home/xyy/COVID19_P')
covid.names <- paste0('c',c(1:7))
healthy.names <- paste0('h',c(1:6))
ob.list <- c(covid.names, healthy.names)


###integation
setwd('/mnt/research/xlab/covid.y/pbmc/addInteg')
nCoV.list = list()
result.matrix = matrix(0,length(ob.list),3)
i = 0
#sample = ob.list[1]
for(sample in ob.list){
  i = i + 1
  print(sample)
  sample.tmp <- paste0('/mnt/research/xlab/covid.y/pbmc/data/imputed/imputed_',sample,'.rds')
  nCoV.data.i <- readRDS(sample.tmp)
  #nCoV.data.i <- nCoV.data.i$estimate 
  
  sample.tmp.seurat <- CreateSeuratObject(counts = nCoV.data.i, min.cells = 10)
  sample.tmp.seurat[['percent.mt']] <- PercentageFeatureSet(sample.tmp.seurat, pattern = "^MT-")
  
 # dpi = 300
  #png(file=paste('ShowImpute/',sample,"_qc.png",sep=''), width = dpi*16, height = dpi*8, units = "px",res = dpi,type='cairo')
  #print(VlnPlot(object = sample.tmp.seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
  #dev.off()
  
  #png(file=paste('ShowImpute/',sample,"_umi-mito.png",sep=''), width = dpi*6, height = dpi*5, units = "px",res = dpi,type='cairo')
  #print(FeatureScatter(object = sample.tmp.seurat, feature1 = "nCount_RNA", feature2 = "percent.mt"))
  #dev.off()
  
 # png(file=paste('ShowImpute/',sample,"_umi-gene.png",sep=''), width = dpi*6, height = dpi*5, units = "px",res = dpi,type='cairo')
  #print(FeatureScatter(object = sample.tmp.seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA"))
  #dev.off()
  
  result.matrix[i,1] = sample
  result.matrix[i,2] = dim(sample.tmp.seurat)[2]
  sample.tmp.seurat <- subset(x = sample.tmp.seurat, subset =  nCount_RNA > 500 ) #& percent.mt < 0.2
  result.matrix[i,3] = dim(sample.tmp.seurat)[2]

  sample.tmp.seurat <- NormalizeData(sample.tmp.seurat, verbose = FALSE)
  sample.tmp.seurat <- FindVariableFeatures(sample.tmp.seurat, selection.method = "vst", nfeatures = 2000,verbose = FALSE)
  nCoV.list[sample] = sample.tmp.seurat
}


result.dataframe = as.data.frame(result.matrix)
colnames(result.dataframe) = c('sample','before','filter')
write.table(result.dataframe,file='ShowImpute/statistics_filter.txt',row.names = FALSE,quote = FALSE,sep='\t')

P.features <- SelectIntegrationFeatures(object.list = nCoV.list, nfeatures = 2000)
P.features = unique(c(gene, P.features))
P.anchors <- FindIntegrationAnchors(object.list = nCoV.list, normalization.method = "LogNormalize", anchor.features = P.features)

B.integrated <- IntegrateData(anchorset = B.anchors, normalization.method = "LogNormalize")

nCoV <- FindIntegrationAnchors(object.list = nCoV.list, dims = 1:50)
nCoV.integrated <- IntegrateData(anchorset = nCoV, dims = 1:50,features.to.integrate = rownames(nCoV))

####add  sample info
samples = read.delim2("meta.txt",header = TRUE, 
            stringsAsFactors = FALSE,check.names = FALSE, sep = "\t")
sample_info = as.data.frame(colnames(nCoV.integrated))
colnames(sample_info) = c('ID')
rownames(sample_info) = sample_info$ID
sample_info$sample = substr(sample_info[,], 1, 2)
sample_info = dplyr::left_join(sample_info, samples)
rownames(sample_info) = sample_info$ID
nCoV.integrated = AddMetaData(object = nCoV.integrated, metadata = sample_info)

save(nCoV.integrated,file = 'nCoV.integrated_09_10_2020.rda')

print('***integrate data with seurat v3***')

load('nCoV.integrated_09_10_2020.rda')

###first generate data and scale data in RNA assay
DefaultAssay(nCoV.integrated) <- "RNA"
nCoV.integrated[['percent.mt']] <- PercentageFeatureSet(nCoV.integrated, pattern = "^MT-")
nCoV.integrated <- NormalizeData(object = nCoV.integrated, normalization.method = "LogNormalize", scale.factor = 1e4)
nCoV.integrated <- FindVariableFeatures(object = nCoV.integrated, selection.method = "vst", nfeatures = 2000,verbose = FALSE)
nCoV.integrated <- ScaleData(nCoV.integrated, verbose = FALSE, vars.to.regress = c("nCount_RNA", "percent.mt"))
print('***first generate data and scale data in RNA assay***')
save(nCoV.integrated,file = 'nCoV.integrated_09_10_2020.rda')

##############################################
##change to integrated assay
##############################################
DefaultAssay(nCoV.integrated) <- "integrated"
# Run the standard workflow for visualization and clustering
nCoV.integrated <- ScaleData(nCoV.integrated, verbose = FALSE, vars.to.regress = c("nCount_RNA", "percent.mito"))
nCoV.integrated <- RunPCA(nCoV.integrated, verbose = FALSE,npcs = 100)
nCoV.integrated <- ProjectDim(object = nCoV.integrated)
png(file="cluster_ElbowPlot.png", width = dpi*10, height = dpi*6, units = "px",res = dpi,type='cairo')
ElbowPlot(object = nCoV.integrated,ndims = 100)
dev.off()
nCoV.integrated <- JackStraw(nCoV.integrated, num.replicate = 100,dims = 60)
nCoV.integrated <- ScoreJackStraw(nCoV.integrated, dims = 1:60)
png(file="cluster_JackStrawPlot.png", width = dpi*10, height = dpi*6, units = "px",res = dpi,type='cairo')
JackStrawPlot(nCoV.integrated, dims = 1:60)
dev.off()
save(nCoV.integrated,file = 'nCoV.integrated_09_10_2020.rda')
print('***Run the standard workflow for visualization and clustering***')

###cluster
nCoV.integrated <- FindNeighbors(object = nCoV.integrated, dims = 1:50)
nCoV.integrated <- FindClusters(object = nCoV.integrated, resolution = 1.2)### 
print('***cluster***')

nCoV.integrated <- RunUMAP(nCoV.integrated, reduction = "pca", dims = 1:50)

save(nCoV.integrated,file = 'nCoV.integrated_09_10_2020.rda')

png(file="umap.png", width = dpi*8, height = dpi*6, units = "px",res = dpi,type='cairo')
DimPlot(object = nCoV.integrated, reduction = 'umap',label = TRUE)
dev.off()

png(file="umap_disease_09_11_2020.png", width = dpi*8, height = dpi*6, units = "px",res = dpi,type='cairo')
 DimPlot(nCoV.integrated, reduction = "umap", group.by = "disease")
dev.off()

png(file="umap_sample_09_11_2020.png", width = dpi*8, height = dpi*6, units = "px",res = dpi,type='cairo')
DimPlot(nCoV.integrated, label = TRUE, split.by = "sample", ncol = 4, reduction = 'umap')
dev.off()

DefaultAssay(nCoV.integrated) <- "RNA"
print('***tsne and umap***')

# find markers for every cluster compared to all remaining cells, report only the positive ones
load('nCoV.integrated_09_10_2020.rda')

nCoV.integrated@misc$markers <- FindAllMarkers(object = nCoV.integrated, assay = 'RNA',only.pos = TRUE, test.use = 'MAST')
write.table(nCoV.integrated@misc$markers,file='marker_MAST.txt',row.names = FALSE,quote = FALSE,sep = '\t')

dpi = 300
png(file="feature_vln_09_12_2020.png", width = dpi*24, height = dpi*5, units = "px",res = dpi,type='cairo')
VlnPlot(object = nCoV.integrated, features = c("nFeature_RNA", "nCount_RNA"))
dev.off()

save(nCoV.integrated,file = 'nCoV.integrated_09_10_2020.rda')

hc.markers = read.delim2("marker_MAST.txt",header = TRUE, stringsAsFactors = FALSE,check.names = FALSE, sep = "\t")
hc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC) -> top10
tt1 = DoHeatmap(object = subset(nCoV.integrated, downsample = 500), features = top10$gene) + NoLegend()
ggplot2::ggsave(file="marker_heatmap_MAST.pdf",plot = tt1,device = 'pdf',width = 20, height = 16, units = "in",dpi = dpi,limitsize = FALSE)
print('***find markers***')

top_gene = hc.markers %>% group_by(cluster) %>% top_n(n=100, wt=avg_logFC)
write.csv(cbind(top_gene$gene, top_gene$cluster), file = 'genes_list_COVID_PBMC_saver_09_11_2020.csv')

res =Get_mean(Epi)
write.csv(res, 'COVID_PBMC_saver_mean_sd_09_11_2020.csv')

###SingleR annotation
ref <- HumanPrimaryCellAtlasData()
dat <- GetAssayData(nCoV.integrated[['RNA']], slot = 'data')
nCoV.integrated_clusters = nCoV.integrated$seurat_clusters
singler.pred <- SingleR(test = dat, ref = ref, labels = ref$label.main,
                        method = "cluster", clusters = nCoV.integrated_clusters, de.method="wilcox")
save(singler_cluster,file = 'singler.pred.rda')
write.csv(singler.pred,'singler_pred_labels.csv')

###plot
dpi = 300
png(file=paste('singler.pred.orig.ident.png',sep=''), width = dpi*16, height = dpi*8, units = "px",res = dpi,type='cairo')
plotScoreHeatmap(singler.pred, clusters = nCoV.integrated$seurat_clusters)
dev.pff()

png(file=paste('singler.pred.seurat_clusters.png',sep=''), width = dpi*16, height = dpi*8, units = "px",res = dpi,type='cairo')
plotScoreHeatmap(singler.pred, clusters = nCoV.integrated@meta.data$seurat_clusters)
dev.pff()

###???
singler.results <- merge(data.frame(cell = rownames(singler.pred), singler = singler.pred$labels), 
                         data.frame(cell = rownames(nCoV.integrated@meta.data), 
                                    cluster = nCoV.integrated@meta.data$seurat_clusters), 
                         by = "cell", 
                         all.y = FALSE)
singler.results$cell <- NULL
singler.results$count <- 1
singler.results <- aggregate(count ~ ., singler.results, FUN = sum)
singler.final <- singler.results %>% group_by(cluster) %>% top_n(n = 1, wt = count)

###use marker???
setwd('/mnt/research/xlab/covid.y/pbmc/addInteg/DotPlot')
DefaultAssay(nCoV.integrated) = 'RNA'

############ 
# using original count data
load("/mnt/research/xlab/covid.y/pbmc/original/covid_combined.rda")
gene_id = match(rownames(nCoV.integrated), rownames(covid_combined))
id = match(colnames(nCoV.integrated), colnames(covid_combined))
A = GetAssayData(covid_combined[['RNA']], slot = 'counts')
A1 = GetAssayData(covid_combined[['RNA']], slot = 'data')

B = subset(nCoV.integrated, cells = colnames(nCoV.integrated)[!is.na(id)])
id = match(colnames(B), colnames(covid_combined))
A = A[gene_id, id]
A1 = A1[gene_id, id]

B@assays$RNA@counts = A
B@assays$RNA@data = A1

DefaultAssay(B) = 'RNA'
png(file="Covid_P_dotplot_09_12_2020.png", width = dpi*16, height = dpi*8, units = "px",res = dpi,type='cairo')
DotPlot(B, features =  gene) +  RotatedAxis()
dev.off()

png(file="p5.png", width = dpi*16, height = dpi*8, units = "px",res = dpi,type='cairo')
DotPlot(Seurat:::subset.Seurat(nCoV.integrated, 
                               idents = unlist(lapply(singler.final[singler.final$singler=="Monocyte","cluster"],
                                                      as.character))), 
        features = c("CD14", "LYZ", "FCGR3A", "IL1B", "IDO1", "FCER1A", "FLT3", "IL3RA", "NRP1"))
dev.off()
png(file="p6.png", width = dpi*16, height = dpi*8, units = "px",res = dpi,type='cairo')
DotPlot(Seurat:::subset.Seurat(nCoV.integrated, 
                               idents = unlist(lapply(singler.final[singler.final$singler=="B_cell","cluster"],
                                                      as.character))), 
        features = c("MME", "CD22", "FCER2", "CD38", "CD44", "CD27", "SDC1"))

dev.off()
#png(file="p7.png", width = dpi*16, height = dpi*8, units = "px",res = dpi,type='cairo')
#DotPlot(Seurat:::subset.Seurat(nCoV.integrated, idents = c(covid_myeloid.idents)), features = unique(c("CD14", "LYZ", "FCGR3B", "FCGR3A", "CLC", "ELANE", "LTF", "MPO", "CTSG", "IDO1", "FCER1A", "FLT3", "IL3RA", "NRP1", "MME", "CD22", "FCER2", "CD44", "CD27", "SDC1", "CD4", "CD8A", "ITGAL", "SELL", "GZMB", "CD3E", "CD3G", "CD3D")), group.by = "cell.type.fine") + ggpubr::rotate_x_text()
#dev.off()
ggsave("p.pdf", path = "/mnt/research/xlab/covid.y/pbmc/addInteg/Dotplot/", height = 5, width = 10)

