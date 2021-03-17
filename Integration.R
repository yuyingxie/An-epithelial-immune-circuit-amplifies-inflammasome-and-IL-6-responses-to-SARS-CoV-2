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
source('COVID_Function.R')

library(phateR) 
library(dplyr)   # define funciton %>%
library(hdf5r)

T_SCT.features <- SelectIntegrationFeatures(object.list = T.list, nfeatures = 2000)
T_SCT.features = unique(c(gene, T_SCT.features))
for(i in 1:length(T.list)){
    T_SCT.features = intersect(T_SCT.features, rownames(T.list[[i]]))
}

T.list <- PrepSCTIntegration(object.list = T.list, anchor.features = T_SCT.features)

reference_dataset = 22

T.anchors <- FindIntegrationAnchors(object.list = T.list, normalization.method = "SCT", 
    anchor.features = T_SCT.features, reference = reference_dataset)

T.integrated <- IntegrateData(anchorset = T.anchors, normalization.method = "SCT")

T.integrated <- RunPCA(object = T.integrated, verbose = FALSE, npcs = 100)
T.integrated <- RunUMAP(object = T.integrated, dims = 1:30)

###cluster
T.integrated <- FindNeighbors(object = T.integrated, dims = 1:30)
T.integrated <- FindClusters(object = T.integrated, resolution = 0.8)### 

saveRDS(T.integrated,file = 'T_SCT_integrated_11_26_2020.rds')

old_name= c('L_N1', 'L_N2', 'L_N3', 'L_N4', 'A_N1',
            'A_N2', 'A_N3', 'A_N4', 'A_N5', 'W_N1',
            'W_N2', 'W_N3', 'W_N4', 'W_N5', 'W_N6',
            'L_C2', 'L_C4', 'L_C5', 'L_C7_2', 'L_C8',
            'A_C2', 'A_C3', 'A_C7', 'W_C7', 'L_C1',
           'L_C3_1', 'L_C3_2', 'L_C6_1', 'L_C6_2', 'L_C7_1',
            'A_C1', 'A_C4', 'A_C5', 'A_C6', 'W_C1_1',
            'W_C1_2', 'W_C2', 'W_C3', 'W_C4', 'W_C6' ) 

new_name= c('L_H1', 'L_H2', 'L_H3', 'L_H4', 'A_H1',
            'A_H2', 'A_H3', 'A_H4', 'A_H5', 'W_H1',
            'W_H2', 'W_H3', 'W_H4', 'W_H5', 'W_H6',
            'L_M2', 'L_M4', 'L_M5', 'L_M7_2', 'L_M8',
            'A_M2', 'A_M3', 'A_M7', 'W_M7', 'L_S1',
            'L_S3_1', 'L_S3_2', 'L_S6_1', 'L_S6_2', 'L_S7_1',
            'A_S1', 'A_S4', 'A_S5', 'A_S6', 'W_S1_1',
            'W_S1_2', 'W_S2', 'W_S3', 'W_S4', 'W_S6' ) 
for(i in 1:40){
    id = T.integrated$Sample == old_name[i]
    T.integrated$Sample[id] = new_name[i]
}

DefaultAssay(T.integrated) = 'integrated'

pdf(file="T_SCT_Exp_umap.pdf", width = 8, height = 6)
DimPlot(object = T.integrated, reduction = 'umap',group.by = 'Experiment')
dev.off()

pdf(file="T_SCT_Exp_umap_split.pdf", width = 8, height = 6)
DimPlot(T.integrated, label = TRUE, split.by = "Experiment", ncol = 4, reduction = 'umap')
dev.off()

pdf(file="T_SCT_umap.pdf", width = 8, height = 6)
DimPlot(object = T.integrated, reduction = 'umap',label = TRUE)
dev.off()

pdf(file="umap_disease.pdf", width = 8, height = 6)
 DimPlot(T.integrated, reduction = "umap", group.by = "Disease")
dev.off()

pdf(file="umap_disease.pdf", width = 12, height = 4)
 DimPlot(T.integrated, reduction = "umap", split.by = "Disease")
dev.off()

pdf(file="umap_sample.pdf", width = 18, height = 24)
DimPlot(T.integrated, label = TRUE, split.by = "Patient", ncol = 8, reduction = 'umap')
dev.off()

pdf(file="umap_sample_split.pdf", width = 18, height = 24)
DimPlot(T.integrated, label = TRUE, split.by = "Sample", ncol = 8, reduction = 'umap')
dev.off()

T.integrated = ScaleData(T.integrated, vars.to.regress=  c("nCount_RNA", "percent.mt"))

#########################################################################
##########  Annotation 
##################################################################
library(SingleR)
DefaultAssay(T.integrated) = 'RNA'
#ImmGe = ImmGenData() # This is for mouse 
ref <- HumanPrimaryCellAtlasData()
logdata <- GetAssayData(T.integrated[['RNA']], slot = 'data')
T_clusters = T.integrated$seurat_clusters

singler_cluster <- SingleR(test = logdata, ref = ref, labels = ref$label.main,  
                           method = "cluster", clusters = T_clusters, de.method="wilcox")
###################################
## Find Marker
##################################
DefaultAssay(T.integrated) <- "RNA"

T.integrated@misc$markers <- FindAllMarkers(object = T.integrated, assay = 'RNA', only.pos = TRUE, test.use = 'wilcox')
write.table(T.integrated@misc$markers,file='marker_wilcox.txt',row.names = FALSE,quote = FALSE,sep = '\t')

dpi = 300
png(file="feature_vln_09_12_2020.png", width = dpi*24, height = dpi*5, units = "px",res = dpi,type='cairo')
VlnPlot(object = T.integrated, features = c("nFeature_RNA", "nCount_RNA"))
dev.off()

hc.markers = read.delim2("marker_wilcox.txt",header = TRUE, stringsAsFactors = FALSE,check.names = FALSE, sep = "\t")
hc.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC) -> top50
tt1 = DoHeatmap(object = subset(T.integrated, downsample = 500), features = top50$gene) + NoLegend()
ggplot2::ggsave(file="marker_heatmap_MAST.pdf", plot = tt1, device = 'pdf',width = 20, height = 16, units = "in",dpi = dpi,limitsize = FALSE)
print('***find markers***')

top_gene = hc.markers %>% group_by(cluster) %>% top_n(n=100, wt=avg_logFC)
write.csv(cbind(top_gene$gene, top_gene$cluster), file = 'genes_list_COVID_PBMC_saver_09_11_2020.csv')

res =Get_mean(T.integrated)

##################################
##Subset
######################################
filter_gene_list = readRDS( 'filter_gene_list_11_25_2020.rds')
T_SCT.features = readRDS('filter_SCT_gene_list_11_25_2020.rds')

gene_sub = unique(filter_gene_list, T_SCT.features)

for(i in 1:length(T.list)){
    gene_sub = intersect(gene_sub, rownames(T.list[[i]]))
} 


T_sub.list = list
for(i in 1:length(T.list)){

}

