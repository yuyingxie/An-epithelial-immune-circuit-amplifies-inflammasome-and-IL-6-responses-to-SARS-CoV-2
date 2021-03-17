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

source('additionalFunctions_repo.R')
library(phateR) 
library(dplyr)   # define funciton %>%
library(hdf5r)

source('additionalFunctions_repo.R')

#   Sample GEO ID        Sample name    Age   Gender   Disease
# 1  C1_1  GSM4557327 	covid_555_1     65    M        Severe
# 2  C1_2  GSM4557328 	covid_555_2     65    M        Severe
# 3  C2    GSM4557329 	covid_556       45    M        Severe
# 4  C3    GSM4557330 	covid_557       35    M        Severe
# 5  C4    GSM4557331 	covid_558       35    M        Severe

# 6  C5    GSM4557332 	covid_559       55    M        Severe
# 7  C6    560          covid_560       80    M        Severe
# 8  C7    GSM4557333 	covid_561       25    M        Moderate
# 9  N1    GSM4557334 	HIP002          49    F        Healthy
# 10 N2    GSM4557335 	HIP015          49    M        Healthy

# 11 N3    GSM4557336 	HIP023          36    F        Healthy
# 12 N4    GSCM4557337 	HIP043          49    M        Healthy
# 13 N5    GSM4557338 	HIP044          48    M        Healthy
# 14 N6    GSM4557339 	HIP045          37    M        Healthy


path = "Data/"
cm.list = paste0(path, list.files(pattern = "*.matrices.rds", path = path))
cm.files <- lapply(cm.list, readRDS)
names(cm.files) <- sub(path,"",
		sub("\\_cell.counts.matrices.rds", "", cm.list))
cm.pp <- mapply(EpicPreHS, cm.files, orig.ident = names(cm.files), SIMPLIFY = F)
rm(cm.files)
sp1_1 = cm.pp$GSM4557327_555_1$emat
sp1_2 = cm.pp$GSM4557328_555_2$emat
sp2 = cm.pp$GSM4557329_556$emat
sp3 = cm.pp$GSM4557330_557$emat
sp4 = cm.pp$GSM4557331_558$emat
sp5 = cm.pp$GSM4557332_559$emat
sp6 = cm.pp$`560`$emat
sp7 = cm.pp$GSM4557333_561$emat
hp1 = cm.pp$GSM4557334_HIP002$emat
hp2 = cm.pp$GSM4557335_HIP015$emat
hp3 = cm.pp$GSM4557336_HIP023$emat
hp4 = cm.pp$GSM4557337_HIP043$emat
hp5 = cm.pp$GSM4557338_HIP044$emat
hp6 = cm.pp$GSM4557339_HIP045$emat


disease_type = c('Severe', 'Severe', 'Severe', 'Severe', 'Severe',
                 'Severe', 'Severe', 'Moderate', 'Healthy', 'Healthy',
              'Healthy', 'Healthy', 'Healthy', 'Healthy')  

Patient_id = c('W_C1', 'W_C1', 'W_C2', 'W_C3', 'W_C4',
               'W_C5', 'W_C6', 'W_C7', 'W_N1', 'W_N2',
               'W_N3', 'W_N4', 'W_N5', 'W_N6')

Sample_id = c('W_C1_1', 'W_C1_2', 'W_C2', 'W_C3', 'W_C4',
              'W_C5',   'W_C6',   'W_C7', 'W_N1', 'W_N2',
              'W_N3',   'W_N4', 'W_N5', 'W_N6')

Gender_tmp = c('M', 'M', 'M', 'M', 'M', 
               'M', 'M', 'M', 'F', 'M',
               'F', 'M', 'M', 'M')

Age_tmp = c(65, 65, 45, 35, 35,
            55, 80, 25, 49, 49,
            36, 49, 48, 37)

p_id = c(2:7, 1, 8:14)

#sample = ob.list[1]
result.matrix = matrix(0, 14,3)
i = 0
dpi = 300
W.list = list()
nFeature.low = 200 
nFeature.hi = 10000
mt.hi = 20
p_id 
i = 0
for(sample_tmp in p_id){
    i = i + 1
    mat = cm.pp[[sample_tmp]]$emat
    n = dim(mat)[2]
    sample.tmp.seurat <- CreateSeuratObject(counts = mat, min.cells = 10)
    sample.tmp.seurat[['percent.mt']] <- PercentageFeatureSet(sample.tmp.seurat, pattern = "^MT-")

    meta = data.frame(Disease = rep(disease_type[i], n), Age = rep(Age_tmp[i], n), Gender = rep(Gender_tmp[i], n),
                     Sample = rep(Sample_id[i], n), Patient = rep(Patient_id[i], n), row.names = colnames(sample.tmp.seurat),
                     Experiment = rep('W', n))
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
                        nFeature_RNA < nFeature.hi & percent.mt < mt.hi) 
  result.matrix[i,3] = dim(sample.tmp.seurat)[2]
  sample.tmp.seurat <- NormalizeData(sample.tmp.seurat, verbose = FALSE)
  sample.tmp.seurat <- FindVariableFeatures(sample.tmp.seurat, selection.method = "vst", nfeatures = 4000,verbose = FALSE)
  W.list[Sample_id[i]] = sample.tmp.seurat
  print(Sample_id[i])
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
plot(W.list[[i]]$nCount_RNA, (W.list[[i]]$nFeature_RNA))
abline(v = 2000, col = 'red')
abline(h = 500, col = 'red')

for(i in 1:14){
    id = W.list[[i]]$nCount_RNA > 2000 & W.list[[i]]$nFeature_RNA < 500
    W.list[[i]] = subset(x = W.list[[i]], cells = colnames(W.list[[i]])[!id]  )
    print(i)
}

i = 1
plot(W.list[[i]]$nCount_RNA, (W.list[[i]]$nFeature_RNA))
abline(v = 2000, col = 'red')
abline(h = 500, col = 'red')

for(i in 1:14){
    sample_tmp = Sample_id[i]
    print(i)
    sample.tmp.seurat = W.list[[i]]
    result.matrix[i,1] = sample_tmp
    result.matrix[i,2] = dim(sample.tmp.seurat)[2]
    result.matrix[i,3] = dim(sample.tmp.seurat)[2]
    sample.tmp.seurat <- NormalizeData(sample.tmp.seurat, verbose = FALSE)
    sample.tmp.seurat <- FindVariableFeatures(sample.tmp.seurat, selection.method = "vst", nfeatures = 4000,verbose = FALSE)
    W.list[sample_tmp] = sample.tmp.seurat
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

result.dataframe = as.data.frame(result.matrix)
colnames(result.dataframe) = c('sample','before','filter')
write.table(result.dataframe, file='statistics_filter.txt',row.names = FALSE,quote = FALSE,sep='\t')

gene2[duplicated(gene2)]

gene2[is.na(match(gene2, rownames(W.list$W_C2)))]

P.features <- SelectIntegrationFeatures(object.list = W.list, nfeatures = 6000)
saveRDS(P.features, 'W_gene6000_11_25_2020.rds')
saveRDS(W.list, file = 'W_list.rds')





for(i in 1:14){
    sample.tmp.seurat = W.list[[i]]
    sample_tmp = Sample_id[i]
    print(sample_tmp)
    sample.tmp.seurat <- SCTransform(sample.tmp.seurat, vars.to.regress = "percent.mt", verbose = FALSE)
    W.list[sample_tmp] = sample.tmp.seurat
}
