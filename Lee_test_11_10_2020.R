library(Matrix)
library(readr)
library(viridis)
library(gridExtra)
source( "/mnt/ufs18/home-133/xyy/COVID19_B/seuratWrapper.R")
source('../COVID_Function.R')
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
library(pscl)
library(boot)
#The suffix of each barcode sequence represents patient information.
#'-1' - 'Sample1' - 'nCoV 1 scRNA-seq'
#'-2' - 'Sample2' - 'nCoV 2 scRNA-seq'
#'-3' - 'Sample3' - 'Flu 1 scRNA-seq'
#'-4' - 'Sample4' - 'Flu 2 scRNA-seq'
#'-5' - 'Sample5' - 'Normal 1 scRNA-seq'
#'-6' - 'Sample6' - 'Flu 3 scRNA-seq'
#'-7' - 'Sample7' - 'Flu 4 scRNA-seq'
#'-8' - 'Sample8' - 'Flu 5 scRNA-seq'
#'-9' - 'Sample9' - 'nCoV 3 scRNA-seq'
#'-10' - 'Sample10' - 'nCoV 4 scRNA-seq'
#'-11' - 'Sample11' - 'nCoV 5 scRNA-seq'
#'-12' - 'Sample12' - 'nCoV 6 scRNA-seq'
#'-13' - 'Sample13' - 'Normal 2 scRNA-seq'
#'-14' - 'Sample14' - 'Normal 3 scRNA-seq'
#'-15' - 'Sample15' - 'nCoV 7 scRNA-seq'
#'-16' - 'Sample16' - 'nCoV 8 scRNA-seq'
#'-17' - 'Sample17' - 'nCoV 9 scRNA-seq'
#'-18' - 'Sample18' - 'nCoV 10 scRNA-seq'
#'-19' - 'Sample19' - 'Normal 4 scRNA-seq'
#'-20' - 'Sample20' - 'nCoV 11 scRNA-seq'

matrix_dir = "/mnt/ufs18/home-133/xyy/COVID19_P/Lee/"
barcode.path <- paste0(matrix_dir, "GSE149689_barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "GSE149689_features.tsv.gz")
matrix.path <- paste0(matrix_dir, "GSE149689_matrix.mtx.gz")
mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path, header = FALSE, stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path, header = FALSE, stringsAsFactors = FALSE)
colnames(mat) = barcode.names$V1
rownames(mat) = feature.names$V2

temp =  substr(colnames(mat), 17, 19)
Disease = Gender = Patient = Sample = Batch = temp
Age = 1:length(temp)
temp_id = paste0('-', 1:20)
disease_type = c('Severe', 'Moderate', 'Flu', 'Flu', 'Healthy',
                 'Flu',    'Flu',  'Flu', 'Severe', 'Severe',
                 'Moderate', 'Moderate', 'Healthy', 'Healthy', 'Severe',
                 'Severe', 'Severe', 'Moderate', 'Healthy', 'Moderate') 

Patient_id = c('L_C1', 'L_C2', 'L_F1', 'L_F2', 'L_N1',
               'L_F3', 'L_F4', 'L_F5', 'L_C3', 'L_C3',
               'L_C4', 'L_C5', 'L_N2', 'L_N3', 'L_C6',
               'L_C6', 'L_C7', 'L_C7', 'L_N4', 'L_C8')

Sample_id = c('L_C1', 'L_C2', 'L_F1', 'L_F2', 'L_N1',
              'L_F3', 'L_F4', 'L_F5', 'L_C3_1', 'L_C3_2',
            'L_C4', 'L_C5', 'L_N2', 'L_N3', 'L_C6_1',
                'L_C6_2', 'L_C7_1', 'L_C7_2', 'L_N4', 'L_C8')

Gender_tmp = c('M', 'F', 'M', 'F', 'F', 
        'M', 'F', 'M', 'F', 'F',
        'M', 'F', 'F', 'F', 'F',
        'F', 'M', 'M', 'M', 'F')
        
Age_tmp = c(63, 82, 68, 75, 63,
        70, 56, 78, 67, 67,
        46, 38, 54, 67, 73,
        73, 61, 61, 63, 62)
Batch_temp = c(1, 1, 1, 1, 2,
                2, 2, 2, 3, 3,
               3, 3, 3, 3, 4,
               4, 4, 4, 4, 4)

P_id = c('L_C1', 'L_C2',  'L_N1',  'L_C3_1', 'L_C3_2',
          'L_C4', 'L_C5', 'L_N2', 'L_N3', 'L_C6_1',
          'L_C6_2', 'L_C7_1', 'L_C7_2', 'L_N4', 'L_C8')

for(i in 1:20){
    id = (temp == temp_id[i])
    Disease[id]  = disease_type[i]
    Age[id] = Age_tmp[i]
    Patient[id] = Patient_id[i]
    Gender[id] = Gender_tmp[i]
    Sample[id] = Sample_id[i]
    Batch[id] = Batch_temp[i]
}

Lee.seurat <- CreateSeuratObject(counts = mat, min.cells = 10)
Lee.seurat[['percent.mt']] <- PercentageFeatureSet(Lee.seurat, pattern = "^MT-")

meta = data.frame(Disease = Disease, Age = Age, Gender = Gender, Sample = Sample,
                 Patient = Patient, Batch = Batch, row.names = colnames(Lee.seurat),Experiment = rep('Lee', dim(mat)[2] ))

Lee.seurat = AddMetaData(Lee.seurat, metadata = meta)
temp_id = Lee.seurat@meta.data$Disease == 'Flu'

Lee = subset(Lee.seurat, cells = colnames(Lee.seurat)[!temp_id])
rm(Lee.seurat)
Lee.list = list()
result.matrix = matrix(0,15,3)
i = 0
#sample = ob.list[1]
nFeature.low = 200 
nFeature.hi=10000
mt.hi=20

for(sample in P_id){
  i = i + 1
  print(sample)
  temp_id = Lee@meta.data$Sample == sample
  sample.tmp.seurat <- subset(Lee, cells = colnames(Lee)[temp_id]) 
  dpi = 300
  png(file=paste0('QC/', sample, "_qc.png"), width = dpi*16, height = dpi*8, units = "px",res = dpi,type='cairo')
  print(VlnPlot(object = sample.tmp.seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
  dev.off()
   print(1)
  png(file=paste0('QC/', sample,"_umi-mito.png"), width = dpi*6, height = dpi*5, units = "px",res = dpi,type='cairo')
  print(FeatureScatter(object = sample.tmp.seurat, feature1 = "nCount_RNA", feature2 = "percent.mt"))
  dev.off()
  print(2)
  png(file=paste0('QC/', sample,"_umi-gene.png"), width = dpi*6, height = dpi*5, units = "px",res = dpi,type='cairo')
  print(FeatureScatter(object = sample.tmp.seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA"))
  dev.off()
  print(3)
  result.matrix[i, 1] = sample
  result.matrix[i, 2] = dim(sample.tmp.seurat)[2]
  sample.tmp.seurat <- subset(x = sample.tmp.seurat, subset =  nCount_RNA > 500  &  nFeature_RNA > nFeature.low & 
                        nFeature_RNA < nFeature.hi & percent.mt < mt.hi) #& percent.mt < 0.2
  result.matrix[i,3] = dim(sample.tmp.seurat)[2]

  sample.tmp.seurat <- NormalizeData(sample.tmp.seurat, verbose = FALSE)
  sample.tmp.seurat <- FindVariableFeatures(sample.tmp.seurat, selection.method = "vst", nfeatures = 4000,verbose = FALSE)
  Lee.list[sample] = sample.tmp.seurat
}

result.dataframe = as.data.frame(result.matrix)
colnames(result.dataframe) = c('sample','before','filter')
write.table(result.dataframe,file='statistics_filter.txt',row.names = FALSE,quote = FALSE,sep='\t')
rm(Lee)

gene2[duplicated(gene2)]

gene2[is.na(match(gene2, rownames(Lee.list$L_C8)))]

P.features <- SelectIntegrationFeatures(object.list = Lee.list, nfeatures = 6000)
saveRDS(P.features, 'Lee_gene6000_11_25_2020.rds')
saveRDS(Lee.list, file = 'L_list.rds')


match(gene, rownames(Lee.list$L_C8))
P.features = unique(c(gene, P.features))
P.features = intersect(P.features, rownames(Lee.list$C8))

P.anchors <- FindIntegrationAnchors(object.list = Lee.list, normalization.method = "LogNormalize", anchor.features = P.features)
Lee.integrated <- IntegrateData(anchorset = P.anchors, normalization.method = "LogNormalize")

save(Lee.integrated,file = 'nCoV.integrated_11_24_2020.rda')

DefaultAssay(Lee.integrated) <- 'integrated'
Lee.integrated[['percent.mt']] <- PercentageFeatureSet(Lee.integrated, pattern = "^MT-")
Lee.integrated <- NormalizeData(object = Lee.integrated, normalization.method = "LogNormalize", scale.factor = 1e4)
Lee.integrated <- ScaleData(Lee.integrated, verbose = FALSE, vars.to.regress = c("nCount_RNA", "percent.mt"))

save(Lee.integrated,file = 'nCoV.integrated_09_10_2020.rda')

