# TODO: Add comment
# 
# Author: xyy
###############################################################################
library(Matrix)
library(Matrix.utils)
library(plyr)
library(dplyr)
library(Seurat)
library(sctransform)
library(igraph)
library(factoextra)
library(ComplexHeatmap)
library(circlize)
require(Hmisc)
require(openxlsx)
require(ggplot2)
library(ggpubr)
require(cowplot)
library(data.table)
library(RColorBrewer)
library(rowr)
library(SingleR)
library(scater)
library(pheatmap)
library(nichenetr)
library(tidyverse)

library(SAVER)
library(readr)
library(viridis)
library(gridExtra)
setwd('/mnt/ufs18/home-133/xyy/COVID19_P')
source("seuratWrapper.R")
source('additionalFunctions_repo.R')
library(phateR) 
library(dplyr)   # define funciton %>%
library(hdf5r)
#setwd("/Users/xyy/Dropbox (Personal)/single_cell/result")


setwd('C:/Users/xyy/Dropbox (Personal)/Covid19 data/GSE150728')
source('additionalFunctions_repo.R')

# GSM4557327 	PBMCs from COVID-19 sample covid_555_1     sp1_1
# GSM4557328 	PBMCs from COVID-19 sample covid_555_2     sp1_2
# GSM4557329 	PBMCs from COVID-19 sample covid_556       sp2
# GSM4557330 	PBMCs from COVID-19 sample covid_557       sp3
# GSM4557331 	PBMCs from COVID-19 sample covid_558       sp4
# GSM4557332 	PBMCs from COVID-19 sample covid_559       sp5
# 560           PBMCs from COVID-19 sample covid_560       sp6
# GSM4557333 	PBMCs from COVID-19 sample covid_561       sp7
# GSM4557334 	PBMCs from healthy control HIP002          hp1 
# GSM4557335 	PBMCs from healthy control HIP015          hp2
# GSM4557336 	PBMCs from healthy control HIP023          hp3
# GSCM4557337 	PBMCs from healthy control HIP043          hp4
# GSM4557338 	PBMCs from healthy control HIP044          hp5
# GSM4557339 	PBMCs from healthy control HIP045          hp6


path = "Raw/"
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

rm(cm.pp)
#Merge
covid_combined.emat <- mergeCM(cm.pp, type = "emat")
#covid_combined.nmat <- mergeCM(cm.pp, type = "nmat")

# filter out CD45 negative and  KRT18 positive
sp1_1 <-CreateSeuratObject(sp1_1, project = "sp1_1", min.cells = 3, min.features = 200)
sp1_2 <-CreateSeuratObject(sp1_2, project = "sp1_2", min.cells = 3, min.features = 200)
sp2 <-CreateSeuratObject(sp2, project = "sp2", min.cells = 3, min.features = 200)
sp3 <-CreateSeuratObject(sp3, project = "sp3", min.cells = 3, min.features = 200)
sp4 <-CreateSeuratObject(sp4, project = "sp4", min.cells = 3, min.features = 200)
sp5 <-CreateSeuratObject(sp5, project = "sp5", min.cells = 3, min.features = 200)
sp6 <-CreateSeuratObject(sp6, project = "sp6", min.cells = 3, min.features = 200)
sp7 <-CreateSeuratObject(sp7, project = "sp7", min.cells = 3, min.features = 200)

hp1 <-CreateSeuratObject(hp1, project = "hp1", min.cells = 3, min.features = 200)
hp2 <-CreateSeuratObject(hp2, project = "hp2", min.cells = 3, min.features = 200)
hp3 <-CreateSeuratObject(hp3, project = "hp3", min.cells = 3, min.features = 200)
hp4 <-CreateSeuratObject(hp4, project = "hp4", min.cells = 3, min.features = 200)
hp5 <-CreateSeuratObject(hp5, project = "hp5", min.cells = 3, min.features = 200)
hp6 <-CreateSeuratObject(hp6, project = "hp6", min.cells = 3, min.features = 200)

#####QC
sp1_1[["percent.mt"]] <- PercentageFeatureSet(sp1_1, pattern = "^MT-") 
sp1_2[["percent.mt"]] <- PercentageFeatureSet(sp1_2, pattern = "^MT-") 
sp2[["percent.mt"]] <- PercentageFeatureSet(sp2, pattern = "^MT-") 
sp3[["percent.mt"]] <- PercentageFeatureSet(sp3, pattern = "^MT-") 
sp4[["percent.mt"]] <- PercentageFeatureSet(sp4, pattern = "^MT-") 
sp5[["percent.mt"]] <- PercentageFeatureSet(sp5, pattern = "^MT-") 
sp6[["percent.mt"]] <- PercentageFeatureSet(sp6, pattern = "^MT-") 
sp7[["percent.mt"]] <- PercentageFeatureSet(sp7, pattern = "^MT-") 
hp1[["percent.mt"]] <- PercentageFeatureSet(hp1, pattern = "^MT-") 
hp2[["percent.mt"]] <- PercentageFeatureSet(hp2, pattern = "^MT-") 
hp3[["percent.mt"]] <- PercentageFeatureSet(hp3, pattern = "^MT-") 
hp4[["percent.mt"]] <- PercentageFeatureSet(hp4, pattern = "^MT-") 
hp5[["percent.mt"]] <- PercentageFeatureSet(hp5, pattern = "^MT-") 
hp6[["percent.mt"]] <- PercentageFeatureSet(hp6, pattern = "^MT-") 

VlnPlot(sp1_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(sp1_2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(sp2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(sp3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(sp4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(sp5, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(sp6, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(sp7, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(hp1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(hp2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(hp3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(hp4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(hp5, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(hp6, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Normalization : lognormalizaton 
sp1_1 = seuratWrap2(sp1_1, nFeature.low = 200, nFeature.hi = 4000,mt.hi = 15)
sp1_2 = seuratWrap2(sp1_2, nFeature.low = 200, nFeature.hi = 4000,mt.hi = 10)
sp2 = seuratWrap2(sp2, nFeature.low = 200, nFeature.hi = 4000,mt.hi = 20)
sp3 = seuratWrap2(sp3, nFeature.low = 200, nFeature.hi = 5000,mt.hi = 15)
sp4 = seuratWrap2(sp4, nFeature.low = 200, nFeature.hi = 5000,mt.hi = 15)
sp5 = seuratWrap2(sp5, nFeature.low = 200, nFeature.hi = 4000,mt.hi = 15)
sp6 = seuratWrap2(sp6, nFeature.low = 200, nFeature.hi = 4000,mt.hi = 20)
sp7 = seuratWrap2(sp7, nFeature.low = 200, nFeature.hi = 4000,mt.hi = 15)

hp1 = seuratWrap2(hp1, nFeature.low = 200, nFeature.hi = 5200,mt.hi = 20)
hp2 = seuratWrap2(hp2, nFeature.low = 200, nFeature.hi = 5500,mt.hi = 20)
hp3 = seuratWrap2(hp3, nFeature.low = 200, nFeature.hi = 4000,mt.hi = 20)
hp4 = seuratWrap2(hp4, nFeature.low = 200, nFeature.hi = 4000,mt.hi = 20)
hp5 = seuratWrap2(hp5, nFeature.low = 200, nFeature.hi = 4000,mt.hi = 20)
hp6 = seuratWrap2(hp6, nFeature.low = 200, nFeature.hi = 5000,mt.hi = 20)

# filter out epithiel cell by 'PTPRC' negative and 'KRT18' positive
sp1_1 = subset(sp1_1, PTPRC > 0 )
sp1_1_count = sp1_1@assays$RNA@counts
saveRDS(sp1_1_count, file = 'sp1_1_count.rds')
sp1_2 = subset(sp1_2, PTPRC > 0 )
sp1_2_count  = sp1_2@assays$RNA@counts
saveRDS(sp1_2_count, file = 'sp1_2_count.rds')
sp2 = subset(sp2, PTPRC > 0 )
sp2_count  = sp2@assays$RNA@counts
saveRDS(sp2_count, file = 'sp2_count.rds')
sp3 = subset(sp3, PTPRC > 0 )
sp3_count  = sp3@assays$RNA@counts
saveRDS(sp3_count, file = 'sp3_count.rds')
sp4 = subset(sp4, PTPRC > 0 )
sp4_count  = sp4@assays$RNA@counts
saveRDS(sp4_count, file = 'sp4_count.rds')
sp5 = subset(sp5, PTPRC > 0 )
sp5_count  = sp5@assays$RNA@counts
saveRDS(sp5_count, file = 'sp5_count.rds')
sp6 = subset(sp6, PTPRC > 0 )
sp6_count  = sp6@assays$RNA@counts
saveRDS(sp6_count, file = 'sp6_count.rds')
sp7 = subset(sp7, PTPRC > 0 )
sp7_count  = sp7@assays$RNA@counts
saveRDS(sp7_count, file = 'sp7_count.rds')

hp1 = subset(hp1, PTPRC > 0 )
hp1_count  = hp1@assays$RNA@counts
saveRDS(hp1_count, file = 'hp1_count.rds')
hp2 = subset(hp2, PTPRC > 0 )
hp2_count  = hp2@assays$RNA@counts
saveRDS(hp2_count, file = 'hp2_count.rds')
hp3 = subset(hp3, PTPRC > 0 )
hp3_count  = hp3@assays$RNA@counts
saveRDS(hp3_count, file = 'hp3_count.rds')
hp4 = subset(hp4, PTPRC > 0 )
hp4_count  = hp4@assays$RNA@counts
saveRDS(hp4_count, file = 'hp4_count.rds')
hp5 = subset(hp5, PTPRC > 0 )
hp5_count  = hp5@assays$RNA@counts
saveRDS(hp5_count, file = 'hp5_count.rds')
hp6 = subset(hp6, PTPRC > 0 )
hp6_count  = hp6@assays$RNA@counts
saveRDS(hp6_count, file = 'hp6_count.rds')



sp1_1 = readRDS('Counts/sp1_1_count.rds')
sp1_2 = readRDS('Counts/sp1_2_count.rds')
sp2 = readRDS('Counts/sp2_count.rds')
sp3 = readRDS('Counts/sp3_count.rds')
sp4 = readRDS('Counts/sp4_count.rds')
sp5 = readRDS('Counts/sp5_count.rds')
sp6 = readRDS('Counts/sp6_count.rds')
sp7 = readRDS('Counts/sp7_count.rds')

hp1 = readRDS('Counts/hp1_count.rds')
hp2 = readRDS('Counts/hp2_count.rds')
hp3 = readRDS('Counts/hp3_count.rds')
hp4 = readRDS('Counts/hp4_count.rds')
hp5 = readRDS('Counts/hp5_count.rds')
hp6 = readRDS('Counts/hp6_count.rds')



sp1_1 <-CreateSeuratObject(sp1_1, project = "sp1_1", min.cells = 3, min.features = 200)
sp1_2 <-CreateSeuratObject(sp1_2, project = "s1_2", min.cells = 3, min.features = 200)
sp2 <-CreateSeuratObject(sp2, project = "sp2", min.cells = 3, min.features = 200)
sp3 <-CreateSeuratObject(sp3, project = "sp3", min.cells = 3, min.features = 200)
sp4 <-CreateSeuratObject(sp4, project = "sp4", min.cells = 3, min.features = 200)
sp5 <-CreateSeuratObject(sp5, project = "sp5", min.cells = 3, min.features = 200)
sp6 <-CreateSeuratObject(sp6, project = "sp6", min.cells = 3, min.features = 200)
sp7 <-CreateSeuratObject(sp7, project = "sp7", min.cells = 3, min.features = 200)
hp1 <-CreateSeuratObject(hp1, project = "hp1", min.cells = 3, min.features = 200)
hp2 <-CreateSeuratObject(hp2, project = "hp2", min.cells = 3, min.features = 200)
hp3 <-CreateSeuratObject(hp3, project = "hp3", min.cells = 3, min.features = 200)
hp4 <-CreateSeuratObject(hp4, project = "hp4", min.cells = 3, min.features = 200)
hp5 <-CreateSeuratObject(hp5, project = "hp5", min.cells = 3, min.features = 200)
hp6 <-CreateSeuratObject(hp6, project = "hp6", min.cells = 3, min.features = 200)


sp1_1[["percent.mt"]] <- PercentageFeatureSet(sp1_1, pattern = "^MT-") 
sp1_2[["percent.mt"]] <- PercentageFeatureSet(sp1_2, pattern = "^MT-") 
sp2[["percent.mt"]] <- PercentageFeatureSet(sp2, pattern = "^MT-") 
sp3[["percent.mt"]] <- PercentageFeatureSet(sp3, pattern = "^MT-") 
sp4[["percent.mt"]] <- PercentageFeatureSet(sp4, pattern = "^MT-") 
sp5[["percent.mt"]] <- PercentageFeatureSet(sp5, pattern = "^MT-") 
sp6[["percent.mt"]] <- PercentageFeatureSet(sp6, pattern = "^MT-") 
sp7[["percent.mt"]] <- PercentageFeatureSet(sp7, pattern = "^MT-") 
hp1[["percent.mt"]] <- PercentageFeatureSet(hp1, pattern = "^MT-") 
hp2[["percent.mt"]] <- PercentageFeatureSet(hp2, pattern = "^MT-") 
hp3[["percent.mt"]] <- PercentageFeatureSet(hp3, pattern = "^MT-") 
hp4[["percent.mt"]] <- PercentageFeatureSet(hp4, pattern = "^MT-") 
hp5[["percent.mt"]] <- PercentageFeatureSet(hp5, pattern = "^MT-") 
hp6[["percent.mt"]] <- PercentageFeatureSet(hp6, pattern = "^MT-") 


sp1_1$disease = 'COVID'
sp1_2$disease = 'COVID'
sp2$disease = 'COVID'
sp3$disease = 'COVID'
sp4$disease = 'COVID'
sp5$disease = 'COVID'
sp6$disease = 'COVID'
sp7$disease = 'COVID'
hp1$disease = 'Healthy'
hp2$disease = 'Healthy'
hp3$disease = 'Healthy'
hp4$disease = 'Healthy'
hp5$disease = 'Healthy'
hp6$disease = 'Healthy'


sp1_1$id = 'sp1'
sp1_2$id = 'sp1'
sp2$id = 'sp2'
sp3$id = 'sp3'
sp4$id = 'sp4'
sp5$id = 'sp5'
sp6$id = 'sp6'
sp7$id = 'sp7'
hp1$id = 'hp1'
hp2$id = 'hp2'
hp3$id = 'hp3'
hp4$id = 'hp4'
hp5$id = 'hp5'
hp6$id = 'hp6'



##########################################################################
####  Integration 

anchor_dim = 30

covid_p_anchors <- FindIntegrationAnchors(object.list = list(sp1_1, sp1_2, sp2, sp3, sp4, sp5, sp6, sp7,
				hp1, hp2, hp3, hp4, hp5, hp6), dims = 1:anchor_dim) # may change 30 to 50
#sox2.anchors <- FindIntegrationAnchors(object.list = list(sox, wk2, wk4), dims = 1:20)
covid.combined <- IntegrateData(anchorset = covid.anchors, dims = 1:anchor_dim)
DefaultAssay(covid.combined) <- "integrated"




covid_P_S_anchors1 <- FindIntegrationAnchors(object.list = list(sp1_1, sp1_2, sp2, sp3),  
		dims = 1:anchor_dim) # may change 30 to 50
save(covid_P_S_anchors1, file = 'covid_P_S_anchors1.rda')

covid_P_S_combined1 <- IntegrateData(anchorset = covid_P_S_anchors1, dims = 1:anchor_dim)


covid_P_S_anchors <- FindIntegrationAnchors(object.list = list(sp1_1, sp1_2, sp2, sp3), sp4, sp5, sp6, sp7), 
dims = 1:anchor_dim, k.filter) # may change 30 to 50
save(covid_P_S_anchors, file = 'covid_P_S_anchors.rda')


covid_P_H_anchors <- FindIntegrationAnchors(object.list = list(hp1, hp2, hp3, hp4, hp5, hp6), 
		dims = 1:anchor_dim, k.filter = 150) # may change 30 to 50
save(covid_P_H_anchors, file = 'covid_P_H_anchors.rda')

rm(sp1_1, sp1_2, sp2, sp3, sp4, sp5, sp6, sp7, hp1, hp2, hp3, hp4, hp5, hp6)

covid_P_S_combined <- IntegrateData(anchorset = covid_P_S_anchors, dims = 1:anchor_dim)

covid_P_H_combined <- IntegrateData(anchorset = covid_P_H_anchors, dims = 1:anchor_dim)











