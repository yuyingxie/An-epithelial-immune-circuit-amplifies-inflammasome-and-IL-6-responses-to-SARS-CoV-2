library(Seurat)
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
library(phateR) 
library(dplyr)   # define funciton %>%
library(hdf5r)
library(RColorBrewer)


set.seed(1)
source('COVID_Function.R')

T.integrated = readRDS('T_SCT_integrated_11_26_2020.rds')


darkcols <- brewer.pal(12, "Paired")
bar.col = c(brewer.pal(12, 'Set3'))
hist(discoveries, col = darkcols)

cluster_per = matrix(0, 40, 30)
rownames(cluster_per) = unique(T.integrated$Sample)
colnames(cluster_per) = 1:30

Disease = list()
Sample_id = unique(T.integrated$Sample)
for(i in 1:40){
    temp_1 =  T.integrated$Sample == Sample_id[i]
    Disease[Sample_id[i]] = unique(T.integrated$Disease[temp_1])
    n = sum(temp_1)
    for(j in 1:30){
        temp_2 = T.integrated$seurat_clusters == (j - 1)
        cluster_per[i, j] = sum((temp_1 + temp_2) == 2)/n
        cat('(', i, j, ')', '\n')
    }
}
 
D = unlist(Disease)

id = c((1:40)[D == 'Healthy'], (1:40)[D == 'Moderate'], (1:40)[D == 'Severe'])

dat = cluster_per[id, ]

pdf("Cluster_bar_11_28_2020.pdf", width = 20, height = 5)
opar <- par(lwd = 0.0003)
par(mar = c(5, 3, 4.1, 12), xpd = TRUE)
barplot(height = t(dat), col = bar.col, width = 2, axes = F, xlab = "",  space = 0.1,  border = NA, las = 2) 
axis(2)
#axis(1)
legend(x = 96, y = 0.95, 
		legend = (colnames(dat)), #in order from top to bottom
		fill = bar.col, # 6:1 reorders so legend order matches graph
		title = "Cluster", ncol = 1, cex = 0.75, bty = "n")
dev.off()


id = c((1:40)[D == 'Healthy'], (1:40)[D == 'Moderate'], (1:40)[D == 'Severe'])

dat = cluster_per[id, ]

pdf("Cluster_bar_horiz_11_28_2020.pdf", width = 20, height = 15)
opar <- par(lwd = 0.0003)
par(mar = c(5, 3, 4.1, 12), xpd = TRUE)
barplot(height = t(dat), col = bar.col, width = 2, axes = F, xlab = "",  space = 0.1,  border = NA, las = 2, horiz = T) 
#axis(2)
axis(1)
legend(x = 1.05, y = 90, 
		legend = (colnames(dat)), #in order from top to bottom
		fill = bar.col, # 6:1 reorders so legend order matches graph
		title = "Cluster", ncol = 1, cex = 0.75, bty = "n")
dev.off()

cluster_m = matrix(0, 3, 30)
group_name = c('Healthy', 'Moderate', 'Severe')
for(i in 1:3){
     tmp_id =  (1:40)[D == group_name[i]]
     tmp = cluster_per[tmp_id, ]
     cluster_m[i, ] = apply(tmp, 2, mean)

}


pdf("Cluster_bar_mean_horiz_11_29_2020.pdf", width = 20, height = 7)
opar <- par(lwd = 0.0003)
par(mar = c(5, 6, 4.1, 12), xpd = TRUE)
barplot(height = t(cluster_m), col = bar.col, width = 2, axes = F, xlab = "",  space = 0.1,
      border = NA, las = 2, horiz = T,  names.arg=c("Healthy", "Moderate", "Severe")) 
#axis(2)
axis(1)
legend(x = 1.05, y = 7, 
		legend = (colnames(dat)), #in order from top to bottom
		fill = bar.col, # 6:1 reorders so legend order matches graph
		title = "Cluster", ncol = 1, cex = 0.75, bty = "n")
dev.off()


res = rep(0, 30)
for(i in 1:30){
    tmp = data.frame(clus = cluster_per[, i], Dis = D)
    A = kruskal.test(clus ~ Dis, data = tmp)
    res[i] = A$p.value
}

# pvlue
# [1] 0.294 0.012 0.118 0.035 0.004,  0.071 0.001 0.158 0.002 0.611 
#     0.045 0.231 0.117 0.002 0.006,  0.940 0.009 0.117 0.060 0.331
#     0.367 0.055 0.014 0.002 0.001   0.291 0.000 0.334 0.055 0.293

# 4, 5, 7, 9, 14, 15, 17, 24, 25, 27, 

sig_id = c(4, 5, 7, 9, 14, 15, 17, 24, 25, 27)


for(i in sig_id){
    pdf(paste('T_Cluster', i, '_percentage_11_29_2020.pdf'),
     width = 6, height = 6)
    dat_tmp = cluster_per[, i]
    dat1 = data.frame(cluster_percentage = dat_tmp, Disease = D)
    print(ggplot(data = dat1, aes(y = cluster_percentage, fill = Disease)) + geom_boxplot() +
    ggtitle(paste("Cluster", i)) + ylab( 'Percentage'))
    dev.off()
    cat(i, '\n')
}

