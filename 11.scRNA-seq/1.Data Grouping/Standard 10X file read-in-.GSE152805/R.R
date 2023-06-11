library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(SingleR)
library(CCA)
library(clustree)
library(cowplot)
library(monocle)
library(tidyverse)
library(SCpubr)
library(UCell)
library(irGSEA)
library(GSVA)
library(GSEABase)
library(harmony)
library(plyr)



data_dir <- paste0(getwd(),"/data")       
###Download the GSE152805 standard 10X file from 
###the website (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152805)

samples=list.files(data_dir)
dir=file.path(data_dir,samples)
afdata <- Read10X(data.dir = dir)

#创建seurat对象
# min.cell：每个feature至少在多少个细胞中表达
# min.features：每个细胞中至少有多少个feature被检测到
#nFeature_RNA是每个细胞中检测到的基因数量
#nCount_RNA是细胞内检测到的分子总数
#nFeature_RNA过低，表示该细胞可能已死/将死或是空液滴。高nCount_RNA和/或nFeature_RNA表明“细胞”实际上可以是两个或多个细胞。
#结合线粒体基因（percent.mt）与核糖体基因（percent.rb）除去异常值，即可除去大多数双峰/死细胞/空液滴，因此它们过滤是常见的预处理步骤
af <- CreateSeuratObject(counts = afdata, 
                         project = "SeuratObject", 
                         min.cells = 3,
                         min.features = 200)
afidens=mapvalues(Idents(af), from = levels(Idents(af)), to = samples)
Idents(af)=afidens
af$Type=Idents(af)
af[["percent.mt"]] <- PercentageFeatureSet(af, pattern = "^MT-")
#af[["percent.rb"]] <- PercentageFeatureSet(af, pattern = "^RP")
mask1 <- af$nCount_RNA >= 1000 & af$nCount_RNA<= 35000
mask2 <- af$nFeature_RNA >= 200 & af$nFeature_RNA <= 10000
mask3 <- af$percent.mt <= 10
#mask4<-af$percent.rb<= 20
mask <-mask1 & mask2 & mask3
af <- af[, mask]
saveRDS(af,"af.rds")

