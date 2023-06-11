
#install.packages("celldex")
#devtools::install_github("sqjin/CellChat")
#BiocManager::install("celldex")
library(limma)
library(NMF)
library(ggalluvial)
library(svglite)
library(CellChat)
library(celldex)
library(Seurat)
library(tidyverse)
library(Matrix)
library(stringr)
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
library(randomcoloR)
library(CellChat)
library(ggpubr)
##############################################
##############################################
###################################数据标准化
af=readRDS("af.rds")
afgs="c-FRGs.txt"

# QC指标使用小提琴图可视化,ncol为图片排列列数
pdf(file = "01.vlnplot.pdf",width = 8,height = 5)
VlnPlot(af, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)+scale_fill_manual(values = c("#C77CFF","#7CAE00","#00BFC4","#F8766D"))
dev.off()

# 指标之间的相关性
plot1 <- FeatureScatter(af, feature1 = "nCount_RNA", feature2 = "percent.mt")+ RotatedAxis()
#plot2 <- FeatureScatter(af, feature1 = "nCount_RNA", feature2 = "percent.rb")+ RotatedAxis()
plot3 <- FeatureScatter(af, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")+ RotatedAxis()
#组图
pdf(file = "01.corqc.pdf",width =12,height = 5)
plot1+plot3+plot_layout(ncol = 2)      #plot_layout，patchwork函数，指定一行有几个图片
dev.off()

##标准化,使用LogNormalize方法
af <- NormalizeData(af, normalization.method = "LogNormalize", scale.factor = 10000)

## 鉴定高变基因
# 高变基因：在一些细胞中表达高，另一些细胞中表达低的基因
# 变异指标： mean-variance relationship
# 返回2000个高变基因，用于下游如PCA降维分析。
af <- FindVariableFeatures(af, selection.method = "vst", nfeatures = 2000)

# 提取前10的高变基因
top10 <- head(VariableFeatures(af), 10)
top10

# 展示高变基因
plot1 <- VariableFeaturePlot(af)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

pdf(file = "01.topgene.pdf",width =7,height = 6)
plot2                   #plot_layout，patchwork函数，指定一行有几个图片
dev.off()

## 归一化
# 归一化处理：每一个基因在所有细胞中的均值变为0，方差标为1，对于降维来说是必需步骤
# 归一化后的值保存在：af[["RNA"]]@scale.data
af <- ScaleData(af)

# 可以选择全部基因归一化
all.genes <- rownames(af)
af <- ScaleData(af, features = all.genes)

##########0.3.part3 降维(绘制原始分布)##########################
# PCA降维，用前面1500个高变基因，可以使用features改变用于降维的基因集
af <- Seurat::RunPCA(af, features = VariableFeatures(object = af))
af <- Seurat::RunTSNE(af,dims = 1:20)
pdf(file = "02.rawtsne.pdf",width =7.5,height = 5.5)
DimPlot(af, reduction = "tsne",pt.size = 1)+theme_classic()+theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),legend.position = "right") #top为图列位置最上方，除此之外还有right、left、bottom(意思同英文)
dev.off()
pdf(file = "02.rawpca.pdf",width =7.5,height = 5.5)
DimPlot(af, reduction = "pca",pt.size = 1)+theme_classic()+theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),legend.position = "right")
dev.off()
colaa=distinctColorPalette(100)
pdf(file = "02.raw.tsne.split.pdf",width =8,height =5)
do_DimPlot(sample = af,
           plot.title = "",
           reduction = "tsne",
           legend.position = "bottom",
           dims = c(1,2),split.by = "Type",pt.size =0.5
) #选择展示的主成分，这边是PC2与PC1
dev.off()

#####################################（选做 harmony 去批次与降维）################################################################################################
af <- RunHarmony(af, group.by.vars = "Type")
###########################################################################################################
#####################################################################
#######################################0.3.part3 矫正后结果可视化
#################################################################
# PCA降维，用前面1500个高变基因，可以使用features改变用于降维的基因集
pdf(file = "03.harmony.pdf",width =7.5,height = 5.5)
DimPlot(af, reduction = "harmony",pt.size = 1)+theme_classic()+theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),legend.position = "right")
dev.off()
af <- Seurat::RunTSNE(af,dims = 1:20,reduction ='harmony')
pdf(file = "03.tsne.pdf",width =7.5,height = 5.5)
DimPlot(af, reduction = "tsne",pt.size = 1)+theme_classic()+theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),legend.position = "right")
dev.off()
collist=c(ggsci::pal_nejm()(8))
names(collist)=names(table(af$Type))
pdf(file = "03.tsne.split.pdf",width =12,height = 7.5)
do_DimPlot(sample = af,
           plot.title = "",
           reduction = "tsne",
           legend.position = "bottom",
           dims = c(1,2),split.by = "Type",pt.size =0.5
) #选择展示的主成分，这边是PC2与PC1
dev.off()

#########################################################################################################################################
collist=c(ggsci::pal_nejm()(8))
names(collist)=names(table(af$Type))
# 前两个PC特征基因可视化
VizDimLoadings(af, dims = 1:2, reduction = "pca")
#热图可视化前15个PC
pdf(file = "04.pc_heatmap.pdf",width =7.5,height = 9)
DimHeatmap(af, dims = 1:20, cells = 1000, balanced = TRUE)
dev.off()
##确定使用PC个数
# each PC essentially representing a ‘metafeature’
af <- JackStraw(af, num.replicate = 100)
af <- ScoreJackStraw(af, dims = 1:20)
pdf(file = "04.jackstrawplot.pdf",width =7.5,height = 5.5)
JackStrawPlot(af, dims = 1:20)
dev.off()
pdf(file = "04.ElbowPlot.pdf",width =5,height = 4)
ElbowPlot(af,ndims = 30,reduction = "harmony")
dev.off()

#选择PC
afPC=12
##对细胞聚类
# 首先基于PCA空间构建一个基于欧氏距离的KNN图
#af <- FindNeighbors(af, dims = 1:15)
#设置不同的分辨率，观察分群效果，dim为PCA选择的主成分数
af=FindNeighbors(af, dims = 1:afPC, reduction = "harmony")
for (res in c(0.01, 0.05, 0.1, 0.2, 0.3, 0.6,0.8,1,1.2,1.5,2,2.5,3)) {
  af=FindClusters(af, graph.name = "RNA_snn", resolution = res, algorithm = 1)}
apply(af@meta.data[,grep("RNA_snn_res",colnames(af@meta.data))],2,table)

p2_tree=clustree(af@meta.data, prefix = "RNA_snn_res.")
pdf(file = "04.clustertree.pdf",width =12,height =10)
p2_tree
dev.off()
# 聚类并最优化
# resolution参数：值越大，细胞分群数越多，根据前面进行选择
# 0.4-1.2 typically returns good results for single-cell datasets of around 3K cells
# Optimal resolution often increases for larger datasets. 
#选择分辨率进行降维
af=FindNeighbors(af, dims = 1:afPC, reduction = "harmony")
af <- FindClusters(af, resolution = 0.6) #关键！

# 查看聚类数ID
head(Idents(af), 5)

# 查看每个类别多少个细胞
head(af@meta.data)
table(af@meta.data$seurat_clusters)
# 鉴定各个细胞集群的标志基因only.pos：只保留上调差异表达的基因
af.markers <- FindAllMarkers(af, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(af.markers,file = "05.cluster_markers.csv")

## 将细胞在低维空间可视化UMAP/tSNE
af <- RunUMAP(af, dims = 1:afPC, reduction = "harmony")
af <- RunTSNE(af, dims = 1:afPC, reduction = "harmony")

# 可视化UMAP/tSNE
pdf(file = "05-cluster.UMAP.pdf",width =7,height = 5.5)
DimPlot(af, reduction = "umap", label = T, label.size = 3.5,pt.size = 1)+theme_classic()+theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),legend.position = "right")
dev.off()
pdf(file = "05-cluster.TSEN.pdf",width =7,height = 5.5)
DimPlot(af, reduction = "tsne", label = T, label.size = 3.5,pt.size = 1)+theme_classic()+theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),legend.position = "right")
dev.off()


#################
# AddModuleScore计算基因集评分
afgenes=read.table("SignGenes.txt",header = F,sep = "\t")[,1]
afgss=as.character(read.table(afgs,header = F,sep = "\t")[,1])
af <- AddModuleScore(
  object = af,
  features =list(afgss[afgss %in% rownames(af)]) ,
  ctrl = 100, #默认值是100
  name = gsub(".txt","",afgs)
)
colnames(af@meta.data)[length(colnames(af@meta.data))] <- gsub(".txt","",afgs) 
colsa = distinctColorPalette(100) #随机颜色

#####################显著性分析
#group diff
source(file = "vnplot.R")
gene_sig <- gsub(".txt","",afgs)
#comparisons <- list(names(table(af$Type)))
#afvp(af,gene_signature = gene_sig, file_name = "06-group.GS_VlnPlot", test_sign = comparisons,pta=0.1,cols=colsa,label="p.format",group="Type",widplot=6,heiplot=5.5) #
#报错添加点
#comparisons <- list(names(table(af$Type)))
#afvp(af,gene_signature = intersect(afgenes,rownames(af)), file_name = "06-group.hub_VlnPlot", test_sign = comparisons,pta=0.1,cols=colsa,label="p.signif",group="Type",widplot=15,heiplot=12,ak=0.9) 
#cell diff
comparisons <- list()
comp=combn(names(table(af$seurat_clusters)),2)
names(table(af$seurat_clusters))
for(j in 1:ncol(comp)){comparisons[[j]]<-comp[,j]}
#afvp(af,gene_signature = gene_sig, file_name = "06-cluster.GS_VlnPlot", test_sign = comparisons,pta=0.1,cols=colsa,label="p.signif",group="seurat_clusters",widplot=20,heiplot=8,ak=0.01,split = "Type")
#cell diff
pdf(file = "06-cluster.hub_VlnPlot.pdf",width =10,height = 6)
VlnPlot(af, features = afgenes,group.by = "seurat_clusters", stack=TRUE,cols = colsa, slot = "data")+ NoLegend()   
dev.off()



################人工注释###########################################################################################################
#实用网址：
#https://www.thermofisher.cn/cn/zh/home/life-science/cell-analysis/cell-analysis-learning-center/immunology-at-work.html
#http://xteam.xbio.top/CellMarker/
#https://www.jianshu.com/p/15dddefc7038
#https://toppgene.cchmc.org/)
af <- FindClusters(af, resolution = 0.6)
genes <- list("EC" = c("PLVAP","PECAM115","VWF16"),
              "SMC"=c("ACTA2"),
              "Macrophage" = c("CD68", "CD163","CD14"),
              "T cells" = c("CD3D", "CD2", "CD4", "CD5", "CD7", "GZMK", "CCL5", "IL7R", "IL32", "NKG7"),
              "Mast cells" = c("TPSAP1","TPSAB1","CPA3","HPGDS","MS4A2","GATA2"),
              "Fibroblasts" = c("PRG4","HAS1","HTRA1","CXCL12"),
              "NK cells" = c("NCAM1","CD56","IFNγ"),
              "NKT cells" = c("ZBTB16","PLZF"),
              "B cells"=c("MS4A1","CD79A","CD19", "CD20", "CD22","RALGPS2","MZB1","IER5","CD37"),
              "DC" = c("CD11b","CD13","CD33", "CD80", "CD83")
)
pdf(file = "08.ann_cluster_marker.pdf",width =20,height = 15)
do_DotPlot(sample = af,features = genes,dot.scale = 12,colors.use = c("yellow","red"),legend.length = 50,
           legend.framewidth = 2, font.size =12)
dev.off()
#人工注释
table(af@active.ident)
ann.ids <- c("Fibroblasts",  #cluster0
             "Fibroblasts",  #cluster1
             "Fibroblasts",      #以下按顺序操作
             "Fibroblasts",
             "Fibroblasts",
             "Macrophage",
             "Fibroblasts",
             "Fibroblasts",
             "DC",
             "EC",
             "SMC",
             "T cells",
             "Mast cells",
             "Fibroblasts")
length(ann.ids)

afidens=mapvalues(Idents(af), from = levels(Idents(af)), to = ann.ids)
Idents(af)=afidens
af$cellType=Idents(af)
##saveRDS(af,"af1.rds")
#########人工注释后结果可视化
# 可视化UMAP/tSNE
colors=c('#2874C5','#b42e20','#ebc03e','#377b4c',
         '#7bc7cd','#5d84a4','#4619CC')
pdf(file = "08-ann.scRNA.UMAP.pdf",width =7.5,height = 5.5)
DimPlot(af, reduction = "umap", label = T, label.size = 3.5,pt.size = 1,cols = colors)+theme_classic()+theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),legend.position = "right")
dev.off()
pdf(file = "08-ann.scRNA.TSEN.pdf",width =7.5,height = 5.5)
DimPlot(af, reduction = "tsne", label = T, label.size = 3.5,pt.size = 1,cols = colors)+theme_classic()+theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),legend.position = "right")
dev.off()
################人工注释###########################################################################################################
###########################################
#########################
# only.pos：只保留上调差异表达的基因
af.markers <- FindAllMarkers(af, only.pos = F, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(af.markers,file = "08.cell_markers.csv")
# get top 10 genes
top5af.markers <- af.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)

# plot
colors2=c('#2874C5','#b42e20','#ebc03e','#377b4c',
         '#7bc7cd','#5d84a4','#4619CC')
pdf(file = "09-cell_marker.hetmap-1.pdf",width =15,height = 10)
DoHeatmap(af,features = top5af.markers$gene,
          group.colors = colors2) +
  scale_fill_gradient2(low = '#2874C5',mid = 'white',high = '#b42e20',
                       name = 'Z-score')
dev.off()


#################
#11SignGenes分布图
colaa=distinctColorPalette(100)
afgenes=read.table("SignGenes.txt",header = F,sep = "\t")[,1]
pdf(file = "09-cell_FeaturePlot.pdf",width =12,height = 10)
FeaturePlot(af, features = afgenes, cols = c("grey", "#d42e20"),min.cutoff = 0.1, max.cutoff = 1,ncol=4,pt.size = 0.5, slot = "counts")    #min.cutoff与max.cutoff修改截断以更好可视化结果，通过颜色强调基因的分布
dev.off()

colors2=c('#2874C5','#2874C5','#b42e20','#ebc03e','#377b4c','#2874C5','#2874C5',
         '#7bc7cd','#2874C5','#5d84a4','#4619CC')
pdf(file = "09-cell_VlnPlot.pdf",width =11,height = 7)
VlnPlot(af, features = afgenes,group.by = "cellType", stack=TRUE,cols = colors2, slot = "counts")+ NoLegend()   
dev.off()
#剩下29个c-FRGs基因
colaa=distinctColorPalette(100)
afgenes1=read.table("Nonsign-c-FRGs.txt",header = F,sep = "\t")[,1]
pdf(file = "09-cell_FeaturePlot1.pdf",width =15,height = 12)
FeaturePlot(af, features = afgenes1, cols = c("grey", "#2874C5"),min.cutoff = 0.1, max.cutoff = 1,ncol=6,pt.size = 0.5, slot = "counts")    #min.cutoff与max.cutoff修改截断以更好可视化结果，通过颜色强调基因的分布
dev.off()
pdf(file = "09-cell_VlnPlot1.pdf",width =20,height = 11)
VlnPlot(af, features = afgenes1,group.by = "cellType", stack=TRUE,cols = colaa, slot = "counts")+ NoLegend()   
dev.off()
#####################显著性分析
#group diff
source(file = "vnplot.R")
gene_sig <- gsub(".txt","",afgs)
#cell diff
comparisons <- list()
comp=combn(names(table(af$cellType)),2)
names(table(af$cellType))
for(j in 1:ncol(comp)){comparisons[[j]]<-comp[,j]}
afvp(af,gene_signature = gene_sig, file_name = "09-cell.GS.stat_VlnPlot", test_sign = comparisons,pta=0.1,cols=colsa,label="p.signif",group="cellType",widplot=10,heiplot=10,ak=0.9)
#cell diff
comparisons <- list()
comp=combn(names(table(af$cellType)),2)
names(table(af$cellType))
for(j in 1:ncol(comp)){comparisons[[j]]<-comp[,j]}
afvp(af,gene_signature = intersect(afgenes,rownames(af)), file_name = "09-cell.hub.stat_VlnPlot", test_sign = comparisons,pta=0.1,cols=colsa,label="p.signif",group="cellType",widplot=20,heiplot=15,ak=3)

#####################特定基因的细胞亚群比例
geneselect="PTGS2"
cellselect="Fibroblasts"
af$geneType=ifelse(as.numeric(as.matrix(af@assays$RNA@scale.data)[geneselect,])>median(sort(as.numeric(as.matrix(af[,which(af$cellType %in% c(cellselect))]@assays$RNA@scale.data)[geneselect,]))),paste0("High ",geneselect," ",cellselect),paste0("Low ",geneselect," ",cellselect))
Cellratio <- prop.table(table( af[,which(af$cellType %in% c(cellselect))]$geneType,af[,which(af$cellType %in% c(cellselect))]$Type), margin = 2)#计算各组样本不同细胞群比例
Cellratio <- as.data.frame(Cellratio)
colnames(Cellratio)[1]="Celltype"
colourCount = length(unique(Cellratio$Celltype))
colaa=distinctColorPalette(100)
ggplot(Cellratio) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Celltype),stat = "identity",width = 0.7,linewidth = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='Type',y = 'Ratio')+
  coord_flip()+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),legend.position = "right")+   # 图例："left" 左, "right" 右,  "bottom" 下, "top" 上
  scale_fill_manual(values=colaa)
ggsave("10-cell_geneselect_ration.pdf",width = 6,height = 3.5)

afvp(af,gene_signature = geneselect, file_name = "06-cell.geneselect._VlnPlot", test_sign = comparisons,pta=0.1,cols=colsa,label="p.signif",group="Type",widplot=12,heiplot=8,ak=0.9,split = "cellType")

afc=af[,which(af$cellType %in% c(cellselect))]
Idents(afc)=afc$geneType
# 可视化UMAP/tSNE
pdf(file = "10-sg.scRNA.UMAP.pdf",width =7.5,height = 5.5)
DimPlot(afc, reduction = "umap", label = T, label.size = 3.5,pt.size = 1)+theme_classic()+theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),legend.position = "right")
dev.off()
pdf(file = "10-sg.scRNA.TSEN.pdf",width =7.5,height = 5.5)
DimPlot(afc, reduction = "tsne", label = T, label.size = 3.5,pt.size = 1)+theme_classic()+theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),legend.position = "right")
dev.off()

##############################################################
######################使用irGSEA进行基因集分析：https://github.com/chuiqin/irGSEA#########################################
library(irGSEA)
library(GSVA)
library(GSEABase)
#irGSEA调用msigdbr中的基因集进行后续分析，注意哈，更改category参数即可指定基因集进行分析，支持的基因集有如下:
#"【category"】	 【description】
#"H"	        hallmark gene sets
#"C1"	        positional gene sets
#"C2"	        curated gene sets
#"C3"	        motif gene sets
#"C4"	        computational gene sets
#"C5"	        GO gene sets
#"C6"	        oncogenic signatures
#"C7"	        immunologic signatures
af.final <- irGSEA.score(object = af[,which(af$cellType %in% c(cellselect))], 
                         assay = "RNA", 
                         slot = "data", 
                         seeds = 123, 
                         ncores = 1,
                         min.cells = 3, 
                         min.feature = 0,
                         custom = F, 
                         msigdb = T, 
                         species = "Homo sapiens", 
                         category = "H",  
                         geneid = "symbol",
                         method = c("AUCell", "UCell", "singscore", 
                                    "ssgsea"),
                         kcdf = 'Gaussian')
result.dge <- irGSEA.integrate(object = af.final, 
                               group.by = "geneType",
                               method = c("AUCell","UCell","singscore",
                                          "ssgsea"))
irGSEA.heatmap.plot <- irGSEA.heatmap(object = result.dge, 
                                      method = "RRA",
                                      top = 50)
pdf(file = "11-hallmark.heatmap.pdf",width =10,height = 8)
irGSEA.heatmap.plot
dev.off()


##############################################
#######################
#######高低细胞组的差异分析
af.markers <- FindAllMarkers(afc, only.pos = F, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(af.markers,file = "11.geneGroup_Diff.csv")
###enrichment
library(clusterProfiler)
library(org.Hs.eg.db)
ids=bitr(af.markers$gene,'SYMBOL','ENTREZID','org.Hs.eg.db') ## 将SYMBOL转成ENTREZID
af.markers=merge(af.markers,ids,by.x='gene',by.y='SYMBOL')
View(af.markers)
## 函数split()可以按照分组因子，把向量，矩阵和数据框进行适当的分组。
## 它的返回值是一个列表，代表分组变量每个水平的观测。
gcSample=split(af.markers$ENTREZID, af.markers$cluster) 
## KEGG,12,15
xx <- compareCluster(gcSample,
                     fun = "enrichKEGG",
                     organism = "hsa",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05
)
write.csv(as.data.frame(xx),file = "11.geneGroup_KEGG.csv")
p <- dotplot(xx)
pdf(file = "11-geneGroup_KEGG.pdf",width =8,height = 8)
p +scale_y_discrete(labels=function(x) stringr::str_wrap(x, width=60))+ theme(axis.text.x = element_text(
  angle = 45,
  vjust = 0.5, hjust = 0.5
))
dev.off()
## GO
xx <- compareCluster(gcSample,
                     fun = "enrichGO",
                     OrgDb = "org.Hs.eg.db",
                     #ont = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.05
)
write.csv(as.data.frame(xx),file = "11.geneGroup_GO.csv")
p <- dotplot(xx)
pdf(file = "11-geneGroup_GO.pdf",width =8,height = 8)
p+scale_y_discrete(labels=function(x) stringr::str_wrap(x, width=60)) + theme(axis.text.x = element_text(
  angle = 45,
  vjust = 0.5, hjust = 0.5
))
dev.off()

# get top 10 genes
top5af.markers <- af.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)

comparisons <- list()
comp=combn(names(table(afc$geneType)),2)
names(table(afc$geneType))
for(j in 1:ncol(comp)){comparisons[[j]]<-comp[,j]}
afvp(af=afc,gene_signature = top5af.markers$gene, file_name = "11-geneGroup_VlnPlot", test_sign = comparisons,pta=0.1,cols=colsa,label="p.signif",group="geneType",widplot=16,heiplot=12,ak=0.9)

##############################################
#######################
#######高低细胞组的拟时序分析
logFCfilter=1           
adjPvalFilter=0.05
af.markers=af.markers[(abs(as.numeric(as.vector(af.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(af.markers$p_val_adj))<adjPvalFilter),]
monocle.matrix=as.matrix(afc@assays$RNA@counts, 'sparseMatrix')
afmetadata=afc@meta.data
monocle.sample=afmetadata[,8,drop=F]
monocle.geneAnn=data.frame(gene_short_name = row.names(monocle.matrix), row.names = row.names(monocle.matrix))
monocle.geneAnn$gene_kk_name=monocle.geneAnn$gene_short_name
monocle.markers=af.markers

#将seurat对象转化为monocle输入格式
data <- as(as.matrix(monocle.matrix), 'sparseMatrix')
pd<-new("AnnotatedDataFrame", data = monocle.sample)
fd<-new("AnnotatedDataFrame", data = monocle.geneAnn)
cds <- newCellDataSet(data, phenoData = pd, featureData = fd)
names(pData(cds))[names(pData(cds))=="Size_Factor"]="Cluster"
ssss=as.data.frame(Idents(afc))
pData(cds)[,"geneType"]=paste0(ssss$`Idents(afc)`)
pData(cds)$Type=afc$Type
source("order_cells.R")
source("beam.R")
library(igraph)
library(SingleCellExperiment)
library(ggsci)
#开始细胞轨迹分析
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
#monocle选择高变基因
disp_table <- dispersionTable(cds)
disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
cds <- setOrderingFilter(cds, disp.genes)
#plot_ordering_genes(cds)
cds <- reduceDimension(cds, max_components = 2, reduction_method = 'DDRTree',auto_param_selection = F)
cds <- orderCells(cds)
#树枝的细胞轨迹图
pdf(file = "12.cds_geneType.pdf",width =5,height = 5)
m1=plot_cell_trajectory(cds,color_by = "geneType", cell_size = 0.5)#+facet_wrap(~cell_type2,ncol=5)  #适当利用分面
m1
dev.off()
#时间的细胞轨迹图
pdf(file = "12.cds_time.pdf",width =5,height = 5)
m2=plot_cell_trajectory(cds,color_by = "Pseudotime", cell_size = 0.5) #+facet_wrap(~cell_type2,ncol=5)  #适当利用分面
m2
dev.off()
#细胞名称的细胞轨迹图
pdf(file = "12.cds_state.pdf",width =6.5,height = 7)
m3=plot_cell_trajectory(cds,color_by = "State", cell_size = 0.5) #+facet_wrap(~celltype,ncol=3) #适当利用分面
m3
dev.off()

pdf(file = "12.cds_GS.pdf",width =8,height = 7)
pData(cds)[,gene_sig] = afc@meta.data[,gene_sig]
plot_cell_trajectory(cds, color_by = gene_sig)  + scale_color_gsea()
dev.off()

pdf(file = "12.cds_geneselect.pdf",width =8,height = 7)
pData(cds)[,geneselect] = afc@assays$RNA@scale.data[geneselect,]
plot_cell_trajectory(cds, color_by = geneselect)  + scale_color_gsea()
dev.off()

#这里用的是disp.genes，https://www.jianshu.com/p/9995cd707002
BEAM_res <- BEAM(cds[disp.genes,], branch_point = 2, cores = 2) 
#BEAM_res <- BEAM(cds, branch_point = 1, cores = 2)  #也可以对所有基因基因进行排序
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
head(BEAM_res)
write.csv(BEAM_res, "12.BEAM_res.csv", row.names = F)


pdf(file = "12.BEAM_cluster.pdf",width =14,height = 20)
plot_genes_branched_heatmap(cds[row.names(subset(BEAM_res,
                                                 qval < 1e-6)),],
                            branch_point = 1, #绘制的是哪个分支
                            num_clusters = 4, #分成几个cluster，根据需要调整
                            cores = 2,
                            use_gene_short_name = T,
                            show_rownames = T)#展示行名
dev.off()

a=plot_genes_branched_heatmap(cds[row.names(subset(BEAM_res, qval < 1e-6)),],
                              branch_point = 1, #绘制的是哪个分支
                              num_clusters = 4, #分成几个cluster，根据需要调整
                              cores = 2,
                              use_gene_short_name = T,
                              show_rownames = T, return_heatmap = T)

#保存基因及其对应的聚类类型
clusters <- a$annotation_row
clusters <- data.frame(clusters)
write.csv(clusters,file = "12.BEAM_cluster.csv",quote = F,row.names = T,col.names = T)
##############################细胞通讯的比较分析
##############################
aflist=SplitObject(af,split.by = "Type") 
#创建CellChat 对象：
#用户可以从数据矩阵、Seurat 或SingleCellExperiment对象创建新的 CellChat 对象。如果输入是 Seurat 或SingleCellExperiment对象，则默认情况下将使用对象中的meta data，用户必须提供该数据来定义细胞分组。例如，group.by=Seurat 对象中默认的细胞标识（我们在此处使用前面放置的细胞类型注释“celltype”）。
#https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/CellChat-vignette.html
#创建cellchat对象
af1=readRDS("af1.rds")
af1
cellchat <- createCellChat(object = af1, group.by = 'cellType', assay = "RNA")
cellchat@DB <- CellChatDB.human
ppi = PPI.human

#对表达数据进行预处理
#CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation costfuture::plan("multiprocess", workers = 2) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat) ##识别每种细胞中过表达的基因
cellchat <- identifyOverExpressedInteractions(cellchat)  #识别基因的相互作用
cellchat <- projectData(cellchat, ppi)

#计算细胞通讯的概率
cellchat <- computeCommunProb(cellchat)
#过滤掉小于10个细胞的细胞通讯
cellchat <- filterCommunication(cellchat, min.cells = 10)
#输出细胞间的通讯关系
df.net=subsetCommunication(cellchat)
write.table(file="COMM02.Comm.network.xls", df.net, sep="\t", row.names=F, quote=F)
#在信号通路的水平进一步推测胞间的通讯, 推断通路水平的互作网络
cellchat <- computeCommunProbPathway(cellchat)
#对计算结果汇总整合，展示整体细胞通讯状态
cellchat <- aggregateNet(cellchat)

#输出细胞通讯的图形(互作数量的图形) #输出细胞通讯的图形(互作强度的图形)
colors=c('#2874C5','#b42e20','#ebc03e','#377b4c',
         '#7bc7cd','#5d84a4','#4619CC')
groupSize <- as.numeric(table(cellchat@idents))
pdf('12.netVisual_circle-Number.pdf')
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions", edge.label.cex = 1.5,color.use = colors,edge.width.max = 15, arrow.size = 0.4,arrow.width = 2,vertex.label.cex=1.5)
dev.off()
pdf('12.netVisual_circle-strength.pdf')
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength",edge.label.cex = 1.5,color.use = colors,edge.width.max = 15, arrow.size = 0.4,arrow.width = 2,vertex.label.cex=1.5)
dev.off()

#分细胞类型展示(把单个细胞类型提取出来,观察这个细胞与其他细胞的通讯)
pdf(file="12.singleCell.pdf", width=10, height=6)
weight_mat <- cellchat@net$weight
par(mfrow = c(2,4), mgp=c(0,0,0), xpd=TRUE)
for (cel in unique(cellchat@idents)){
  cir_mat <- matrix(0, nrow = nrow(weight_mat), ncol = ncol(weight_mat), dimnames = dimnames(weight_mat))
  cir_mat[cel, ] <- weight_mat[cel, ]
  netVisual_circle(cir_mat, vertex.weight= groupSize, weight.scale= T,edge.weight.max = max(weight_mat), vertex.label.cex=1.5,title.name=cel,edge.label.cex =1,color.use = colors,label.edge = FALSE,edge.width.max = 15, arrow.size = 0.4,arrow.width = 2)
}
dev.off()

#绘制受体配体对的气泡图
pdf(file="12.bubble.pdf", width=10, height=10)
netVisual_bubble(cellchat, remove.isolate = FALSE, angle.x = 45)
dev.off()

pdf(file="12.bubble1.pdf", width=5, height=15)
netVisual_bubble(cellchat, sources.use = 6,targets.use = c(1,2,3,4,5),remove.isolate = FALSE, angle.x = 45)
dev.off()

#指定基因集
cellchatgenes=read.table("SignGenes.txt",header = F,sep = "\t")[,1]
def.hub=df.net[((df.net$ligand %in% cellchatgenes) | (df.net$receptor %in% cellchatgenes)),]
write.table(file="COMM07.Comm.hubNetwork.xls", def.hub, sep="\t", row.names=F, quote=F)

#通路水平的可视化
cellchat@netP$pathways     #展示所有相关通路的名称
pathways.show="EGF"       #选择需要展示的通路(可修改)
#通路的细胞通讯图
pdf(file=paste0("COMM08.", pathways.show , ".circle.pdf"), width=8, height=6)
circle=netVisual_aggregate(cellchat, signaling=pathways.show, layout="circle")#, vertex.size = groupSize
print(circle)
dev.off()
#使用层次图展示通路的细胞通讯
pdf(file=paste0("COMM09.", pathways.show , ".hierarchy.pdf"), width=12, height=6)
hierarchy=netVisual_aggregate(cellchat, signaling=pathways.show, layout="hierarchy",  vertex.receiver=seq(1,4))#, vertex.size = groupSize
print(hierarchy)
dev.off()
#细胞通讯的热图
pdf(file=paste0("COMM10.", pathways.show , ".heatmap.pdf"), width=8, height=6)
heatmap=netVisual_heatmap(cellchat, signaling=pathways.show, color.heatmap = "Reds", measure= 'weight')	
print(heatmap)
dev.off()
#细胞作用类型分析
pdf(file=paste0("COMM11.", pathways.show , ".netAnalysis.pdf"), width=6, height=5)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") 
netAnalysis=netAnalysis_signalingRole_network(cellchat, signaling =pathways.show, width = 8, height = 5, font.size = 12)
print(netAnalysis)
dev.off()

#查看哪些配体和受体对在通路中起作用(配体受体对的贡献程度)
pdf(file=paste0("COMM12.", pathways.show , ".contribution.pdf"), width=8, height=6)
contribution=netAnalysis_contribution(cellchat, signaling= pathways.show)
print(contribution)
dev.off()
#查看通路中互作基因的表达水平
pdf(file=paste0("COMM13.", pathways.show , ".geneExp.pdf"), width=8, height=6)
geneExp=plotGeneExpression(cellchat, signaling=pathways.show)
print(geneExp)
dev.off()

#配体受体对水平的细胞通讯展示
pairLR <- extractEnrichedLR(cellchat, signaling=pathways.show, geneLR.return=FALSE)
pdf(file=paste0("COMM14.", pathways.show , ".pairLR.pdf"), width=9, height=8)
pairCircos=netVisual_individual(cellchat, signaling=pathways.show, pairLR.use=pairLR[1] , layout="circle" )
print(pairCircos)
dev.off()
#对通路中的受体配体对进行循环, 以和弦图形式展示
for(i in 1:nrow(pairLR)){
  pdf(file=paste0("COMM15.", pairLR[i,], ".pairLR.pdf"), width=8, height=6)
  pairChord=netVisual_individual(cellchat, signaling=pathways.show, pairLR.use=pairLR[i,], layout="chord" )
  print(pairChord)
  dev.off()
}



