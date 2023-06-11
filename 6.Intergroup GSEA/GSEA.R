
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")
#BiocManager::install("org.Hs.eg.db")
#BiocManager::install("DOSE")
#BiocManager::install("clusterProfiler")
#BiocManager::install("enrichplot")


#引用包
library(limma)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)

diffFile="all.txt"            #差异分析的结果文件
gmtFile="h.all.v2022.1.Hs.symbols.gmt.txt"     #免疫基因集文件
#setwd("C:\\biowolf\\wgcnaDiagnostic\\10.GSEA")      #设置工作目录

#读取文件,并对输入文件进行整理
rt=read.table(diffFile, header=T, sep="\t", check.names=F)
rt=rt[order(rt$logFC, decreasing=T),]
logFC=as.vector(rt[,2])
names(logFC)=as.vector(rt[,1])

#读入基因集文件
gmt=read.gmt(gmtFile)

#GSEA富集分析
kk=GSEA(logFC, TERM2GENE=gmt, pvalueCutoff=1)
kkTab=as.data.frame(kk)
kkTab=kkTab[kkTab$p.adjust<0.05,]
write.table(kkTab,file="11.GSEA-50symbol.result.txt",sep="\t",quote=F,row.names = F)

#绘制实验组富集的图形
termNum=5     #展示前5个基因集
kkUp=kkTab[kkTab$NES>0,]
if(nrow(kkUp)>=termNum){
  showTerm=row.names(kkUp)[1:termNum]       #需要展示的基因集名称
  gseaplot=gseaplot2(kk, showTerm, base_size=8, title="Enriched in treat group")
  pdf(file="11.GSEA-50symbol.treat.pdf", width=8, height=6)
  print(gseaplot)
  dev.off()
}

#绘制对照组富集的图形
termNum=8     #展示前5个基因集
kkDown=kkTab[kkTab$NES<0,]
if(nrow(kkDown)>=termNum){
  showTerm=row.names(kkDown)[1:termNum]       #需要展示的基因集名称
  gseaplot=gseaplot2(kk, showTerm, base_size=8, title="Enriched in control group")
  pdf(file="11.GSEA-50symbol.con.pdf", width=8, height=6)
  print(gseaplot)
  dev.off()
}


