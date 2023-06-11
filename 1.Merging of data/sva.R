
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("sva")


#引用包
library(limma)
library(sva)
outFile="merge.txt"       #输出文件
#setwd("C:\\biowolf\\wgcnaDiagnostic\\05.sva")      #设置工作目录

#获取目录下所有".txt"结尾的文件
files=dir()
files=grep("txt$", files, value=T)
geneList=list()

#读取所有txt文件中的基因信息，保存到geneList
for(file in files){
  if(file==outFile){next}
  rt=read.table(file, header=T, sep="\t", check.names=F)      #读取输入文件
  geneNames=as.vector(rt[,1])      #提取基因名称
  uniqGene=unique(geneNames)       #基因取unique
  header=unlist(strsplit(file, "\\.|\\-"))
  geneList[[header[1]]]=uniqGene
}

#获取交集基因
interGenes=Reduce(intersect, geneList)

#数据合并
allTab=data.frame()
batchType=c()
for(i in 1:length(files)){
  inputFile=files[i]
  header=unlist(strsplit(inputFile, "\\.|\\-"))
  #读取输入文件，并对输入文件进行整理
  rt=read.table(inputFile, header=T, sep="\t", check.names=F)
  rt=as.matrix(rt)
  rownames(rt)=rt[,1]
  exp=rt[,2:ncol(rt)]
  dimnames=list(rownames(exp),colnames(exp))
  data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
  rt=avereps(data)
  
  #对数值大的数据取log2
  qx=as.numeric(quantile(rt, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
  LogC=( (qx[5]>100) || ( (qx[6]-qx[1])>50 && qx[2]>0) )
  if(LogC){
    rt[rt<0]=0
    rt=log2(rt+1)}
  rt=normalizeBetweenArrays(rt)
  
  #数据合并
  if(i==1){
    allTab=rt[interGenes,]
  }else{
    allTab=cbind(allTab, rt[interGenes,])
  }
  batchType=c(batchType, rep(i,ncol(rt)))
}

#对数据进行矫正，输出矫正后的结果
outTab=ComBat(allTab, batchType, par.prior=TRUE)
outTab=rbind(geneNames=colnames(outTab), outTab)
write.table(outTab, file="merge.txt", sep="\t", quote=F, col.names=F)

