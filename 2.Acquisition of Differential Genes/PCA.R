

#install.packages("ggplot2")


library(ggplot2)        #引用包
#setwd("C:\\Users\\lexb4\\Desktop\\subgroup\\06.PCA")    #设置工作目录

bioPCA=function(inputFile=null, outFile=null){
  #读取输入文件,提取数据
  rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
  data=t(rt)
  Project=gsub("(.*?)\\_.*", "\\1", rownames(data))
  rownames(data)=gsub("(.*?)\\_(.*?)", "\\2", rownames(data))
  
  #PCA分析
  data.pca=prcomp(data)
  pcaPredict=predict(data.pca)
  PCA=data.frame(PC1=pcaPredict[,1], PC2=pcaPredict[,2], Type=Project)
  
  #定义颜色
  bioCol=c("#33FF33","#7700FF","#FF0000","#0066FF","#FF9900","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121")
  bioCol=bioCol[1:length(levels(factor(Project)))]
  
  #绘制PCA图
  pdf(file=outFile, height=5, width=6)       #保存输入出文件
  p=ggplot(data = PCA, aes(PC1, PC2)) + geom_point(aes(color = Type)) +
    scale_colour_manual(name="",  values=bioCol)+
    theme_bw()+
    theme(plot.margin=unit(rep(1.5,4),'lines'))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(p)
  dev.off()
}

#批次矫正前的PCA图
bioPCA(inputFile="merge.preNorm.txt", outFile="PCA.preNorm.pdf")
#批次矫正后的PCA图
bioPCA(inputFile="merge.normalzie.txt", outFile="PCA.normalzie.pdf")


