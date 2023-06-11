

#install.packages("corrplot")


library(corrplot)      #引用包
inputFile="diffGeneExp-cFRG.txt"     #差异基因的表达文件
#setwd("C:\\biowolf\\geoFRG\\08.corrplot")      #设置工作目录

#读取输入文件
data=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)

#删除正常样品
group=gsub("(.*)\\_(.*)", "\\2", colnames(data))
rt=data[,group=="Treat",drop=F]

#得到相关性矩阵
rt=t(rt)
M=cor(rt)
#进行相关性检验, 得到显著性pvalue
res1=cor.mtest(rt, conf.level = .95)

#绘制相关性图形
pdf(file="08.corpot.pdf", width=15, height=15)
corrplot(M,
         type = "upper",       #图形展示在右上方
         method = "circle",    #图形以圆圈的形式展示
         col=colorRampPalette(c('blue', 'white', 'red'),alpha = TRUE)(100), tl.pos="lt",   #颜色的设置
         p.mat=res1$p, insig="label_sig", sig.level = c(.001, .01, .05), pch.cex = 0.85)   #加上显著性标记
corrplot(M, type="lower", add=TRUE, method="number",col=colorRampPalette(c('blue', 'white', 'red'), alpha = TRUE)(100), tl.pos = "n", cl.pos="n", diag=FALSE, number.cex = 0.3) 
dev.off()


