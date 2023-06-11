

#install.packages("corrplot")
#install.packages("circlize")


#引用包
library(corrplot)
library(circlize)

inputFile="diffGeneExp-cFRG.txt"    #输入文件
#setwd("C:\\biowolf\\geoCRG\\10.cor")     #设置工作目录

#读取输入文件
data=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)

#去除对照组样品
group=gsub("(.*)\\_(.*)", "\\2", colnames(data))
data=data[,group=="Treat",drop=F]
rt=t(data)

#计算基因间相关系数
cor1=cor(rt)

#设置图形颜色
col = c(rgb(1,0,0,seq(1,0,length=32)),rgb(0,1,0,seq(0,1,length=32)))
cor1[cor1==1]=0
c1 = ifelse(c(cor1)>=0,rgb(1,0,0,abs(cor1)),rgb(0,1,0,abs(cor1)))
col1 = matrix(c1,nc=ncol(rt))

#绘制圈图
pdf(file="08.circos.pdf", width=10, height=10)
par(mar=c(2,2,2,4))
circos.par(gap.degree=c(3,rep(2, nrow(cor1)-1)), start.degree = 180)
chordDiagram(cor1, grid.col=rainbow(ncol(rt)), col=col1, transparency = 0.5, symmetric = T)
par(xpd=T)
#绘制图例
colorlegend(col, vertical = T,labels=c(1,0,-1),xlim=c(1.1,1.3),ylim=c(-0.4,0.4))
dev.off()
circos.clear()

#绘制相关性图形
pdf(file="08.corrplot.pdf", width=12, height=12)
corrplot(cor1,
         method = "pie",
         order = "hclust",
         type = "upper",
         tl.cex=1,
         col=colorRampPalette(c("green", "white", "red"))(50)
)
dev.off()


