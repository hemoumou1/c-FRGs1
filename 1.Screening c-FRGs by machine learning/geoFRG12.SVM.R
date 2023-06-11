

#install.packages("e1071")


#引用包
#set.seed(12345)
library(e1071)

inputFile="diffGeneExp.txt"     #输入文件
#setwd("C:\\biowolf\\geoFRG\\12.SVM")      #设置工作目录
source("geoFRG12.msvmRFE.R")

#读取输入文件
data=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
data=as.data.frame(t(data))
#获取样品的分组信息
group=gsub("(.*)\\_(.*)", "\\2", row.names(data))
data=cbind(group, data)
data$group=factor(data$group, levels=c("Control","Treat"))

#构建机器学习-支持向量机递归特征消除算法(SVM-RFE)
svmRFE(data, k=10, halve.above=10)
nfold=10
geneNum=nrow(data)
folds=rep(1:nfold, len=geneNum)[sample(geneNum)]
folds=lapply(1:nfold, function(x) which(folds == x))
results=lapply(folds, svmRFE.wrap, data, k=10, halve.above=10)

#对特征基因的重要性进行排序
top.features=WriteFeatures(results, data, save=F)
#输出排序的结果
write.table(top.features, file="feature_svm.txt", sep="\t", quote=F,row.names=F)

#交叉验证
featsweep=lapply(1:40, FeatSweep.wrap, results, data)

#获取交叉验证的误差
no.info=min(prop.table(table(data[,1])))
errors=sapply(featsweep, function(x) ifelse(is.null(x), NA, x$error))

#绘制交叉验证误差的图形
pdf(file="errors.pdf", width=5, height=5)
PlotErrors(errors, no.info=no.info)
dev.off()

#绘制交叉验证准确性的图形
pdf(file="accuracy.pdf", width=5, height=5)
Plotaccuracy(1-errors, no.info=no.info)
dev.off()

#输出SVM的特征基因
featureGenes=top.features[1:which.min(errors),1,drop=F]
write.table(file="SVM-RFE.gene.txt", featureGenes, sep="\t", quote=F, row.names=F, col.names=F)



