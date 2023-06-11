
#install.packages("caret")
#install.packages("DALEX")
#install.packages("ggplot2")
#install.packages("randomForest")
#install.packages("kernlab")
#install.packages("pROC")
#install.packages("xgboost")


#引用包
library(caret)
library(DALEX)
library(ggplot2)
library(randomForest)
library(kernlab)
library(xgboost)
library(pROC)

set.seed(123)      #设置种子
inputFile="test.normalize.txt"         #表达数据文件
geneFile="interGenes.txt"      #基因列表文件
#setwd("C:\\biowolf\\geoCRG\\25.testROC")      #设置工作目录

#读取表达数据文件
data=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
row.names(data)=gsub("-", "_", row.names(data))

#读取基因列表文件, 提取疾病特征基因的表达量
geneRT=read.table(geneFile, header=T, sep="\t", check.names=F)
data=data[as.vector(geneRT[,1]),]

#获取样品分组信息
data=t(data)
group=gsub("(.*)\\_(.*)", "\\2", row.names(data))
data=as.data.frame(data)
data$Type=group

#对数据进行分组
inTrain<-createDataPartition(y=data$Type, p=0.7, list=F)
train<-data[inTrain,]
test<-data[-inTrain,]

#选择模型
control=trainControl(method="repeatedcv", number=8, savePredictions=TRUE)
if(geneFile=="interGenes.txt"){
  #随机森林树模型
  model=train(Type ~ ., data = train, method='rf', trControl = control)
}else if(geneFile=="interGenes.txt"){

  
}

#绘制ROC曲线
yTest=ifelse(test$Type=="Control", 0, 1)
pred1=predict(model, newdata=test, type="prob")
roc1=roc(yTest, as.numeric(pred1[,2]))
ci1=ci.auc(roc1, method="bootstrap")
ciVec=as.numeric(ci1)
pdf(file="ROC.pdf", width=5, height=5)
plot(roc1, print.auc=T, legacy.axes=T, main="", col="red")
text(0.39, 0.43, paste0("95% CI: ",sprintf("%.03f",ciVec[1]),"-",sprintf("%.03f",ciVec[3])), col="red")
dev.off()


######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056
######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

#install.packages("glmnet")
#install.packages("pROC")


#引用包
library(glmnet)
library(pROC)

expFile="test.normalize.txt"      #表达数据文件
geneFile="interGenes.txt"      #交集基因的列表文件
#setwd("C:\\biowolf\\geoFRG\\14.ROC")    #设置工作目录

#读取输入文件
rt=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)

#获取样品的分组信息
y=gsub("(.*)\\_(.*)", "\\2", colnames(rt))
y=ifelse(y=="Control", 0, 1)  #

#读取基因的列表文件
geneRT=read.table(geneFile, header=F, sep="\t", check.names=F)

#对交集基因进行循环，绘制ROC曲线
bioCol=rainbow(nrow(geneRT), s=0.9, v=0.9)    #定义图形的颜色
aucText=c()
k=0
for(x in as.vector(geneRT[,1])){
  k=k+1
  #绘制ROC曲线
  roc1=roc(y, as.numeric(rt[x,]))     #得到ROC曲线的参数
  if(k==1){
    pdf(file="ROC.genes-test.pdf", width=5, height=4.75)
    plot(roc1, print.auc=F, col=bioCol[k], legacy.axes=T, main="")
    aucText=c(aucText, paste0(x,", AUC=",sprintf("%.3f",roc1$auc[1])))
  }else{
    plot(roc1, print.auc=F, col=bioCol[k], legacy.axes=T, main="", add=TRUE)
    aucText=c(aucText, paste0(x,", AUC=",sprintf("%.3f",roc1$auc[1])))
  }
}
#绘制图例，得到ROC曲线下的面积
legend("bottomright", aucText, lwd=2, bty="n", col=bioCol[1:(ncol(rt)-1)])
dev.off()

#构建逻辑模型
rt=rt[as.vector(geneRT[,1]),]
rt=as.data.frame(t(rt))
logit=glm(y ~ ., family=binomial(link='logit'), data=rt)
pred=predict(logit, newx=rt)     #得到模型的打分

