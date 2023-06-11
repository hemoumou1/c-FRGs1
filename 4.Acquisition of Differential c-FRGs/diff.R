

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("pheatmap")
#install.packages("reshape2")
#install.packages("ggpubr")


#引用包
library(limma)
library(pheatmap)
library(reshape2)
library(ggpubr)

expFile="FRGEXP-0.5.txt"      #表达数据文件
#setwd("C:\\biowolf\\geoCRG\\07.diff")      #设置工作目录

#读取表达数据文件
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
exp=data

#提取样品的分组信息
Type=gsub("(.*)\\_(.*)", "\\2", colnames(data))

#基因差异分析
sigVec=c()
sigGeneVec=c()
for(i in row.names(data)){
  test=wilcox.test(data[i,] ~ Type)
  pvalue=test$p.value
  Sig=ifelse(pvalue<0.001,"***",ifelse(pvalue<0.01,"**",ifelse(pvalue<0.05,"*","")))
  if(pvalue<0.05){
    sigVec=c(sigVec, paste0(i, Sig))
    sigGeneVec=c(sigGeneVec, i)}
}
#输出差异基因的表达量
data=data[sigGeneVec,]
outTab=rbind(ID=colnames(data), data)
write.table(outTab, file="diffGeneExp-cFRGs.txt", sep="\t", quote=F, col.names=F)
row.names(data)=sigVec

#对差异基因进行可视化，绘制热图
names(Type)=colnames(data)
Type=as.data.frame(Type)
pdf(file="07.heatmap.pdf", width=10, height=6)
pheatmap(data,
         annotation=Type,
         color = colorRampPalette(c(rep("blue",2), "white", rep("red",2)))(100),
         cluster_cols =F,
         cluster_rows =T,
         scale="row",
         show_colnames=F,
         show_rownames=T,
         fontsize=7,
         fontsize_row=7,
         fontsize_col=7)
dev.off()

#把表达数据转换成ggplot2输入文件
exp=as.data.frame(t(exp))
exp=cbind(exp, Type=Type)
data=melt(exp, id.vars=c("Type"))
colnames(data)=c("Type", "Gene", "Expression")

#绘制箱线图
p=ggboxplot(data, x="Gene", y="Expression", color = "Type", 
            xlab="",
            ylab="Gene expression",
            legend.title="Type",
            palette = c("blue", "red"),
            add="point",
            width=0.8)
p=p+rotate_x_text(60)
p1=p+stat_compare_means(aes(group=Type),
                        method="wilcox.test",
                        symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                        label = "p.signif")

#输出箱线图
pdf(file="07.boxplot.pdf", width=20, height=10)
print(p1)
dev.off()


