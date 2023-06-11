

library(limma)
library(reshape2)
library(ggplot2)
library(ggpubr)

#主要参数
expFile="ssgseaOut.txt"       #ssgsea结果
typeFile="sample.txt"          #样本属性文件
C="C"                     #正常控制组命名
P="P"                     #疾病实验组名称
Ccol="green"               #正常控制组的颜色
Pcol="red"               #疾病实验组的颜色
afmethod="t.test"        #t检验："t.test"（小样本量）,秩和检验:"wilcox.test"（大样本量）

rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=t(data)

type=read.table(typeFile, sep="\t", header=F, check.names=F, row.names=1)
colnames(type)="data"
sameSample=intersect(row.names(data),row.names(type))
rt1=cbind(data[sameSample,],type[sameSample,])
colnames(rt1)[ncol(rt1)]="data"
rt1=as.data.frame(rt1)
rt1[,1:(ncol(rt1)-1)]=lapply(rt1[,1:(ncol(rt1)-1)],as.numeric)
rt1=melt(rt1,id.vars=c("data"))
colnames(rt1)=c("data","Gene","Expression")
group=levels(factor(rt1$data))
rt1$data=factor(rt1$data, levels=c(C,P))
comp=combn(group,2)
my_comparisons=list()
for(j in 1:ncol(comp)){my_comparisons[[j]]<-comp[,j]}

boxplot=ggboxplot(rt1, x="Gene", y="Expression", fill ="data",
				  xlab="",
				  ylab="Score",
				  legend.title="Type",
				  width=0.8,
				  palette = c(Ccol,Pcol) )+
				  rotate_x_text(50)+
	stat_compare_means(aes(group=data),
	method=afmethod,
	symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 0.2,1), symbols=c("***", "**", "*", "#","ns")), label="p.signif")+
  theme(axis.text=element_text(size=7,face = "bold"))+rotate()

pdf(file="h.all.diff.pdf",  width=7, height=12)
print(boxplot)
dev.off()


