

#setwd("C:\\Users\\77632\\Desktop\\cuproptosis_ferroptosis\\3_step3")
library(limma)
cur<-read.table("cuproptosis_exp.txt",header = T,sep = "\t",check.names = F)
cur=as.matrix(cur)
rownames(cur)=cur[,1]
curExp=cur[,2:ncol(cur)]
curExp1=matrix(as.numeric(as.matrix(curExp)),nrow=nrow(curExp),dimnames=list(rownames(curExp),colnames(curExp)))
curExp1=avereps(curExp1)
curExp1=curExp1[rowMeans(curExp1)>0,]
#st<- which(substr(colnames(curExp1),14,15) == '11')
#curtumor=curExp1[,-st]

####

immue<-read.table("normalize.txt",header = T,sep = "\t",check.names = F)
immue=as.matrix(immue)
rownames(immue)=immue[,1]
immueExp=immue[,2:ncol(immue)]
immueExp1=matrix(as.numeric(as.matrix(immueExp)),nrow=nrow(immueExp),dimnames=list(rownames(immueExp),colnames(immueExp)))
immueExp1=avereps(immueExp1)
immueExp1=immueExp1[rowMeans(immueExp1)>0,]
#st<- which(substr(colnames(immueExp1),14,15) == '11')
#immuetumor=immueExp1[,-st]

####关注微信公众号生信狂人团队
###遇到代码报错等不懂的问题可以添加微信scikuangren进行答疑
###作者邮箱：sxkrteam@shengxinkuangren.com
ndf=data.frame()
for(a in row.names(immueExp1)){
  for(b in row.names(curExp1)){
    v1=as.numeric(immueExp1[a,])
    v2=as.numeric(curExp1[b,])
    mycor=cor.test(v1,v2)
    cor=mycor$estimate
    p=mycor$p.value
    if((cor>0.5) & (p<0.05)){
      ndf=rbind(ndf,cbind(CRG=b,FRG=a,cor,p,Regulation="postive"))
    }
    if((cor<-0.5) & (p<0.05)){
      ndf=rbind(ndf,cbind(CRG=b,FRG=a,cor,p,Regulation="negative"))
    }
  }
}

write.table(ndf,"CRG_FRG-0.5.txt",sep="\t",quote=F,row.names=F)

type1=data.frame(id=unique(as.vector(ndf[,"FRG"])), Type="FRG")
type2=data.frame(id=unique(as.vector(ndf[,"CRG"])), Type="CRG")
type3=rbind(type1, type2)
write.table(type3, "type.txt", sep="\t", quote=F, row.names=F)

ICGid=unique(as.vector(ndf[,"FRG"]))
ICGEXP=immueExp1[ICGid,]
ICGEXP2=rbind(id=colnames(ICGEXP), ICGEXP)
write.table(ICGEXP2,"FRGEXP-0.5.txt",sep="\t",quote=F,col.names=F)



