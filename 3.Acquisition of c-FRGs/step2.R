

#setwd("C:\\Users\\77632\\Desktop\\cuproptosis_ferroptosis\\1_step1")
cupr=read.table("ferroptosis.txt",header = T,sep = "\t",check.names = F)
mrna=read.table("normalize.txt",header = T,sep = "\t",check.names=F)
cuproptosis_gene=merge(cupr,mrna,by="id")
write.table(cuproptosis_gene,"ferroptosis_exp.txt",quote = F,row.names = F,sep = "\t")



