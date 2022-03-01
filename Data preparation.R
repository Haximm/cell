###################################0.data preparation###################################
library(limma)               
ReadFile="rowread.txt"        
NSFile="nsnames.txt"             
AWFile="awnames.txt"          
PUFile="punames.txt"           

#Read the rowread file and reorganize
rt=read.table(ReadFile,sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
#Read sample information
sample1=read.table(NSFile,sep="\t",header=F,check.names=F)
sample2=read.table(AWFile,sep="\t",header=F,check.names=F)
sample3=read.table(PUFile,sep="\t",header=F,check.names=F)
sampleName1=gsub("^ | $", "", as.vector(sample1[,1]))
sampleName2=gsub("^ | $", "", as.vector(sample2[,1]))
sampleName3=gsub("^ | $", "", as.vector(sample3[,1]))
NSData=data[,sampleName1]
AWData=data[,sampleName2]
PUData=data[,sampleName3]
#output NS samples expression levels
NSData=rbind(id=colnames(NSData), NSData)
write.table(NSData, file="NS.txt",sep="\t",quote=F,col.names=F)
#output AW,PU samples expression levels
AWData=rbind(id=colnames(AWData), AWData)
write.table(AWData, file="AW.txt",sep="\t",quote=F,col.names=F)
PUData=rbind(id=colnames(PUData), PUData)
write.table(PUData, file="PU.txt",sep="\t",quote=F,col.names=F)
