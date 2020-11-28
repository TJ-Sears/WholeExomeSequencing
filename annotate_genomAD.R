#!/usr/bin/Rscript

args = commandArgs(trailingOnly=TRUE)

print(args[1:3])

genomAD<-read.table(args[2])
variants<-read.csv(args[1],header = F,sep = "\t")

AFs<-genomAD[,5:8]
AFs<-apply(AFs,1,paste,collapse=" ")

AFs<-substring(AFs, regexpr("AF=", AFs) + 3)
AFs<-substring(AFs,1,11)
AFs[AFs=="N=0 rf_tp_p"]<-"0"
AFs<-as.data.frame(as.numeric(AFs))

AFindex<-cbind(genomAD[,1:4],AFs)
AFcollapsed<-apply(AFindex[1:4],1,paste,collapse="|")
AFcollapsed<-gsub(" ", "", AFcollapsed, fixed = TRUE)
AFcollapsed<-unlist(AFcollapsed)

variantsCollapsed<-apply(variants[c(1,2,4,5)],1,paste,collapse="|")
variantsCollapsed<-gsub(" ", "", variantsCollapsed, fixed = TRUE)
variantsCollapsed<-unlist(variantsCollapsed)

AFVector<-c()
for (i in 1:nrow(variants)){
  if (variantsCollapsed[i]%in%AFcollapsed) {
    index<-which(AFcollapsed%in%variantsCollapsed[i])
    AFVector[i]<-AFindex[index,5]
    } else {
    AFVector[i]<-0
  }
}
AFVector<-as.data.frame(AFVector)
colnames(AFVector)<-"POP_AF"
variants<-cbind(variants,AFVector)

write.csv(variants,args[3])