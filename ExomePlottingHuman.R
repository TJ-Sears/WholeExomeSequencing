library(pheatmap)
library(splitstackshape)
library(RColorBrewer)
library(viridis)

setwd("/Users/tsears/code/variants/hiseq_consensus_11_25/")
files<-list.files(pattern = "results"
                  ,recursive = T,full.names = F,include.dirs = F,)


names<-substr(files,1,5)
all_samples <- lapply(files,read.table,header=T,fill=T,sep=",")
names(all_samples) <- names

for (i in names){
  all_samples[[i]]<-cSplit(all_samples[[i]],25,sep = " ")
}

for (i in names){
  AF<-all_samples[[i]]$V27/all_samples[[i]]$V28
  all_samples[[i]]<-cbind(all_samples[[i]],AF)
}

gene_list<-COSMIC #can be any list of genes, I used a custom curated list of COSMIC genes

#or, make long fancy gene list to cover every possible gene
gene_list<-c()
for (i in names){
  gene_list<-c(gene_list,all_samples[[i]]$V7)
}

gene_list<-gene_list[(!duplicated(gene_list))]


#Tally samples by severity of mutation

df <- data.frame()
afs <- read.csv(text="AF,gene,sampleName")
af_counter<-0

for (j in gene_list){
  counter = 0
  SNVector <- c()
  for (i in names){
    counter=counter+1
    
    if (j %in% all_samples[[i]]$V7){
      index <- which(all_samples[[i]]$V7 == j)
      storage <- all_samples[[i]][index,]
      storage <- storage[which.max(storage$V24_09),]
      if (storage$V26>=30 & storage$V28>=30 & storage$AF>=0.05 & storage$POP_AF < 0.01 & storage$V30 == -1 
          & storage$V6=="exonic" & storage$V25==0 & storage$V9!="."){
        
        #Store AF data for later plotting
        af_counter=af_counter+1
        afs[af_counter,"AF"]<-storage$AF
        afs[af_counter,"gene"]<-storage$V7
        afs[af_counter,"sampleName"]<-i
        
        if (storage$V9=="synonymous SNV"){
          SNVector[counter]<-1
        }
        if (storage$V9=="nonsynonymous SNV"){
          SNVector[counter]<-2
        }
        if (storage$V9=="stopgain"){
          SNVector[counter]<-3
        } 
      } else {
        SNVector[counter]<-0
      }
    } else {
    SNVector[counter]<-0
  }
  }
  df<-rbind(df,SNVector)
}

rownames(df)<-gene_list
colnames(df)<-names
df<-df[!rowSums(df)==0,]
df<-df[!rowSums(df)==1,]
df<-df[order(rowSums(df),decreasing = T),]
df<-df[1:40,]

cols<-c(brewer.pal(4,name = "Dark2"))
cols[1]<-"grey90"

pheatmap(df,color = cols,
         legend_breaks = c(0,1,2,3),legend_labels = c("No Mutation or No Description","Synonymous","Nonsynonymous","Stopgain")
         ,border_color = "grey76",
         annotation_names_row = F,
         annotation_legend = T,
         annotation_names_col = F,cluster_cols = F,cluster_rows = F
         ,fontsize_col = 10
         ,filename = "HumanTestLong.pdf"
         ,fontsize_row = 1,width = 10)


## Plot histograms of AFs for each patient
#use afs from earlier

cols<-c(brewer.pal(7,name="Dark2"))

library(ggplot2)
ggplot(afs, aes(x=sampleName,y=AF, colour = sampleName)) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.7)) 
  + scale_color_manual(values=cols)







