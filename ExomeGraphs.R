#This script produces a heatmap of mutation impact from running VEP on mutect2 results

library(pheatmap)
library(splitstackshape)

## read in all files into one object which can be referenced by sample name

setwd("/Users/tsears/code/variants/Variant_11_13/")
files<-list.files(pattern = "results"
                  ,recursive = T,full.names = F,include.dirs = F,)

metadata<-read.csv("../../MouseExomeResults/exome_metadatax.csv")

names<-substr(files,1,5)
all_samples_tabs <- lapply(files,read.table,header=F,fill=T)
names(all_samples_tabs) <- names

#split up VCFs by delimiter

for (i in metadata$Label){
  all_samples_tabs[[i]]<-cSplit(all_samples_tabs[[i]],11,sep = ":")
  all_samples_tabs[[i]]<-cSplit(all_samples_tabs[[i]],8,sep = "|")
}

#Add custom gene list

gene_list<-scan(text="Notch1 Akt1 Ikzf1 Pten Trp53 Phlda3 Rnf213 Ret Elf4 Chd4 Fkbp9 Tet2 
Wwtr1 Vav1 Sh3gl1 Asxl2 Setbp1 Etv5 Smarcb1 Mpl Top1 Ctnnb1 Cntnap2 Muc16 Lrp1b Med12 Ddx3x Crtc1",what="")

#Set NA depth reads to 0 in case of indels
for (i in metadata$Label){
  all_samples_tabs[[i]]$V19[is.na(all_samples_tabs[[i]]$V19)]<-0
  all_samples_tabs[[i]]$V21[is.na(all_samples_tabs[[i]]$V21)]<-0
}


#Tally samples by severity of mutation
df <- data.frame()

for (i in gene_list) {
  counter = 0
  SNVector <- c()
  for (j in (metadata$Label)) {
    counter <- counter + 1
    
    if ((i %in% all_samples_tabs[[j]]$V8_05))   
      {
      index <- which(all_samples_tabs[[j]]$V8_05 == i)
      storage <- all_samples_tabs[[j]][index, ]
      
      if ("MODERATE" %in% storage$V8_04) {
        mod_storage <- storage[which(storage$V8_04 == "MODERATE"),]
        mod_storage <- mod_storage[which.max(mod_storage$V11_03),]
        if ((mod_storage$V11_03 >= 0.05 & mod_storage$V19 < 2 & mod_storage$V21 >= 2) | (nchar(mod_storage$V4)!=nchar(mod_storage$V5)) & mod_storage$V13==-1) {
          SNVector[counter] <- 1
          } else SNVector[counter] <- 0
      } 
      
      if ("HIGH" %in% storage$V8_04) {
        hi_storage <- storage[which(storage$V8_04 == "HIGH"),]
        hi_storage <- hi_storage[which.max(hi_storage$V11_03),]
        if ((hi_storage$V11_03 >= 0.05 & hi_storage$V19 < 2 & hi_storage$V21 >= 2) | (nchar(hi_storage$V4)!=nchar(hi_storage$V5)) & hi_storage$V13==-1) {
          SNVector[counter] <- 2
          } 
      } else if (is.na(SNVector[counter])|is.null(SNVector[counter])) {
        SNVector[counter] <- 0
      }
    } else SNVector[counter] <- 0 
  } 
  df <- rbind(df, SNVector) 
}

#Clean up dataframe by removing genes with <2 SNVs and ordering by frequency. The exeption being Akt1 and Phlda3.
akt_holder<-df["Akt1",]
rownames(df)<-gene_list
colnames(df)<-metadata$Label
df<-df[!rowSums(df)==0,]
df<-df[!rowSums(df)==1,]
df<-df[order(rowSums(df),decreasing = T),]
phlda3<-as.data.frame(c(rep(0,ncol(df))))
colnames(phlda3)<-"Phlda3"
df<-rbind(df,akt_holder,(phlda3$Phlda3))
rownames(df)[length(rownames(df))]<-"Phlda3"

# Time to plot heatmap!
library(pheatmap)
library(RColorBrewer)
library(viridis)

levels<-as.data.frame(metadata$Mouse.genotype)
rownames(levels)<-metadata$Label
colnames(levels)<-"Group"
Gene<-c("Oncogene","Oncogene","Tumor Suppressor Gene","Tumor Suppressor Gene","Tumor Suppressor Gene","Tumor Suppressor Gene")
names(Gene)<-gene_list
Gene<-as.data.frame(Gene)

pheatmap(df,color = c("gray98","grey62","black"),
         legend_breaks = c(0,1,2),legend_labels = c("Low Impact or N/A","Moderate Impact","High Impact")
         ,border_color = "grey78",
         #annotation_row = Gene,
         annotation_names_row = F,
         annotation_legend = T,
         annotation_col = levels,
         annotation_names_col = F,cluster_cols = F,cluster_rows = F
         ,fontsize_col = 10
         ,filename = "COSMIC_V2.pdf"
         ,fontsize_row = 10,width = 12.5)

## Repeating the tallying process for calculating mean AF of a mutation

#Checking average AF for a gene of interest

df<-data.frame()
gene_list<-c("Notch1")
for (i in gene_list){
  
  counter=0
  SNVector<-c()
  
  for (j in (metadata$Label)){
    counter<-counter+1
    
    if ( (i %in% all_samples_tabs[[j]]$V8_05))   {
      index<-which(all_samples_tabs[[j]]$V8_05==i)
      storage<-all_samples_tabs[[j]][index,]
      SNVector[counter]<-max(storage$V11_03)
    } 
    else SNVector[counter]<-0
}
  names(SNVector)<-metadata$Label
}

## Add Sample ID to each sample

#Read in with just tab sep for clarity

setwd("/Users/tsears/code/variants/Variant_11_13/")
files<-list.files(pattern = "results"
                  ,recursive = T,full.names = F,include.dirs = F,)


names<-substr(files,1,5)
all_samples_tabs <- lapply(files,read.table,header=F,fill=T)
names(all_samples_tabs) <- names


for (i in metadata$Label){
  sample_ID<-as.data.frame(c(rep(i,nrow(all_samples_tabs[[i]]))))
  colnames(sample_ID)<-"Sample_ID"
  all_samples_tabs[[i]]<-cbind(all_samples_tabs[[i]],sample_ID)
}

full<-data.frame()
for (i in metadata$Label){
  full<-rbind(full,all_samples_tabs[[i]]) 
}
colnames(full)<-c("CHROM","Position","V3","REF","ALT","V4","Mutect_Filter","Mutect_Info"
                  ,"INFO","INFO2","INFO3","RepeatCHROM","Start","End","Gene","RepeatScore","Strand"
                  ,"Repeat","GermlineReads","GermlineDepth","TumorReads","TumorDepth","Sample_ID")
write.csv(full,"full_variant_list_draft2.csv",row.names = F)





