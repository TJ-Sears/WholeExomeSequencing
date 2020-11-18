# This script produces a heatmap for CNVkit results

setwd("/Users/tsears/code/variants/CNV_11_17/")
files<-list.files(pattern = "genemetrics"
                  ,recursive = T,full.names = F,include.dirs = F,)

metadata<-read.csv("../../MouseExomeResults/exome_metadatax.csv")

names<-paste("x",substr(files,1,4),sep = "")
all_samples_tabs <- lapply(files,read.table,header=T,fill=T)
names(all_samples_tabs) <- names

gene_list<-scan(text="Notch1 Akt1 Ikzf1 Pten Trp53 Phlda3 Rnf213 Ret Elf4 Chd4 Fkbp9 Tet2 
Wwtr1 Vav1 Sh3gl1 Asxl2 Setbp1 Etv5 Smarcb1 Mpl Top1 Ctnnb1 Cntnap2 Muc16 Lrp1b Med12 Ddx3x Crtc1",what="")

df <- data.frame()

for (i in gene_list) {
  counter = 0
  SNVector <- c()
  for (j in (metadata$Label)) {
    counter <- counter + 1
    
    if ((i %in% all_samples_tabs[[j]]$gene))   
    {
      index <- which(all_samples_tabs[[j]]$gene == i)
      storage <- all_samples_tabs[[j]][index, ]
      SNVector[counter] <- storage$log2[1]
    }
    else SNVector[counter] <- 0
  }
  df <- rbind(df, SNVector)
}

rownames(df)<-gene_list
colnames(df)<-metadata$Label

library(pheatmap)
library(RColorBrewer)
library(viridis)

levels<-as.data.frame(metadata$Mouse.genotype)
rownames(levels)<-metadata$Label
colnames(levels)<-"Group"
#Gene<-c("Oncogene","Oncogene","Tumor Suppressor Gene","Tumor Suppressor Gene","Tumor Suppressor Gene","Tumor Suppressor Gene")
#names(Gene)<-gene_list
#Gene<-as.data.frame(Gene)

cols <- brewer.pal(n = 8,name = "RdBu")
cols <- rev(cols)
cols<- cols[1:4]
cols<- c(cols,"gray98","darksalmon")

pheatmap(df
         ,color = cols
         ,breaks = c(-5,-3,-2,-1,-0.01,0,1)
         ,border_color = "grey78",
         main = "CNV Cosmic Genes Log2 Fold Change",
         annotation_names_row = F,
         annotation_legend = T,
         annotation_col = levels,
         annotation_names_col = F,cluster_cols = F,cluster_rows = F
         ,fontsize_col = 10
         ,filename = "COSMIC_CNV.pdf"
         ,fontsize_row = 10,width = 10)

