## Project: β-diketone accumulation in response to drought stress is weakened in modern bread wheat varieties (Triticum aestivum L.)
## Purpose: R codes for transcriptome downstream analysis
## Created By: Peng Gao, Aswini Kuruparan, Eliana Gonzales-Vigil,Raju Soolanayakanahally
## Created Date: 2024.05.30

##  'All output from Salmon has been transferred to folder ./Salmon under the current folder'

setwd("./")
#packages installation
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("tximport")
BiocManager::install("DESeq2")

#Importing transcript abundance with tximport
library(tximport)

#SampleInfo.csv including salmon quansi folder names in "Name" column
samples = read.csv("SampleInfo.csv")

#tx2gene_wheat.csv including 2 column, 1st are transcript names, 2nd are corresponding gene names, information from Biomart wheat database
tx2gene = read.csv("tx2gene_wheat.csv")

#list of Salmon output
files <- file.path("./Salmon", samples$Name, "quant.sf")

#change sample name from sequencing name
names(files) =paste0(samples$Sample)

#check file existence
all(file.exists(files)) 

#generate counts of each sample at gene level
txi = tximport(files, type = "salmon", tx2gene =tx2gene)

#transform tximport output to DESeq2
library(DESeq2)
colData=read.csv("SampleInfo.csv")
dds <- DESeqDataSetFromTximport(txi, colData = colData, design = ~ Cultivar)

#remove low-expressed genes to get the list of expressed genes
filter <- rowSums(counts(dds, normalized=T) >= 10) >= 3
dds <- dds[filter,]

#DESeq2 cooking
dds <- DESeq(dds)

#output normalized counts and raw counts of expressed genes
write.csv(counts(dds,normalized=T),"Normalized_counts_expressed_genes.csv")
write.csv(counts(dds),"Raw_counts_expressed_genes.csv")

#PCA
vsd <- vst(dds,blind=F)
png("PCA.png", height=100, width=130, units='mm', res = 300)
plotPCA(vsd,"Cultivar")
dev.off()

#DEG analysis
res1 <- results(dds, contrast =c("Cultivar", "Magnet", "Concord")) 
res2 <- results(dds, contrast =c("Cultivar", "Starbuck", "Concord")) 
res3 <- results(dds, contrast =c("Cultivar", "Tradition", "Concord"))  
res4 <- results(dds, contrast =c("Cultivar", "Starbuck", "Magnet")) 
res5 <- results(dds, contrast =c("Cultivar", "Tradition", "Magnet"))
res6 <- results(dds, contrast =c("Cultivar", "Tradition", "Starbuck"))
write.csv(res1,"DEG_Magnet_vs_Concord.csv")
write.csv(res2,"DEG_Starbuck_vs_Concord.csv")
write.csv(res3,"DEG_Tradition_vs_Concord.csv")
write.csv(res4,"DEG_Starbuck_vs_Magnet.csv")
write.csv(res5,"DEG_Tradition_vs_Magnet.csv")
write.csv(res6,"DEG_Tradition_vs_Starbuck.csv")
