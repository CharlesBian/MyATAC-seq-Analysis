##############################################################################
# TODO: calcuate overlap between enhancers and tf
#	By Jialiang Huang 20150313
###############################################################################

#source("http://bioconductor.org/biocLite.R")
#options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/")
#BiocInstaller::biocLite("Rmisc")


library(RColorBrewer)
library(Rmisc)
library(ggplot2)
library(gridExtra)
library(cowplot)
library(data.table)
library(dplyr)
library(lattice)
library(tidyr)
library(DESeq2)
library(DESeq)

#########################################################
### in/out files
#########################################################
rm(list = ls())
setwd("C:/Users/77214/Desktop/ATAC/DEseq");

getwd()
path_in = "C:/Users/77214/Desktop/ATAC/after_seq/";

#file_deg = '/cluster/huanglab/zyang/wenqing/ExpressionEnrichmentMultipleMapping/Data/20181028_E48/RoadmapRNAseq48_degs_ENID.txt'

path_ins = list.files(path = path_in, pattern = '.*rp10m.bed');
length(path_ins)
print(path_ins)
path_out = "C:/Users/77214/Desktop/ATAC/DEseq/";

#file_out1 = 'A_UP.txt'
#file_out2 = 'A_down.txt'
#file_out3 = 'A_share.txt'

#file_out4 = 'S_UP.txt'
#file_out5 = 'S_down.txt'
#file_out6 = 'S_share.txt'
# check the file number
#a = grep("20181113_E013",path_ins);
#b = grep("20181113_E004",a);
#print(a)
#print(total)


# read gene lists
#creat empty matrix
matrix <- data.frame(matrix(ncol = 3+length(path_ins), nrow = 94526))

total = c("Chr","Start","End");
for (x in 1:length(path_ins)) {
    num1 = gsub('_cvg_rp10m.bed','',path_ins[x]);
    total = append(total, num1);
  }

colnames(matrix) = total



#########################################################
### Reading in data 
#########################################################
i=1
for (i in 1:length(path_ins)) {
  data= read.delim(paste0(path_in,path_ins[i]),skip=0, header=FALSE,sep = "\t", check.names = FALSE);
  num1 = gsub('_cvg_rp10m.bed','',path_ins[i]);
  matrix[,num1]  = data[,9]
  
}
matrix[,1]  = data[,1]
matrix[,2]  = data[,2]
matrix[,3]  = data[,3]

#### write file 
file_out = 'All.txt'
write.table(matrix,file=paste(path_out,file_out,sep=''),sep='\t',quote=F,row.names = FALSE,col.names = TRUE);
################################################################################################################################

######round#####
matrix2 = matrix
matrix[,4:6] = round(matrix[,4:6],0);


######Deseq2 analysis #######

#####contion 1 #######
type <- factor(c(rep("CTRL",1), rep("AOAA",1)))
AOAA = matrix[,c(4,5)]
AOAA <- AOAA %>% select(2,1)
AOAA = AOAA + 1
AOAA <- as.matrix(AOAA)
colData <- data.frame(row.names=colnames(AOAA), type)
cds <- newCountDataSet(AOAA,type)
bbs <-  estimateSizeFactors(cds)
cds <- estimateDispersions(bbs, method='blind' ,sharingMode="fit-only", fitType = 'local')
res <- nbinomTest(cds,"CTRL","AOAA")
res2 = cbind(matrix[,1:3],res)
res2 = res2[,-4]
head(res2)
length(res2)
plotMA(res2)

UP <-subset(res2,log2FoldChange > 1 & pval < 0.05, select=c(1,2,3,8,9,10))#####pval < 0.05 & 
Down <-subset(res2,log2FoldChange < -1 & pval < 0.05,select=c(1,2,3,8,9,10))#####pval < 0.05 & 
DE_ALL <- subset(res2, log2FoldChange > -1 & log2FoldChange < 1 & pval < 0.05,select=c(1,2,3,8,9,10))
write.csv(UP,file=paste0(path_out,"AOAA_up.csv"),sep='\t',quote=F,row.names = FALSE)
write.csv(Down,file=paste0(path_out,'AOAA_down.csv'),sep='\t',quote=F,row.names = FALSE)

write.csv(DE_ALL,file=paste0(path_out,"AOAA_DE_ALL.csv"),sep='\t',quote=F,row.names = FALSE)

pdf("AOAA.pdf")
res2$change <- as.factor(ifelse(res2$pval < 0.05 & abs(res2$log2FoldChange) > 1,ifelse(res2$log2FoldChange > 1,'UP','DOWN'),'NOT')) 

p <- ggplot(res2, aes(x = log2FoldChange, y = -log10(pval), color = change))+geom_point(alpha=0.8, size = 1)+scale_color_manual(values =c("blue","black","red"))+ geom_hline(yintercept = -log10(0.05),lty=4,lwd=0.6,alpha=0.8)+ geom_vline(xintercept = c(1,-1),lty=4,lwd=0.6,alpha=0.8)+xlim(-9,9) + labs(title="Volcanoplot")+ theme_bw(base_size = 18)+theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())
p
dev.off()


#####contion 2 #######
type <- factor(c(rep("CTRL",1), rep("HS",1)))
HS= matrix[,c(5,6)]
HS = HS + 1
HS <- as.matrix(HS)
colData <- data.frame(row.names=colnames(HS), type)
cds <- newCountDataSet(HS,type)
bbs <-  estimateSizeFactors(cds)
cds <- estimateDispersions(bbs, method='blind' ,sharingMode="fit-only", fitType = 'local')
res <- nbinomTest(cds,"CTRL","HS")
res2 = cbind(matrix[,1:3],res)
res2 = res2[,-4]
head(res2)
length(res2)
plotMA(res2)

UP <-subset(res2,log2FoldChange > 1 & pval < 0.05, select=c(1,2,3,8,9,10))#####pval < 0.05 & 
Down <-subset(res2,log2FoldChange < -1 & pval < 0.05,select=c(1,2,3,8,9,10))#####pval < 0.05 & 
DE_ALL <- subset(res2, log2FoldChange > -1 & log2FoldChange < 1 & pval < 0.05,select=c(1,2,3,8,9,10))
write.csv(UP,file=paste0(path_out,"HS_up.csv"),sep='\t',quote=F,row.names = FALSE)
write.csv(Down,file=paste0(path_out,'HS_down.csv'),sep='\t',quote=F,row.names = FALSE)

write.csv(DE_ALL,file=paste0(path_out,"HS_DE_ALL.csv"),sep='\t',quote=F,row.names = FALSE)

pdf("HS.pdf")
res2$change <- as.factor(ifelse(res2$pval < 0.05 & abs(res2$log2FoldChange) > 1,ifelse(res2$log2FoldChange > 1,'UP','DOWN'),'NOT')) 

p <- ggplot(res2, aes(x = log2FoldChange, y = -log10(pval), color = change))+geom_point(alpha=0.8, size = 1)+scale_color_manual(values =c("blue","black","red"))+ geom_hline(yintercept = -log10(0.05),lty=4,lwd=0.6,alpha=0.8)+ geom_vline(xintercept = c(1,-1),lty=4,lwd=0.6,alpha=0.8)+xlim(-9,9) + labs(title="Volcanoplot")+ theme_bw(base_size = 18)+theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())
p
dev.off()


########condition2#######
type <- factor(c(rep("Control",1), rep("shLSD1_2",1)))
shLSD1_2 = matrix[,c(4,6)]
shLSD1_2 = shLSD1_2 + 1
shLSD1_2 <- as.matrix(shLSD1_2)
colData <- data.frame(row.names=colnames(shLSD1_2), type)
cds <- newCountDataSet(shLSD1_2,type)
bbs <-  estimateSizeFactors(cds)
cds <- estimateDispersions(bbs, method='blind' ,sharingMode="fit-only", fitType = 'local')

res <- nbinomTest(cds,"Control","shLSD1_2")
res2 = cbind(matrix[,1:3],res)
res2 = res2[,-4]
head(res2)
length(res2)
plotMA(res2)

UP <-subset(res2,log2FoldChange > 1 & pval < 0.05,select=c(1,2,3,8,9,10))#####pval < 0.05 & 
Down <-subset(res2,log2FoldChange < -1 & pval < 0.05,select=c(1,2,3,8,9,10))#####pval < 0.05 & 
share <- subset(res2, log2FoldChange > -1 & log2FoldChange < 1 & pval < 0.05, select=c(1,2,3,8,9,10))
write.table(UP,file=paste0(path_out,"shLSD1_2_up.bed"),sep='\t',quote=F,row.names = FALSE,col.names = TRUE)
write.table(Down,file=paste0(path_out,'shLSD1_2_down.bed'),sep='\t',quote=F,row.names = FALSE,col.names = TRUE)
write.table(share,file=paste0(path_out,"shLSD1_2_share.bed"),sep='\t',quote=F,row.names = FALSE,col.names = TRUE)
pdf("shLSD1_2.pdf")
res2$change <- as.factor(ifelse(res2$pval < 0.05 & abs(res2$log2FoldChange) > 1,ifelse(res2$log2FoldChange > 1,'UP','DOWN'),'NOT')) 

p <- ggplot(res2, aes(x = log2FoldChange, y = -log10(pval), color = change))+geom_point(alpha=0.8, size = 1)+scale_color_manual(values =c("blue","black","red"))+ geom_hline(yintercept = -log10(0.05),lty=4,lwd=0.6,alpha=0.8)+ geom_vline(xintercept = c(1,-1),lty=4,lwd=0.6,alpha=0.8)+xlim(-9,9) + labs(title="Volcanoplot")+ theme_bw(base_size = 18)+theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())
p
dev.off()
########condition3#######
type <- factor(c(rep("Control",1), rep("shLSD1_E",1)))
shLSD1_E = matrix[,c(4,7)]
shLSD1_E = shLSD1_E + 1
shLSD1_E <- as.matrix(shLSD1_E)
colData <- data.frame(row.names=colnames(shLSD1_E), type)
cds <- newCountDataSet(shLSD1_E,type)
bbs <-  estimateSizeFactors(cds)
cds <- estimateDispersions(bbs, method='blind' ,sharingMode="fit-only", fitType = 'local')
res <- nbinomTest(cds,"Control","shLSD1_E")
res2 = cbind(matrix[,1:3],res)
res2 = res2[,-4]
head(res2)
length(res2)
plotMA(res2)

UP <-subset(res2,log2FoldChange > 1 & pval < 0.05,select=c(1,2,3,8,9,10))#####pval < 0.05 & 
Down <-subset(res2,log2FoldChange < -1 & pval < 0.05,select=c(1,2,3,8,9,10))#####pval < 0.05 & 
share <- subset(res2, log2FoldChange > -1 & log2FoldChange < 1 & pval < 0.05, select=c(1,2,3,8,9,10))
write.table(UP,file=paste0(path_out,"shLSD1_E_up.bed"),sep='\t',quote=F,row.names = FALSE,col.names = TRUE)
write.table(Down,file=paste0(path_out,'shLSD1_E_down.bed'),sep='\t',quote=F,row.names = FALSE,col.names = TRUE)

write.table(share,file=paste0(path_out,"shLSD1_E_share.bed"),sep='\t',quote=F,row.names = FALSE,col.names = TRUE)
pdf("shLSD1_E.pdf")
res2$change <- as.factor(ifelse(res2$pval < 0.05 & abs(res2$log2FoldChange) > 1,ifelse(res2$log2FoldChange > 1,'UP','DOWN'),'NOT')) 

p <- ggplot(res2, aes(x = log2FoldChange, y = -log10(pval), color = change))+geom_point(alpha=0.8, size = 1)+scale_color_manual(values =c("blue","black","red"))+ geom_hline(yintercept = -log10(0.05),lty=4,lwd=0.6,alpha=0.8)+ geom_vline(xintercept = c(1,-1),lty=4,lwd=0.6,alpha=0.8)+xlim(-9,9) + labs(title="Volcanoplot")+ theme_bw(base_size = 18)+theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())
p
dev.off()
###### two replicates ######
condition <- factor(c(rep("ATACseq-EryD",2), rep("ATACseq-EryP",2)))
matrix4 = matrix[,c(4,5,6,7)]
matrix4 = matrix4 + 1
matrix4 <- as.matrix(matrix4)

colData <- data.frame(row.names=colnames(matrix4), condition)

dds <- DESeqDataSetFromMatrix(matrix4, colData, design= ~ condition)
dds <- DESeq(dds)
resultsNames(dds)

res <- results(dds, contrast=c("condition","ATACseq-EryD","ATACseq-EryP"))
res <- res[order(res$pvalue),]
res2 = as.data.frame(res)
res2 = cbind(matrix[,1:3],res2)

head(res2)
summary(res2)


####save all result ####

write.table(res2,file=paste0(path_out2,file_out5),sep='\t',quote=F,row.names = FALSE,col.names = TRUE)

####save diff result ####
diff_gene_deseq2 <-subset(res2,padj < 0.05 & (log2FoldChange > 1 | log2FoldChange < -1))
diff_D_peak <- subset(diff_gene_deseq2, log2FoldChange < -1 ,select=c(1,2,3,5,9))
diff_P_peak <- subset(diff_gene_deseq2, log2FoldChange >  1 ,select=c(1,2,3,5,9))

dim(diff_gene_deseq2)
head(diff_gene_deseq2)

write.table(diff_gene_deseq2,file=paste0(path_out2,file_out6),sep='\t',quote=F,row.names = FALSE,col.names = TRUE)
write.table(diff_D_peak,file=paste0(path_out2,file_out7),sep='\t',quote=F,row.names = FALSE,col.names = FALSE)
write.table(diff_P_peak,file=paste0(path_out2,file_out8),sep='\t',quote=F,row.names = FALSE,col.names = FALSE)