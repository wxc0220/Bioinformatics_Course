library(tidyverse)
setwd("/home/WXC/Endseq/R_pipeline/peak.counts/")

DKO.1 <- read_tsv(file = "./DKO-1.count_table",col_names = F)
#DKO.2 <- read_tsv(file = "./DKO-2.count_table",col_names = F)
#DKO.3 <- read_tsv(file = "./DKO-3.count_table",col_names = F)
DKO.4 <- read_tsv(file = "./DKO-4.count_table",col_names = F)
DKO.5 <- read_tsv(file = "./DKO-5.count_table",col_names = F)
colnames(DKO.1) <- c("chrom","start","end","DKO.1")
#colnames(DKO.2) <- c("chrom","start","end","DKO.2")
#colnames(DKO.3) <- c("chrom","start","end","DKO.3")
colnames(DKO.4) <- c("chrom","start","end","DKO.2")
colnames(DKO.5) <- c("chrom","start","end","DKO.3")

WT.1 <- read_tsv(file = "./WT-1.count_table",col_names = F)
WT.2 <- read_tsv(file = "./WT-2.count_table",col_names = F)
#WT.3 <- read_tsv(file = "./WT-3.count_table",col_names = F)
#WT.4 <- read_tsv(file = "./WT-4.count_table",col_names = F)
WT.5 <- read_tsv(file = "./WT-5.count_table",col_names = F)
colnames(WT.1) <- c("chrom","start","end","WT.1")
colnames(WT.2) <- c("chrom","start","end","WT.2")
#colnames(WT.3) <- c("chrom","start","end","WT.3")
#colnames(WT.4) <- c("chrom","start","end","WT.4")
colnames(WT.5) <- c("chrom","start","end","WT.3")

count_df <- DKO.1
#count_df <- cbind(count_df,DKO.2[,4])
#count_df <- cbind(count_df,DKO.3[,4])
count_df <- cbind(count_df,DKO.4[,4])
count_df <- cbind(count_df,DKO.5[,4])

count_df <- cbind(count_df,WT.1[,4])
count_df <- cbind(count_df,WT.2[,4])
#count_df <- cbind(count_df,WT.3[,4])
#count_df <- cbind(count_df,WT.4[,4])
count_df <- cbind(count_df,WT.5[,4])


row.names(count_df) <- paste0(count_df$chrom, ":", count_df$start, "-", count_df$end)

# 差异分析 & PCA
colnames(count_df)
samples=data.frame(
  sample=c("DKO.1","DKO.2","DKO.3",
           "WT.1","WT.2","WT.3"),
  Group = c("DKO","DKO","DKO",
            "WT","WT","WT"))
rownames(samples)=samples$sample
samples   


# 将 featurecount 转化为 Matrix
counts.DEseq=as.matrix(count_df[rownames(samples)])

# DEseq读取
library(DESeq2)
samples$Group=factor(samples$Group) 
dds = DESeqDataSetFromMatrix(countData = counts.DEseq, colData=samples, design = ~ Group)
# 差异分析
dds = DESeq(dds)
### PCA 查看样本间关系 ###
vsd=vst(dds,blind=F,nsub = sum(rowMeans(counts(dds, normalized=TRUE)) > 5 ))

plotPCA(vsd, intgroup=c('Group'))
PCAreturn<- plotPCA(vsd,intgroup=c('Group'),returnData = TRUE)

PCAreturn <- PCAreturn[,c(1,2)]
PCAreturn <- merge(PCAreturn,samples,by =0)
ggplot(PCAreturn) +
  geom_point(aes(PC1,PC2,color=Group),size=3) +
  labs(x="PC1",y="PC2",face='bold') +
  ggrepel::geom_label_repel(aes(PC1,PC2,label = Row.names),color="black",size =2.5)+
  scale_color_manual(values=c("#CC0000","#2F5597")) +
  theme_light(base_size = 16)


# CMP 转换 #####################################################################
setwd("/home/Method/zhufang/Brca1-P53.Endseq/bed.quality_control/")
merge_df <- NULL
case_info.vec = c("DKO-1","DKO-2","DKO-3",
                  "WT-1","WT-2","WT-3")

for(case in case_info.vec){
  print(case)
  total_df_filename <- sprintf("./BAM.idxstats/%s.tsv",case)  
  total_df <- read_tsv(file = total_df_filename,col_names = F)
  colnames(total_df) <- c("chrom","length","map_read_count","unmap_read_count")
  total_df$case_info = case
  merge_df<-bind_rows(merge_df,total_df)
}
merge_df

colnames(count_df[,-c(1:3)])
table(merge_df$case_info)
merge_df$case_info <- factor(merge_df$case_info)

total_count.vec = filter(merge_df,chrom != "Y",chrom != "*")%>%
  group_by(case_info)%>%
  summarise(total_count = sum(map_read_count))%>%
  pull(total_count)


colnames(count_df[,-c(1:3)])
CMP_df <- as.matrix(count_df[,-c(1:3)]) 
boxplot(CMP_df, outline = F,notch = T ,las = 2,
        col=c(rep('#CC0000',3),rep('#2F5597',3)))

CMP_df =  (CMP_df / total_count.vec) * 1e6
boxplot(CMP_df, outline = F,notch = T ,las = 2,
        col=c(rep('#CC0000',3),rep('#2F5597',3)))

CMP_df <- as.data.frame(CMP_df)

save(count_df, CMP_df, file = "Endseq.count.CMP.Rdata")
