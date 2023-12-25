library(DESeq2)

#setwd("/home/Method/RNF168/End_seq.6.rep/01.exper1/")
#load("./Endseq.count.CMP.Rdata")


colnames(count_df)
# 设置样本 对照
rownames(samples)=samples$sample
samples$Group <- factor(samples$Group,levels = c("WT","DKO"))   


# PCA #########################################################################
  counts.DEseq=as.matrix(count_df[rownames(samples)])
  samples$Group=factor(samples$Group) 
  dds = DESeqDataSetFromMatrix(countData = counts.DEseq, colData=samples, design = ~ Group)
  dds = DESeq(dds)
  #vsd=vst(dds,blind=F,nsub = sum(rowMeans(counts(dds, normalized=TRUE)) > 5 ))
  vsd=vst(dds,blind=F)

  plotPCA(vsd, intgroup=c('Group'))
  PCAreturn<- plotPCA(vsd,intgroup=c('Group'),returnData = TRUE)
  PCAreturn <- PCAreturn[,c(1,2)]
  PCAreturn <- merge(PCAreturn,samples,by =0)
  ggplot(PCAreturn,aes(PC1,PC2,color=Group)) +
    geom_point(size=3) +
    labs(x="PC1",y="PC2",face='bold') +
    ggrepel::geom_text_repel(aes(label = Row.names),color="black",size =3)+
    scale_color_manual(values=c("#CC0000","#2F5597")) +
    theme_light(base_size = 16)


# 差异分析 #####################################################################
resultsNames(dds)

KOvsWT <- results(dds, name="Group_DKO_vs_WT")
summary(KOvsWT)
head(KOvsWT)
KOvsWT=as.data.frame(KOvsWT)
KOvsWT$padj[is.na(KOvsWT$padj)] = 1

plot(KOvsWT$log2FoldChange,-log(KOvsWT$pvalue))
abline(h=2,lty=2)
abline(v=1,lty=2)
abline(v=-1,lty=2)
head(KOvsWT)
write.csv2(KOvsWT,file = "./diffEndseq.csv")

