# 导入相关 R包
  
  # 安装 DESeq2 
    # install.packages("BiocManager")
    # library(BiocManager)
    # BiocManager::install("DESeq2")  
  # 导入  
    library(DESeq2)
    library(ggplot2)
  
  # GENEID 转化包
    # BiocManager::install("clusterProfiler")
    library(clusterProfiler)
    # BiocManager::install("org.Hs.eg.db")
    # BiocManager::install("org.Mm.eg.db")
    library(org.Hs.eg.db)
    # library(org.Mm.eg.db)
    
    # install.packages("ggrepel")
    library(ggrepel)

# 设置工作目录
  setwd('/home/biotrainee.workspace/RNAseq/')
  getwd()

# 读取数据 #####################################################################
  counts=read.csv('./featurecount.txt',
                  sep = '\t',
                  row.names = 'Geneid',
                  comment.char = '#',
                  header = T,
                  check.names = F)

# 删除前5列
  counts = counts[,-c(1:4)]
# 保留行相加 大于 10的行
  counts=counts[rowSums(counts[,-1])>20,]
  counts


# 样本对照关系 抽出 VDR.DBD.DEL ==========================================================================
  colnames(counts)
  samples <- data.frame(Sample = c("CA-CRT-TLS",
                                   "CA-high-TLS-1","CA-high-TLS-2","CA-high-TLS-3","CA-high-TLS-4",
                                   "CA-low-TLS-1","CA-low-TLS-2","CA-low-TLS-3","CA-low-TLS-4","CA-low-TLS-5" ,
                                   
                                   "CRT-TLS",
                                   "high-TLS-1","high-TLS-2","high-TLS-3","high-TLS-4",
                                   "low-TLS-1","low-TLS-2","low-TLS-3","low-TLS-4"),
                        
                        Group = c(rep("Cancer",10),
                                  rep("Normal",9)))
  samples$Sample
  row.names(samples) <- samples$Sample
  samples$Group <- factor(samples$Group,levels = c("Normal","Cancer"))
  samples$Group
  samples

#samples <- subset(samples,samples$Group == "Cancer")

# 将 featurecount 转化为 Matrix ================================================
  counts.DEseq <- as.matrix(counts[rownames(samples)])
  dim(counts.DEseq);class(counts.DEseq)

  
# DEseq读取
  dds = DESeqDataSetFromMatrix(countData = counts.DEseq, colData=samples, design = ~ Group)
# 差异分析
  dds = DESeq(dds)
  
  
### PCA 查看样本间关系 ###
  vsd=vst(dds,blind=F)

  plotPCA(vsd, intgroup=c('Group'))

  PCAreturn<- plotPCA(vsd,intgroup=c('Group'),returnData = TRUE)

  ggplot(PCAreturn,aes(PC1,PC2,color=group)) +
    geom_point(size=3) +
    ggsci::scale_color_igv()+
    ggrepel::geom_label_repel(aes(label = name),
                              color="black",size =2.5)+
    labs(x="PC1",y="PC2",face='bold') +
    theme_light(base_size = 16)
  rm(PCAreturn)


################################################################################
###                           DEseq 差异分析                                 ###
################################################################################
  resultsNames(dds)
  Cancer_vs_Normal <- results(dds, name="Group_Cancer_vs_Normal")
  summary(Cancer_vs_Normal)
  Cancer_vs_Normal=as.data.frame(Cancer_vs_Normal)
  
  counts.DEseq <- as.data.frame(counts.DEseq)
  Cancer_vs_Normal <- merge(Cancer_vs_Normal,counts.DEseq,by = 0,Row.names=0)
  Cancer_vs_Normal$padj[is.na(Cancer_vs_Normal$padj)] = 1
  rownames(Cancer_vs_Normal) <- Cancer_vs_Normal$Row.names
  Cancer_vs_Normal <- Cancer_vs_Normal[,-1]

# GeneID -> GeneSymbol 
  Ensembl_ID <- row.names(Cancer_vs_Normal)
  gene_symbol <- bitr(Ensembl_ID, fromType="ENSEMBL", toType=c("SYMBOL", "ENTREZID"), OrgDb="org.Hs.eg.db")
  head(gene_symbol)
  Cancer_vs_Normal <- merge(gene_symbol,Cancer_vs_Normal,by.x="ENSEMBL" ,by.y=0)

# 查看火山图
  plot(Cancer_vs_Normal$log2FoldChange,-log10(Cancer_vs_Normal$padj))
  write.csv(Cancer_vs_Normal,file = 'Cancer_vs_Normal.csv')


# 可视化

  df <- Cancer_vs_Normal
  plot(df$log2FoldChange,-log(df$padj))
  
  P.Value_t = 0.05
  FC = 2
  
  up_genes <- subset(df,df$pvalue < P.Value_t & df$log2FoldChange > FC)
  down_genes <- subset(df,df$pvalue < P.Value_t & df$log2FoldChange < -FC)
  sig_genes <- rbind(up_genes,down_genes)
  
  df$gene_type = ifelse(df$log2FoldChange < -1 & df$pvalue < P.Value_t, 
                        "down",ifelse(df$log2FoldChange > 1 & df$pvalue < P.Value_t, 
                                      "up","ns"))
  
  
  P <- ggplot(df,aes(x = log2FoldChange , y = -log10(pvalue))) +
    geom_point(aes(color = gene_type), alpha = 0.6, shape = 16, size = 1.8) +
    geom_point(data = up_genes, shape = 21, size = 2, fill = "red", colour = "black") +
    geom_point(data = down_genes, shape = 21, size = 2, fill = "steelblue", colour = "black") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_vline(xintercept = c(-1,1), linetype = "dashed") +
    geom_label_repel(data = sig_genes, aes(label = SYMBOL), force = 2, nudge_y = 1)+
    scale_color_manual(values = c("up" = "#ffad73", "down" = "#26b3ff", "ns" = "grey")) +
    scale_x_continuous(breaks = c(seq(-8, 8, 2)), limits = c(-8, 8)) +
    labs(x = "log2(fold change)", y = "-log10(adjusted P-value)", colour = "Expression change") +
    guides(color = guide_legend(override.aes = list(size = 5))) +
    theme_bw() + 
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          plot.background = element_blank(),
          axis.title = element_text(face = "bold", color = "black", size = 10),
          axis.text = element_text(color = "black", size = 9, face = "bold"),
          legend.background = element_blank(),
          legend.title = element_text(face = "bold", color = "black", size = 10),
          legend.text = element_text(face = "bold", color = "black", size = 9),
          legend.spacing.x = unit(0, "cm"),
          legend.position = c(0.12, 0.9))
  P
  ggsave(P,filename = "Cancer_vs_Normal.png",height = 6,width = 7,dpi = 300)
  