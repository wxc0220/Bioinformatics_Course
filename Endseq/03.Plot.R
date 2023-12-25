library(ggplot2)

setwd("/home/Method/RNF168/End_seq-8.16/")
load("./Endseq.Diff.Rdata")

# 可视化 =======================================================================
  plot(Diff.df$In.RNF168KO.log2FoldChange,-log(Diff.df$In.RNF168KO.padj),main="RNF168 KO")
  plot(Diff.df$In.WT.log2FoldChange,-log(Diff.df$In.WT.padj),main="WT")
  
  RNF168.Diff <- subset(Diff.df,Diff.df$In.RNF168KO.pvalue<0.01)
  hist(RNF168.Diff$In.RNF168KO.log2FoldChange,main="Cold vs Ctrl in RNF168 - log2FC",xlab="",breaks = 50)
  #hist(Diff.df$In.RNF168KO.logFC.CMP,main="Cold vs Ctrl in RNF168 - log2FC",xlab="",breaks = 50)
  
  WT.Diff <- subset(Diff.df,Diff.df$In.WT.pvalue<0.01)
  hist(WT.Diff$In.WT.log2FoldChange,main="Cold vs Ctrl in WT - log2FC",xlab="",breaks = 50)
  #hist(Diff.df$In.WT.logFC.CMP,main="Cold vs Ctrl in WT - log2FC",xlab="",breaks = 50)

  ggplot(Diff.df,aes(x=In.WT.log2FoldChange,y=In.RNF168KO.log2FoldChange))+
    geom_bin2d(bins=50)+ scale_fill_continuous(type = 'viridis')+
    theme_bw()+
    theme_classic()+
    theme(plot.title = element_text(hjust = 0.5))+
    geom_hline(yintercept = 0,color="red",linetype="dashed",size= 0.8) +
    geom_vline(xintercept = 0,color="red",linetype="dashed",size= 0.8)+
    labs(title = "End-seq Cold vs Ctrl",x="Cold vs Ctrl in WT - log2FC",y="Cold vs Ctrl in RNF168 KO - log2FC")+
    #stat_poly_eq(aes(label=paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),formula = y~x,parse = T)+
    theme(axis.title = element_text(size = 12,face = "bold"),
          axis.text = element_text(size = 12,face = "bold"),
          plot.title = element_text (size = 12, face = "bold"))

  
  # WT | RNF168 sig =========================================================================================================
  P.Value_t = 0.05
  Diff.df$sig = ifelse(Diff.df$In.WT.pvalue < P.Value_t & 
                             Diff.df$In.RNF168KO.pvalue < P.Value_t, "all.sig",
                      ifelse(Diff.df$In.WT.pvalue < P.Value_t,"WT.sig",
                      ifelse(Diff.df$In.RNF168KO.pvalue < P.Value_t,"RNF168.sig",
                             "no.sig")))
  
  all.sig.Diff <- subset(Diff.df,Diff.df$sig == "all.sig") 
  WT.sig.Diff <- subset(Diff.df,Diff.df$sig == "WT.sig") 
  RNF168.sig.Diff <- subset(Diff.df,Diff.df$sig == "RNF168.sig")
  RNF168.up <- subset(Diff.df,Diff.df$KOvsWT.padj < 0.01 & Diff.df$KOvsWT.log2FoldChange > 0)
  
  ggplot(Diff.df, aes(x=In.WT.log2FoldChange,y=In.RNF168KO.log2FoldChange)) +
    geom_point(colour="gray",size=2,alpha=0.6,na.rm = T)+
    geom_point(size = 1.6, data = WT.sig.Diff,colour="#00A600FF") +
    geom_point(size = 1.6, data = RNF168.sig.Diff,colour="dodgerblue") +
    geom_point(size = 1.6, data = all.sig.Diff,colour="red") +
    geom_hline(yintercept = 0,color="black",linetype="dashed") +
    geom_point(size = 2, data = RNF168.up,colour="black",shape = 1) +
    geom_vline(xintercept = 0,color="black",linetype="dashed")+
    theme_bw()+
    xlab("Cold vs Ctrl in WT - log2FC") + 
    ylab("Cold vs Ctrl in RNF168 KO - log2FC")+ 
    ggtitle("End-seq Cold vs Ctrl")+
    theme(legend.title = element_blank()) +
    theme(axis.title = element_text(size = 12,face = "bold"),
          axis.text = element_text(size = 12,face = "bold"),
          plot.title = element_text (size = 12, face = "bold"))
  
  write.csv(allDiff,file= "./Endseq.diff.csv")

