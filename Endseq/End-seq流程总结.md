# End-seq分析

### 一、数据处理

```shell
#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#
################################################################################################ 
# 测序数据目录
indir=/sibcb1/hurongguilab1/zhumin/ZF/Brca1-P53.Endseq.2/Rawdata;
# 工作目录
workdir=/sibcb1/hurongguilab1/wangxuechen/Endseq;
# 参考基因组路径 ref=/sibcb1/hurongguilab1/zhumin/Ref/GRCm39.hisat2.bulid/mm39.noassembly_patches/mm39;
# 基因组黑名单路径
blacklist=/sibcb1/hurongguilab1/zhumin/Ref/GRCm39.hisat2.bulid/mm10-blacklist.v2.bed;
		 
### Samples List ###############################################################################
samps_list='DKO-1 DKO-2 DKO-3 WT-1  WT-2  WT-3';

# 步骤一 fastp质控
# 提示步骤一开始
echo "<---- Workflow 1 fastp";
# 创建工作目录
mkdir -p $workdir/01_CleanData;
mkdir -p $workdir/01_CleanData/fastp;
			  
for id in $samps_list;
	do
		echo ${id};
		/sibcb1/hurongguilab1/renjin/software/fastp \
		-i ${indir}/$id\_1.fq.gz \
		-I ${indir}/$id\_2.fq.gz \
		-o ${workdir}/01_CleanData/$id\.fastped_1.fastq.gz \
		-O ${workdir}/01_CleanData/$id\.fastped_2.fastq.gz \
		# 生成HTML格式的报告
		-h ${workdir}/01_CleanData/fastp/$id.html \
		# 生成JSON格式的报告
		-j ${workdir}/01_CleanData/fastp/$id.json;
	done
# 提示步骤一完成
echo "fastp is done ---->";

# 步骤二 bowtie2比对
# 提示步骤二开始
echo "<---- workflow bowtie2";
# 创建工作目录
mkdir -p $workdir/02_alignment;
mkdir -p $workdir/03_bam;
for id in $samps_list;
	do
		echo ${id}
		/sibcb/program/install/bowtie2-2.3.1/bowtie2 \
		# 设置线程和参考基因组路径
		-p 24 -x ${ref} \
		-1 ${workdir}/01_CleanData/${id}.fastped_1.fastq.gz \
		-2 ${workdir}/01_CleanData/${id}.fastped_2.fastq.gz |\
		# 对bowtie2的输出进行排序并生成bam
		samtools sort -@ 16 -O bam -o ${workdir}/02_alignment/${id}.bam -;
		# 为生成的bam建立索引
		samtools index ${workdir}/02_alignment/${id}.bam;
		# 生成bam文件的统计信息，并保存
		samtools flagstat ${workdir}/02_alignment/${id}.bam > ${workdir}/02_alignment/${id}.stat;
		# 去重
		/sibcb1/hurongguilab1/zhumin/anaconda3/envs/ATAC/bin/sambamba markdup -r ${workdir}/02_alignment/${id}.bam ${workdir}/02_alignment/${id}.rmdup.bam
		# 过滤线粒体，重新排序
		samtools view -h -f 2 -F 4 -q 30 ${workdir}/02_alignment/${id}.rmdup.bam | grep -v MT | samtools sort -O bam -@ 20 -o - > ${workdir}/03_bam/${id}.rmdup-noChrM.bam
		# 为过滤后的BAM文件建立索引
		samtools index ${workdir}/03_bam/${id}.rmdup-noChrM.bam
	done
# 提示步骤二完成
echo "bowtie2 is done ---->";
																	
# 步骤三 bam转bed
# 提示步骤三开始
echo "<---- bam to bed";
# 创建工作目录
mkdir -p ${workdir}/04_bed;		
for id in $samps_list;
	do
		echo ${id};
		/sibcb/program/bin/bedtools bamtobed -i ${workdir}/03_bam/${id}.rmdup-noChrM.bam > ${workdir}/04_bed/${id}.rmdup-noChrM.bed;
	done
# 提示步骤三结束
echo "bam to bed is done ---->";
																	
# 步骤四 bed转bedgraph转bigwig
fai=/sibcb1/hurongguilab1/zhumin/Ref/GRCm39.hisat2.bulid/mm39.noassembly_patches/Mus_musculus.GRCm39.fa.fai
# 提示步骤四开始
echo "<---- bed --> bedgraph --> bigwig"
# 创建工作目录
mkdir -p $workdir/05_bedgraph-bigwig
for id in $samps_list;
	do
		echo ${id};
		# 通过计算 BED 文件中的行数（reads 数），并将其转换为每百万 reads 数 (RPM)
		rpm=$(wc -l ${workdir}/04_bed/${id}.rmdup-noChrM.bed | cut -f 1 -d ' ' | awk '{printf "%f\n",$1/1000000.0}');
		# 对上一步计算的 RPM 进行归一化处理，将其倒数取倒数
		rpm=$(awk '{printf "%f\n",1/("'$rpm'")}' /sibcb1/hurongguilab1/zhumin/ZF/168-msg_Cold-RT/tst.txt);
		echo ${rpm};																					# 生成bedgraph，并排序					  
		/sibcb/program/bin/bedtools genomecov -bg -scale ${rpm} -i ${workdir}/04_bed/${id}.rmdup-noChrM.bed -g ${fai} \
		| /sibcb/program/bin/bedtools sort -i stdin > ${workdir}/05_bedgraph-bigwig/${id}.bedgraph;
		# 生成bigwig
		bedGraphToBigWig $workdir/05_bedgraph-bigwig/${id}.bedgraph ${fai} $workdir/05_bedgraph-bigwig/${id}.bigwig;
	done
echo "bam --> bedgraph --> bigwig is done ---";

# 步骤五 peakcalling
# 两种选择，用bam或bed文件来做
### 用bam ######################################################################################
#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#
# 设置工作目录
workdir=/sibcb1/hurongguilab1/zhumin/ZF/168-msg_Cold-RT.8-16
# 样本列表		 
samps_list='DKO-1 DKO-2 DKO-3 WT-1 WT-2 WT-3'
# 创建目录，储存峰值文件
mkdir -p ${workdir}/macs2Peaks;

for id in $samps_list;
  do
	echo $id;
	macs2 callpeak -t ${workdir}/03_bam/${id}.rmdup-noChrM.bam \
		# 指定参考基因组，输入为bam，指定输出文件前缀，指定输出文件目录
		-g mm -f BAMPE -n ${id} --outdir ${workdir}/macs2Peaks \
		# 指定 p-value 阈值用于过滤显著性峰值，寻找至少 3 个位点并且每个位点的高度至少是 50
		-q 0.05 -m 3 50 \
		# 不使用默认的建模,在进行峰值调用时对片段进行的偏移,在进行峰值调用时对片段进行的扩展大小
		#--nomodel --shift -100 --extsize 200 \
  done

### 用bed ######################################################################################
#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#

workdir=/sibcb1/hurongguilab1/zhumin/ZF/168-msg_Cold-RT.8-16
		 
samps_list='DKO-1 DKO-2 DKO-3 WT-1 WT-2 WT-3'
mkdir -p ${workdir}/macs2Peaks.bed;

for id in $samps_list;
  do
	echo $id;
	macs2 callpeak -t ${workdir}/04_bed/${id}.rmdup-noChrM.bed \
		-g mm -f BEDPE -n ${id} --outdir ${workdir}/macs2Peaks.bed \
		--nomodel --shift -100 --extsize 200 
  done

# 步骤六 
# 合并所有narrowPeak文件，过滤Y染色体上的峰值
cat *.narrowPeak | grep -v Y > ../../merge.narrowPeak
# 排序
bedtools sort -i  merge.narrowPeak -g \ /home/Method/Reference_Genome/GRCm39.hisat2.bulid/Mus_musculus.GRCm39.dna.primary_assembly.fa.fai > sort.merge.narrowPeak
# 融合
bedtools merge -i sort.merge.narrowPeak > combine.sort.merge.narrowPeak

# 步骤七
# 建立索引
samplelist='DKO-1 DKO-2 DKO-3 WT-1 WT-2 WT-3'
for id in $samplelist;
do 
	echo $id;
	samtools idxstats ../BAM/${id}.rmdup-noChrM.bam > ./BAM.idxstats/${id}.tsv;
done

# 步骤八
# peak counts
samplelist='DKO-1 DKO-2 DKO-3 WT-1 WT-2 WT-3'
mkdir -p peak.counts;
for id in $samplelist;
do
	echo $id;
	bedtools coverage \
	-a combine.sort.merge.narrowPeak \
	-b ../03_bam/${id}.rmdup-noChrM.bam \
	-g /home/Method/Reference_Genome/GRCm39.hisat2.bulid/Mus_musculus.GRCm39.dna.primary_assembly.fa.fai \
	-sorted -counts > ./peak.counts/${id}.count_table;
done
```

### 二、差异分析

* Read merge data

```R
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
```

* Diff

```R
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
```

* Plot

```R
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
```

