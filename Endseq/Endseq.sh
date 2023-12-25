#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#
######################################## 

indir=/sibcb1/hurongguilab1/zhumin/ZF/Brca1-P53.Endseq.2/Rawdata;
workdir=/sibcb1/hurongguilab1/wangxuechen/Endseq;
ref=/sibcb1/hurongguilab1/zhumin/Ref/GRCm39.hisat2.bulid/mm39.noassembly_patches/mm39;
blacklist=/sibcb1/hurongguilab1/zhumin/Ref/GRCm39.hisat2.bulid/mm10-blacklist.v2.bed;
		 
### Samples List ###
samps_list='DKO-1 DKO-2 DKO-3 WT-1  WT-2  WT-3';
		   
echo "<---- Workflow 1 fastp";
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
		-h ${workdir}/01_CleanData/fastp/$id.html \
		-j ${workdir}/01_CleanData/fastp/$id.json;
	done
echo "fastp is done ---->";
											   
echo "<---- workflow bowtie2";
mkdir -p $workdir/02_alignment;
mkdir -p $workdir/03_bam;
for id in $samps_list;
	do
		echo ${id}
		/sibcb/program/install/bowtie2-2.3.1/bowtie2 \
		-p 24 -x ${ref} \
		-1 ${workdir}/01_CleanData/${id}.fastped_1.fastq.gz \
		-2 ${workdir}/01_CleanData/${id}.fastped_2.fastq.gz |\
		samtools sort -@ 16 -O bam -o ${workdir}/02_alignment/${id}.bam -;
		
		samtools index ${workdir}/02_alignment/${id}.bam;
		samtools flagstat ${workdir}/02_alignment/${id}.bam > ${workdir}/02_alignment/${id}.stat;
		/sibcb1/hurongguilab1/zhumin/anaconda3/envs/ATAC/bin/sambamba markdup -r ${workdir}/02_alignment/${id}.bam ${workdir}/02_alignment/${id}.rmdup.bam
		samtools view -h -f 2 -F 4 -q 30 ${workdir}/02_alignment/${id}.rmdup.bam | grep -v MT | samtools sort -O bam -@ 20 -o - > ${workdir}/03_bam/${id}.rmdup-noChrM.bam
		samtools index ${workdir}/03_bam/${id}.rmdup-noChrM.bam
	done
echo "bowtie2 is done ---->";
																				   
																		   
echo "<---- bam to bed";
mkdir -p ${workdir}/04_bed;		
for id in $samps_list;
	do
		echo ${id};
		/sibcb/program/bin/bedtools bamtobed -i ${workdir}/03_bam/${id}.rmdup-noChrM.bam > ${workdir}/04_bed/${id}.rmdup-noChrM.bed;
	done
echo "bam to bed is done ---->";
																													  

fai=/sibcb1/hurongguilab1/zhumin/Ref/GRCm39.hisat2.bulid/mm39.noassembly_patches/Mus_musculus.GRCm39.fa.fai
echo "<---- bed --> bedgraph --> bigwig"
mkdir -p $workdir/05_bedgraph-bigwig
for id in $samps_list;
	do
		echo ${id};
		rpm=$(wc -l ${workdir}/04_bed/${id}.rmdup-noChrM.bed | cut -f 1 -d ' ' | awk '{printf "%f\n",$1/1000000.0}');
		rpm=$(awk '{printf "%f\n",1/("'$rpm'")}' /sibcb1/hurongguilab1/zhumin/ZF/168-msg_Cold-RT/tst.txt);
		echo ${rpm};																																				  
		/sibcb/program/bin/bedtools genomecov -bg -scale ${rpm} -i ${workdir}/04_bed/${id}.rmdup-noChrM.bed -g ${fai} \
		| /sibcb/program/bin/bedtools sort -i stdin > ${workdir}/05_bedgraph-bigwig/${id}.bedgraph;
	
		bedGraphToBigWig $workdir/05_bedgraph-bigwig/${id}.bedgraph ${fai} $workdir/05_bedgraph-bigwig/${id}.bigwig;
	done
echo "bam --> bedgraph --> bigwig is done ---";
