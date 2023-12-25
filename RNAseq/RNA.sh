#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#
########################################


#/sibcb/program/install/hisat2-2.1.0/hisat2-build -p 4 \
# /sibcb1/hurongguilab1/zhumin/YangYi/Ref/Mus_musculus.GRCm39.dna.primary_assembly.fa \
# /sibcb1/hurongguilab1/zhumin/YangYi/Ref/GRCm39;


### Config INFO ####################################################
workdir=/sibcb1/hurongguilab1/wangxuechen/RNAseq;
samps_list='5Y-A-1 5Y-A-2 5Y-A-3';
genome=/sibcb1/hurongguilab1/zhumin/YangYi/Ref/GRCm39;
gtf=/sibcb1/hurongguilab1/zhumin/YangYi/Mus_musculus.GRCm39.108.gtf;
####################################################################


### Clean Data ###
  echo "Workflow 1: fastp";
  mkdir -p ${workdir}/02_CleanData;
  mkdir -p ${workdir}/02_CleanData/fastp;

  for id in $samps_list;
    do
      echo ${id};
      /sibcb1/hurongguilab1/renjin/software/fastp \
                          -i ${workdir}/rawdata/${id}_1.clean.fq.gz \
                          -I ${workdir}/rawdata/${id}_2.clean.fq.gz \
                          -o ${workdir}/02_CleanData/${id}_clean_R1.fq.gz \
                          -O ${workdir}/02_CleanData/${id}_clean_R2.fq.gz \
                          -h ${workdir}/02_CleanData/fastp/$id.html \
                          -j ${workdir}/02_CleanData/fastp/$id.json;
          done
  echo "fastp is done";


### Mapping ###  
  echo "Workflow 2: hisat2";
  mkdir -p ${workdir}/03_hisat2;

  for id in $samps_list;
    do
      echo ${id};
      /sibcb/program/install/hisat2-2.1.0/hisat2 \
        -p 20 --dta -x ${genome} \
        -1 ${workdir}/02_CleanData/${id}_clean_R1.fq.gz \
        -2 ${workdir}/02_CleanData/${id}_clean_R2.fq.gz \
        -S ${workdir}/03_hisat2/${id}.sam;
    done
  echo "hisat2 id done";


### sort Bam ###
  echo "Workflow 3: samtools sort";
  for id in $samps_list;
    do
      echo ${id};
      /sibcb/program/install/samtools-1.4/bin/samtools sort -@ 20 -o ${workdir}/03_hisat2/${samps}.bam ${workdir}/03_hisat2/${samps}.sam;
      /sibcb/program/install/samtools-1.4/bin/samtools view -h -f 2 -F 4 -q 30  ${workdir}/03_hisat2/${samps}.bam > ${workdir}/03_hisat2/${samps}.mapped.bam;
      rm -rf ${workdir}/03_hisat2/${samps}.sam ${workdir}/03_hisat2/${samps}.bam;
    done
  echo "samtools is done";


  echo "Workflow 4: featureCounts";
  /sibcb/program/bin/featureCounts -T 20 -a ${gtf} -o ${workdir}/featurecount.txt ${workdir}/03_hisat2/*.mapped.bam;
  echo "featureCounts is done";

  