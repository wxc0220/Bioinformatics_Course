#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#
######################################## 

workdir=/sibcb1/hurongguilab1/zhumin/ZF/168-msg_Cold-RT.8-16
		 
samps_list='168-cold-1 168-cold-2 168-cold-3 168-cold-4 168-RT-1 168-RT-2 168-RT-3 168-RT-4 msg-cold-1 msg-cold-2 msg-cold-3 msg-cold-4 msg-RT-1 msg-RT-2 msg-RT-3 msg-RT-4'
mkdir -p ${workdir}/macs2Peaks.bed;

for id in $samps_list;
  do
	echo $id;
	macs2 callpeak -t ${workdir}/04_bed/${id}.rmdup-noChrM.bed \
		-g mm -f BEDPE -n ${id} --outdir ${workdir}/macs2Peaks.bed \
		--nomodel --shift -100 --extsize 200 
  done

