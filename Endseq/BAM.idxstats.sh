samplelist='DKO-1 DKO-2 DKO-3 WT-1 WT-2 WT-3'

for id in $samplelist;
do 
	echo $id;
	samtools idxstats ../BAM/${id}.rmdup-noChrM.bam > ./BAM.idxstats/${id}.tsv;
done

