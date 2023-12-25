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
