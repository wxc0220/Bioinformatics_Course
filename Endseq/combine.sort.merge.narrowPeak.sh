
cat *.narrowPeak | grep -v Y > ../../merge.narrowPeak

bedtools sort -i  merge.narrowPeak -g \ /home/Method/Reference_Genome/GRCm39.hisat2.bulid/Mus_musculus.GRCm39.dna.primary_assembly.fa.fai > sort.merge.narrowPeak

bedtools merge -i sort.merge.narrowPeak > combine.sort.merge.narrowPeak