#! /bin/bash

if [ $# -ne 1 ];then
	echo Usage:bash $0 your_sample_name
	echo Note:this script only for pair-end fastq file
	exit 1
fi

#map to genome
echo run $1 hisat
cd /data/yinxx/bioinfo_projects/PMI_smoking/rawdata
time hisat2 -p 30 -x /data/yinxx/bioinfo_projects/PMI_smoking/rawdata/data_genome/Homo_sapiens.GRCh38.dna.chromosome.1_tran -1 ${1}.R1.fastq.gz -2 ${1}.R2.fastq.gz -S /data/yinxx/bioinfo_projects/PMI_smoking/hisat_align/${1}.sam  2>/data/yinxx/bioinfo_projects/PMI_smoking/hisat_align/${1}.align.log
wait
echo done $1 hisat align

# convert sam to bam 
cd /data/yinxx/bioinfo_projects/PMI_smoking/hisat_align
pwd
echo $1 samtools samtobam
time samtools sort -@ 30 -o ${1}.bam ${1}.sam
wait
echo done $1 bam file convert

# run stringtie
echo run $1 stringtie
mkdir /data/yinxx/bioinfo_projects/PMI_smoking/Ballgown_out/${1}
time stringtie -B -e -p 30 -G /data/yinxx/bioinfo_projects/PMI_smoking/rawdata/data_genome/Homo_sapiens.GRCh38.92.chr.gtf -o ../Ballgown_out/${1}/${1}_dta.gtf ${1}.bam -A ../stringtie_abun_out/${1}_abun.txt
wait
echo done $1 stringtie 
