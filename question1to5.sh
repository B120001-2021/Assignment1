#! /bin/bash

#making a directory for this assignment
mkdir assignment1
#working on assignment1 directory for this assignment
cd assignment1

#use fastqc to check the quality of illumina sequencing results
#unzip and integrate the 45 RNA-seq files to one file for output one fastqc report
zcat /localdisk/data/BPSM/AY21/fastq/*.fq.gz > allsample.fq  
#fastqc with 4 threads
fastqc -t 4 allsample.fq
firefox *.html #there would be a report in html after fastqc and can be oepn with firefox

#align the forward and reverse reads with bowtie2
#extract referencegene for reads alignment
zcat /localdisk/data/BPSM/AY21/Tcongo_genome/TriTrypDB-46_TcongolenseIL3000_2019_Genome.fasta.gz > referencegene.fq
bowtie2-build -f referencegene.fq referencegene #build index and name the index as referencegene

#align the forward and reverse reads against reference gene for each sample
for i in /localdisk/data/BPSM/AY21/fastq/*1.fq.gz #wildcard was used here for extract the name of all fq.gz files that needs to be aligned
do
name=$(basename $i _1.fq.gz) #get the name of each variable without suffix
bowtie2 -x referencegene -1 /localdisk/data/BPSM/AY21/fastq/${name}_1.fq.gz -2 /localdisk/data/BPSM/AY21/fastq/${name}_2.fq.gz|samtools view -b -o ${name}.bam
done

rm reference* #remove the indexed referencegene for alignments because we don't need it anymore

#for sort the bam file first to allow following construction of indexed bam file, also makes bedtools run fatser when compare against the reference genotype
for i in *.bam
do
bamname=$(basename $i .bam)
samtools sort ${bamname}.bam -o ${bamname}.sort.bam
done

#construct indexed bam file
for i in *sort.bam
do
sortedname=$(basename $i .sort.bam)
samtools index ${sortedname}.sort.bam ${sortedname}.sort.bam.bai
done

#generate gene counts data for each bam file with provided reference genotype information
bedtools multicov -bams *sort.bam -bed /localdisk/data/BPSM/AY21/TriTrypDB-46_TcongolenseIL3000_2019.bed > countreads

rm *.bam *.bai

#generate a plain text file of average gene counts for each group
awk '{FS="\t";OFS="\t";}{print $4,($6+$9+$12)/3,($15+$17+$19)/3,($16+$18+$20)/3,($7+$10+$13)/3,($8+$11+$14)/3,($21+$24+$27)/3,($30+$32+$34)/3,($31+$33+$35)/3,($22+$25+$28)/3,($23+$26+$29)/3,($36+$39+$42)/3,($45+$47+$49)/3,($46+$48+$50)/3,($37+$40+$43)/3,($38+$41+$44)/3, $5;}' countreads > averagecounts
#add header for averagecounts of each group
awk 'BEGIN{print "Gene\tclone1_0_uninduced\tclone1_24_uninduced\tclone1_48_uninduced\tclone1_24_induced\tclone1_48_induced\tclone2_0_uninduced\tclone2_24_uninduced\tclone2_48_uninduced\tclone2_24_induced\tclone2_48_induced\tWT_0_uninduced\tWT_24_uninduced\tWT_48_uninduced\tWT_24_induced\tWT_48_induced\tDescription"}1' averagecounts >> averagecounts.txt
