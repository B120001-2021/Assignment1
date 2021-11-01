#! /bin/bash

#making a directory for this assignment
mkdir assignment1
#working on assignment1 directory for this assignment
cd assignment1

echo "Hello~ Fastqc report for illumina sequencing results is being generated."
#use fastqc to check the quality of illumina sequencing results
#unzip and integrate the 45 RNA-seq files to one file for output one fastqc report
zcat /localdisk/data/BPSM/AY21/fastq/*.fq.gz > allsample.fq  
#fastqc with 4 threads
fastqc -t 4 allsample.fq

#Assessment of the number and quality of the raw sequence
unzip allsample.fq.zip
cd allsample.fq
echo "The number of analyzed reads: `grep "Total Sequences" fastqc_data.txt`"
echo "The quality of raw RNA-seq: `grep "Sequences flagged as poor quality" fastqc_data.txt`"
echo "Modules that did not pass FastQC are: `grep -v "PASS" summary.txt | cut -f 1,2`" 
echo "Modules that passed FastQC are: `grep "PASS" summary.txt | cut -f 1,2`"

cd assignment1

echo "After viewing the fastqc report, please press control & C to quit."
firefox *.html #there would be a report in html after fastqc and can be oepn with firefox

#interactive prompt for users to choose continue or exist
while true; do
    read -p "Are you satisfied with the quality of sequencing? Please type 'y' or 'n':" answer
    case $answer in
    	[Yy]* ) break;;
        [Nn]* ) echo "Bye ~ ~ ~" ; exit;;
        * ) echo "Please answer yes or no.";;
    esac
done

rm allsample*

echo "Good choice. Hello again~"

#align the forward and reverse reads with bowtie2
#extract referencegene for reads alignment
bowtie2-build -f -q /localdisk/data/BPSM/AY21/Tcongo_genome/TriTrypDB-46_TcongolenseIL3000_2019_Genome.fasta.gz referencegene #build index and name the index as referencegene
echo "Index of Tcongo_genome has been built."

#align the forward and reverse reads against reference gene for each sample
for i in /localdisk/data/BPSM/AY21/fastq/*1.fq.gz #wildcard was used here for extract the name of all fq.gz files that needs to be aligned
do
name=$(basename $i _1.fq.gz) #get the name of each variable without suffix
bowtie2 -x referencegene -1 /localdisk/data/BPSM/AY21/fastq/${name}_1.fq.gz -2 /localdisk/data/BPSM/AY21/fastq/${name}_2.fq.gz|samtools sort -o ${name}.sort.bam
done
echo "Sequence of each sample has been aligned with bowtie2 and sorted by samtools"

rm reference* #remove the indexed referencegene for alignments because we don't need it anymore

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

#For analyzing effect of RNAi in clone1 without treatment
awk '{FS="\t";OFS="\t";} NR>=2{print $1,$2,$12,log(($2+0.000000001)/($12+0.000000001))/log(2);}' averagecounts.txt > strain_unC1toWT_1 #4
awk '{FS="\t";OFS="\t";} NR>=2{print $3,$13,log(($3+0.000000001)/($13+0.000000001))/log(2);}' averagecounts.txt > strain_unC1toWT_2 #7
awk '{FS="\t";OFS="\t";} NR>=2{print $4,$14,log(($4+0.000000001)/($14+0.000000001))/log(2);}' averagecounts.txt > strain_unC1toWT_3 #10
paste strain_unC1toWT_1 strain_unC1toWT_2  strain_unC1toWT_3 > strain_unC1toWT
awk 'BEGIN{print "Gene\tclone1_0_uninduced\tWT_0_uninduced\tLogFC0\tclone1_24_uninduced\tWT_24_uninduced\tLogFC24\tclone1_48_uninduced\tWT_48_uninduced\tLogFC48"}1' strain_unC1toWT >> strain_unC1toWT.txt

#For analyzing effect of RNAi in clone2 without treatment
awk '{FS="\t";OFS="\t";} NR>=2{print $1,$7,$12,log(($7+0.000000001)/($12+0.000000001))/log(2);}' averagecounts.txt > strain_unC2toWT_1 #4
awk '{FS="\t";OFS="\t";} NR>=2{print $8,$13,log(($8+0.000000001)/($13+0.000000001))/log(2);}' averagecounts.txt > strain_unC2toWT_2 #7
awk '{FS="\t";OFS="\t";} NR>=2{print $9,$14,log(($9+0.000000001)/($14+0.000000001))/log(2);}' averagecounts.txt > strain_unC2toWT_3 #10
paste strain_unC2toWT_1 strain_unC2toWT_2  strain_unC2toWT_3 > strain_unC2toWT
awk 'BEGIN{print "Gene\tclone2_0_uninduced\tWT_0_uninduced\tLogFC0\tclone2_24_uninduced\tWT_24_uninduced\tLogFC24\tclone2_48_uninduced\tWT_48_uninduced\tLogFC48"}1' strain_unC2toWT >> strain_unC2toWT.txt

#For analyzing effect of RNAi in clone1 with treatment
awk '{FS="\t";OFS="\t";} NR>=2{print $1,$5,$15,log(($5+0.000000001)/($15+0.000000001))/log(2);}' averagecounts.txt > strain_inC1toWT_1 #4
awk '{FS="\t";OFS="\t";} NR>=2{print $6,$16,log(($6+0.000000001)/($16+0.000000001))/log(2);}' averagecounts.txt > strain_inC1toWT_2 #7
paste strain_inC1toWT_1 strain_inC1toWT_2 > strain_inC1toWT
awk 'BEGIN{print "Gene\tclone1_24_induced\tWT_24_induced\tLogFC24\tclone1_48_induced\tWT_48_induced\tLogFC2"}1' strain_inC1toWT >> strain_inC1toWT.txt

#For analyzing effect of RNAi in clone2 with treatment
awk '{FS="\t";OFS="\t";} NR>=2{print $1,$10,$15,log(($10+0.000000001)/($15+0.000000001))/log(2);}' averagecounts.txt > strain_inC2toWT_1 #4
awk '{FS="\t";OFS="\t";} NR>=2{print $11,$16,log(($11+0.000000001)/($16+0.000000001))/log(2);}' averagecounts.txt > strain_inC2toWT_2 #7
paste strain_inC2toWT_1 strain_inC2toWT_2 > strain_inC2toWT
awk 'BEGIN{print "Gene\tclone2_24_induced\tWT_24_induced\tLogFC24\tclone2_48_induced\tWT_48_induced\tLogFC48"}1' strain_inC2toWT >> strain_inC2toWT.txt

#For analyzing effect of tetracycline in WT
awk '{FS="\t";OFS="\t";} NR>=2{print $1,$15,$13,log(($15+0.000000001)/($13+0.000000001))/log(2);}' averagecounts.txt > treatment_WT_1 #4
awk '{FS="\t";OFS="\t";} NR>=2{print $16,$14,log(($16+0.000000001)/($14+0.000000001))/log(2);}' averagecounts.txt > treatment_WT_2 #7
paste treatment_WT_1 treatment_WT_2 > treatment_WT
awk 'BEGIN{print "Gene\tWT_24_induced\tWT_24_uninduced\tLogFC24\tWT_48_induced\tWT_48_uninduced\tLogFC48"}1' treatment_WT >> treatment_WT.txt

#For analyzing effect of tetracycline in clone1
awk '{FS="\t";OFS="\t";} NR>=2{print $1,$5,$3,log(($5+0.000000001)/($3+0.000000001))/log(2);}' averagecounts.txt > treatment_C1_1 #4
awk '{FS="\t";OFS="\t";} NR>=2{print $6,$4,log(($6+0.000000001)/($4+0.000000001))/log(2);}' averagecounts.txt > treatment_C1_2 #7
paste treatment_C1_1 treatment_C1_2 > treatment_C1
awk 'BEGIN{print "Gene\tC1_24_induced\tC1_24_uninduced\tLogFC24\tC1_48_induced\tC1_48_uninduced\tLogFC48"}1' treatment_C1 >> treatment_C1.txt

#For analyzing effect of tetracycline in clone2
awk '{FS="\t";OFS="\t";} NR>=2{print $1,$10,$8,log(($10+0.000000001)/($8+0.000000001))/log(2);}' averagecounts.txt > treatment_C2_1 #4
awk '{FS="\t";OFS="\t";} NR>=2{print $11,$9,log(($11+0.000000001)/($9+0.000000001))/log(2);}' averagecounts.txt > treatment_C2_2 #7
paste treatment_C2_1 treatment_C2_2 > treatment_C2
awk 'BEGIN{print "Gene\tC2_24_induced\tC2_24_uninduced\tLogFC24\tC2_48_induced\tC2_48_uninduced\tLogFC48"}1' treatment_C2 >> treatment_C2.txt

#For analyzing effect of time in WT without treatment
awk '{FS="\t";OFS="\t";} NR>=2{print $1,$13,$12,log(($13+0.000000001)/($12+0.000000001))/log(2);}' averagecounts.txt > time_unWT_1 #4
awk '{FS="\t";OFS="\t";} NR>=2{print $14,$12,log(($14+0.000000001)/($12+0.000000001))/log(2);}' averagecounts.txt > time_unWT_2 #7
paste time_unWT_1 time_unWT_2 > time_unWT
awk 'BEGIN{print "Gene\tWT_24_uninduced\tWT_0_uninduced\tLogFC240\tWT_48_uninduced\tWT_0_uninduced\tLogFC480"}1' time_unWT >> time_unWT.txt

#For analyzing effect of time in clone1 without treatment
awk '{FS="\t";OFS="\t";} NR>=2{print $1,$3,$2,log(($3+0.000000001)/($2+0.000000001))/log(2);}' averagecounts.txt > time_unC1_1 #4
awk '{FS="\t";OFS="\t";} NR>=2{print $4,$2,log(($4+0.000000001)/($2+0.000000001))/log(2);}' averagecounts.txt > time_unC1_2 #7
paste time_unC1_1 time_unC1_2 > time_unC1
awk 'BEGIN{print "Gene\tC1_24_uninduced\tC1_0_uninduced\tLogFC240\tC1_48_uninduced\tC1_0_uninduced\tLogFC480"}1' time_unC1 >> time_unC1.txt

#For analyzing effect of time in clone1 with treatment
awk '{FS="\t";OFS="\t";} NR>=2{print $1,$5,$2,log(($5+0.000000001)/($2+0.000000001))/log(2);}' averagecounts.txt > time_inC1_1 #4
awk '{FS="\t";OFS="\t";} NR>=2{print $6,$2,log(($6+0.000000001)/($2+0.000000001))/log(2);}' averagecounts.txt > time_inC1_2 #7
paste time_inC1_1 time_inC1_2 > time_inC1
awk 'BEGIN{print "Gene\tC1_24_induced\tC1_0_uninduced\tLogFC240\tC1_48_induced\tC1_0_uninduced\tLogFC480"}1' time_inC1 >> time_inC1.txt

#For analyzing effect of time in clone2 without treatment
awk '{FS="\t";OFS="\t";} NR>=2{print $1,$8,$7,log(($8+0.000000001)/($7+0.000000001))/log(2);}' averagecounts.txt > time_unC2_1 #4
awk '{FS="\t";OFS="\t";} NR>=2{print $9,$7,log(($9+0.000000001)/($7+0.000000001))/log(2);}' averagecounts.txt > time_unC2_2 #7
paste time_unC2_1 time_unC2_2 > time_unC2
awk 'BEGIN{print "Gene\tC1_24_uninduced\tC1_0_uninduced\tLogFC240\tC1_48_uninduced\tC1_0_uninduced\tLogFC480"}1' time_unC2 >> time_unC2.txt

#For analyzing effect of time in clone2 with treatment
awk '{FS="\t";OFS="\t";} NR>=2{print $1,$10,$7,log(($10+0.000000001)/($7+0.000000001))/log(2);}' averagecounts.txt > time_inC2_1 #4
awk '{FS="\t";OFS="\t";} NR>=2{print $11,$7,log(($11+0.000000001)/($7+0.000000001))/log(2);}' averagecounts.txt > time_inC2_2 #7
paste time_inC2_1 time_inC2_2 > time_inC2
awk 'BEGIN{print "Gene\tC2_24_induced\tC2_0_uninduced\tLogFC240\tC2_48_induced\tC2_0_uninduced\tLogFC480"}1' time_inC2 >> time_inC2.txt

mkdir output
cp *.txt output
rm *
cd output
echo "The file of average gene counts and group-wise comparisons are shown in txt format in here."

#additional analyses
echo "The total number of genes that had been assesed: `wc -l averagecounts.txt`" #10313
echo "Number of genes that is signficantly inhibited in clone1 at time 0 is `awk '{FS="\t"; if($4 <= -2){print $0;}}' strain_unC1toWT.txt | wc -l` "
echo "Number of genes that is signficantly inhibited in clone1 at time 24 is `awk '{FS="\t"; if($7 <= -2){print $0;}}' strain_unC1toWT.txt | wc -l` "
echo "Number of genes that is signficantly inhibited in clone1 at time 48 is `awk '{FS="\t"; if($4 <= -2){print $0;}}' strain_unC1toWT.txt | wc -l` "
echo "Please refer the PDF for more information"
#other data was obatined by a same silly way that is not shown in here for saving space on the server.
