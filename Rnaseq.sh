#!/bin/bash

#################################################
# RNA-seq Analysis Pipeline
# Human RNA-seq
#################################################

###############
# 1. Quality Control
###############

fastqc -t 8 fastq/*.fastq.gz -o fastqc_before

###############
# 2. Read Trimming
###############

for sample in SRR1039508 SRR1039510 SRR1039512 SRR1039514 SRR1039516 SRR1039518 SRR1039520 SRR1039522
do
fastp \
-i fastq/${sample}_1.fastq.gz \
-I fastq/${sample}_2.fastq.gz \
-o trimmed_fastq/${sample}_1.trim.fastq \
-O trimmed_fastq/${sample}_2.trim.fastq \
-h fastp_reports/${sample}_fastp.html \
-j fastp_reports/${sample}_fastp.json \
-w 8
done

###############
# 3. FastQC After Trimming
###############

fastqc -t 8 trimmed_fastq/*.trim.fastq -o fastqc_after

###############
# 4. HISAT2 Genome Index
###############

hisat2-build -p 8 \
reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
reference/human_index

###############
# 5. Alignment
###############

for sample in SRR1039508 SRR1039510 SRR1039512 SRR1039514 SRR1039516 SRR1039518 SRR1039520 SRR1039522
do
hisat2 -p 8 \
-x reference/human_index \
-1 trimmed_fastq/${sample}_1.trim.fastq \
-2 trimmed_fastq/${sample}_2.trim.fastq \
-S alignment/${sample}.sam
done

###############
# 6. SAM → BAM
###############

for sample in SRR1039508 SRR1039510 SRR1039512 SRR1039514 SRR1039516 SRR1039518 SRR1039520 SRR1039522
do
samtools view -@ 8 \
-bS alignment/${sample}.sam \
-o bam/${sample}.bam
done

###############
# 7. Sort BAM
###############

for sample in SRR1039508 SRR1039510 SRR1039512 SRR1039514 SRR1039516 SRR1039518 SRR1039520 SRR1039522
do
samtools sort -@ 8 \
bam/${sample}.bam \
-o sorted_bam/${sample}.sorted.bam
done

###############
# 8. BAM Indexing
###############

for sample in SRR1039508 SRR1039510 SRR1039512 SRR1039514 SRR1039516 SRR1039518 SRR1039520 SRR1039522
do
samtools index sorted_bam/${sample}.sorted.bam
done

###############
# 9. Gene Counting
###############

featureCounts \
-T 8 \
-p \
-t exon \
-g gene_id \
-a reference/Homo_sapiens.GRCh38.115.gtf \
-o counts/gene_counts.txt \
sorted_bam/*.sorted.bam
