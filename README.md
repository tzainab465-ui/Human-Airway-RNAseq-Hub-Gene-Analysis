# Human Airway RNA-seq Analysis and Hub Gene Identification

## Project Overview

This project presents a complete bioinformatics workflow for analyzing human airway RNA-seq data and identifying candidate hub genes associated with albumin treatment.

The analysis includes RNA-seq quality control, read trimming, alignment, gene quantification, differential expression analysis, functional enrichment, weighted gene co-expression network analysis (WGCNA), and hub gene identification.

## Workflow

Raw RNA-seq Data
↓
SRA Toolkit
↓
FASTQ Files
↓
FastQC / fastp
↓
Quality Control & Trimming
↓
HISAT2 Alignment
↓
featureCounts
↓
Gene Count Matrix
↓
DESeq2
↓
Differentially Expressed Genes
↓
GO / KEGG Enrichment
↓
VST Transformation
↓
WGCNA
↓
Module–Trait Analysis
↓
Hub Gene Identification
↓
STRING / PPI Analysis

## Dataset

- Organism: Homo sapiens
- Reference genome: GRCh38
- Number of samples: 8
- Experimental comparison: Albumin-treated vs untreated airway samples

## Differential Expression Analysis

DESeq2 was used to identify differentially expressed genes.

Significant DEGs were selected using:

- Adjusted p-value < 0.05
- |log2 Fold Change| > 1

Two significant DEGs were identified after filtering.

## Functional Enrichment

Functional analysis was performed using Gene Ontology (GO) and KEGG pathway enrichment.

The analysis identified pathways related to:

- Circadian rhythm
- Circadian entrainment
- Leukocyte transendothelial migration
- Cell adhesion molecule interactions

## WGCNA

Weighted Gene Co-expression Network Analysis (WGCNA) was performed to identify gene co-expression modules associated with the experimental condition.

Key parameters:

- Soft-thresholding power (β): 11
- Minimum module size: 30
- TOM type: unsigned
- Module merging cut height: 0.25

A total of 26,821 genes were used for WGCNA after preprocessing and filtering.

## Hub Gene Identification

Candidate hub genes were identified using:

- Module Membership (MM) ≥ 0.90
- Gene Significance (GS) ≥ 0.80

The ME24 and ME59 modules were further investigated for hub genes.

## PPI Network Analysis

STRING was used to investigate protein-protein interactions among candidate hub genes.

The resulting networks were evaluated based on:

- Number of nodes
- Number of edges
- Average node degree
- PPI enrichment

## Tools and Software

- SRA Toolkit
- FastQC
- fastp
- HISAT2
- samtools
- featureCounts
- R / RStudio
- DESeq2
- WGCNA
- clusterProfiler
- STRING
- Cytoscape

## Project Structure

```text
Human-Airway-RNAseq-Hub-Gene-Analysis/
│
├── README.md
├── scripts/
├── results/
├── figures/
└── report/
## Objective

The main objective of this project is to analyze human airway transcriptomic data and identify genes and co-expression modules that may be associated with the biological response to albumin treatment.

## Author

**Zainab Tahir**

BS Bioinformatics  
Government College University Faisalabad
