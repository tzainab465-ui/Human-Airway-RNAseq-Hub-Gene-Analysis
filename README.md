# Human Airway RNA-seq Analysis and Hub Gene Identification

**Author:** Zainab Tahir  
**Program:** BS Bioinformatics  
**Institution:** Government College University Faisalabad (GCUF)

## Project Overview

This project presents a bioinformatics workflow for analyzing human airway RNA-seq data and identifying candidate hub genes associated with albumin treatment.

The analysis includes RNA-seq quality control, read trimming, alignment, gene quantification, differential expression analysis, functional enrichment, weighted gene co-expression network analysis (WGCNA), hub gene identification, and protein-protein interaction (PPI) analysis.

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
Module-Trait Analysis  
↓  
Hub Gene Identification  
↓  
STRING / PPI Analysis

## Dataset

- **Organism:** Homo sapiens
- **Reference genome:** GRCh38
- **Number of samples:** 8
- **Experimental comparison:** Albumin-treated vs untreated airway samples

## Differential Expression Analysis

DESeq2 was used to identify differentially expressed genes.

Significant DEGs were selected using:

- Adjusted p-value < 0.05
- |log2 Fold Change| > 1

After filtering, **2 significant DEGs** were identified.

The significant genes were:

- ENSG00000292348
- ENSG00000179094

## Functional Enrichment

Functional enrichment analysis was performed using Gene Ontology (GO) and KEGG pathway analysis.

The analysis identified pathways related to:

- Circadian rhythm
- Circadian entrainment
- Leukocyte transendothelial migration
- Cell adhesion molecule interactions

## WGCNA

Weighted Gene Co-expression Network Analysis (WGCNA) was performed to identify gene co-expression modules associated with the experimental condition.

### Key Parameters

- **Soft-thresholding power (β):** 11
- **Minimum module size:** 30
- **TOM type:** Unsigned
- **Module merging cut height:** 0.25
- **Genes used:** 26,821
- **Samples:** 8

## Hub Gene Identification

Candidate hub genes were identified using the following criteria:

- **|Module Membership (MM)| ≥ 0.90**
- **Gene Significance (GS) ≥ 0.80**

Two modules were selected for further investigation:
| Module | Total Genes | Hub Genes After Filtering |
|--------|------------:|---------------------------:|
| ME24 (darkgrey) | 407 | 11 |
| ME59 (darkseagreen4) | 190 | 20 |
| **Total** | **597** | **31** |

The ME59 filtered result contained **20 genes** satisfying the specified MM and GS criteria.

A total of **31 candidate hub genes** were identified across the ME24 and ME59 modules after applying the filtering criteria of |MM| ≥ 0.90 and GS ≥ 0.80.

## PPI Network Results

Protein-protein interaction (PPI) analysis was performed using STRING to investigate potential interactions among the candidate hub genes.

### ME24 PPI Network

The ME24 network contained:

- **Nodes:** 8
- **Edges:** 1
- **Average node degree:** 0.25
- **Average local clustering coefficient:** 0.25
- **Expected number of edges:** 0
- **PPI enrichment p-value:** 0.183

The ME24 network showed a small number of interactions and was **not significantly enriched for protein-protein interactions** (PPI enrichment p = 0.183).

### ME59 PPI Network

The ME59 network contained:

- **Nodes:** 3
- **Edges:** 0
- **Average node degree:** 0
- **PPI enrichment p-value:** 1.0

No protein-protein interactions were detected among the proteins included in the ME59 STRING network. The network therefore showed **no significant PPI enrichment**.

### Overall PPI Interpretation

The STRING analysis showed limited protein-protein interactions among the investigated candidate hub proteins. The ME24 network contained one observed interaction, whereas the ME59 network contained no observed interactions. Therefore, the candidate hub genes did not form strongly interconnected PPI networks under the STRING analysis conditions used.

## Tools and Software

### Quality Control and Processing

- SRA Toolkit
- FastQC
- fastp
- HISAT2
- samtools
- featureCounts

### Statistical and Network Analysis

- R
- RStudio
- DESeq2
- WGCNA
- clusterProfiler

### Functional and PPI Analysis

- Gene Ontology (GO)
- KEGG
- STRING
- Cytoscape

## Project Structure

```text
Human-Airway-RNAseq-Hub-Gene-Analysis/
│
├── README.md
│
├── scripts/
│
├── results/
│
├── figures/
│
└── report/
