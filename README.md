#  Human Airway RNA-seq Analysis and Hub Gene Identification

**A bioinformatics analysis of human airway RNA-seq data to identify differentially expressed genes, co-expression modules, candidate hub genes, and protein-protein interaction networks associated with albumin treatment.**

**Author:** Zainab Tahir
**Program:** BS Bioinformatics
**Institution:** Government College University Faisalabad (GCUF)

---

##  Project Overview

This project investigates gene expression patterns in human airway samples following albumin treatment using a complete RNA-seq bioinformatics workflow.

The analysis integrates:

**RNA-seq → Quality Control → Trimming → Alignment → Quantification → DESeq2 → GO/KEGG → WGCNA → Hub Genes → STRING/PPI**

The objective was to identify genes and biological pathways potentially associated with the response to albumin treatment.

---

## 🧪 Dataset

| Feature          | Information                  |
| ---------------- | ---------------------------- |
| Organism         | *Homo sapiens*               |
| Reference genome | GRCh38                       |
| Samples          | 8                            |
| Comparison       | Albumin-treated vs untreated |
| Analysis type    | RNA-seq                      |
| Network analysis | WGCNA                        |
| PPI database     | STRING                       |

---

## 🔬 Analysis Workflow

```text
Raw RNA-seq Data
       ↓
SRA Toolkit
       ↓
FASTQ
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
```

---

##  Key Results

| Analysis                  | Main Result |
| ------------------------- | ----------: |
| Genes analyzed by DESeq2  |      78,899 |
| Significant DEGs          |           2 |
| Genes used for WGCNA      |      26,821 |
| WGCNA soft threshold      |      β = 11 |
| ME24 hub genes            |          11 |
| ME59 hub genes            |          20 |
| Total candidate hub genes |          31 |
| ME24 STRING nodes         |           8 |
| ME24 STRING edges         |           1 |
| ME59 STRING nodes         |           3 |
| ME59 STRING edges         |           0 |

---

##  Differential Expression Analysis

DESeq2 was used to identify differentially expressed genes between albumin-treated and untreated airway samples.

Significant DEGs were defined using:

* Adjusted p-value < 0.05
* |log₂ Fold Change| > 1

A total of **2 significant DEGs** were identified:

* `ENSG00000292348`
* `ENSG00000179094`

### Differential Expression Figures

* [Volcano Plot](Volcano_plot.png)
* [MA Plot](MA_plot.png)
* [PCA Plot](PCA_plot.png)
* [Heatmap](Heatmap.png)
* [DEG Summary](DEG_summary_plot.png)

---

## 🧩 Functional Enrichment

GO and KEGG enrichment analyses were performed to investigate the biological functions and pathways associated with the identified genes.

Important pathways included:

* Circadian rhythm
* Circadian entrainment
* Leukocyte transendothelial migration
* Cell adhesion molecule interactions

### Enrichment Results

* [GO Biological Process Results](GO_BP_results.csv)
* [GO BP Dot Plot](GO_BP_dotplot.png)
* [KEGG Dot Plot](Dot%20plot%20KEEG.png)

---

##  WGCNA Analysis

Weighted Gene Co-expression Network Analysis (WGCNA) was performed using **26,821 genes** from 8 samples.

### Parameters

| Parameter               | Value    |
| ----------------------- | -------- |
| Soft-thresholding power | β = 11   |
| Minimum module size     | 30       |
| TOM type                | Unsigned |
| Merge cut height        | 0.25     |
| Samples                 | 8        |
| Genes                   | 26,821   |

Two modules were selected for hub gene analysis:

| Module    | Color         |   Genes | Hub Genes |
| --------- | ------------- | ------: | --------: |
| ME24      | Darkgrey      |     407 |        11 |
| ME59      | Darkseagreen4 |     190 |        20 |
| **Total** |               | **597** |    **31** |

Hub genes were selected using:

**|MM| ≥ 0.90 and GS ≥ 0.80**

### Hub Gene Results

* [ME24 Filtered Hub Genes](HubGenes_ME24_Filtered.csv)
* [ME59 Filtered Hub Genes](HubGenes_ME59_Filtered.csv)

---

##  Protein-Protein Interaction Analysis

STRING was used to investigate potential interactions among the candidate hub proteins.

### ME24

**11 hub genes were submitted → 8 proteins recognized**

* Nodes: **8**
* Edges: **1**
* Average node degree: **0.25**
* PPI enrichment p-value: **0.183**

The ME24 network was not significantly enriched for protein-protein interactions.

### ME59

**20 hub genes were submitted → 3 proteins recognized**

* Nodes: **3**
* Edges: **0**
* Average node degree: **0**
* PPI enrichment p-value: **1.0**

No interactions were detected among the three proteins recognized by STRING.

> **Note:** The PPI results represent only the proteins successfully mapped by STRING. Therefore, the absence of interactions does not indicate that all candidate hub genes lack biological interactions.

### STRING Results

* [ME24 PPI Network](string_hires_image%20ME24.png)
* [ME59 PPI Network](string_hires_image%20ME59.png)
* [ME24 Interactions](string_interactions%20ME24.tsv)
* [ME59 Interactions](string_interactions_short%20ME59.tsv)

---

##  Key Findings

1. Only **2 significant DEGs** were identified using stringent DESeq2 criteria.
2. WGCNA identified two modules associated with the experimental condition.
3. **31 candidate hub genes** were identified from ME24 and ME59.
4. Functional enrichment highlighted pathways involving circadian regulation, immune-cell migration, and cell adhesion.
5. STRING analysis showed **limited PPI connectivity** among the recognized hub proteins.
6. Several hub gene identifiers could not be mapped by STRING, representing an important limitation of the PPI analysis.

---

##  Limitations

* The study included only **8 samples**, which limits statistical power.
* Only two genes reached the predefined DESeq2 significance thresholds.
* Not all hub gene identifiers were successfully recognized by STRING.
* The PPI analysis therefore represents only a subset of the candidate hub genes.
* Experimental validation would be required to confirm the biological relevance of the identified hub genes.

---

##  Tools & Software

### RNA-seq Processing

* SRA Toolkit
* FastQC
* fastp
* HISAT2
* samtools
* featureCounts

### Statistical Analysis

* R
* RStudio
* DESeq2
* WGCNA

### Functional Analysis

* clusterProfiler
* Gene Ontology
* KEGG

### Network Analysis

* STRING
* Cytoscape

---

##  Repository Structure

```text
Human-Airway-RNAseq-Hub-Gene-Analysis/
│
├── README.md
├── scripts/
├── results/
├── figures/
└── report/
```

---

##  Author

**Zainab Tahir**
BS Bioinformatics
Government College University Faisalabad (GCUF)

---

##  Project Status

**Analysis completed:** RNA-seq → DESeq2 → Enrichment → WGCNA → Hub Gene Identification → STRING/PPI

**Status:** Completed
