############################################################
# HUMAN RNA-SEQ ANALYSIS PIPELINE
# Untreated vs Alb
# DESeq2 + DEG + PCA + MA + Volcano + Heatmap + KEGG
############################################################


#############################
# 1. Install packages
#############################

# Run only once

install.packages("BiocManager")

BiocManager::install(c(
  "DESeq2",
  "pheatmap",
  "clusterProfiler",
  "enrichplot",
  "org.Hs.eg.db"
))


#############################
# 2. Load libraries
#############################

library(DESeq2)
library(pheatmap)
library(ggplot2)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(AnnotationDbi)


#############################
# 3. Set working directory
#############################

setwd("D:/RNAseq/Human-rnaseq/counts")


#############################
# 4. Read count file
#############################

counts <- read.delim(
  "gene_counts.txt",
  comment.char="#"
)


# Remove annotation columns

count_matrix <- counts[,7:14]


rownames(count_matrix) <- counts$Geneid


#############################
# 5. Rename samples
#############################

colnames(count_matrix) <- c(
  "SRR1039508",
  "SRR1039510",
  "SRR1039512",
  "SRR1039514",
  "SRR1039516",
  "SRR1039518",
  "SRR1039520",
  "SRR1039522"
)



#############################
# 6. Create sample metadata
#############################

condition <- factor(c(
  "Untreated",
  "Alb",
  "Untreated",
  "Alb",
  "Untreated",
  "Alb",
  "Untreated",
  "Alb"
))


colData <- data.frame(
  row.names = colnames(count_matrix),
  condition
)


colData


#############################
# 7. DESeq2 analysis
#############################

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = colData,
  design = ~ condition
)


dds <- DESeq(dds)


res <- results(
  dds,
  contrast=c(
    "condition",
    "Untreated",
    "Alb"
  )
)


res <- res[order(res$padj),]


head(res)


#############################
# 8. Save DESeq2 results
#############################

write.csv(
  as.data.frame(res),
  "DESeq2_results.csv"
)



#############################
# 9. PCA PLOT
#############################

vsd <- vst(
  dds,
  blind=FALSE
)


pdf("PCA_plot.pdf")

plotPCA(
  vsd,
  intgroup="condition"
)

dev.off()



#############################
# 10. MA PLOT
#############################

pdf("MA_plot.pdf")

plotMA(
  res,
  ylim=c(-5,5)
)

dev.off()



#############################
# 11. Volcano Plot
#############################

volcano_data <- as.data.frame(res)


volcano_data$Significance <- "Not Significant"


volcano_data$Significance[
  volcano_data$padj < 0.05 &
    abs(volcano_data$log2FoldChange)>1
] <- "Significant"



ggplot(
  volcano_data,
  aes(
    x=log2FoldChange,
    y=-log10(padj),
    color=Significance
  )
)+
  
  geom_point()+
  
  theme_classic()+
  
  labs(
    title="Volcano Plot",
    x="log2 Fold Change",
    y="-log10 Adjusted P-value"
  )


ggsave(
  "Volcano_plot.png",
  width=7,
  height=6
)



#############################
# 12. Significant DEGs
#############################

sig <- subset(
  res,
  padj <0.05 &
    abs(log2FoldChange)>1
)


write.csv(
  as.data.frame(sig),
  "Significant_DEGs.csv"
)


nrow(sig)



#############################
# 13. Heatmap
#############################

top_genes <- head(
  rownames(res),
  50
)


heatmap_data <- assay(vsd)[top_genes,]


pheatmap(
  heatmap_data,
  annotation_col=colData,
  filename="Heatmap.png"
)



#############################
# 14. Ensembl ID conversion
#############################

gene_conversion <- bitr(
  rownames(sig),
  fromType="ENSEMBL",
  toType="ENTREZID",
  OrgDb=org.Hs.eg.db
)


gene_conversion


write.csv(
  gene_conversion,
  "Gene_ID_conversion.csv"
)



#############################
# 15. KEGG pathway analysis
#############################

kegg_result <- enrichKEGG(
  gene=gene_conversion$ENTREZID,
  organism="hsa",
  pvalueCutoff=0.05
)


kegg_result



#############################
# 16. KEGG Dot Plot
#############################

pdf(
  "KEGG_dotplot.pdf",
  width=8,
  height=6
)


dotplot(
  kegg_result
)


dev.off()



#############################
# 17. KEGG Bar Plot
#############################

pdf(
  "KEGG_barplot.pdf",
  width=8,
  height=6
)


barplot(
  kegg_result
)


dev.off()



#############################
# 18. Save KEGG table
#############################

write.csv(
  as.data.frame(kegg_result),
  "KEGG_results.csv"
)
#############################
# 15. GO Biological Process Analysis
#############################

GO_BP <- enrichGO(
  gene = gene_conversion$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)


#############################
# 16. GO Dot Plot
#############################

pdf(
  "GO_BP_dotplot.pdf",
  width = 12,
  height = 10
)

dotplot(
  GO_BP,
  showCategory = 10,
  font.size = 10,
  label_format = 40,
  title = "GO Biological Process Enrichment"
)

dev.off()


#############################
# 17. GO Dot Plot (PNG)
#############################

png(
  "GO_BP_dotplot.png",
  width = 3600,
  height = 3000,
  res = 300
)

dotplot(
  GO_BP,
  showCategory = 10,
  font.size = 10,
  label_format = 40,
  title = "GO Biological Process Enrichment"
)

dev.off()


#############################
# 18. Save GO Results
#############################

write.csv(
  as.data.frame(GO_BP),
  "GO_BP_results.csv",
  row.names = FALSE
)
############################################################
# END OF ANALYSIS
############################################################