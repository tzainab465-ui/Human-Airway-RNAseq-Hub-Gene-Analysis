############################################################
# WGCNA ANALYSIS OF RNA-seq DATA
############################################################

#############################
# 1. Load Library
#############################

library(WGCNA)

options(stringsAsFactors = FALSE)
enableWGCNAThreads()

#############################
# 2. Expression Matrix
#############################

expr <- assay(vsd)
datExpr <- t(expr)

#############################
# 3. Quality Control
#############################

gsg <- goodSamplesGenes(datExpr, verbose = 3)

datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]

#############################
# 4. Sample Clustering
#############################

sampleTree <- hclust(dist(datExpr), method = "average")

plot(sampleTree,
     main = "Sample Clustering",
     sub = "",
     xlab = "")

#############################
# 5. Soft Threshold Selection
#############################

powers <- 1:20

sft <- pickSoftThreshold(
  datExpr,
  powerVector = powers,
  verbose = 5
)

softPower <- 11

#############################
# 6. Scale-Free Topology Plot
#############################

plot(
  sft$fitIndices[,1],
  -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
  xlab = "Soft Threshold (Power)",
  ylab = "Scale Free Topology Model Fit",
  type = "n",
  main = "Scale Independence"
)

text(
  sft$fitIndices[,1],
  -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
  labels = powers,
  col = "red"
)

abline(h = 0.90, col = "blue", lty = 2)

#############################
# 7. Mean Connectivity
#############################

plot(
  sft$fitIndices[,1],
  sft$fitIndices[,5],
  xlab = "Soft Threshold (Power)",
  ylab = "Mean Connectivity",
  type = "n",
  main = "Mean Connectivity"
)

text(
  sft$fitIndices[,1],
  sft$fitIndices[,5],
  labels = powers,
  col = "red"
)

#############################
# 8. Network Construction
#############################

# blockwiseModules() internally computes:
# 1. Adjacency matrix
# 2. Topological Overlap Matrix (TOM)
# 3. Dissimilarity (1 - TOM)
# 4. Hierarchical clustering
# 5. Dynamic Tree Cut
# 6. Module eigengenes

net <- blockwiseModules(
  datExpr,
  power = softPower,
  TOMType = "unsigned",
  minModuleSize = 30,
  reassignThreshold = 0,
  mergeCutHeight = 0.25,
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  verbose = 3
)

#############################
# 9. Module Colors
#############################

moduleColors <- labels2colors(net$colors)

table(moduleColors)

#############################
# 10. Gene Dendrogram
#############################

plotDendroAndColors(
  net$dendrograms[[1]],
  moduleColors[net$blockGenes[[1]]],
  "Module Colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05
)

#############################
# 11. Trait Data
#############################

traitData <- data.frame(
  Alb = c(0,1,0,1,0,1,0,1)
)

rownames(traitData) <- rownames(datExpr)

#############################
# 12. Module Eigengenes
#############################

MEs <- net$MEs

#############################
# 13. Module-Trait Relationship
#############################

moduleTraitCor <- cor(
  MEs,
  traitData,
  use = "p"
)

moduleTraitPvalue <- corPvalueStudent(
  moduleTraitCor,
  nSamples = nrow(datExpr)
)

#############################
# 14. Module-Trait Heatmap
#############################

textMatrix <- paste(
  signif(moduleTraitCor,2),
  "\n(",
  signif(moduleTraitPvalue,1),
  ")",
  sep=""
)

dim(textMatrix) <- dim(moduleTraitCor)

labeledHeatmap(
  Matrix = moduleTraitCor,
  xLabels = names(traitData),
  yLabels = names(MEs),
  ySymbols = names(MEs),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix,
  setStdMargins = FALSE,
  cex.text = 0.6,
  zlim = c(-1,1),
  main = "Module-Trait Relationships"
)

#############################
# 15. Gene Significance
#############################

geneTraitSignificance <- as.data.frame(
  cor(
    datExpr,
    traitData$Alb,
    use = "p"
  )
)

#############################
# 16. Module Membership
#############################

geneModuleMembership <- as.data.frame(
  cor(
    datExpr,
    MEs,
    use = "p"
  )
)

#############################
# 17. Hub Genes (ME24)
#############################

moduleGenes24 <- net$colors == 24

hubME24 <- data.frame(
  Gene = colnames(datExpr)[moduleGenes24],
  MM = geneModuleMembership[moduleGenes24,"ME24"],
  GS = geneTraitSignificance[moduleGenes24,1]
)

hubME24 <- hubME24[
  order(abs(hubME24$MM), decreasing = TRUE),
]

head(hubME24,20)

#############################
# 18. Hub Genes (ME59)
#############################

moduleGenes59 <- net$colors == 59

hubME59 <- data.frame(
  Gene = colnames(datExpr)[moduleGenes59],
  MM = geneModuleMembership[moduleGenes59,"ME59"],
  GS = geneTraitSignificance[moduleGenes59,1]
)

hubME59 <- hubME59[
  order(abs(hubME59$MM), decreasing = TRUE),
]

head(hubME59,20)

#############################
# 19. Export Results
#############################

write.csv(
  hubME24,
  "HubGenes_ME24.csv",
  row.names = FALSE
)

write.csv(
  hubME59,
  "HubGenes_ME59.csv",
  row.names = FALSE
)

############################################################
# END OF WGCNA ANALYSIS
############################################################