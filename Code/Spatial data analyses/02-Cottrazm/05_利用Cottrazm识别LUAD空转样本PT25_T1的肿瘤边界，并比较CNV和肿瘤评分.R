######目的：利用Cottrazm识别LUAD空转样本PT25_T1的肿瘤边界，并比较CNV和肿瘤评分
######作者：申奥
######日期：2024-08-01


# 0.环境设置----
rm(list = ls())
options(stringAsFactor = F)

library(Seurat)
library(magrittr)
library(dplyr)
library(Matrix)
library(ggplot2)
library(stringr)
library(RColorBrewer)
library(patchwork)
library(ggtree)
library(BiocGenerics)
library(readr)
library(rtracklayer)
library(infercnv)
library(phylogram)
library(utils)
library(dendextend)
library(assertthat)
library(reticulate)
library(openxlsx)
library(scatterpie)
library(cowplot)
library(stats)
library(quadprog)
library(data.table)
library(Rfast)
library(ggrepel)
library(tibble)
library(clusterProfiler)
library(utils)
library(org.Hs.eg.db)
library(Cottrazm)
library(ggpubr)
library(GSVA)
library(ggsci)
library(RColorBrewer)
library(cols4all)
library(ggalluvial)

source("STCNV_rewrite.R")

# 1.PT15_T1----
## 1.1 识别边界----
InDir <- "1.data/SP/PT25_T1/"
Sample = "PT25_T1"
OutDir <- paste("2.Cottrazm_outputs/", Sample, "/", sep = "")
dir.create(paste0("2.Cottrazm_outputs/", Sample))
TumorST <- STPreProcess(InDir = InDir, Sample = Sample, OutDir = OutDir)

res <- 1.5
TumorST <- STModiCluster(InDir = InDir,
                         Sample = Sample,
                         OutDir = OutDir,
                         TumorST = TumorST,
                         res = res)

STInferCNV <- STCNV(TumorST = TumorST, assay = "Spatial", OutDir = OutDir, num_threads = 48)
TumorST <- STCNVScore(TumorST = TumorST, assay = "Spatial", Sample = Sample, OutDir = OutDir)
TumorSTn <- BoundaryDefine(TumorST = TumorST, MalLabel = c(1,2,3,4,7,8), OutDir = OutDir, Sample = Sample)
TumorST <- BoundaryPlot(TumorSTn = TumorSTn, TumorST = TumorST, OutDir = OutDir, Sample = Sample)

save(TumorST, file = paste0(OutDir, Sample, "_TumorST_Boudary.RData"))

## 1.2 CNV评分比较----
table(TumorST$Location)
# Mal  Bdy nMal 
# 769  309  230
cnv_score <- TumorST@meta.data[,c("cnv_score", "Location")]
class(cnv_score$Location)
class(cnv_score$cnv_score)

comparison <- list(c("Mal", "Bdy"),
                   c("Mal", "nMal"),
                   c("Bdy", "nMal"))
cnvscore_plot <- ggplot(cnv_score, aes(x = Location, y = cnv_score, fill = Location)) +
  geom_boxplot(fill = c("#cb181d", "#1f78b4", "#fdb462"), width = 0.6, alpha = 0.9) +
  theme_pubr() +
  theme(axis.title = element_text(size = 14)) +
  theme(axis.text = element_text(size = 12)) +
  theme(legend.position = "none") +
  ggtitle(Sample) +
  ylab("CNV scores") + xlab(NULL) +
  stat_compare_means(comparisons = comparison, method = "wilcox.test", label = "p.format")
cnvscore_plot
ggsave(cnvscore_plot, filename = paste0(OutDir, Sample, "_CNVscore.pdf"), height = 4, width = 5.5, dpi = 300)