######目的：GSE189357单细胞数据中进行糖酵解和Kla评分
######作者：申奥
######日期：2024-11-01
######Server: yjx@xiyoucloud.top:34022
######Version: R 4.3.3


# 0.环境设置----
rm(list = ls())
options(stringAsFactors = F)

library(Seurat, lib.loc = "~/bioSoft/seurat_v4/")
library(AUCell)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(cols4all)
library(harmony)
library(GSVA)


# 1.读取44脚本中处理好的单细胞数据以及相关基因集----
## 1.1 读取----
scMerge <- readRDS("5.scRNA/GSE189357_final_ref.rds")
table(scMerge$Cell_types)
# Normal epithelials        Tumor cells                  T                 NK                  B                 DC 
# 4601               5717              26427              17264              10067               7748 
# Mast        Macrophages        Endothelial        Fibroblasts          Monocytes             Plasma 
# 7778               7486               4141               3067               2626               2548 
# Proliferative 
# 842

## 1.2 加载基因集----
### H3K18la
load("15.diff/upKla_genes.RData")
upKla <- as.character(upKla_genes)
upKla <- list(upKla)
names(upKla) <- "H3K18la"

### KEGG glycolysis
glyco <- read.table("KEGG_GLYCOLYSIS_GLUCONEOGENESIS.v2024.1.Hs.gmt")
glyco <- glyco[,3:ncol(glyco)]
glyco <- t(glyco)
glyco <- as.character(glyco)
glyco <- list(glyco)
names(glyco) <- "Glycolysis"


# 2.对所有细胞类型都进行AUCell评分----
## 2.1 Glycolysis----
cells_rankings <- AUCell_buildRankings(scMerge@assays$RNA@data, splitByBlocks=TRUE) 

cells_AUC <- AUCell_calcAUC(glyco, cells_rankings, nCores = 10,
                            aucMaxRank = ceiling(0.1 * nrow(cells_rankings))) # 选择基因排序中top10%的基因的rank位置

# set.seed(333)
# cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, assign=TRUE)#选定合适的阈值
# auc_thr = cells_assignment$NEDDylation$aucThr$selected
# auc_thr

#提取结果
geneSet <- "Glycolysis"
AUCell_auc <- as.numeric(getAUC(cells_AUC)[geneSet, ])
#添加至metadata中
scMerge$AUCell_Glycolysis <- AUCell_auc
head(scMerge@meta.data)

## 2.2 Kla----
#cells_rankings <- AUCell_buildRankings(scMerge@assays$RNA@data, splitByBlocks=TRUE) # 与2.1相同

cells_AUC <- AUCell_calcAUC(upKla, cells_rankings, nCores = 10,
                            aucMaxRank = ceiling(0.1 * nrow(cells_rankings))) # 选择基因排序中top10%的基因的rank位置

# set.seed(333)
# cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, assign=TRUE)#选定合适的阈值 
# auc_thr = cells_assignment$NEDDylation$aucThr$selected
# auc_thr

#提取结果
geneSet <- "H3K18la"
AUCell_auc <- as.numeric(getAUC(cells_AUC)[geneSet, ])
#添加至metadata中
scMerge$AUCell_upKla <- AUCell_auc
head(scMerge@meta.data)

## 2.3 可视化比较----
df <- scMerge@meta.data
class(df$Cell_types)
table(df$Cell_types)

zzm60colors <- c('#4b6aa8','#3ca0cf','#c376a7','#ad98c3','#cea5c7',
                 '#53738c','#a5a9b0','#a78982','#696a6c','#92699e',
                 '#d69971','#df5734','#6c408e','#ac6894','#d4c2db',
                 '#537eb7','#83ab8e','#ece399','#405993','#cc7f73',
                 '#b95055','#d5bb72','#bc9a7f','#e0cfda','#d8a0c0',
                 '#e6b884','#b05545','#d69a55','#64a776','#cbdaa9',
                 '#efd2c9','#da6f6d','#ebb1a4','#a44e89','#a9c2cb',
                 '#b85292','#6d6fa0','#8d689d','#c8c7e1','#d25774',
                 '#c49abc','#927c9a','#3674a2','#9f8d89','#72567a',
                 '#63a3b8','#c4daec','#61bada','#b7deea','#e29eaf',
                 '#4490c4','#e6e2a3','#de8b36','#c4612f','#9a70a8',
                 '#76a2be','#408444','#c6adb0','#9d3b62','#2d3462')
mycolors <- c(pal_d3("category20")(20),
              zzm60colors)

ggviolin(df, x="Cell_types", y = "AUCell_upKla", width = 0.6, 
         color = "black",#轮廓颜色
         fill="Cell_types",#填充
         palette = mycolors,
         add = 'mean_sd',
         xlab = F, #不显示x轴的标签
         ylab = "AUCell score of upKla genes",
         bxp.errorbar=T,#显示误差条
         bxp.errorbar.width=0.5, #误差条大小
         size = .6, #箱型图边线的粗细
         outlier.shape=NA, #不显示outlier
         legend = "right") +
  stat_compare_means(method = "kruskal.test", label.x = 2, label.y = 0.12) +
  theme(axis.text.x = element_blank())

ggviolin(df, x="Cell_types", y = "AUCell_Glycolysis", width = 0.8, 
         color = "black",#轮廓颜色
         fill="Cell_types",#填充
         palette = mycolors,
         add = 'mean_sd',
         xlab = F, #不显示x轴的标签
         ylab = "AUCell score of glycolysis",
         bxp.errorbar=T,#显示误差条
         bxp.errorbar.width=0.5, #误差条大小
         size = .6, #箱型图边线的粗细
         outlier.shape=NA, #不显示outlier
         legend = "right") +
  stat_compare_means(method = "kruskal.test", label.x = 2, label.y = 0.35) +
  theme(axis.text.x = element_blank())


# 3.只在上皮细胞中进行AUCell评分----
rm(df, scMerge, cells_AUC, cells_rankings, AUCell_auc, geneSet)
## 3.1 加载44脚本中运行完的sc.epi结果----
sc.epi <- readRDS("5.scRNA/GSE189357_epi.rds")

## 3.2 Glycolysis AUCell评分----
cells_rankings <- AUCell_buildRankings(sc.epi@assays$RNA@data, splitByBlocks=TRUE) 

cells_AUC <- AUCell_calcAUC(glyco, cells_rankings, nCores = 10,
                            aucMaxRank = ceiling(0.1 * nrow(cells_rankings))) # 选择基因排序中top10%的基因的rank位置

# set.seed(333)
# cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, assign=TRUE) #选定合适的阈值
# auc_thr = cells_assignment$NEDDylation$aucThr$selected
# auc_thr

#提取结果
geneSet <- "Glycolysis"
AUCell_auc <- as.numeric(getAUC(cells_AUC)[geneSet, ])
#添加至metadata中
sc.epi$AUCell_Glycolysis <- AUCell_auc
head(sc.epi@meta.data)

## 3.3 Kla AUCell评分----
cells_rankings <- AUCell_buildRankings(sc.epi@assays$RNA@data, splitByBlocks=TRUE) # 与2.1相同

cells_AUC <- AUCell_calcAUC(upKla, cells_rankings, nCores = 10,
                            aucMaxRank = ceiling(0.1 * nrow(cells_rankings))) # 选择基因排序中top10%的基因的rank位置

# set.seed(333)
# cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, assign=TRUE)#选定合适的阈值 
# auc_thr = cells_assignment$NEDDylation$aucThr$selected
# auc_thr

#提取结果
geneSet <- "H3K18la"
AUCell_auc <- as.numeric(getAUC(cells_AUC)[geneSet, ])
#添加至metadata中
sc.epi$AUCell_upKla <- AUCell_auc
head(sc.epi@meta.data)

## 3.4 可视化比较AUCell评分----
df <- sc.epi@meta.data
class(df$TumorCell_pred)
table(df$TumorCell_pred)
# Normal epithelials        Tumor cells 
# 4601               5717
my_comparisons <- list(c("Tumor cells", "Normal epithelials"))
aucell_upKla_plot <- ggviolin(df, x="TumorCell_pred", y = "AUCell_upKla", width = 0.8, 
                              color = "black",#轮廓颜色
                              fill="TumorCell_pred",#填充
                              palette = mycolors,
                              add = 'mean_sd',
                              xlab = F, #不显示x轴的标签
                              ylab = "AUCell score of upKla genes",
                              bxp.errorbar=T,#显示误差条
                              bxp.errorbar.width=0.5, #误差条大小
                              size = .6, #箱型图边线的粗细
                              outlier.shape=NA, #不显示outlier
                              legend = "none") +
  scale_y_continuous(
    breaks = seq(0, 0.1, by = 0.02),   # 设置 y 轴的刻度
    labels = scales::label_number(accuracy = 0.02) # 保留小数点后两位
  ) +
  coord_cartesian(ylim = c(0, 0.1)) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label.y = 0.095) +
  theme(axis.title.x = element_blank())
aucell_upKla_plot
ggsave(filename = "5.scRNA/GSE189357_aucell_upKla_plot.pdf", aucell_upKla_plot,
       height = 4, width = 4, dpi = 300)

aucell_glyco_plot <- ggviolin(df, x="TumorCell_pred", y = "AUCell_Glycolysis", width = 0.8, 
                              color = "black",#轮廓颜色
                              fill="TumorCell_pred",#填充
                              palette = mycolors,
                              add = 'mean_sd',
                              xlab = F, #不显示x轴的标签
                              ylab = "AUCell score of glycolysis",
                              bxp.errorbar=T,#显示误差条
                              bxp.errorbar.width=0.5, #误差条大小
                              size = .6, #箱型图边线的粗细
                              outlier.shape=NA, #不显示outlier
                              legend = "none") +
  scale_y_continuous(
    breaks = seq(0, 0.3, by = 0.1),   # 设置 y 轴的刻度
    labels = scales::label_number(accuracy = 0.02) # 保留小数点后两位
  ) +
  coord_cartesian(ylim = c(0, 0.3)) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label.y = 0.28) +
  theme(axis.title.x = element_blank())
aucell_glyco_plot
ggsave(filename = "5.scRNA/GSE189357_aucell_glyco_plot.pdf", aucell_glyco_plot,
       height = 4, width = 4, dpi = 300)

#### UMAP图
dat<- data.frame(sc.epi@meta.data, 
                 sc.epi@reductions$umap@cell.embeddings,
                 sc.epi@active.ident)

class_avg <- dat %>%
  group_by(TumorCell_pred) %>%
  summarise(
    UMAP_1 = median(UMAP_1),
    UMAP_2 = median(UMAP_2)
  )
#也可以是细胞类型，可以自己设置
aucell_glyco_umap <- ggplot(dat, aes(UMAP_1, UMAP_2))  +
  geom_point(aes(colour  = AUCell_Glycolysis), size = .6) + 
  viridis::scale_color_viridis(option="A") +
  # ggrepel::geom_label_repel(aes(label = Cell_Type),
  #                           data = class_avg,
  #                           label.size = 0,
  #                           segment.color = NA)+
  theme(legend.position = "none") + 
  theme_classic() + 
  theme(panel.grid = element_blank()) +
  theme(axis.line = element_line(linewidth = .5),
        axis.text = element_text(size = 12)) +
  tidydr::theme_dr(xlength = 0.2, 
                   ylength = 0.2,
                   arrow = arrow(length = unit(0.1, "inches"),type = "closed")) + 
  ggtitle("Glycolysis score") +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2, hjust = 0.01),
        legend.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.6, face = "bold", size = 14))
aucell_glyco_umap
ggsave(filename = "5.scRNA/GSE189357_aucell_glyco_UMAPlot.pdf", aucell_glyco_umap,
       height = 4.5, width = 6, dpi = 300)

aucell_upKla_umap <- ggplot(dat, aes(UMAP_1, UMAP_2))  +
  geom_point(aes(colour  = AUCell_upKla), size = .6) + 
  viridis::scale_color_viridis(option="A") +
  # ggrepel::geom_label_repel(aes(label = Cell_Type),
  #                           data = class_avg,
  #                           label.size = 0,
  #                           segment.color = NA)+
  theme(legend.position = "none") + 
  theme_classic() + 
  theme(panel.grid = element_blank()) +
  theme(axis.line = element_line(linewidth = .5),
        axis.text = element_text(size = 12)) +
  tidydr::theme_dr(xlength = 0.2, 
                   ylength = 0.2,
                   arrow = arrow(length = unit(0.1, "inches"),type = "closed")) + 
  ggtitle("upKla score") +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2, hjust = 0.01),
        legend.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.6, face = "bold", size = 14))
aucell_upKla_umap
ggsave(filename = "5.scRNA/GSE189357_aucell_upKla_UMAPlot.pdf", aucell_upKla_umap,
       height = 4.5, width = 6, dpi = 300)


## 3.5 Glycolysis ssGSEA评分----
ssgsea_glyco <- gsva(as.matrix(sc.epi@assays$RNA$counts),
                     glyco,
                     method = "ssgsea", kcdf = "Poisson", abs.ranking = T)
ssgsea_glyco <- as.data.frame(t(ssgsea_glyco))

identical(rownames(ssgsea_glyco), rownames(sc.epi@meta.data))
sc.epi$ssGSEA_Glycolysis <- ssgsea_glyco

df <- sc.epi@meta.data
class(df$TumorCell_pred)
table(df$TumorCell_pred)
my_comparisons <- list(c("Tumor cells", "Normal epithelials"))

ggviolin(df, x="TumorCell_pred", y = "ssGSEA_Glycolysis", width = 0.8, 
         color = "black",#轮廓颜色
         fill="TumorCell_pred",#填充
         palette = "npg",
         add = 'mean_sd',
         xlab = F, #不显示x轴的标签
         ylab = "ssGSEA score of glycolysis",
         bxp.errorbar=T,#显示误差条
         bxp.errorbar.width=0.5, #误差条大小
         size = .6, #箱型图边线的粗细
         outlier.shape=NA, #不显示outlier
         legend = "right") +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label.x = 2, label.y = 1.2) +
  theme(axis.text.x = element_blank())

## 3.6 Kla ssGSEA评分----
ssgsea_upKla <- gsva(as.matrix(sc.epi@assays$RNA$counts),
                     upKla,
                     method = "ssgsea", kcdf = "Poisson", abs.ranking = T)
ssgsea_upKla <- as.data.frame(t(ssgsea_upKla))

identical(rownames(ssgsea_upKla), rownames(sc.epi@meta.data))
sc.epi$ssGSEA_upKla <- ssgsea_upKla

df <- sc.epi@meta.data
class(df$TumorCell_pred)
table(df$TumorCell_pred)
my_comparisons <- list(c("Tumor cells", "Normal epithelials"))

ggviolin(df, x="TumorCell_pred", y = "ssGSEA_upKla", width = 0.8, 
         color = "black",#轮廓颜色
         fill="TumorCell_pred",#填充
         palette = "npg",
         add = 'mean_sd',
         xlab = F, #不显示x轴的标签
         ylab = "ssGSEA score of upKla",
         bxp.errorbar=T,#显示误差条
         bxp.errorbar.width=0.5, #误差条大小
         size = .6, #箱型图边线的粗细
         outlier.shape=NA, #不显示outlier
         legend = "right") +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label.x = 2, label.y = 0.6) +
  theme(axis.text.x = element_blank())


# 4.上皮细胞CytoTRACE评分和GIS评分----
## 4.1 CytoTRACE----
#devtools::install_github("digitalcytometry/cytotrace2", subdir = "cytotrace2_r") #installing
library(CytoTRACE2) #loading

Idents(sc.epi) <- "TumorCell_pred"
cytotrace_score <- cytotrace2(sc.epi, species = "human", is_seurat = T, seed = 123)

df <- cytotrace_score@meta.data
class(df$TumorCell_pred)
table(df$TumorCell_pred)
my_comparisons <- list(c("Tumor cells", "Normal epithelials"))

ggviolin(df, x="TumorCell_pred", y = "CytoTRACE2_Score", width = 0.8, 
         color = "black",#轮廓颜色
         fill="TumorCell_pred",#填充
         palette = "npg",
         add = 'mean_sd',
         xlab = F, #不显示x轴的标签
         ylab = "CytoTRACE2 score",
         bxp.errorbar=T,#显示误差条
         bxp.errorbar.width=0.5, #误差条大小
         size = .6, #箱型图边线的粗细
         outlier.shape=NA, #不显示outlier
         legend = "right") +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label.x = 2, label.y = 1) +
  theme(axis.text.x = element_blank())

cor.test(cytotrace_score$CytoTRACE2_Score, cytotrace_score$ssGSEA_upKla)
cor.test(cytotrace_score$CytoTRACE2_Score, cytotrace_score$AUCell_upKla)
cor.test(cytotrace_score$CytoTRACE2_Score, cytotrace_score$ssGSEA_Glycolysis)
cor.test(cytotrace_score$CytoTRACE2_Score, cytotrace_score$AUCell_Glycolysis)
cor.test(cytotrace_score$AUCell_upKla, cytotrace_score$ssGSEA_upKla)

## 4.2 基因组不稳定性评分----
library(genomicInstability)
library(clusterProfiler)
library(org.Hs.eg.db)

count_matrix <- GetAssayData(sc.epi, layer = "counts")

gene.df <- bitr(rownames(count_matrix),
                fromType = "SYMBOL",
                toType = c("ENTREZID"),
                OrgDb = org.Hs.eg.db)
count_matrix <- count_matrix[gene.df$SYMBOL,]
rownames(count_matrix) <- gene.df$ENTREZID

cnv <- inferCNV(as.matrix(count_matrix))
cnv <- genomicInstabilityScore(cnv)

gis_score <- data.frame(row.names = names(cnv[['gis']]), GIS_score = cnv[['gis']])

identical(rownames(gis_score), rownames(sc.epi@meta.data))
sc.epi <- AddMetaData(sc.epi, metadata = gis_score)

df <- sc.epi@meta.data
class(df$TumorCell_pred)
table(df$TumorCell_pred)
my_comparisons <- list(c("Tumor cells", "Normal epithelials"))

GIscore_plot <- ggviolin(df, x="TumorCell_pred", y = "GIS_score", width = 0.8, 
                         color = "black",#轮廓颜色
                         fill="TumorCell_pred",#填充
                         palette = mycolors,
                         add = 'mean_sd',
                         #xlab = F, #不显示x轴的标签
                         ylab = "Genomic instability score", 
                         bxp.errorbar=T,#显示误差条
                         bxp.errorbar.width=0.5, #误差条大小
                         size = .6, #箱型图边线的粗细
                         outlier.shape=NA, #不显示outlier
                         legend = "none") +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label.x = 2, label.y = 1) +
  theme(axis.title.x = element_blank())
GIscore_plot
ggsave(filename = "5.scRNA/GSE189357_GIscore_plot.pdf", GIscore_plot,
       height = 4, width = 4, dpi = 300)

identical(rownames(sc.epi@meta.data), rownames(cytotrace_score@meta.data))
cor.test(sc.epi$GIS_score, sc.epi$ssGSEA_upKla)
cor.test(sc.epi$GIS_score, sc.epi$AUCell_upKla)
cor.test(sc.epi$GIS_score, sc.epi$AUCell_Glycolysis)
cor.test(sc.epi$GIS_score, sc.epi$ssGSEA_Glycolysis)
cor.test(sc.epi$GIS_score, cytotrace_score$CytoTRACE2_Score)


# 5.比较不同进展阶段的评分----
class(df$Sample_group)
table(df$Sample_group)
# AIS  IAC  MIA 
# 2821 6269 1228
df$Sample_group <- factor(df$Sample_group, levels = c("AIS", "MIA", "IAC"))
df$Sample_group2 <- ifelse(df$Sample_group == "AIS", "in situ", "invasive")
df$Sample_group2 <- factor(df$Sample_group2, levels = c("in situ", "invasive"))
table(df$Sample_group2)
table(df$Sample_group2, df$TumorCell_pred)
#               Normal epithelials Tumor cells
# in situ                2372         449
# invasive               2229        5268

my_comparisons <- list(c("MIA", "AIS"),
                       c("IAC", "AIS"),
                       c("IAC", "MIA"))
my_comparisons2 <- list(c("invasive", "in situ"))

ggviolin(df, x="Sample_group", y = "AUCell_Glycolysis", width = 0.8, 
         color = "black",#轮廓颜色
         fill="Sample_group",#填充
         palette = "npg",
         add = 'mean_sd',
         xlab = F, #不显示x轴的标签
         ylab = "AUCell score of glycolysis",
         bxp.errorbar=T,#显示误差条
         bxp.errorbar.width=0.5, #误差条大小
         size = .6, #箱型图边线的粗细
         outlier.shape=NA, #不显示outlier
         legend = "right") +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") +
  theme(axis.text.x = element_blank())

ggviolin(df, x="Sample_group", y = "AUCell_upKla", width = 0.8, 
         color = "black",#轮廓颜色
         fill="Sample_group",#填充
         palette = "npg",
         add = 'mean_sd',
         xlab = F, #不显示x轴的标签
         ylab = "AUCell score of upKla",
         bxp.errorbar=T,#显示误差条
         bxp.errorbar.width=0.5, #误差条大小
         size = .6, #箱型图边线的粗细
         outlier.shape=NA, #不显示outlier
         legend = "right") +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") +
  theme(axis.text.x = element_blank())

glyco_vlnplot <- ggviolin(df, x="Sample_group2", y = "AUCell_Glycolysis", width = 0.8, 
                          color = "black",#轮廓颜色
                          fill="Sample_group2",#填充
                          palette = c("#ADA4CB", "#4B3C97"),
                          add = 'mean_sd',
                          xlab = F, #不显示x轴的标签
                          ylab = "AUCell score of glycolysis",
                          bxp.errorbar=T,#显示误差条
                          bxp.errorbar.width=0.5, #误差条大小
                          size = .6, #箱型图边线的粗细
                          outlier.shape=NA, #不显示outlier
                          legend = "none") +
  scale_y_continuous(
    breaks = seq(0, 0.3, by = 0.1),   # 设置 y 轴的刻度
    labels = scales::label_number(accuracy = 0.01) # 保留小数点后两位
  ) +
  coord_cartesian(ylim = c(0, 0.3)) +
  stat_compare_means(comparisons = my_comparisons2, method = "wilcox.test", label.y = 0.28) +
  theme(axis.title.x = element_blank())
glyco_vlnplot
ggsave(filename = "5.scRNA/GSE189357_aucell_glyco_vlnPlot.pdf", glyco_vlnplot,
       height = 4, width = 4, dpi = 300)

upkla_vlnplot <- ggviolin(df, x="Sample_group2", y = "AUCell_upKla", width = 0.8, 
                          color = "black",#轮廓颜色
                          fill="Sample_group2",#填充
                          palette = c("#ADA4CB", "#4B3C97"),
                          add = 'mean_sd',
                          xlab = F, #不显示x轴的标签
                          ylab = "AUCell score of upKla genes",
                          bxp.errorbar=T,#显示误差条
                          bxp.errorbar.width=0.5, #误差条大小
                          size = .6, #箱型图边线的粗细
                          outlier.shape=NA, #不显示outlier
                          legend = "none") +
  scale_y_continuous(
    breaks = seq(0, 0.1, by = 0.02),   # 设置 y 轴的刻度
    labels = scales::label_number(accuracy = 0.02) # 保留小数点后两位
  ) +
  coord_cartesian(ylim = c(0, 0.1)) +
  stat_compare_means(comparisons = my_comparisons2, method = "wilcox.test", label.y = 0.095) +
  theme(axis.title.x = element_blank())
upkla_vlnplot
ggsave(filename = "5.scRNA/GSE189357_aucell_upKla_vlnPlot.pdf", upkla_vlnplot,
       height = 4, width = 4, dpi = 300)


# 6.保存结果----
saveRDS(sc.epi, file = "5.scRNA/GSE189357_scEpi.rds")