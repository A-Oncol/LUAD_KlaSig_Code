######目的：根据SpaCET定义的区域重新比较Cottrazm的CNV评分
######作者：申奥
######日期：2024-10-24
######Server: yjx@xiyoucloud.top:34022
######Version: R 4.3.3


# 0.环境设置----
rm(list = ls())
options(stringAsFactors = F)

library(ggplot2)
library(ggpubr)
library(ggsci)
library(ggsignif)
library(ggdist)


# 1.P10----
## 1.1 加载CNV评分----
load("2.Cottrazm_outputs/PT10_T1/PT10_T1_TumorST_Boudary.RData")
cnv_score <- data.frame(barcodes = rownames(TumorST@meta.data),
                        cnv_score = TumorST@meta.data$cnv_score)

## 1.2 加载SpaCET定义的区域----
load("13.SpaCET/P10_T1/P10_T1_SpaCET_obj.RData")
location <- data.frame(barcodes = SpaCET_obj@input[["spotCoordinates"]][["barcode"]],
                       location = as.character(SpaCET_obj@results[["CCI"]][["interface"]]))

## 1.3 合并画图----
identical(cnv_score$barcodes, location$barcodes)
merge_df <- inner_join(location, cnv_score, by = "barcodes")
table(merge_df$location)
# Interface    Stroma     Tumor 
# 464       809       927
class(merge_df$location)
merge_df$location <- factor(merge_df$location, levels = c("Stroma", "Interface", "Tumor"))

comparison <- list(c("Stroma", "Interface"),
                   c("Stroma", "Tumor"),
                   c("Interface", "Tumor"))
cnvscore_plot <- ggplot(merge_df, aes(x = location, y = cnv_score, fill = location)) +
  geom_boxplot(fill = c("#fdb462", "#1f78b4", "#cb181d"), width = 0.6, alpha = 0.9) +
  theme_pubr() +
  theme(axis.title = element_text(size = 14)) +
  theme(axis.text = element_text(size = 12)) +
  theme(legend.position = "none") +
  ylab("CNV scores") + xlab(NULL) +
  stat_compare_means(comparisons = comparison, method = "wilcox.test", label = "p.format")
cnvscore_plot
ggsave(cnvscore_plot, filename = "13.SpaCET/P10_T1/P10_CNVscore.pdf", height = 4, width = 4, dpi = 300)


# 2.P15----
rm(list = ls())
## 2.1 加载CNV评分----
load("2.Cottrazm_outputs/PT15_T1/PT15_T1_TumorST_Boudary.RData")
cnv_score <- data.frame(barcodes = rownames(TumorST@meta.data),
                        cnv_score = TumorST@meta.data$cnv_score)

## 2.2 加载SpaCET定义的区域----
load("13.SpaCET/P15_T1/P15_T1_SpaCET_obj.RData")
location <- data.frame(barcodes = SpaCET_obj@input[["spotCoordinates"]][["barcode"]],
                       location = as.character(SpaCET_obj@results[["CCI"]][["interface"]]))

## 2.3 合并画图----
identical(cnv_score$barcodes, location$barcodes)
merge_df <- inner_join(location, cnv_score, by = "barcodes")
table(merge_df$location)
# Interface    Stroma     Tumor 
# 897       287      1440 
class(merge_df$location)
merge_df$location <- factor(merge_df$location, levels = c("Stroma", "Interface", "Tumor"))

comparison <- list(c("Stroma", "Interface"),
                   c("Stroma", "Tumor"),
                   c("Interface", "Tumor"))
cnvscore_plot <- ggplot(merge_df, aes(x = location, y = cnv_score, fill = location)) +
  geom_boxplot(fill = c("#fdb462", "#1f78b4", "#cb181d"), width = 0.6, alpha = 0.9) +
  theme_pubr() +
  theme(axis.title = element_text(size = 14)) +
  theme(axis.text = element_text(size = 12)) +
  theme(legend.position = "none") +
  ylab("CNV scores") + xlab(NULL) +
  stat_compare_means(comparisons = comparison, method = "wilcox.test", label = "p.format")
cnvscore_plot
ggsave(cnvscore_plot, filename = "13.SpaCET/P15_T1/P15_CNVscore.pdf", height = 4, width = 4, dpi = 300)


# 3.P16----
rm(list = ls())
## 3.1 加载CNV评分----
load("2.Cottrazm_outputs/PT16_T1/PT16_T1_TumorST_Boudary.RData")
cnv_score <- data.frame(barcodes = rownames(TumorST@meta.data),
                        cnv_score = TumorST@meta.data$cnv_score)

## 3.2 加载SpaCET定义的区域----
load("13.SpaCET/P16_T1/P16_T1_SpaCET_obj.RData")
location <- data.frame(barcodes = SpaCET_obj@input[["spotCoordinates"]][["barcode"]],
                       location = as.character(SpaCET_obj@results[["CCI"]][["interface"]]))

## 3.3 合并画图----
identical(cnv_score$barcodes, location$barcodes)
merge_df <- inner_join(location, cnv_score, by = "barcodes")
table(merge_df$location)
# Interface    Stroma     Tumor 
# 424      1135       703
class(merge_df$location)
merge_df$location <- factor(merge_df$location, levels = c("Stroma", "Interface", "Tumor"))

comparison <- list(c("Stroma", "Interface"),
                   c("Stroma", "Tumor"),
                   c("Interface", "Tumor"))
cnvscore_plot <- ggplot(merge_df, aes(x = location, y = cnv_score, fill = location)) +
  geom_boxplot(fill = c("#fdb462", "#1f78b4", "#cb181d"), width = 0.6, alpha = 0.9) +
  theme_pubr() +
  theme(axis.title = element_text(size = 14)) +
  theme(axis.text = element_text(size = 12)) +
  theme(legend.position = "none") +
  ylab("CNV scores") + xlab(NULL) +
  stat_compare_means(comparisons = comparison, method = "wilcox.test", label = "p.format")
cnvscore_plot
ggsave(cnvscore_plot, filename = "13.SpaCET/P16_T1/P16_CNVscore.pdf", height = 4, width = 4, dpi = 300)


# 4.P24----
rm(list = ls())
## 4.1 加载CNV评分----
load("2.Cottrazm_outputs/PT24_T1/PT24_T1_TumorST_Boudary.RData")
cnv_score <- data.frame(barcodes = rownames(TumorST@meta.data),
                        cnv_score = TumorST@meta.data$cnv_score)

## 4.2 加载SpaCET定义的区域----
load("13.SpaCET/P24_T1/P24_T1_SpaCET_obj.RData")
location <- data.frame(barcodes = SpaCET_obj@input[["spotCoordinates"]][["barcode"]],
                       location = as.character(SpaCET_obj@results[["CCI"]][["interface"]]))

## 4.3 合并画图----
identical(cnv_score$barcodes, location$barcodes)
merge_df <- inner_join(location, cnv_score, by = "barcodes")
table(merge_df$location)
# Interface    Stroma     Tumor 
# 1279       483      1664
class(merge_df$location)
merge_df$location <- factor(merge_df$location, levels = c("Stroma", "Interface", "Tumor"))

comparison <- list(c("Stroma", "Interface"),
                   c("Stroma", "Tumor"),
                   c("Interface", "Tumor"))
cnvscore_plot <- ggplot(merge_df, aes(x = location, y = cnv_score, fill = location)) +
  geom_boxplot(fill = c("#fdb462", "#1f78b4", "#cb181d"), width = 0.6, alpha = 0.9) +
  theme_pubr() +
  theme(axis.title = element_text(size = 14)) +
  theme(axis.text = element_text(size = 12)) +
  theme(legend.position = "none") +
  ylab("CNV scores") + xlab(NULL) +
  stat_compare_means(comparisons = comparison, method = "wilcox.test", label = "p.format")
cnvscore_plot
ggsave(cnvscore_plot, filename = "13.SpaCET/P24_T1/P24_CNVscore.pdf", height = 4, width = 4, dpi = 300)


# 5.P25----
rm(list = ls())
## 5.1 加载CNV评分----
load("2.Cottrazm_outputs/PT25_T1/PT25_T1_TumorST_Boudary.RData")
cnv_score <- data.frame(barcodes = rownames(TumorST@meta.data),
                        cnv_score = TumorST@meta.data$cnv_score)

## 5.2 加载SpaCET定义的区域----
load("13.SpaCET/P25_T1/P25_T1_SpaCET_obj.RData")
location <- data.frame(barcodes = SpaCET_obj@input[["spotCoordinates"]][["barcode"]],
                       location = as.character(SpaCET_obj@results[["CCI"]][["interface"]]))

## 5.3 合并画图----
identical(cnv_score$barcodes, location$barcodes)
merge_df <- inner_join(location, cnv_score, by = "barcodes")
table(merge_df$location)
# Interface    Stroma     Tumor 
# 290       311       707
class(merge_df$location)
merge_df$location <- factor(merge_df$location, levels = c("Stroma", "Interface", "Tumor"))

comparison <- list(c("Stroma", "Interface"),
                   c("Stroma", "Tumor"),
                   c("Interface", "Tumor"))
cnvscore_plot <- ggplot(merge_df, aes(x = location, y = cnv_score, fill = location)) +
  geom_boxplot(fill = c("#fdb462", "#1f78b4", "#cb181d"), width = 0.6, alpha = 0.9) +
  theme_pubr() +
  theme(axis.title = element_text(size = 14)) +
  theme(axis.text = element_text(size = 12)) +
  theme(legend.position = "none") +
  ylab("CNV scores") + xlab(NULL) +
  stat_compare_means(comparisons = comparison, method = "wilcox.test", label = "p.format")
cnvscore_plot
ggsave(cnvscore_plot, filename = "13.SpaCET/P25_T1/P25_CNVscore.pdf", height = 4, width = 4, dpi = 300)
