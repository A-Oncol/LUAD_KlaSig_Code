######目的：空转基因集评分可视化
######作者：申奥
######日期：2024-10-25
######Server: yjx@xiyoucloud.top:34022
######Version: R 4.3.3


# 0.环境设置----
rm(list = ls())
options(stringAsFactors = F)

library(Seurat)
library(GSVA)
library(tidyverse)


# 1.加载基因集----
load("15.diff/upKla_genes.RData")
upKla <- list(upKla_genes)
names(upKla) <- "upKla"


# 2.P10----
## 2.1 读取----
st_obj = Load10X_Spatial(data.dir = "1.data/SP/PT10_T1/", slice = "PT10_T1")
st_obj <- SCTransform(st_obj, assay = "Spatial")

## 2.2 upKla评分----
expr <- as.matrix(st_obj@assays$SCT@counts)
upKla_score <- gsva(expr,
                    upKla,
                    kcdf = "Poisson",
                    method = "ssgsea",
                    abs.ranking = T, parallel.sz = 48)
upKla_score <- as.data.frame(t(upKla_score))
save(upKla_score, file = "14.Location_heterogeneity/metabolism/P10_T1_upKla_ssGSEA_score.RData")

identical(rownames(st_obj@meta.data), rownames(upKla_score))
st_obj@meta.data$upKla_score <- upKla_score$upKla

## 2.3 把之前scMetabolism计算的糖酵解评分也添加到meta.data中----
load("14.Location_heterogeneity/metabolism/P10_T1_KEGG_metabolism_score.RData")
glycolysis_score <- as.data.frame(t(meta.score["Glycolysis / Gluconeogenesis",]))
colnames(glycolysis_score) <- "Glycolysis_score"
rownames(glycolysis_score) <- str_replace_all(rownames(glycolysis_score), "\\.", "-")
save(glycolysis_score, file = "14.Location_heterogeneity/metabolism/P10_T1_glycolysis_ssGSEA_score.RData")

identical(rownames(st_obj@meta.data), rownames(glycolysis_score))
st_obj@meta.data$glycolysis_score <- glycolysis_score$Glycolysis_score
save(st_obj, file = "14.Location_heterogeneity/metabolism/P10_T1_seuObj_glyco_upKla.RData")

## 2.4 画图----
p1 <- SpatialFeaturePlot(st_obj, features = "glycolysis_score") +
  theme(legend.position = "right",
        legend.title = element_blank(),
        legend.key.size = unit(0.4, "cm")) +
  ggtitle("P10")
p1

p2 <- SpatialFeaturePlot(st_obj, features = "upKla_score") +
  theme(legend.position = "right",
        legend.title = element_blank(),
        legend.key.size = unit(0.4, "cm")) +
  ggtitle("P10")
p2


# 3.P15----
rm(expr, glycolysis_score, meta.score, st_obj, upKla_score)
## 3.1 读取----
st_obj = Load10X_Spatial(data.dir = "1.data/SP/PT15_T1/", slice = "PT15_T1")
st_obj <- SCTransform(st_obj, assay = "Spatial")

## 3.2 upKla评分----
expr <- as.matrix(st_obj@assays$SCT@counts)
upKla_score <- gsva(expr,
                    upKla,
                    kcdf = "Poisson",
                    method = "ssgsea",
                    abs.ranking = T, parallel.sz = 48)
upKla_score <- as.data.frame(t(upKla_score))
save(upKla_score, file = "14.Location_heterogeneity/metabolism/P15_T1_upKla_ssGSEA_score.RData")

identical(rownames(st_obj@meta.data), rownames(upKla_score))
st_obj@meta.data$upKla_score <- upKla_score$upKla

## 3.3 把之前scMetabolism计算的糖酵解评分也添加到meta.data中----
load("14.Location_heterogeneity/metabolism/P15_T1_KEGG_metabolism_score.RData")
glycolysis_score <- as.data.frame(t(meta.score["Glycolysis / Gluconeogenesis",]))
colnames(glycolysis_score) <- "Glycolysis_score"
rownames(glycolysis_score) <- str_replace_all(rownames(glycolysis_score), "\\.", "-")
save(glycolysis_score, file = "14.Location_heterogeneity/metabolism/P15_T1_glycolysis_ssGSEA_score.RData")

identical(rownames(st_obj@meta.data), rownames(glycolysis_score))
st_obj@meta.data$glycolysis_score <- glycolysis_score$Glycolysis_score
save(st_obj, file = "14.Location_heterogeneity/metabolism/P15_T1_seuObj_glyco_upKla.RData")

## 3.4 画图----
p3 <- SpatialFeaturePlot(st_obj, features = "glycolysis_score") +
  theme(legend.position = "right",
        legend.title = element_blank(),
        legend.key.size = unit(0.4, "cm")) +
  ggtitle("P15")
p3

p4 <- SpatialFeaturePlot(st_obj, features = "upKla_score") +
  theme(legend.position = "right",
        legend.title = element_blank(),
        legend.key.size = unit(0.4, "cm")) +
  ggtitle("P15")
p4


# 4.P16----
rm(expr, glycolysis_score, meta.score, st_obj, upKla_score)
## 4.1 读取----
st_obj = Load10X_Spatial(data.dir = "1.data/SP/PT16_T1/", slice = "PT16_T1")
st_obj <- SCTransform(st_obj, assay = "Spatial")

## 4.2 upKla评分----
expr <- as.matrix(st_obj@assays$SCT@counts)
upKla_score <- gsva(expr,
                    upKla,
                    kcdf = "Poisson",
                    method = "ssgsea",
                    abs.ranking = T, parallel.sz = 48)
upKla_score <- as.data.frame(t(upKla_score))
save(upKla_score, file = "14.Location_heterogeneity/metabolism/P16_T1_upKla_ssGSEA_score.RData")

identical(rownames(st_obj@meta.data), rownames(upKla_score))
st_obj@meta.data$upKla_score <- upKla_score$upKla

## 4.3 把之前scMetabolism计算的糖酵解评分也添加到meta.data中----
load("14.Location_heterogeneity/metabolism/P16_T1_KEGG_metabolism_score.RData")
glycolysis_score <- as.data.frame(t(meta.score["Glycolysis / Gluconeogenesis",]))
colnames(glycolysis_score) <- "Glycolysis_score"
rownames(glycolysis_score) <- str_replace_all(rownames(glycolysis_score), "\\.", "-")
save(glycolysis_score, file = "14.Location_heterogeneity/metabolism/P16_T1_glycolysis_ssGSEA_score.RData")

identical(rownames(st_obj@meta.data), rownames(glycolysis_score))
st_obj@meta.data$glycolysis_score <- glycolysis_score$Glycolysis_score
save(st_obj, file = "14.Location_heterogeneity/metabolism/P16_T1_seuObj_glyco_upKla.RData")

## 4.4 画图----
p5 <- SpatialFeaturePlot(st_obj, features = "glycolysis_score") +
  theme(legend.position = "right",
        legend.title = element_blank(),
        legend.key.size = unit(0.4, "cm")) +
  ggtitle("P16")
p5

p6 <- SpatialFeaturePlot(st_obj, features = "upKla_score") +
  theme(legend.position = "right",
        legend.title = element_blank(),
        legend.key.size = unit(0.4, "cm")) +
  ggtitle("P16")
p6


# 5.P24----
rm(expr, glycolysis_score, meta.score, st_obj, upKla_score)
## 5.1 读取----
st_obj = Load10X_Spatial(data.dir = "1.data/SP/PT24_T1/", slice = "PT24_T1")
st_obj <- SCTransform(st_obj, assay = "Spatial")

## 5.2 upKla评分----
expr <- as.matrix(st_obj@assays$SCT@counts)
upKla_score <- gsva(expr,
                    upKla,
                    kcdf = "Poisson",
                    method = "ssgsea",
                    abs.ranking = T, parallel.sz = 48)
upKla_score <- as.data.frame(t(upKla_score))
save(upKla_score, file = "14.Location_heterogeneity/metabolism/P24_T1_upKla_ssGSEA_score.RData")

identical(rownames(st_obj@meta.data), rownames(upKla_score))
st_obj@meta.data$upKla_score <- upKla_score$upKla

## 5.3 把之前scMetabolism计算的糖酵解评分也添加到meta.data中----
load("14.Location_heterogeneity/metabolism/P24_T1_KEGG_metabolism_score.RData")
glycolysis_score <- as.data.frame(t(meta.score["Glycolysis / Gluconeogenesis",]))
colnames(glycolysis_score) <- "Glycolysis_score"
rownames(glycolysis_score) <- str_replace_all(rownames(glycolysis_score), "\\.", "-")
save(glycolysis_score, file = "14.Location_heterogeneity/metabolism/P24_T1_glycolysis_ssGSEA_score.RData")

identical(rownames(st_obj@meta.data), rownames(glycolysis_score))
st_obj@meta.data$glycolysis_score <- glycolysis_score$Glycolysis_score
save(st_obj, file = "14.Location_heterogeneity/metabolism/P24_T1_seuObj_glyco_upKla.RData")

## 5.4 画图----
p7 <- SpatialFeaturePlot(st_obj, features = "glycolysis_score") +
  theme(legend.position = "right",
        legend.title = element_blank(),
        legend.key.size = unit(0.4, "cm")) +
  ggtitle("P24")
p7

p8 <- SpatialFeaturePlot(st_obj, features = "upKla_score") +
  theme(legend.position = "right",
        legend.title = element_blank(),
        legend.key.size = unit(0.4, "cm")) +
  ggtitle("P24")
p8


# 6.P25----
rm(expr, glycolysis_score, meta.score, st_obj, upKla_score)
## 6.1 读取----
st_obj = Load10X_Spatial(data.dir = "1.data/SP/PT25_T1/", slice = "PT25_T1")
st_obj <- SCTransform(st_obj, assay = "Spatial")

## 6.2 upKla评分----
expr <- as.matrix(st_obj@assays$SCT@counts)
upKla_score <- gsva(expr,
                    upKla,
                    kcdf = "Poisson",
                    method = "ssgsea",
                    abs.ranking = T, parallel.sz = 48)
upKla_score <- as.data.frame(t(upKla_score))
save(upKla_score, file = "14.Location_heterogeneity/metabolism/P25_T1_upKla_ssGSEA_score.RData")

identical(rownames(st_obj@meta.data), rownames(upKla_score))
st_obj@meta.data$upKla_score <- upKla_score$upKla

## 6.3 把之前scMetabolism计算的糖酵解评分也添加到meta.data中----
load("14.Location_heterogeneity/metabolism/P25_T1_KEGG_metabolism_score.RData")
glycolysis_score <- as.data.frame(t(meta.score["Glycolysis / Gluconeogenesis",]))
colnames(glycolysis_score) <- "Glycolysis_score"
rownames(glycolysis_score) <- str_replace_all(rownames(glycolysis_score), "\\.", "-")
save(glycolysis_score, file = "14.Location_heterogeneity/metabolism/P25_T1_glycolysis_ssGSEA_score.RData")

identical(rownames(st_obj@meta.data), rownames(glycolysis_score))
st_obj@meta.data$glycolysis_score <- glycolysis_score$Glycolysis_score
save(st_obj, file = "14.Location_heterogeneity/metabolism/P25_T1_seuObj_glyco_upKla.RData")

## 6.4 画图----
p9 <- SpatialFeaturePlot(st_obj, features = "glycolysis_score") +
  theme(legend.position = "right",
        legend.title = element_blank(),
        legend.key.size = unit(0.4, "cm")) +
  ggtitle("P25")
p9

p10 <- SpatialFeaturePlot(st_obj, features = "upKla_score") +
  theme(legend.position = "right",
        legend.title = element_blank(),
        legend.key.size = unit(0.4, "cm")) +
  ggtitle("P25")
p10


# 7.拼图----
library(patchwork)

## 7.1 glycolysis----
glyco_plot <- p1 + p3 + p5 + p7 + p9 + plot_layout(ncol = 2)
glyco_plot
ggsave(filename = "14.Location_heterogeneity/metabolism/glycolysis_score.pdf", glyco_plot,
       height = 4, width = 4, dpi = 300)

## 7.2 upKla----
upkla_plot <- p2 + p4 + p6 + p8 + p10 + plot_layout(ncol = 5)
upkla_plot
ggsave(filename = "14.Location_heterogeneity/metabolism/upKla_score.pdf", upkla_plot,
       height = 4, width = 12, dpi = 300)

## 7.3 保存小图----
save(p1, p3, p5, p7, p9, file = "14.Location_heterogeneity/metabolism/glycolysis_plots.RData")
save(p2, p4, p6, p8, p10, file = "14.Location_heterogeneity/metabolism/upKla_plots.RData")


# 8.分析glycolysis与upKla的相关性----
rm(list = ls())

## 8.1 P10----
load("14.Location_heterogeneity/metabolism/P10_T1_glycolysis_ssGSEA_score.RData")
load("14.Location_heterogeneity/metabolism/P10_T1_upKla_ssGSEA_score.RData")
load("14.Location_heterogeneity/P10_T1_LocInfo.RData")

identical(rownames(loc_df), rownames(glycolysis_score))
identical(rownames(loc_df), rownames(upKla_score))

P10_df <- cbind(loc_df, glycolysis_score, upKla_score)

### All
glyco_kla_cor <- cor(P10_df$Glycolysis_score, P10_df$upKla, method = "pearson")
glyco_kla_cor_P <- cor.test(P10_df$Glycolysis_score, P10_df$upKla, method = "pearson")
glyco_kla_cor_Pvalue <- glyco_kla_cor_P$p.value

### Tumor
P10_tumor_df <- P10_df %>% 
  filter(Location == "Tumor")

tumor_glyco_kla_cor <- cor(P10_tumor_df$Glycolysis_score, P10_tumor_df$upKla, method = "pearson")
tumor_glyco_kla_cor_P <- cor.test(P10_tumor_df$Glycolysis_score, P10_tumor_df$upKla, method = "pearson")
tumor_glyco_kla_cor_Pvalue <- tumor_glyco_kla_cor_P$p.value

### Interface
P10_interface_df <- P10_df %>% 
  filter(Location == "Interface")

interface_glyco_kla_cor <- cor(P10_interface_df$Glycolysis_score, P10_interface_df$upKla, method = "pearson")
interface_glyco_kla_cor_P <- cor.test(P10_interface_df$Glycolysis_score, P10_interface_df$upKla, method = "pearson")
interface_glyco_kla_cor_Pvalue <- interface_glyco_kla_cor_P$p.value

### Stroma
P10_stroma_df <- P10_df %>% 
  filter(Location == "Stroma")

stroma_glyco_kla_cor <- cor(P10_stroma_df$Glycolysis_score, P10_stroma_df$upKla, method = "pearson")
stroma_glyco_kla_cor_P <- cor.test(P10_stroma_df$Glycolysis_score, P10_stroma_df$upKla, method = "pearson")
stroma_glyco_kla_cor_Pvalue <- stroma_glyco_kla_cor_P$p.value

### 汇总
P10_corR <- data.frame(ID = "P10",
                       nMal = stroma_glyco_kla_cor,
                       Bdy = interface_glyco_kla_cor,
                       Mal = tumor_glyco_kla_cor,
                       All = glyco_kla_cor)

P10_corP <- data.frame(ID = "P10",
                       nMal = stroma_glyco_kla_cor_Pvalue,
                       Bdy = interface_glyco_kla_cor_Pvalue,
                       Mal = tumor_glyco_kla_cor_Pvalue,
                       All = glyco_kla_cor_Pvalue)

## 8.2 P15----
rm(list = setdiff(ls(), c("P10_corR", "P10_corP")))

load("14.Location_heterogeneity/metabolism/P15_T1_glycolysis_ssGSEA_score.RData")
load("14.Location_heterogeneity/metabolism/P15_T1_upKla_ssGSEA_score.RData")
load("14.Location_heterogeneity/P15_T1_LocInfo.RData")

identical(rownames(loc_df), rownames(glycolysis_score))
identical(rownames(loc_df), rownames(upKla_score))

P15_df <- cbind(loc_df, glycolysis_score, upKla_score)

### All
glyco_kla_cor <- cor(P15_df$Glycolysis_score, P15_df$upKla, method = "pearson")
glyco_kla_cor_P <- cor.test(P15_df$Glycolysis_score, P15_df$upKla, method = "pearson")
glyco_kla_cor_Pvalue <- glyco_kla_cor_P$p.value

### Tumor
P15_tumor_df <- P15_df %>% 
  filter(Location == "Tumor")

tumor_glyco_kla_cor <- cor(P15_tumor_df$Glycolysis_score, P15_tumor_df$upKla, method = "pearson")
tumor_glyco_kla_cor_P <- cor.test(P15_tumor_df$Glycolysis_score, P15_tumor_df$upKla, method = "pearson")
tumor_glyco_kla_cor_Pvalue <- tumor_glyco_kla_cor_P$p.value

### Interface
P15_interface_df <- P15_df %>% 
  filter(Location == "Interface")

interface_glyco_kla_cor <- cor(P15_interface_df$Glycolysis_score, P15_interface_df$upKla, method = "pearson")
interface_glyco_kla_cor_P <- cor.test(P15_interface_df$Glycolysis_score, P15_interface_df$upKla, method = "pearson")
interface_glyco_kla_cor_Pvalue <- interface_glyco_kla_cor_P$p.value

### Stroma
P15_stroma_df <- P15_df %>% 
  filter(Location == "Stroma")

stroma_glyco_kla_cor <- cor(P15_stroma_df$Glycolysis_score, P15_stroma_df$upKla, method = "pearson")
stroma_glyco_kla_cor_P <- cor.test(P15_stroma_df$Glycolysis_score, P15_stroma_df$upKla, method = "pearson")
stroma_glyco_kla_cor_Pvalue <- stroma_glyco_kla_cor_P$p.value

### 汇总
P15_corR <- data.frame(ID = "P15",
                       nMal = stroma_glyco_kla_cor,
                       Bdy = interface_glyco_kla_cor,
                       Mal = tumor_glyco_kla_cor,
                       All = glyco_kla_cor)

P15_corP <- data.frame(ID = "P15",
                       nMal = stroma_glyco_kla_cor_Pvalue,
                       Bdy = interface_glyco_kla_cor_Pvalue,
                       Mal = tumor_glyco_kla_cor_Pvalue,
                       All = glyco_kla_cor_Pvalue)

## 8.3 P16----
rm(list = setdiff(ls(), c("P10_corR", "P10_corP", "P15_corR", "P15_corP")))

load("14.Location_heterogeneity/metabolism/P16_T1_glycolysis_ssGSEA_score.RData")
load("14.Location_heterogeneity/metabolism/P16_T1_upKla_ssGSEA_score.RData")
load("14.Location_heterogeneity/P16_T1_LocInfo.RData")

identical(rownames(loc_df), rownames(glycolysis_score))
identical(rownames(loc_df), rownames(upKla_score))

P16_df <- cbind(loc_df, glycolysis_score, upKla_score)

### All
glyco_kla_cor <- cor(P16_df$Glycolysis_score, P16_df$upKla, method = "pearson")
glyco_kla_cor_P <- cor.test(P16_df$Glycolysis_score, P16_df$upKla, method = "pearson")
glyco_kla_cor_Pvalue <- glyco_kla_cor_P$p.value

### Tumor
P16_tumor_df <- P16_df %>% 
  filter(Location == "Tumor")

tumor_glyco_kla_cor <- cor(P16_tumor_df$Glycolysis_score, P16_tumor_df$upKla, method = "pearson")
tumor_glyco_kla_cor_P <- cor.test(P16_tumor_df$Glycolysis_score, P16_tumor_df$upKla, method = "pearson")
tumor_glyco_kla_cor_Pvalue <- tumor_glyco_kla_cor_P$p.value

### Interface
P16_interface_df <- P16_df %>% 
  filter(Location == "Interface")

interface_glyco_kla_cor <- cor(P16_interface_df$Glycolysis_score, P16_interface_df$upKla, method = "pearson")
interface_glyco_kla_cor_P <- cor.test(P16_interface_df$Glycolysis_score, P16_interface_df$upKla, method = "pearson")
interface_glyco_kla_cor_Pvalue <- interface_glyco_kla_cor_P$p.value

### Stroma
P16_stroma_df <- P16_df %>% 
  filter(Location == "Stroma")

stroma_glyco_kla_cor <- cor(P16_stroma_df$Glycolysis_score, P16_stroma_df$upKla, method = "pearson")
stroma_glyco_kla_cor_P <- cor.test(P16_stroma_df$Glycolysis_score, P16_stroma_df$upKla, method = "pearson")
stroma_glyco_kla_cor_Pvalue <- stroma_glyco_kla_cor_P$p.value

### 汇总
P16_corR <- data.frame(ID = "P16",
                       nMal = stroma_glyco_kla_cor,
                       Bdy = interface_glyco_kla_cor,
                       Mal = tumor_glyco_kla_cor,
                       All = glyco_kla_cor)

P16_corP <- data.frame(ID = "P16",
                       nMal = stroma_glyco_kla_cor_Pvalue,
                       Bdy = interface_glyco_kla_cor_Pvalue,
                       Mal = tumor_glyco_kla_cor_Pvalue,
                       All = glyco_kla_cor_Pvalue)

## 8.4 P24----
rm(list = setdiff(ls(), c("P10_corR", "P10_corP", "P15_corR", "P15_corP", "P16_corR", "P16_corP")))

load("14.Location_heterogeneity/metabolism/P24_T1_glycolysis_ssGSEA_score.RData")
load("14.Location_heterogeneity/metabolism/P24_T1_upKla_ssGSEA_score.RData")
load("14.Location_heterogeneity/P24_T1_LocInfo.RData")

identical(rownames(loc_df), rownames(glycolysis_score))
identical(rownames(loc_df), rownames(upKla_score))

P24_df <- cbind(loc_df, glycolysis_score, upKla_score)

### All
glyco_kla_cor <- cor(P24_df$Glycolysis_score, P24_df$upKla, method = "pearson")
glyco_kla_cor_P <- cor.test(P24_df$Glycolysis_score, P24_df$upKla, method = "pearson")
glyco_kla_cor_Pvalue <- glyco_kla_cor_P$p.value

### Tumor
P24_tumor_df <- P24_df %>% 
  filter(Location == "Tumor")

tumor_glyco_kla_cor <- cor(P24_tumor_df$Glycolysis_score, P24_tumor_df$upKla, method = "pearson")
tumor_glyco_kla_cor_P <- cor.test(P24_tumor_df$Glycolysis_score, P24_tumor_df$upKla, method = "pearson")
tumor_glyco_kla_cor_Pvalue <- tumor_glyco_kla_cor_P$p.value

### Interface
P24_interface_df <- P24_df %>% 
  filter(Location == "Interface")

interface_glyco_kla_cor <- cor(P24_interface_df$Glycolysis_score, P24_interface_df$upKla, method = "pearson")
interface_glyco_kla_cor_P <- cor.test(P24_interface_df$Glycolysis_score, P24_interface_df$upKla, method = "pearson")
interface_glyco_kla_cor_Pvalue <- interface_glyco_kla_cor_P$p.value

### Stroma
P24_stroma_df <- P24_df %>% 
  filter(Location == "Stroma")

stroma_glyco_kla_cor <- cor(P24_stroma_df$Glycolysis_score, P24_stroma_df$upKla, method = "pearson")
stroma_glyco_kla_cor_P <- cor.test(P24_stroma_df$Glycolysis_score, P24_stroma_df$upKla, method = "pearson")
stroma_glyco_kla_cor_Pvalue <- stroma_glyco_kla_cor_P$p.value

### 汇总
P24_corR <- data.frame(ID = "P24",
                       nMal = stroma_glyco_kla_cor,
                       Bdy = interface_glyco_kla_cor,
                       Mal = tumor_glyco_kla_cor,
                       All = glyco_kla_cor)

P24_corP <- data.frame(ID = "P24",
                       nMal = stroma_glyco_kla_cor_Pvalue,
                       Bdy = interface_glyco_kla_cor_Pvalue,
                       Mal = tumor_glyco_kla_cor_Pvalue,
                       All = glyco_kla_cor_Pvalue)

## 8.5 P25----
rm(list = setdiff(ls(), c("P10_corR", "P10_corP", "P15_corR", "P15_corP", "P16_corR", "P16_corP",
                          "P24_corR", "P24_corP")))

load("14.Location_heterogeneity/metabolism/P25_T1_glycolysis_ssGSEA_score.RData")
load("14.Location_heterogeneity/metabolism/P25_T1_upKla_ssGSEA_score.RData")
load("14.Location_heterogeneity/P25_T1_LocInfo.RData")

identical(rownames(loc_df), rownames(glycolysis_score))
identical(rownames(loc_df), rownames(upKla_score))

P25_df <- cbind(loc_df, glycolysis_score, upKla_score)

### All
glyco_kla_cor <- cor(P25_df$Glycolysis_score, P25_df$upKla, method = "pearson")
glyco_kla_cor_P <- cor.test(P25_df$Glycolysis_score, P25_df$upKla, method = "pearson")
glyco_kla_cor_Pvalue <- glyco_kla_cor_P$p.value

### Tumor
P25_tumor_df <- P25_df %>% 
  filter(Location == "Tumor")

tumor_glyco_kla_cor <- cor(P25_tumor_df$Glycolysis_score, P25_tumor_df$upKla, method = "pearson")
tumor_glyco_kla_cor_P <- cor.test(P25_tumor_df$Glycolysis_score, P25_tumor_df$upKla, method = "pearson")
tumor_glyco_kla_cor_Pvalue <- tumor_glyco_kla_cor_P$p.value

### Interface
P25_interface_df <- P25_df %>% 
  filter(Location == "Interface")

interface_glyco_kla_cor <- cor(P25_interface_df$Glycolysis_score, P25_interface_df$upKla, method = "pearson")
interface_glyco_kla_cor_P <- cor.test(P25_interface_df$Glycolysis_score, P25_interface_df$upKla, method = "pearson")
interface_glyco_kla_cor_Pvalue <- interface_glyco_kla_cor_P$p.value

### Stroma
P25_stroma_df <- P25_df %>% 
  filter(Location == "Stroma")

stroma_glyco_kla_cor <- cor(P25_stroma_df$Glycolysis_score, P25_stroma_df$upKla, method = "pearson")
stroma_glyco_kla_cor_P <- cor.test(P25_stroma_df$Glycolysis_score, P25_stroma_df$upKla, method = "pearson")
stroma_glyco_kla_cor_Pvalue <- stroma_glyco_kla_cor_P$p.value

### 汇总
P25_corR <- data.frame(ID = "P25",
                       nMal = stroma_glyco_kla_cor,
                       Bdy = interface_glyco_kla_cor,
                       Mal = tumor_glyco_kla_cor,
                       All = glyco_kla_cor)

P25_corP <- data.frame(ID = "P25",
                       nMal = stroma_glyco_kla_cor_Pvalue,
                       Bdy = interface_glyco_kla_cor_Pvalue,
                       Mal = tumor_glyco_kla_cor_Pvalue,
                       All = glyco_kla_cor_Pvalue)

## 8.6 汇总5例相关性分析结果----
merge_R <- rbind(P10_corR, P15_corR, P16_corR, P24_corR, P25_corR) %>% 
  column_to_rownames("ID")
save(merge_R, file = "14.Location_heterogeneity/metabolism/AllPatients_kla_glyco_corR.RData")

merge_P <- rbind(P10_corP, P15_corP, P16_corP, P24_corP, P25_corP) %>% 
  column_to_rownames("ID")
save(merge_P, file = "14.Location_heterogeneity/metabolism/AllPatients_kla_glyco_corP.RData")

### 相关性热图可视化----
library(pheatmap)
library(grid)

merge_P2 <- ifelse(merge_P < 0.0001, "****",
                   ifelse(merge_P < 0.001, "***",
                          ifelse(merge_P < 0.01, "**",
                                 ifelse(merge_P < 0.05, "*", ""))))
bk <- c(seq(-2,-0.1,by=0.01),seq(0,2,by=0.01))

cor_heatmap <- pheatmap(
  merge_R, # 相关性系数矩阵，转置后的结果
  scale = "row", # 不对数据进行缩放
  border_color = "white", # 设置单元格边框颜色为白色
  number_color = "black", # 显示数字的颜色为白色
  fontsize_number = 10, # 数字字体大小设置为14
  fontsize_row = 8, # 行标签字体大小设置为8
  fontsize_col = 9, # 列标签字体大小设置为9
  cellwidth = 18, # 单元格宽度设置为15
  cellheight = 18, # 单元格高度设置为15
  cluster_rows = F, # 对行进行聚类
  cluster_cols = F, # 对列进行聚类
  #color = mycolor, # 使用自定义的颜色调色板
  color = c(colorRampPalette(colors = c("#094687","#DDEAF3"))(length(bk)/2),colorRampPalette(colors = c("#DDEAF3","#79C9C7"))(length(bk)/2)),
  #legend_breaks=seq(-2,2,0.5),
  breaks = bk,
  gaps_col = 3,
  angle_col = "45", # 列标签旋转角度设置为45度
  display_numbers = merge_P2, # 在单元格中显示显著性标记
  show_rownames = TRUE # 显示行标签
)

ggsave(filename = "14.Location_heterogeneity/metabolism/AllPatients_kla_glyco_cor_heatmap.pdf",
       cor_heatmap, height = 4, width = 4, dpi = 300)
