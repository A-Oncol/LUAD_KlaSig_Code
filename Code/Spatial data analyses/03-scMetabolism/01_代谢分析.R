######目的：代谢分析
######作者：申奥
######日期：2024-09-28
######Server: yjx@xiyoucloud.top:34022
######Version: R 4.3.3


# 0.环境设置----
rm(list = ls())
options(stringAsFactors = F)

library(scMetabolism)
library(Seurat)
packageVersion('Seurat') # ‘5.1.0’
library(ggplot2)
library(rsvd)
library(ggpubr)
library(cols4all)


# 1.P10----
## 1.1 读取空转数据，并进行SCTransfrom和代谢评估----
st_obj = Load10X_Spatial(data.dir = "1.data/SP/PT10_T1/", slice = "P10_T1")
st_obj <- SCTransform(st_obj, assay = "Spatial")

countexp.Seurat <- st_obj
countexp.Seurat@assays$RNA <- countexp.Seurat@assays$SCT

countexp.Seurat <- sc.metabolism.Seurat(obj = countexp.Seurat, method = "ssGSEA", 
                                        imputation = F, ncores = 48, metabolism.type = "KEGG")

meta.score <- countexp.Seurat@assays[["METABOLISM"]][["score"]]
save(meta.score, file = "14.Location_heterogeneity/metabolism/P10_T1_KEGG_metabolism_score.RData")

## 1.2 加载SpaCET推断的空间区域信息----
load("13.SpaCET/P10_T1/P10_T1_SpaCET_obj.RData")
locations <- t(SpaCET_obj@results[["CCI"]][["interface"]])
barcodes <- SpaCET_obj@input[["spotCoordinates"]][["barcode"]]

loc_df <- data.frame(Barcodes = barcodes,
                     Location = locations)
rownames(loc_df) <- loc_df$Barcodes
colnames(loc_df)[2] <- "Location"
table(loc_df$Location)
save(loc_df, file = "14.Location_heterogeneity/P10_T1_LocInfo.RData")

## 1.3 合并空间区域信息到代谢评分的meta.data中----
identical(rownames(countexp.Seurat@meta.data), loc_df$Barcodes)

countexp.Seurat$Location <- loc_df$Location
class(countexp.Seurat$Location)
table(countexp.Seurat$Location)
countexp.Seurat$Location <- ifelse(countexp.Seurat$Location == "Stroma", "nMal",
                                   ifelse(countexp.Seurat$Location == "Tumor", "Mal", "Bdy"))
table(countexp.Seurat$Location)
# Bdy  Mal nMal 
# 464  927  809
countexp.Seurat$Location <- factor(countexp.Seurat$Location, levels = c("nMal", "Bdy", "Mal"))
save(countexp.Seurat, file = "14.Location_heterogeneity/metabolism/P10_T1_seuObj_runscMeta.RData")

## 1.4 画图----
input.pathway <- rownames(meta.score)
p1 <- DotPlot.metabolism(obj = countexp.Seurat,
                         pathway = input.pathway, phenotype = "Location", norm = "y") +
  scale_x_discrete(limits = c("nMal", "Bdy", "Mal")) +
  ylab(NULL) + xlab("P10") +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 12, family = "Arial", color = "black"),
        axis.text.x = element_text(size = 11, family = "Arial", color = "black"),
        axis.title.x = element_text(size = 14))
p1
save(p1, file = "14.Location_heterogeneity/metabolism/P10_T1_dotplot.RData")


# 2.P15----
rm(list = ls())
## 2.1 读取空转数据，并进行SCTransfrom和代谢评估----
st_obj = Load10X_Spatial(data.dir = "1.data/SP/PT15_T1/", slice = "P15_T1")
st_obj <- SCTransform(st_obj, assay = "Spatial")

countexp.Seurat <- st_obj
countexp.Seurat@assays$RNA <- countexp.Seurat@assays$SCT

countexp.Seurat <- sc.metabolism.Seurat(obj = countexp.Seurat, method = "ssGSEA", 
                                        imputation = F, ncores = 48, metabolism.type = "KEGG")

meta.score <- countexp.Seurat@assays[["METABOLISM"]][["score"]]
save(meta.score, file = "14.Location_heterogeneity/metabolism/P15_T1_KEGG_metabolism_score.RData")

## 2.2 加载SpaCET推断的空间区域信息----
load("13.SpaCET/P15_T1/P15_T1_SpaCET_obj.RData")
locations <- t(SpaCET_obj@results[["CCI"]][["interface"]])
barcodes <- SpaCET_obj@input[["spotCoordinates"]][["barcode"]]

loc_df <- data.frame(Barcodes = barcodes,
                     Location = locations)
rownames(loc_df) <- loc_df$Barcodes
colnames(loc_df)[2] <- "Location"
table(loc_df$Location)
save(loc_df, file = "14.Location_heterogeneity/P15_T1_LocInfo.RData")

## 2.3 合并空间区域信息到代谢评分的meta.data中----
identical(rownames(countexp.Seurat@meta.data), loc_df$Barcodes)

countexp.Seurat$Location <- loc_df$Location
class(countexp.Seurat$Location)
table(countexp.Seurat$Location)
countexp.Seurat$Location <- ifelse(countexp.Seurat$Location == "Stroma", "nMal",
                                   ifelse(countexp.Seurat$Location == "Tumor", "Mal", "Bdy"))
table(countexp.Seurat$Location)
# Bdy  Mal nMal 
# 897 1440  287
countexp.Seurat$Location <- factor(countexp.Seurat$Location, levels = c("nMal", "Bdy", "Mal"))
save(countexp.Seurat, file = "14.Location_heterogeneity/metabolism/P15_T1_seuObj_runscMeta.RData")

## 2.4 画图----
input.pathway <- rownames(meta.score)
p2 <- DotPlot.metabolism(obj = countexp.Seurat,
                         pathway = input.pathway, phenotype = "Location", norm = "y") +
  scale_x_discrete(limits = c("nMal", "Bdy", "Mal")) +
  ylab(NULL) + xlab("P15") +
  theme(legend.position = "none",
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 11, family = "Arial", color = "black"),
        axis.title.x = element_text(size = 14))
p2
save(p2, file = "14.Location_heterogeneity/metabolism/P15_T1_dotplot.RData")


# 3.P16----
rm(list = ls())
## 3.1 读取空转数据，并进行SCTransfrom和代谢评估----
st_obj = Load10X_Spatial(data.dir = "1.data/SP/PT16_T1/", slice = "P16_T1")
st_obj <- SCTransform(st_obj, assay = "Spatial")

countexp.Seurat <- st_obj
countexp.Seurat@assays$RNA <- countexp.Seurat@assays$SCT

countexp.Seurat <- sc.metabolism.Seurat(obj = countexp.Seurat, method = "ssGSEA", 
                                        imputation = F, ncores = 48, metabolism.type = "KEGG")

meta.score <- countexp.Seurat@assays[["METABOLISM"]][["score"]]
save(meta.score, file = "14.Location_heterogeneity/metabolism/P16_T1_KEGG_metabolism_score.RData")

## 3.2 加载SpaCET推断的空间区域信息----
load("13.SpaCET/P16_T1/P16_T1_SpaCET_obj.RData")
locations <- t(SpaCET_obj@results[["CCI"]][["interface"]])
barcodes <- SpaCET_obj@input[["spotCoordinates"]][["barcode"]]

loc_df <- data.frame(Barcodes = barcodes,
                     Location = locations)
rownames(loc_df) <- loc_df$Barcodes
colnames(loc_df)[2] <- "Location"
table(loc_df$Location)
save(loc_df, file = "14.Location_heterogeneity/P16_T1_LocInfo.RData")

## 3.3 合并空间区域信息到代谢评分的meta.data中----
identical(rownames(countexp.Seurat@meta.data), loc_df$Barcodes)

countexp.Seurat$Location <- loc_df$Location
class(countexp.Seurat$Location)
table(countexp.Seurat$Location)
countexp.Seurat$Location <- ifelse(countexp.Seurat$Location == "Stroma", "nMal",
                                   ifelse(countexp.Seurat$Location == "Tumor", "Mal", "Bdy"))
table(countexp.Seurat$Location)
# Bdy  Mal nMal 
# 424  703 1135
countexp.Seurat$Location <- factor(countexp.Seurat$Location, levels = c("nMal", "Bdy", "Mal"))
save(countexp.Seurat, file = "14.Location_heterogeneity/metabolism/P16_T1_seuObj_runscMeta.RData")

## 3.4 画图----
input.pathway <- rownames(meta.score)
p3 <- DotPlot.metabolism(obj = countexp.Seurat,
                         pathway = input.pathway, phenotype = "Location", norm = "y") +
  scale_x_discrete(limits = c("nMal", "Bdy", "Mal")) +
  ylab(NULL) + xlab("P16") +
  theme(legend.position = "none",
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 11, family = "Arial", color = "black"),
        axis.title.x = element_text(size = 14))
p3
save(p3, file = "14.Location_heterogeneity/metabolism/P16_T1_dotplot.RData")


# 4.P24----
rm(list = ls())
## 4.1 读取空转数据，并进行SCTransfrom和代谢评估----
st_obj = Load10X_Spatial(data.dir = "1.data/SP/PT24_T1/", slice = "P24_T1")
st_obj <- SCTransform(st_obj, assay = "Spatial")

countexp.Seurat <- st_obj
countexp.Seurat@assays$RNA <- countexp.Seurat@assays$SCT

countexp.Seurat <- sc.metabolism.Seurat(obj = countexp.Seurat, method = "ssGSEA", 
                                        imputation = F, ncores = 48, metabolism.type = "KEGG")

meta.score <- countexp.Seurat@assays[["METABOLISM"]][["score"]]
save(meta.score, file = "14.Location_heterogeneity/metabolism/P24_T1_KEGG_metabolism_score.RData")

## 4.2 加载SpaCET推断的空间区域信息----
load("13.SpaCET/P24_T1/P24_T1_SpaCET_obj.RData")
locations <- t(SpaCET_obj@results[["CCI"]][["interface"]])
barcodes <- SpaCET_obj@input[["spotCoordinates"]][["barcode"]]

loc_df <- data.frame(Barcodes = barcodes,
                     Location = locations)
rownames(loc_df) <- loc_df$Barcodes
colnames(loc_df)[2] <- "Location"
table(loc_df$Location)
save(loc_df, file = "14.Location_heterogeneity/P24_T1_LocInfo.RData")

## 4.3 合并空间区域信息到代谢评分的meta.data中----
identical(rownames(countexp.Seurat@meta.data), loc_df$Barcodes)

countexp.Seurat$Location <- loc_df$Location
class(countexp.Seurat$Location)
table(countexp.Seurat$Location)
countexp.Seurat$Location <- ifelse(countexp.Seurat$Location == "Stroma", "nMal",
                                   ifelse(countexp.Seurat$Location == "Tumor", "Mal", "Bdy"))
table(countexp.Seurat$Location)
# Bdy  Mal nMal 
# 1279 1664  483
countexp.Seurat$Location <- factor(countexp.Seurat$Location, levels = c("nMal", "Bdy", "Mal"))
save(countexp.Seurat, file = "14.Location_heterogeneity/metabolism/P24_T1_seuObj_runscMeta.RData")

## 4.4 画图----
input.pathway <- rownames(meta.score)
p4 <- DotPlot.metabolism(obj = countexp.Seurat,
                         pathway = input.pathway, phenotype = "Location", norm = "y") +
  scale_x_discrete(limits = c("nMal", "Bdy", "Mal")) +
  ylab(NULL) + xlab("P24") +
  theme(legend.position = "none",
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 11, family = "Arial", color = "black"),
        axis.title.x = element_text(size = 14))
p4
save(p4, file = "14.Location_heterogeneity/metabolism/P24_T1_dotplot.RData")


# 5.P25----
rm(list = ls())
## 5.1 读取空转数据，并进行SCTransfrom和代谢评估----
st_obj = Load10X_Spatial(data.dir = "1.data/SP/PT25_T1/", slice = "P25_T1")
st_obj <- SCTransform(st_obj, assay = "Spatial")

countexp.Seurat <- st_obj
countexp.Seurat@assays$RNA <- countexp.Seurat@assays$SCT

countexp.Seurat <- sc.metabolism.Seurat(obj = countexp.Seurat, method = "ssGSEA", 
                                        imputation = F, ncores = 48, metabolism.type = "KEGG")

meta.score <- countexp.Seurat@assays[["METABOLISM"]][["score"]]
save(meta.score, file = "14.Location_heterogeneity/metabolism/P25_T1_KEGG_metabolism_score.RData")

## 5.2 加载SpaCET推断的空间区域信息----
load("13.SpaCET/P25_T1/P25_T1_SpaCET_obj.RData")
locations <- t(SpaCET_obj@results[["CCI"]][["interface"]])
barcodes <- SpaCET_obj@input[["spotCoordinates"]][["barcode"]]

loc_df <- data.frame(Barcodes = barcodes,
                     Location = locations)
rownames(loc_df) <- loc_df$Barcodes
colnames(loc_df)[2] <- "Location"
table(loc_df$Location)
save(loc_df, file = "14.Location_heterogeneity/P25_T1_LocInfo.RData")

## 5.3 合并空间区域信息到代谢评分的meta.data中----
identical(rownames(countexp.Seurat@meta.data), loc_df$Barcodes)

countexp.Seurat$Location <- loc_df$Location
class(countexp.Seurat$Location)
table(countexp.Seurat$Location)
countexp.Seurat$Location <- ifelse(countexp.Seurat$Location == "Stroma", "nMal",
                                   ifelse(countexp.Seurat$Location == "Tumor", "Mal", "Bdy"))
table(countexp.Seurat$Location)
# Bdy  Mal nMal 
# 290  707  31
countexp.Seurat$Location <- factor(countexp.Seurat$Location, levels = c("nMal", "Bdy", "Mal"))
save(countexp.Seurat, file = "14.Location_heterogeneity/metabolism/P25_T1_seuObj_runscMeta.RData")

## 5.4 画图----
input.pathway <- rownames(meta.score)
p5 <- DotPlot.metabolism(obj = countexp.Seurat,
                         pathway = input.pathway, phenotype = "Location", norm = "y") +
  scale_x_discrete(limits = c("nMal", "Bdy", "Mal")) +
  ylab(NULL) + xlab("P25") +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(size = 11, family = "Arial", color = "black"),
        axis.title.x = element_text(size = 14))
p5
save(p5, file = "14.Location_heterogeneity/metabolism/P25_T1_dotplot.RData")


# 6.dotplot拼图----
rm(list = ls())
library(patchwork)
load("14.Location_heterogeneity/metabolism/P10_T1_dotplot.RData")
load("14.Location_heterogeneity/metabolism/P15_T1_dotplot.RData")
load("14.Location_heterogeneity/metabolism/P16_T1_dotplot.RData")
load("14.Location_heterogeneity/metabolism/P24_T1_dotplot.RData")
load("14.Location_heterogeneity/metabolism/P25_T1_dotplot.RData")

p <- p1 + p2 + p3 + p4 + p5 + plot_layout(ncol = 5)
p # width: 1200, height: 1400