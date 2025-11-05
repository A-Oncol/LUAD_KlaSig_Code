######目的：semla分析以Bdy为中心的距离
######作者：申奥
######日期：2024-12-02
######Server: yjx@xiyoucloud.top:34022
######Version: R 4.4.1


# 0.环境设置----
rm(list = ls())
options(stringAsFactors = F)

library(semla)
library(tibble)
library(ggplot2)
library(patchwork)
library(scico)
library(tidyr)
library(dplyr)


# 1.P10----
## 1.1 semla读取数据----
samples <- "1.data/SP/PT10_T1/filtered_feature_bc_matrix.h5"
imgs <- "1.data/SP/PT10_T1/spatial/tissue_hires_image.png"
spotfiles <- "1.data/SP/PT10_T1/spatial/tissue_positions_list.csv"
json <- "1.data/SP/PT10_T1/spatial/scalefactors_json.json"

infoTable <- tibble(samples, imgs, spotfiles, json)

p10_se <- ReadVisiumData(infoTable = infoTable)

## 1.2 添加SpaCET识别的区域信息到SPATA2对象的meta.data中----
# load("2.Cottrazm_outputs/PT10_T1/PT10_T1_TumorST_Boudary.RData") # 利用Cottrazm识别的边界信息
load("14.Location_heterogeneity/P10_T1_LocInfo.RData") # 36代码1.2

identical(rownames(loc_df), rownames(p10_se@meta.data))
p10_se$Location <- loc_df$Location

p10_se <- LoadImages(p10_se)
MapLabels(p10_se, column_name = "Location", 
          image_use = "raw", pt_alpha = 0.6, pt_size = 2) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right") &
  guides(fill = guide_legend(override.aes = list(size = 3), ncol = 2))

table(p10_se$Location)
p10_se$isTumor <- ifelse(p10_se$Location == "Tumor", "Tumor", "non-Tumor")
MapLabels(p10_se, column_name = "isTumor", override_plot_dims = TRUE, 
          image_use = "raw", drop_na = TRUE, pt_size = 2) +
  theme(legend.position = "right") &
  guides(fill = guide_legend(override.aes = list(size = 5), ncol = 2))

## 1.3 计算距离----
p10_se <- RadialDistance(p10_se, column_name = "isTumor", selected_groups = "Tumor")
MapFeatures(p10_se, features = "r_dist_Tumor", center_zero = TRUE, pt_size = 2, 
            colors = RColorBrewer::brewer.pal(n = 11, name = "RdBu") |> rev(),
            override_plot_dims = TRUE)

### 转换成um
#se$r_dist_Tumor_um <- (100/273)*se$r_dist_3
p10_se <- RadialDistance(p10_se, column_name = "isTumor", 
                         selected_groups = "Tumor", convert_to_microns = TRUE)
MapFeatures(p10_se, features = "r_dist_Tumor", center_zero = TRUE, pt_size = 1.5, 
            colors = RColorBrewer::brewer.pal(n = 11, name = "RdBu") |> rev(),
            override_plot_dims = TRUE) # 高5 宽5 

## 1.4 添加糖酵解评分和组蛋白乳酸化评分到meta.data中----
load("14.Location_heterogeneity/metabolism/P10_T1_seuObj_glyco_upKla.RData") # 42脚本2.3
identical(rownames(st_obj@meta.data), rownames(p10_se@meta.data))

p10_se$Glycolysis_score <- st_obj$glycolysis_score
p10_se$upKla_score <- st_obj$upKla_score

save(p10_se, file = "19.semla/P10_semlaObj.RData")

## 1.5 画图----
score_df <- p10_se@meta.data[,c("r_dist_Tumor", "Glycolysis_score", "upKla_score")]

ggplot(score_df, aes(r_dist_Tumor, Glycolysis_score)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs")) +
  geom_vline(aes(xintercept = 0), linetype = "dashed", color = "gray") +
  theme_pubr()

ggplot(score_df, aes(r_dist_Tumor, upKla_score)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs")) +
  geom_vline(aes(xintercept = 0), linetype = "dashed", color = "gray") +
  theme_pubr()


# 2.P15----
rm(list = ls())
## 2.1 semla读取数据----
samples <- "1.data/SP/PT15_T1/filtered_feature_bc_matrix.h5"
imgs <- "1.data/SP/PT15_T1/spatial/tissue_hires_image.png"
spotfiles <- "1.data/SP/PT15_T1/spatial/tissue_positions_list.csv"
json <- "1.data/SP/PT15_T1/spatial/scalefactors_json.json"

infoTable <- tibble(samples, imgs, spotfiles, json)

p15_se <- ReadVisiumData(infoTable = infoTable)

## 2.2 添加SpaCET识别的区域信息到SPATA2对象的meta.data中----
# load("2.Cottrazm_outputs/PT10_T1/PT10_T1_TumorST_Boudary.RData") # 利用Cottrazm识别的边界信息
load("14.Location_heterogeneity/P15_T1_LocInfo.RData") # 36代码2.2

identical(rownames(loc_df), rownames(p15_se@meta.data))
p15_se$Location <- loc_df$Location

p15_se <- LoadImages(p15_se)
MapLabels(p15_se, column_name = "Location", 
          image_use = "raw", pt_alpha = 0.6, pt_size = 2) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right") &
  guides(fill = guide_legend(override.aes = list(size = 3), ncol = 2))

table(p15_se$Location)
p15_se$isTumor <- ifelse(p15_se$Location == "Tumor", "Tumor", "non-Tumor")
MapLabels(p15_se, column_name = "isTumor", override_plot_dims = TRUE, 
          image_use = "raw", drop_na = TRUE, pt_size = 2) +
  theme(legend.position = "right") &
  guides(fill = guide_legend(override.aes = list(size = 5), ncol = 2))

## 2.3 计算距离----
p15_se <- RadialDistance(p15_se, column_name = "isTumor", selected_groups = "Tumor")
MapFeatures(p15_se, features = "r_dist_Tumor", center_zero = TRUE, pt_size = 2, 
            colors = RColorBrewer::brewer.pal(n = 11, name = "RdBu") |> rev(),
            override_plot_dims = TRUE)

### 转换成um
#se$r_dist_Tumor_um <- (100/273)*se$r_dist_3
p15_se <- RadialDistance(p15_se, column_name = "isTumor", 
                         selected_groups = "Tumor", convert_to_microns = TRUE)
MapFeatures(p15_se, features = "r_dist_Tumor", center_zero = TRUE, pt_size = 1.5, 
            colors = RColorBrewer::brewer.pal(n = 11, name = "RdBu") |> rev(),
            override_plot_dims = TRUE)

## 2.4 添加糖酵解评分和组蛋白乳酸化评分到meta.data中----
load("14.Location_heterogeneity/metabolism/P15_T1_seuObj_glyco_upKla.RData") # 42脚本3.3
identical(rownames(st_obj@meta.data), rownames(p15_se@meta.data))

p15_se$Glycolysis_score <- st_obj$glycolysis_score
p15_se$upKla_score <- st_obj$upKla_score

save(p15_se, file = "19.semla/P15_semlaObj.RData")

## 2.5 画图----
score_df <- p15_se@meta.data[,c("r_dist_Tumor", "Glycolysis_score", "upKla_score")]

ggplot(score_df, aes(r_dist_Tumor, Glycolysis_score)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs")) +
  geom_vline(aes(xintercept = 0), linetype = "dashed", color = "gray") +
  theme_pubr()

ggplot(score_df, aes(r_dist_Tumor, upKla_score)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs")) +
  geom_vline(aes(xintercept = 0), linetype = "dashed", color = "gray") +
  theme_pubr()


# 3.P16----
rm(list = ls())
## 3.1 semla读取数据----
samples <- "1.data/SP/PT16_T1/filtered_feature_bc_matrix.h5"
imgs <- "1.data/SP/PT16_T1/spatial/tissue_hires_image.png"
spotfiles <- "1.data/SP/PT16_T1/spatial/tissue_positions_list.csv"
json <- "1.data/SP/PT16_T1/spatial/scalefactors_json.json"

infoTable <- tibble(samples, imgs, spotfiles, json)

p16_se <- ReadVisiumData(infoTable = infoTable)

## 3.2 添加SpaCET识别的区域信息到SPATA2对象的meta.data中----
# load("2.Cottrazm_outputs/PT10_T1/PT10_T1_TumorST_Boudary.RData") # 利用Cottrazm识别的边界信息
load("14.Location_heterogeneity/P16_T1_LocInfo.RData") # 36代码3.2

identical(rownames(loc_df), rownames(p16_se@meta.data))
p16_se$Location <- loc_df$Location

p16_se <- LoadImages(p16_se)
MapLabels(p16_se, column_name = "Location", 
          image_use = "raw", pt_alpha = 0.6, pt_size = 2) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right") &
  guides(fill = guide_legend(override.aes = list(size = 3), ncol = 2))

table(p16_se$Location)
p16_se$isTumor <- ifelse(p16_se$Location == "Tumor", "Tumor", "non-Tumor")
MapLabels(p16_se, column_name = "isTumor", override_plot_dims = TRUE, 
          image_use = "raw", drop_na = TRUE, pt_size = 2) +
  theme(legend.position = "right") &
  guides(fill = guide_legend(override.aes = list(size = 5), ncol = 2))

## 3.3 计算距离----
p16_se <- RadialDistance(p16_se, column_name = "isTumor", selected_groups = "Tumor")
MapFeatures(p16_se, features = "r_dist_Tumor", center_zero = TRUE, pt_size = 2, 
            colors = RColorBrewer::brewer.pal(n = 11, name = "RdBu") |> rev(),
            override_plot_dims = TRUE)

### 转换成um
#se$r_dist_Tumor_um <- (100/273)*se$r_dist_3
p16_se <- RadialDistance(p16_se, column_name = "isTumor", 
                         selected_groups = "Tumor", convert_to_microns = TRUE)
MapFeatures(p16_se, features = "r_dist_Tumor", center_zero = TRUE, pt_size = 1.5, 
            colors = RColorBrewer::brewer.pal(n = 11, name = "RdBu") |> rev(),
            override_plot_dims = TRUE) # 高5 宽5 

## 3.4 添加糖酵解评分和组蛋白乳酸化评分到meta.data中----
load("14.Location_heterogeneity/metabolism/P16_T1_seuObj_glyco_upKla.RData") # 42脚本4.3
identical(rownames(st_obj@meta.data), rownames(p16_se@meta.data))

p16_se$Glycolysis_score <- st_obj$glycolysis_score
p16_se$upKla_score <- st_obj$upKla_score

save(p16_se, file = "19.semla/P16_semlaObj.RData")

## 3.5 画图----
score_df <- p16_se@meta.data[,c("r_dist_Tumor", "Glycolysis_score", "upKla_score")]

ggplot(score_df, aes(r_dist_Tumor, Glycolysis_score)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs")) +
  geom_vline(aes(xintercept = 0), linetype = "dashed", color = "gray") +
  theme_pubr()

ggplot(score_df, aes(r_dist_Tumor, upKla_score)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs")) +
  geom_vline(aes(xintercept = 0), linetype = "dashed", color = "gray") +
  theme_pubr()


# 4.P24----
rm(list = ls())
## 4.1 semla读取数据----
samples <- "1.data/SP/PT24_T1/filtered_feature_bc_matrix.h5"
imgs <- "1.data/SP/PT24_T1/spatial/tissue_hires_image.png"
spotfiles <- "1.data/SP/PT24_T1/spatial/tissue_positions_list.csv"
json <- "1.data/SP/PT24_T1/spatial/scalefactors_json.json"

infoTable <- tibble(samples, imgs, spotfiles, json)

p24_se <- ReadVisiumData(infoTable = infoTable)

## 4.2 添加SpaCET识别的区域信息到SPATA2对象的meta.data中----
# load("2.Cottrazm_outputs/PT10_T1/PT10_T1_TumorST_Boudary.RData") # 利用Cottrazm识别的边界信息
load("14.Location_heterogeneity/P24_T1_LocInfo.RData") # 36代码4.2

identical(rownames(loc_df), rownames(p24_se@meta.data))
p24_se$Location <- loc_df$Location

p24_se <- LoadImages(p24_se)
MapLabels(p24_se, column_name = "Location", 
          image_use = "raw", pt_alpha = 0.6, pt_size = 2) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right") &
  guides(fill = guide_legend(override.aes = list(size = 3), ncol = 2))

table(p24_se$Location)
p24_se$isTumor <- ifelse(p24_se$Location == "Tumor", "Tumor", "non-Tumor")
MapLabels(p24_se, column_name = "isTumor", override_plot_dims = TRUE, 
          image_use = "raw", drop_na = TRUE, pt_size = 2) +
  theme(legend.position = "right") &
  guides(fill = guide_legend(override.aes = list(size = 5), ncol = 2))

## 4.3 计算距离----
p24_se <- RadialDistance(p24_se, column_name = "isTumor", selected_groups = "Tumor")
MapFeatures(p24_se, features = "r_dist_Tumor", center_zero = TRUE, pt_size = 2, 
            colors = RColorBrewer::brewer.pal(n = 11, name = "RdBu") |> rev(),
            override_plot_dims = TRUE)

### 转换成um
#se$r_dist_Tumor_um <- (100/273)*se$r_dist_3
p24_se <- RadialDistance(p24_se, column_name = "isTumor", 
                         selected_groups = "Tumor", convert_to_microns = TRUE)
MapFeatures(p24_se, features = "r_dist_Tumor", center_zero = TRUE, pt_size = 1.5, 
            colors = RColorBrewer::brewer.pal(n = 11, name = "RdBu") |> rev(),
            override_plot_dims = TRUE) # 高5 宽5 

## 4.4 添加糖酵解评分和组蛋白乳酸化评分到meta.data中----
load("14.Location_heterogeneity/metabolism/P24_T1_seuObj_glyco_upKla.RData") # 42脚本5.3
identical(rownames(st_obj@meta.data), rownames(p24_se@meta.data))

p24_se$Glycolysis_score <- st_obj$glycolysis_score
p24_se$upKla_score <- st_obj$upKla_score

save(p24_se, file = "19.semla/P24_semlaObj.RData")

## 4.5 画图----
score_df <- p24_se@meta.data[,c("r_dist_Tumor", "Glycolysis_score", "upKla_score")]

ggplot(score_df, aes(r_dist_Tumor, Glycolysis_score)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs")) +
  geom_vline(aes(xintercept = 0), linetype = "dashed", color = "gray") +
  theme_pubr()

ggplot(score_df, aes(r_dist_Tumor, upKla_score)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs")) +
  geom_vline(aes(xintercept = 0), linetype = "dashed", color = "gray") +
  theme_pubr()


# 5.P25----
rm(list = ls())
## 5.1 semla读取数据----
samples <- "1.data/SP/PT25_T1/filtered_feature_bc_matrix.h5"
imgs <- "1.data/SP/PT25_T1/spatial/tissue_hires_image.png"
spotfiles <- "1.data/SP/PT25_T1/spatial/tissue_positions_list.csv"
json <- "1.data/SP/PT25_T1/spatial/scalefactors_json.json"

infoTable <- tibble(samples, imgs, spotfiles, json)

p25_se <- ReadVisiumData(infoTable = infoTable)

## 5.2 添加SpaCET识别的区域信息到SPATA2对象的meta.data中----
# load("2.Cottrazm_outputs/PT10_T1/PT10_T1_TumorST_Boudary.RData") # 利用Cottrazm识别的边界信息
load("14.Location_heterogeneity/P25_T1_LocInfo.RData") # 36代码5.2

identical(rownames(loc_df), rownames(p25_se@meta.data))
p25_se$Location <- loc_df$Location

p25_se <- LoadImages(p25_se)
MapLabels(p25_se, column_name = "Location", 
          image_use = "raw", pt_alpha = 0.6, pt_size = 2) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right") &
  guides(fill = guide_legend(override.aes = list(size = 3), ncol = 2))

table(p25_se$Location)
p25_se$isTumor <- ifelse(p25_se$Location == "Tumor", "Tumor", "non-Tumor")
MapLabels(p25_se, column_name = "isTumor", override_plot_dims = TRUE, 
          image_use = "raw", drop_na = TRUE, pt_size = 2) +
  theme(legend.position = "right") &
  guides(fill = guide_legend(override.aes = list(size = 5), ncol = 2))

## 5.3 计算距离----
p25_se <- RadialDistance(p25_se, column_name = "isTumor", selected_groups = "Tumor")
MapFeatures(p25_se, features = "r_dist_Tumor", center_zero = TRUE, pt_size = 2, 
            colors = RColorBrewer::brewer.pal(n = 11, name = "RdBu") |> rev(),
            override_plot_dims = TRUE)

### 转换成um
#se$r_dist_Tumor_um <- (100/273)*se$r_dist_3
p25_se <- RadialDistance(p25_se, column_name = "isTumor", 
                         selected_groups = "Tumor", convert_to_microns = TRUE)
MapFeatures(p25_se, features = "r_dist_Tumor", center_zero = TRUE, pt_size = 1.5, 
            colors = RColorBrewer::brewer.pal(n = 11, name = "RdBu") |> rev(),
            override_plot_dims = TRUE) # 高5 宽5 

## 5.4 添加糖酵解评分和组蛋白乳酸化评分到meta.data中----
load("14.Location_heterogeneity/metabolism/P25_T1_seuObj_glyco_upKla.RData") # 42脚本6.3
identical(rownames(st_obj@meta.data), rownames(p25_se@meta.data))

p25_se$Glycolysis_score <- st_obj$glycolysis_score
p25_se$upKla_score <- st_obj$upKla_score

save(p25_se, file = "19.semla/P25_semlaObj.RData")

## 5.5 画图----
score_df <- p25_se@meta.data[,c("r_dist_Tumor", "Glycolysis_score", "upKla_score")]

ggplot(score_df, aes(r_dist_Tumor, Glycolysis_score)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs")) +
  geom_vline(aes(xintercept = 0), linetype = "dashed", color = "gray") +
  theme_pubr()

ggplot(score_df, aes(r_dist_Tumor, upKla_score)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs")) +
  geom_vline(aes(xintercept = 0), linetype = "dashed", color = "gray") +
  theme_pubr()


# 6.综合绘制5例的糖酵解和upKla评分随距离的变化----
rm(list = ls())

load("19.semla/P10_semlaObj.RData")
load("19.semla/P15_semlaObj.RData")
load("19.semla/P16_semlaObj.RData")
load("19.semla/P24_semlaObj.RData")
load("19.semla/P25_semlaObj.RData")

score_df1 <- p10_se@meta.data[,c("r_dist_Tumor", "Glycolysis_score", "upKla_score")]
score_df1$Patient <- "P10"
score_df2 <- p15_se@meta.data[,c("r_dist_Tumor", "Glycolysis_score", "upKla_score")]
score_df2$Patient <- "P15"
score_df3 <- p16_se@meta.data[,c("r_dist_Tumor", "Glycolysis_score", "upKla_score")]
score_df3$Patient <- "P16"
score_df4 <- p24_se@meta.data[,c("r_dist_Tumor", "Glycolysis_score", "upKla_score")]
score_df4$Patient <- "P24"
score_df5 <- p25_se@meta.data[,c("r_dist_Tumor", "Glycolysis_score", "upKla_score")]
score_df5$Patient <- "P25"

merge_score <- rbind(score_df1, score_df2, score_df3, score_df4, score_df5)

df <- merge_score %>% 
  pivot_longer(all_of(sel_genes), names_to = "variable", values_to = "value")


p1 <- ggplot(merge_score, aes(r_dist_Tumor, Glycolysis_score, color = Patient)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE) +
  geom_vline(aes(xintercept = 0), linetype = "dashed", color = "gray") +
  scale_x_continuous(
    limits = c(-1000, 2000),   # 设置横轴范围
    breaks = seq(-1000, 2000, by = 1000)  # 设置横轴步长
  ) +
  theme_pubr() +
  xlab("Distance (μm)")
p1
ggsave(filename = "19.semla/all_glycolysis_distance.pdf", p1, height = 4, width = 5, dpi = 300,
       device = cairo_pdf)

p2 <- ggplot(merge_score, aes(r_dist_Tumor, upKla_score, color = Patient)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE) +
  geom_vline(aes(xintercept = 0), linetype = "dashed", color = "gray") +
  scale_x_continuous(
    limits = c(-1000, 2000),   # 设置横轴范围
    breaks = seq(-1000, 2000, by = 1000)  # 设置横轴步长
  ) +
  theme_pubr() +
  xlab("Distance (μm)")
p2
ggsave(filename = "19.semla/all_upKla_distance.pdf", p2, height = 4, width = 5, dpi = 300,
       device = cairo_pdf)
