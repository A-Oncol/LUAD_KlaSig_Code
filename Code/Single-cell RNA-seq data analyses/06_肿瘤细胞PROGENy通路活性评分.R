######目的：肿瘤细胞PROGENy通路活性评分
######作者：申奥
######日期：2024-12-27
######Server: yjx@xiyoucloud.top:34022
######Version: R 4.4.2


# 0.环境设置----
rm(list = ls())
options(stringAsFactors = F)

library(tidyverse)
library(Seurat)
packageVersion("Seurat") # ‘5.1.0’
library(progeny)


# 1.加载肿瘤细胞seurat对象----
sc.tumor <- readRDS("5.scRNA/GSE189357_tumorCell.rds") # 47脚本中对肿瘤细胞进行了评分分组
table(Idents(sc.tumor))
table(sc.tumor$TumorCell_group)

CellsClusters <- data.frame(Cell = rownames(sc.tumor@meta.data),
                            CellType = sc.tumor$TumorCell_group,
                            stringsAsFactors = F)


# 2.计算通路活性----
sc.tumor2 <- progeny(sc.tumor, scale = F, organism = "Human", top = 500, perm = 1, return_assay = T)
sc.tumor2 <- Seurat::ScaleData(sc.tumor2, assay = "progeny")
# scale_score <- progeny(sc.tumor, scale = T, organism = "Human", top = 500, perm = 1, return_assay = F) # 这行代码与上面两行效果一致，只不过这行代码只返回一个评分矩阵
qs::qsave(sc.tumor2, file = "27.progeny/scTumor_calProgenyScore.qs")


# 3.可视化----
library(pheatmap)

scale_score <- as.data.frame(t(sc.tumor2@assays[["progeny"]]@scale.data)) %>% 
  rownames_to_column("Cell") %>% 
  gather(Pathway, Activity, -Cell)

merge_score <- inner_join(scale_score, CellsClusters, by = "Cell")

summarized_progeny_scores <- merge_score %>% 
  group_by(Pathway, CellType) %>% 
  summarise(avg = mean(Activity), std = sd(Activity))

summarized_progeny_scores_df <- summarized_progeny_scores %>% 
  dplyr::select(-std) %>% 
  spread(Pathway, avg) %>% 
  data.frame(row.names = 1, check.names = F, stringsAsFactors = F)

paletteLength <- 100
myColors <- colorRampPalette(c("#75c7f2", "#ffffff", "#997ab9"))(paletteLength)

progenyBreaks = c(seq(min(summarized_progeny_scores_df), 0, 
                      length.out = ceiling(paletteLength/2) + 1),
                  seq(max(summarized_progeny_scores_df)/paletteLength, 
                      max(summarized_progeny_scores_df), 
                      length.out = floor(paletteLength/2)))

summarized_progeny_scores_df <- as.data.frame(t(summarized_progeny_scores_df))
qs::qsave(summarized_progeny_scores_df, file = "27.progeny/summarized_progeny_scores_df.qs")

p1 <- pheatmap(summarized_progeny_scores_df[,c(2,1)], 
               color= myColors, breaks = progenyBreaks, 
               angle_col = 45, scale = "none",
               cellheight = 12, cellwidth = 15,
               treeheight_col = 0,  border_color = NA)
p1
ggsave(filename = "27.progeny/progeny_pheatmap.pdf", p1, height = 4, width = 4, dpi = 300)
