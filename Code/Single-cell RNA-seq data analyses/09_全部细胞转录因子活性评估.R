######目的：全部细胞转录因子活性评估
######作者：申奥
######日期：2025-01-08
######Server: yjx@xiyoucloud.top:34022
######Version: R 4.4.2


# 0.环境设置----
rm(list = ls())
options(stringAsFactors = F)

library(tidyverse)
library(Seurat, lib.loc = "~/bioSoft/seurat_v4/")
packageVersion("Seurat") # ‘4.4.0’
library(decoupleR)


# 1.加载全部细胞的分组----
sc_obj <- readRDS("/home/yjx/projects/Rprojects/Bdy_LUAD/5.scRNA/GSE189357_final_ref.rds") # 44脚本中的所有细胞
table(sc_obj$Cell_types)
# Normal epithelials        Tumor cells                  T                 NK                  B                 DC 
# 4601               5717              26427              17264              10067               7748 
# Mast        Macrophages        Endothelial        Fibroblasts          Monocytes             Plasma 
# 7778               7486               4141               3067               2626               2548 
# Proliferative 
# 842
sc.tumor <- readRDS("/home/yjx/projects/Rprojects/Bdy_LUAD/5.scRNA/GSE189357_tumorCell.rds") # 47脚本中对肿瘤细胞进行了评分分组
table(sc.tumor$TumorCell_group)
# High-upKla_TumorCells  Low-upKla_TumorCells 
# 2859                  2858 
Idents(sc.tumor) <- sc.tumor$TumorCell_group
Idents(sc_obj, cells = colnames(sc.tumor)) <- Idents(sc.tumor)
sc_obj$Cell_types <- Idents(sc_obj)
table(sc_obj$Cell_types)
# Low-upKla_TumorCells High-upKla_TumorCells    Normal epithelials                     T                    NK 
# 2858                  2859                  4601                 26427                 17264 
# B                    DC                  Mast           Macrophages           Endothelial 
# 10067                  7748                  7778                  7486                  4141 
# Fibroblasts             Monocytes                Plasma         Proliferative 
# 3067                  2626                  2548                   842 

table(Idents(sc_obj))
DimPlot(sc_obj)


# 2.TF分析----
## 2.1 CollecTRI network----
plan("multisession", workers = 4)
net <- get_collectri(organism = 'human', split_complexes = FALSE)
net

## 2.2 Activity inference with Univariate Linear Model (ULM)----
# Extract the normalized log-transformed counts
mat <- as.matrix(sc_obj@assays$RNA@data)
# Run ulm
acts <- decoupleR::run_ulm(mat = mat, 
                           net = net, 
                           .source = 'source', 
                           .target = 'target',
                           .mor='mor', 
                           minsize = 5)

acts
qs::qsave(acts, file = "/home/yjx/projects/Rprojects/Bdy_LUAD/31.decoupleR/acts.qs") ## 后台跑：nohup Rscript 66_全部细胞转录因子活性评估.R &

## 2.3 Visualization----
acts <- qs::qread("/home/yjx/projects/Rprojects/Bdy_LUAD/31.decoupleR/acts.qs")

# Extract ulm and store it in tfsulm in pbmc
sc_obj[['tfsulm']] <- acts %>%
  pivot_wider(id_cols = 'source', names_from = 'condition',
              values_from = 'score') %>%
  column_to_rownames('source') %>%
  Seurat::CreateAssayObject(.)

# Change assay
DefaultAssay(object = sc_obj) <- "tfsulm"

# Scale the data
sc_obj <- ScaleData(sc_obj)
sc_obj@assays$tfsulm@data <- sc_obj@assays$tfsulm@scale.data

DimPlot(sc_obj, reduction = "umap", label = TRUE, pt.size = 0.5) + 
  NoLegend() + ggtitle('Cell types')

## 2.4 Exploration----
n_tfs <- 25

# Extract activities from object as a long dataframe
df <- t(as.matrix(sc_obj@assays$tfsulm@data)) %>%
  as.data.frame() %>%
  dplyr::mutate(cluster = Seurat::Idents(sc_obj)) %>%
  tidyr::pivot_longer(cols = -cluster, 
                      names_to = "source", 
                      values_to = "score") %>%
  dplyr::group_by(cluster, source) %>%
  dplyr::summarise(mean = mean(score))

# Get top tfs with more variable means across clusters
tfs <- df %>%
  dplyr::group_by(source) %>%
  dplyr::summarise(std = stats::sd(mean)) %>%
  dplyr::arrange(-abs(std)) %>%
  head(n_tfs) %>%
  dplyr::pull(source)

# Subset long data frame to top tfs and transform to wide matrix
top_acts_mat <- df %>%
  dplyr::filter(source %in% tfs) %>%
  tidyr::pivot_wider(id_cols = 'cluster', 
                     names_from = 'source',
                     values_from = 'mean') %>%
  tibble::column_to_rownames('cluster') %>%
  as.matrix()

# Choose color palette
palette_length = 100
my_color = colorRampPalette(c("#75c7f2", "white","#997ab9"))(palette_length)

range(top_acts_mat)
# -2.540718  3.400334
my_breaks <- c(seq(-4, 0, length.out = ceiling(palette_length/2) + 1),
               seq(0.05, 4, length.out = floor(palette_length/2)))

# Plot
p1 <- pheatmap::pheatmap(top_acts_mat, 
                         #border_color = "gray90", 
                         color = my_color, 
                         breaks = my_breaks,
                         cellwidth = 15,
                         cellheight = 15,
                         treeheight_row = 20,
                         treeheight_col = 20, angle_col = 45) 
p1
ggplot2::ggsave(filename = "31.decoupleR/top25_TF_heatmap.pdf", p1, height = 4, width = 8, dpi = 300)
