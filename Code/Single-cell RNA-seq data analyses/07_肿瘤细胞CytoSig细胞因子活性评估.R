######目的：肿瘤细胞CytoSig细胞因子活性评估
######作者：申奥
######日期：2024-12-30
######Server: yjx@xiyoucloud.top:34022
######Version: R 4.4.2


# 0.环境设置----
rm(list = ls())
options(stringAsFactors = F)

library(tidyverse)
library(Seurat, lib.loc = "~/bioSoft/seurat_v4/")
packageVersion("Seurat") # ‘4.4.0’
library(DropletUtils)
library(pheatmap)


# 1.CytoSig输入文件准备----
# sc.tumor <- readRDS("5.scRNA/GSE189357_tumorCell.rds") # 47脚本中对肿瘤细胞进行了评分分组
# table(Idents(sc.tumor))
# table(sc.tumor$TumorCell_group)

# DropletUtils::write10xCounts(x = sc.tumor@assays$RNA@counts, path = "28.cytosig/scTumor_for_cytosig", version = "3") # 不用这种

# sce_average <- AverageExpression(sc.tumor, group.by = "TumorCell_group")
# sce_average_exp <- sce_average[['RNA']]
# write.csv(sce_average_exp, file = "28.cytosig/sce_average_exp.csv")

sc_obj <- readRDS("5.scRNA/GSE189357_final_ref.rds") # 44脚本中的所有细胞
table(sc_obj$Cell_types)
# Normal epithelials        Tumor cells                  T                 NK                  B                 DC 
# 4601               5717              26427              17264              10067               7748 
# Mast        Macrophages        Endothelial        Fibroblasts          Monocytes             Plasma 
# 7778               7486               4141               3067               2626               2548 
# Proliferative 
# 842
sc.tumor <- readRDS("5.scRNA/GSE189357_tumorCell.rds") # 47脚本中对肿瘤细胞进行了评分分组
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
sce_average <- AverageExpression(sc_obj, group.by = "Cell_types")
sce_average_exp <- sce_average[['RNA']]
write.csv(sce_average_exp, file = "28.cytosig/sce_average_exp_allCells.csv")


# 2.跑run_cytosig.sh----
# nohup bash run_cytosig.sh &


# 3.CytoSig输出分析----
cytosig_zscore <- read.table("28.cytosig/sce_average_exp_allCells_output.Zscore", header = T, sep = "\t")
colnames(cytosig_zscore) <- gsub("\\.", "-", colnames(cytosig_zscore))
cytosig_pvalue <- read.table("28.cytosig/sce_average_exp_allCells_output.Pvalue", header = T, sep = "\t")
sig_cytosig_pvalue <- cytosig_pvalue %>% 
  mutate(Cytokines = rownames(.)) %>% 
  filter(High.upKla_TumorCells < 0.05 & Low.upKla_TumorCells < 0.05 & Normal.epithelials < 0.05) %>% 
  pull(Cytokines)

bk <- c(seq(-1,-0.1,by=0.01),seq(0,1,by=0.01))

p1 <- pheatmap(cytosig_zscore[sig_cytosig_pvalue,c(3,1,2)], cellheight = 10, cellwidth = 15, 
               cluster_cols = F, cluster_rows = T, angle_col = 45,
               scale = "row",
               color = c(colorRampPalette(colors = c("#75c7f2","#ffffff"))(length(bk)/2),
                         colorRampPalette(colors = c("#ffffff","#997ab9"))(length(bk)/2)))
p1
ggsave(filename = "28.cytosig/cytosig_pheatmap.pdf", p1, height = 4, width = 4, dpi = 300)
