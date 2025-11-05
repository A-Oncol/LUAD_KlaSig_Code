######目的：单细胞虚拟敲减
######作者：申奥
######日期：2024-12-31
######Server: yjx@xiyoucloud.top:34022
######Version: R 4.4.2


# 0.环境设置----
rm(list = ls())
options(stringAsFactors = F)

library(Seurat, lib.loc = "~/bioSoft/seurat_v4/") # ‘4.4.0’
library(scTenifoldKnk)
library(qs)
library(tidyverse)
library(clusterProfiler)
library(msigdbr)
library(enrichplot)
library(org.Hs.eg.db)
library(openxlsx)


# 1.读取Seurat对象----
sc_obj <- readRDS("~/projects/Rprojects/Bdy_LUAD/5.scRNA/GSE189357_final_ref.rds") # 44脚本中的所有细胞
table(sc_obj$Cell_types)
# Normal epithelials        Tumor cells                  T                 NK                  B                 DC 
# 4601               5717              26427              17264              10067               7748 
# Mast        Macrophages        Endothelial        Fibroblasts          Monocytes             Plasma 
# 7778               7486               4141               3067               2626               2548 
# Proliferative 
# 842

sc.tumor <- readRDS("~/projects/Rprojects/Bdy_LUAD/5.scRNA/GSE189357_tumorCell.rds") # 47脚本中对肿瘤细胞进行了评分分组
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

ExportGeneExpr <- sc.tumor@assays$RNA@counts
dim(ExportGeneExpr)
# 23451  5717
# ExportGeneExpr[1:5,1:5]


# 2.scTenifoldKnk----
# KO_PLK1 <- scTenifoldKnk(countMatrix = ExportGeneExpr, qc = F, gKO = "PLK1", nCores = 8)
# qs::qsave(KO_PLK1, file = "~/projects/Rprojects/Bdy_LUAD/30.scTenifoldKnk/KO_PLK1_inTumorCells.qs")

KO_MIF <- scTenifoldKnk(countMatrix = ExportGeneExpr, qc = F, gKO = "MIF", nCores = 8)
qs::qsave(KO_MIF, file = "~/projects/Rprojects/Bdy_LUAD/30.scTenifoldKnk/KO_MIF_inTumorCells.qs")


# 3.分析敲减后的数据----
## 3.1 读取并筛选有意义的基因----
KO_MIF <- qs::qread("30.scTenifoldKnk/KO_MIF_inTumorCells.qs")
KO_MIF_regulation <- write.table(KO_MIF[["diffRegulation"]], file = "30.scTenifoldKnk/KO_MIF_inTumorCells.txt", quote = F, sep = "\t")

diff_res <- KO_MIF[["diffRegulation"]]
sig_diff <- diff_res %>% 
  dplyr::filter(p.adj < 0.05) %>% 
  dplyr::pull(gene)

## 3.2 富集分析----
### 3.2.1 KEGG----
com_gene <- sig_diff

gene_KEGG <- bitr(geneID = com_gene,
                  fromType = "SYMBOL",
                  toType = "ENTREZID",
                  OrgDb = "org.Hs.eg.db",
                  drop = T)
KEGG <- enrichKEGG(gene = gene_KEGG$ENTREZID,
                   keyType = "kegg",
                   organism = "hsa", #小鼠：mmu
                   pvalueCutoff  = 0.05,
                   pAdjustMethod  = "BH",
                   qvalueCutoff  = 0.2)
write.table(KEGG@result, file = "30.scTenifoldKnk/KO_MIF_KEGG.txt", sep = "\t", row.names = F, quote = F)
write.xlsx(KEGG@result, file = "30.scTenifoldKnk/KO_MIF_KEGG.xlsx")
qs::qsave(KEGG, file = "30.scTenifoldKnk/KO_MIF_KEGG.qs")

### 3.2.2 GO-MF富集----
GO_MF <- enrichGO(gene = com_gene,
                  OrgDb = "org.Hs.eg.db", #小鼠："org.Mm.eg.db"
                  keyType = "SYMBOL",
                  ont = "MF",
                  pvalueCutoff = 0.05,
                  pAdjustMethod = "BH",
                  qvalueCutoff = 0.2,
                  minGSSize = 10,
                  maxGSSize = 500,
                  readable = FALSE,
                  pool = FALSE)
write.table(GO_MF@result, file = "30.scTenifoldKnk/KO_MIF_GO_MF.txt", sep = "\t", row.names = F, col.names = T, quote = F)
write.xlsx(GO_MF@result, file = "30.scTenifoldKnk/KO_MIF_GO_MF.xlsx")
qs::qsave(GO_MF, file = "30.scTenifoldKnk/KO_MIF_GO_MF.qs")

### 3.2.3 GO-BP富集----
GO_BP <- enrichGO(gene = com_gene,
                  OrgDb = "org.Hs.eg.db",
                  keyType = "SYMBOL",
                  ont = "BP",
                  pvalueCutoff = 0.05,
                  pAdjustMethod = "BH",
                  qvalueCutoff = 0.2,
                  minGSSize = 10,
                  maxGSSize = 500,
                  readable = FALSE,
                  pool = FALSE)
write.table(GO_BP@result, file = "30.scTenifoldKnk/KO_MIF_GO_BP.txt", sep = "\t", row.names = F, col.names = T, quote = F)
write.xlsx(GO_BP@result, file = "30.scTenifoldKnk/KO_MIF_GO_BP.xlsx")
qs::qsave(GO_BP, file = "30.scTenifoldKnk/KO_MIF_GO_BP.qs")

### 3.2.4 GO-CC富集----
GO_CC <- enrichGO(gene = com_gene,
                  OrgDb = "org.Hs.eg.db",
                  keyType = "SYMBOL",
                  ont = "CC",
                  pvalueCutoff = 0.05,
                  pAdjustMethod = "BH",
                  qvalueCutoff = 0.2,
                  minGSSize = 10,
                  maxGSSize = 500,
                  readable = FALSE,
                  pool = FALSE)
write.table(GO_CC@result, file = "30.scTenifoldKnk/KO_MIF_GO_CC.txt", sep = "\t", row.names = F, col.names = T, quote = F)
write.xlsx(GO_CC@result, file = "30.scTenifoldKnk/KO_MIF_GO_CC.xlsx")
qs::qsave(GO_CC, file = "30.scTenifoldKnk/KO_MIF_GO_CC.qs")

### 3.2.5 GSEA富集----
gsea_df <- diff_res[,c("gene", "FC")]
gsea_df$log2FC <- log2(gsea_df$FC)
colnames(gsea_df)[1] <- "SYMBOL"

gene2entrezID <- bitr(gsea_df$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

merged_df <- inner_join(gsea_df, gene2entrezID, by = "SYMBOL")
merged_df2 <- merged_df[!is.infinite(merged_df$log2FC),]

geneList <- merged_df2$log2FC
names(geneList) <- merged_df2$ENTREZID
head(geneList)
geneList <- na.omit(geneList)
geneList <- sort(geneList, decreasing = T)
head(geneList)

gsea_kegg <- gseKEGG(geneList,
                     organism = "hsa",
                     pvalueCutoff = 0.2,
                     pAdjustMethod = "BH",
                     minGSSize = 10,
                     maxGSSize = 500,
                     keyType = "kegg",
                     seed = 123)

gsea_go <- gseGO(geneList,
                 OrgDb = "org.Hs.eg.db",
                 ont = "ALL",
                 pvalueCutoff = 0.05,
                 keyType = "ENTREZID")

msigdbr_species() # 人：10 小鼠：13
a <- msigdbr_collections()

### 以人+Hallmark为例
Hs_hallmark_msigdbr <- msigdbr(species = "Homo sapiens", category = "H")
#Hs_reactome_msigdbr <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME")
colnames(Hs_hallmark_msigdbr)
Hs_hallmark_df <- as.data.frame(Hs_hallmark_msigdbr[,c('gs_name','entrez_gene','gene_symbol')])

### GSEA分析
hallmark_msig <- GSEA(geneList, TERM2GENE = Hs_hallmark_df[,c(1,2)], eps = 0)
head(hallmark_msig, 20)
hallmark_results <- hallmark_msig@result

# save(hallmark_msig, hallmark_results, file = "enrichment/GSEA_Hallmark_results.RData")
# write.csv(hallmark_results, "enrichment/GSEA_Hallmark_results.csv", row.names = F, quote = F)

### 3.2.6 可视化KEGG----
sig_KEGG <- KEGG@result %>% 
  filter(pvalue < 0.05)
colnames(sig_KEGG)[1] <- "Category"
colnames(sig_KEGG)[7] <- "P-value"
sig_KEGG$Description <- factor(sig_KEGG$Description, levels = sig_KEGG$Description)

p <- ggplot(sig_KEGG, aes(x = Description, y = Count, fill = `P-value`)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Count, y = Count + 1), size = 4) +
  coord_flip() +
  scale_fill_gradient(low = "#08306b", high = "#d4e4f4") +
  labs(x = '', y = "Gene Count") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "right",
        axis.ticks.y = element_blank(),
        axis.text = element_text(color = "black", size = 11))
p
ggsave(filename = "30.scTenifoldKnk/KO_MIF_KEGG_plot.pdf", p, height = 6, width = 6, dpi = 300)
