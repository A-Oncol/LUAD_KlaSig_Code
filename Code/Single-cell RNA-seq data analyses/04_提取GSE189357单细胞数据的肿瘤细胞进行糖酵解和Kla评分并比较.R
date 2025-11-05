######目的：提取GSE189357单细胞数据的肿瘤细胞进行糖酵解和Kla评分并比较
######作者：申奥
######日期：2024-11-06
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
library(monocle)


# 1.读取45中最后保存的Epi对象，并提取肿瘤细胞----
sc.epi <- readRDS("5.scRNA/GSE189357_scEpi.rds")
table(sc.epi$TumorCell_pred)

sc.tumor <- subset(sc.epi, TumorCell_pred == "Tumor cells")


# 2.走seurat流程----
sc.tumor <- NormalizeData(sc.tumor, normalization.method = "LogNormalize", scale.factor = 10000)
sc.tumor <- FindVariableFeatures(sc.tumor, selection.method = "vst", nfeatures = 2000)
sc.tumor <- ScaleData(sc.tumor)


# 3.评分----
## 3.1 加载基因集----
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

## 3.2 Glycolysis AUCell评分----
cells_rankings <- AUCell_buildRankings(sc.tumor@assays$RNA@data, splitByBlocks=TRUE) 

cells_AUC <- AUCell_calcAUC(glyco, cells_rankings, nCores = 10,
                            aucMaxRank = ceiling(0.1 * nrow(cells_rankings))) # 选择基因排序中top10%的基因的rank位置

set.seed(333)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, assign=TRUE) #选定合适的阈值
auc_thr = cells_assignment$Glycolysis$aucThr$selected
auc_thr

#提取结果
geneSet <- "Glycolysis"
AUCell_auc <- as.numeric(getAUC(cells_AUC)[geneSet, ])
#添加至metadata中
sc.tumor$AUCell_Glycolysis <- AUCell_auc
head(sc.tumor@meta.data)

## 3.3 Kla AUCell评分----
cells_rankings <- AUCell_buildRankings(sc.tumor@assays$RNA@data, splitByBlocks=TRUE) # 与2.1相同

cells_AUC <- AUCell_calcAUC(upKla, cells_rankings, nCores = 10,
                            aucMaxRank = ceiling(0.1 * nrow(cells_rankings))) # 选择基因排序中top10%的基因的rank位置

set.seed(333)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, assign=TRUE)#选定合适的阈值
auc_thr = cells_assignment$H3K18la$aucThr$selected
auc_thr

#提取结果
geneSet <- "H3K18la"
AUCell_auc <- as.numeric(getAUC(cells_AUC)[geneSet, ])
#添加至metadata中
sc.tumor$AUCell_upKla <- AUCell_auc
head(sc.tumor@meta.data)

### 按upKla评分把肿瘤细胞进行分组
sc.tumor$TumorCell_group <- ifelse(sc.tumor$AUCell_upKla >= median(sc.tumor$AUCell_upKla), 
                                   "High-upKla_TumorCells", "Low-upKla_TumorCells")
table(sc.tumor$TumorCell_group)
# High-upKla_TumorCells  Low-upKla_TumorCells 
# 2859                  2858
saveRDS(sc.tumor, file = "5.scRNA/GSE189357_tumorCell.rds")

## 3.4 可视化比较AUCell评分----
df <- sc.tumor@meta.data
class(df$Sample_group)
table(df$Sample_group)
# AIS  IAC  MIA 
# 449 4667  601
df$Sample_group2 <- ifelse(df$Sample_group == "AIS", "in situ", "invasive")
table(df$Sample_group2)
# in situ invasive 
# 449     5268
df$Sample_group2 <- factor(df$Sample_group2, levels = c("in situ", "invasive"))

my_comparisons <- list(c("invasive", "in situ"))
mycolors <- c("#79C9C7", "#6c92b8")

tumor_aucell_upKla_plot <- ggviolin(df, x="Sample_group2", y = "AUCell_upKla", width = 0.8, 
                              color = "black",#轮廓颜色
                              fill="Sample_group2",#填充
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
tumor_aucell_upKla_plot
ggsave(filename = "5.scRNA/GSE189357_tumor_aucell_upKla_plot.pdf", tumor_aucell_upKla_plot,
       height = 4, width = 4, dpi = 300)

tumor_aucell_glyco_plot <- ggviolin(df, x="Sample_group2", y = "AUCell_Glycolysis", width = 0.8, 
                              color = "black",#轮廓颜色
                              fill="Sample_group2",#填充
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
tumor_aucell_glyco_plot
ggsave(filename = "5.scRNA/GSE189357_tumor_aucell_glyco_plot.pdf", tumor_aucell_glyco_plot,
       height = 4, width = 4, dpi = 300)


# 4.拟时序分析----
#参考教程：https://mp.weixin.qq.com/s/Q8R2PGTUGzbRTjOpM2ovkg
## 4.1 数据预处理----
data <- as.matrix(sc.tumor@assays$RNA@counts)
data <- as(data, 'sparseMatrix')

pd <- new('AnnotatedDataFrame', data = sc.tumor@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)

## 4.2 创建Monocle2对象----
mycds <- newCellDataSet(data,
                        phenoData = pd,
                        featureData = fd,
                        expressionFamily = negbinomial.size())
### 预处理：计算数据的size factors和dispersions离散度
mycds <- estimateSizeFactors(mycds) #估计细胞的大小因子，以校正不同细胞的RNA测序深度
mycds <- estimateDispersions(mycds, cores = 4, relative_expr = TRUE) # 估计基因的离散度，以确定哪些基因的表达变化最大

### 质控：去除死细胞、空液滴、碎细胞液滴、双液滴，提高计算准确度、降低预算时间
mycds <- detectGenes(mycds, min_expr = 0.5) #过滤基因表达量小于等于阈值0.5的基因
head(fData(mycds))
#expressed_genes <- row.names(subset(fData(mycds), num_cells_expressed > nrow(sc.epi@meta.data) * 0.01)) #至少要在1%的细胞中表达
expressed_genes2 <- row.names(subset(fData(mycds), num_cells_expressed >= 10)) #至少要在10个细胞中表达
#mycds2 <- mycds[expressed_genes,]
# 查看筛选后的基因个数（用于后面基因的筛选）
#length(expressed_genes)
length(expressed_genes2)
length(row.names(sc.tumor))

## 4.3 构建轨迹----
### 4.3.1 挑选细胞间高度离散的基因----
disp_table <- dispersionTable(mycds)
#### 高度离散基因的筛选标准，可根据数据情况设置mean_expression的值
disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id #筛选出平均表达量不低于0.1，且离散度高于拟合离散度1倍的高变基因
mycds <- setOrderingFilter(mycds[expressed_genes2,], disp.genes)
nrow(mycds)
plot_ordering_genes(mycds) #深色点即为之后拟时序分析与可视化会用到的点，浅色点不会被涉及

### 4.3.2 细胞降维----
# 选取的基因数目为每个细胞的维度，基于默认的'DDRTree'方法进行数据降维
mycds <- reduceDimension(mycds, max_components = 2, method = 'DDRTree') #max_components:降维后的组件数，数据被降维到2个主要组件以用来可视化或进行进一步的分析

### 4.3.3 细胞排序----
# 对细胞进行排序，由于排序无法区分起点和终点，若分析所得时序与实际相反，根据“reverse”参数进行调整，默认reverse=F
mycds <- orderCells(mycds)

### 4.3.4 轨迹可视化----
#### 按细胞类型展示
plot_cell_trajectory(mycds,
                     color_by = "TumorCell_group", #与"seurat_clusters"效果基本等同
                     cell_size = 0.8, # 点大小
                     cell_link_size = 0.5) # 线粗细

plot_cell_trajectory(mycds, 
                     color_by = "TumorCell_group", cell_size = 0.8) + 
  facet_wrap("~TumorCell_group", ncol = 3)

#### 按state和拟时序展示
plot_cell_trajectory(mycds, color_by = "State")
plot_cell_trajectory(mycds, color_by = "Pseudotime", cell_size = 0.6) + 
  scale_color_gsea()
plot_cell_trajectory(mycds, color_by = "AUCell_Glycolysis", cell_size = 0.6) + 
  scale_color_gsea()
plot_cell_trajectory(mycds, color_by = "AUCell_upKla", cell_size = 0.6) + 
  scale_color_gsea()
