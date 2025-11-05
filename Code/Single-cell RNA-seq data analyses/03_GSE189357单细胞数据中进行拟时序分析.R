######目的：GSE189357单细胞数据中进行拟时序分析
######作者：申奥
######日期：2024-11-05
######Server: yjx@xiyoucloud.top:34022
######Version: R 4.3.3


# 0.环境设置----
rm(list = ls())
options(stringAsFactors = F)

library(Seurat, lib.loc = "~/bioSoft/seurat_v4/") #‘4.4.0’
#devtools::install_github("cole-trapnell-lab/monocle-release@develop")
library(monocle) #‘2.34.0’
library(slingshot) #‘2.14.0’


# 1.读取45中最后保存的Epi对象----
sc.epi <- readRDS("5.scRNA/GSE189357_scEpi.rds")


# 2.monocle2----
#参考教程：https://mp.weixin.qq.com/s/Q8R2PGTUGzbRTjOpM2ovkg
## 2.1 数据预处理----
data <- as.matrix(sc.epi@assays$RNA@counts)
data <- as(data, 'sparseMatrix')

pd <- new('AnnotatedDataFrame', data = sc.epi@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)

## 2.2 创建Monocle2对象----
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
length(expressed_genes)
length(expressed_genes2)
length(row.names(sc.epi))

## 2.3 构建轨迹----
### 2.3.1 挑选细胞间高度离散的基因----
disp_table <- dispersionTable(mycds)
#### 高度离散基因的筛选标准，可根据数据情况设置mean_expression的值
disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id #筛选出平均表达量不低于0.1，且离散度高于拟合离散度1倍的高变基因
mycds <- setOrderingFilter(mycds[expressed_genes2,], disp.genes)
nrow(mycds)
plot_ordering_genes(mycds) #深色点即为之后拟时序分析与可视化会用到的点，浅色点不会被涉及

### 2.3.2 细胞降维----
# 选取的基因数目为每个细胞的维度，基于默认的'DDRTree'方法进行数据降维
mycds <- reduceDimension(mycds, max_components = 2, method = 'DDRTree') #max_components:降维后的组件数，数据被降维到2个主要组件以用来可视化或进行进一步的分析

### 2.3.3 细胞排序----
# 对细胞进行排序，由于排序无法区分起点和终点，若分析所得时序与实际相反，根据“reverse”参数进行调整，默认reverse=F
mycds <- orderCells(mycds)

### 2.3.4 轨迹可视化----
#### 按细胞类型展示
plot_cell_trajectory(mycds,
                     color_by = "TumorCell_pred", #与"seurat_clusters"效果基本等同
                     cell_size = 0.8, # 点大小
                     cell_link_size = 0.5) # 线粗细

plot_cell_trajectory(mycds, 
                     color_by = "TumorCell_pred", cell_size = 0.8) + 
  facet_wrap("~TumorCell_pred", ncol = 3)

#### 按state和拟时序展示
plot_cell_trajectory(mycds, color_by = "State")
plot_cell_trajectory(mycds, color_by = "Pseudotime", cell_size = 0.6) + 
  scale_color_gsea()
plot_cell_trajectory(mycds, color_by = "AUCell_Glycolysis", cell_size = 0.6) + 
  scale_color_gsea()
plot_cell_trajectory(mycds, color_by = "AUCell_upKla", cell_size = 0.6) + 
  scale_color_gsea()


#### 沿时间轴的细胞密度图
library(ggpubr) 

df <- pData(mycds)  
# pData(cds)取出的是cds对象中cds@phenoData@data的内容

ggplot(df, 
       aes(Pseudotime, 
           colour = TumorCell_pred, 
           fill = TumorCell_pred)) +   
  geom_density( # 绘制密度图
    bw = 0.5,  # 设定带宽，影响密度曲线的平滑程度
    size = 0.8, # 密度曲线的粗细
    alpha = 0.5) +  # 填充颜色的透明度
  theme_classic2() # 取消网格背景


# 3.slingshot----
library(RColorBrewer)
sim <- as.SingleCellExperiment(sc.epi)

sim <- slingshot(data = sim, 
                 clusterLabels = 'TumorCell_pred', 
                 reducedDim = 'UMAP',
                 start.clus = "Normal epithelials", # 可指定起点亚群
                 end.clus = NULL # 可指定终点亚群
)

# 可以通过“type”参数查看基于聚类的最小生成树最初是如何估计谱系结构的
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plot(reducedDims(sim)$UMAP, col = colors, pch=16, asp=1)
lines(SlingshotDataSet(sim), lwd=2, col = 'black', type = 'lineages')
