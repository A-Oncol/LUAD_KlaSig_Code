######目的：GSE189357单细胞数据处理
######作者：申奥
######日期：2024-10-31
######Server: yjx@xiyoucloud.top:34022
######Version: R 4.3.3


# 0.环境设置----
rm(list = ls())
options(stringAsFactors = F)

library(data.table)
library(Seurat, lib.loc = "~/bioSoft/seurat_v4/")
library(harmony)
library(ggplot2)
library(cols4all)
library(ggsci)
library(TCfinder)


# 1.读取GSE189357单细胞数据并进行整合和分组----
samples <- list.dirs("1.data/GSE189357/", full.names = F)
samples <- samples[-1]

sceList <- lapply(samples, function(samplei){
  print(samplei)
  sce <- CreateSeuratObject(counts = Read10X(paste0("1.data/GSE189357/", samplei)),
                            project = samplei,
                            min.cells = 5, min.features = 300)
  
  return(sce)
})
names(sceList) =  samples
which(sapply(sceList, is.null))

sce.all <- merge(sceList[[1]], 
                 y = sceList[-1] ,
                 add.cell.ids = samples) 
sce.all$Samples <- sce.all$orig.ident
sce.all$Patient = stringr::str_split(sce.all$Samples,'[_]', simplify = T)[,2]
table(sce.all$Patient)

sce.all$Sample_group <- sce.all$Patient
sce.all$Sample_group = gsub("TD5|TD7|TD8", "AIS", sce.all$Sample_group) 
sce.all$Sample_group = gsub("TD3|TD4|TD6", "MIA", sce.all$Sample_group) 
sce.all$Sample_group = gsub("TD1|TD2|TD9", "IAC", sce.all$Sample_group) 
table(sce.all$Sample_group)


# 2.Seurat流程----
sc <- sce.all
sc[["percent.mt"]] <- PercentageFeatureSet(sc, pattern = "^MT-")
sc[["percent.rp"]] <- PercentageFeatureSet(sc, pattern = "^RP[SL]")
sc[["percent.hb"]] <- PercentageFeatureSet(sc, pattern = "^HB[^(P)]")
VlnPlot(sc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rp", "percent.hb"), ncol = 3)

scMerge <- subset(sc, subset = nFeature_RNA > 300 & nFeature_RNA < 10000 & percent.mt < 10 &
                    nCount_RNA < quantile(nCount_RNA, 0.95) & nCount_RNA > 1000 &  # 每个细胞测序的UMI count含量大于1000，且剔除最大的前5%的细胞
                    percent.hb < 10) 

scMerge <- NormalizeData(scMerge, normalization.method = "LogNormalize", scale.factor = 10000)
scMerge <- FindVariableFeatures(scMerge, selection.method = "vst", nfeatures = 2000)
scMerge <- ScaleData(scMerge)
scMerge <- RunPCA(scMerge, features = VariableFeatures(object = scMerge))
ElbowPlot(scMerge, ndims = 50)

scMerge <- RunHarmony(scMerge, group.by.vars = "Samples")

dim.usage <- 20
scMerge <- FindNeighbors(scMerge, dims = 1:dim.usage, reduction = "harmony") %>%
  FindClusters(resolution = 0.3) %>%
  RunUMAP(dims = 1:dim.usage, reduction = "harmony")
table(scMerge$RNA_snn_res.0.3) # 选用0.3

zzm60colors <- c('#4b6aa8','#3ca0cf','#c376a7','#ad98c3','#cea5c7',
                 '#53738c','#a5a9b0','#a78982','#696a6c','#92699e',
                 '#d69971','#df5734','#6c408e','#ac6894','#d4c2db',
                 '#537eb7','#83ab8e','#ece399','#405993','#cc7f73',
                 '#b95055','#d5bb72','#bc9a7f','#e0cfda','#d8a0c0',
                 '#e6b884','#b05545','#d69a55','#64a776','#cbdaa9',
                 '#efd2c9','#da6f6d','#ebb1a4','#a44e89','#a9c2cb',
                 '#b85292','#6d6fa0','#8d689d','#c8c7e1','#d25774',
                 '#c49abc','#927c9a','#3674a2','#9f8d89','#72567a',
                 '#63a3b8','#c4daec','#61bada','#b7deea','#e29eaf',
                 '#4490c4','#e6e2a3','#de8b36','#c4612f','#9a70a8',
                 '#76a2be','#408444','#c6adb0','#9d3b62','#2d3462')
mycolors <- c(pal_d3("category20")(20),
              zzm60colors)

p1 <- DimPlot(scMerge, reduction = "umap", pt.size = 0.8, label = T, label.box = F,
              cols = mycolors, group.by = "RNA_snn_res.0.3") +
  tidydr::theme_dr(xlength = 0.2, 
                   ylength = 0.2,
                   arrow = arrow(length = unit(0.1, "inches"),type = "closed")) + 
  ggtitle("Clusters") +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2, hjust = 0.01),
        legend.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.6, face = "bold", size = 14))
p1
ggsave(filename = '5.scRNA/GSE189357_Clusters_UMAP.pdf', p1, height = 4.5, width = 5, dpi = 300)

p2 <- DimPlot(scMerge, reduction = "umap", pt.size = 0.6, label = F, label.box = F,
              cols = mycolors, group.by = "Patient") +
  tidydr::theme_dr(xlength = 0.2, 
                   ylength = 0.2,
                   arrow = arrow(length = unit(0.1, "inches"),type = "closed")) + 
  ggtitle("Samples") +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2, hjust = 0.01),
        legend.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.6, face = "bold", size = 14))
p2
ggsave(filename = '5.scRNA/GSE189357_Samples_UMAP.pdf', p2, height = 4.5, width = 5, dpi = 300)

p3 <- DimPlot(scMerge, reduction = "umap", pt.size = 0.6, label = F, label.box = F,
              cols = mycolors, group.by = "Sample_group") +
  tidydr::theme_dr(xlength = 0.2, 
                   ylength = 0.2,
                   arrow = arrow(length = unit(0.1, "inches"),type = "closed")) + 
  ggtitle("Sample group") +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2, hjust = 0.01),
        legend.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.6, face = "bold", size = 14))
p3
ggsave(filename = '5.scRNA/GSE189357_SampleGroup_UMAP.pdf', p3, height = 4.5, width = 5, dpi = 300)


genes_to_check = c('EPCAM','KRT19','CDH1',  #上皮
                   'PECAM1' , #基质
                   'VWF',  'CDH5', 'PLVAP', #内皮
                   'COL1A1', 'COL1A2', 'DCN', 'ACTA2', #成纤维
                   'PTPRC', # immune
                   'CD3D', 'CD3E', 'CD8A', 'CD4','CD2', #T
                   'MS4A1','CD19', 'CD79A',  #B
                   'JCHAIN', 'MZB1', # plasma
                   'NKG7', "GNLY", "KLRD1",  #NK
                   'CSF1R', 'CD68', #髓系
                   'C1QC','C1QB','LYZ',  #巨噬
                   'CD14', 'S100A8', 'S100A9', # mono
                   'CD86', 'CD80', 'XCR1', 'LILRA4', 'HLA-DQA1',
                   'TPSAB1','TPSB2', "KIT",#肥大
                   'MKI67', 'TOP2A', 'STMN1', 'PCNA' #增殖
)

p4 <- DotPlot(scMerge, features = genes_to_check) +
  scale_color_continuous_c4a_seq("wes.zissou1") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
#coord_flip()
p4
ggsave(filename = '5.scRNA/GSE189357_MarkersAnnot_dotplot.pdf', p4, height = 4.5, width = 12, dpi = 300)

# degs <- FindAllMarkers(scMerge, 
#                        min.pct = 0.25, # 基因表达在25%以上的细胞中
#                        logfc.threshold = 0.25, # Log fold change阈值
#                        only.pos = T) 
# deg12 <- degs[degs$cluster == 12,]

celltype <- c('T',
              'NK', 'B', 'DC', 'Epithelial', 'Mast',
              'Macrophages', 'Endothelial', 'Fibroblasts', 'Monocytes', 'Plasma',
              'Epithelial', 'Proliferative', 'Epithelial', 'Epithelial', 'Macrophages',
              'Mast', 'T')
Idents(scMerge) <- scMerge@meta.data$RNA_snn_res.0.3
names(celltype) <- levels(scMerge)
scMerge <- RenameIdents(scMerge, celltype)
scMerge@meta.data$Cell_Type <- Idents(scMerge)
Idents(scMerge) <- scMerge@meta.data$Cell_Type
table(scMerge@meta.data$Cell_Type)
saveRDS(scMerge, file = "5.scRNA/GSE189357_allSamples_allCells_annot.rds")

p5 <- DimPlot(scMerge, reduction = "umap", pt.size = 0.6, label = F, label.box = F,
              cols = mycolors) +
  tidydr::theme_dr(xlength = 0.2, 
                   ylength = 0.2,
                   arrow = arrow(length = unit(0.1, "inches"),type = "closed")) + 
  ggtitle("Major cell types") +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2, hjust = 0.01),
        legend.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.6, face = "bold", size = 14))
p5
ggsave(filename = '5.scRNA/GSE189357_CellTypes_UMAP.pdf', p5, height = 4.5, width = 5.5, dpi = 300)

save(p1,p2,p3,p4,p5, file = "5.scRNA/GSE189357_SeuratPlots.RData")


# 3.提取Epi，推断肿瘤细胞----
table(scMerge$Cell_Type)
# T            NK             B            DC    Epithelial          Mast   Macrophages 
# 26427         17264         10067          7748         10318          7778          7486 
# Endothelial   Fibroblasts     Monocytes        Plasma Proliferative 
# 4141          3067          2626          2548           842
sc.epi <- subset(scMerge, Cell_Type %in% "Epithelial")
dim(sc.epi)

sc.epi <- NormalizeData(sc.epi, normalization.method = "LogNormalize", scale.factor = 10000)
sc.epi <- FindVariableFeatures(sc.epi, selection.method = "vst", nfeatures = 2000)
sc.epi <- ScaleData(sc.epi)
sc.epi <- RunPCA(sc.epi, features = VariableFeatures(object = sc.epi))
ElbowPlot(sc.epi, ndims = 50)
#DimHeatmap(sc.epi, dims = 30, cells = 500, balanced = TRUE)

sc.epi <- RunHarmony(sc.epi, group.by.vars = 'Samples')

dim.usage <- 20
sc.epi <- FindNeighbors(sc.epi, dims = 1:dim.usage, reduction = "harmony") %>%
  FindClusters(resolution = 0.3) %>%
  RunUMAP(dims = 1:dim.usage, reduction = "harmony")
table(sc.epi$RNA_snn_res.0.3)

DimPlot(sc.epi, reduction = "umap", pt.size = 0.6, label = T, label.box = F,
        cols = mycolors)
DimPlot(sc.epi, reduction = "umap", pt.size = 0.6, label = F, label.box = F,
        cols = mycolors, group.by = "Samples")
DimPlot(sc.epi, reduction = "umap", pt.size = 0.6, label = F, label.box = F,
        cols = mycolors, group.by = "Sample_group")

## 3.1 TCfinder----
library(reticulate)
reticulate::use_python("/home/yjx/anaconda3/envs/TCfinder/bin/python")
reticulate::py_config()

expr_data <- sc.epi@assays$RNA@counts
head(expr_data)[1:4, 1:8]

norm_data <- data_normalized(expr_data = expr_data,
                             method = "method", genome = "hg38")
head(norm_data)[1:4, 1:8]
colnames(norm_data) <- colnames(expr_data)
score_data <- pathway_score(expr_data = norm_data, normalized = T,
                            method = "method", genome = "hg38")

TCfinder_predict_results <- predict_cell(path_score = score_data)
table(TCfinder_predict_results$cell_type)

save(TCfinder_predict_results, file = "5.scRNA/TCfinder_GSE189357/TCfinder_predict_results.RData")
save(score_data, file = "5.scRNA/TCfinder_GSE189357/score_data.RData")
save(norm_data, file = "5.scRNA/TCfinder_GSE189357/norm_data.RData")

## 3.2 copyKat----
library(copykat)

copykat.res <- copykat(
  rawmat = sc.epi@assays$RNA@counts,    # 原始表达矩阵，TPM数据也能使用
  id.type = "S",           # 基因名称的类型，S代表symbol
  cell.line = "no",
  ngene.chr = 3,           # 过滤细胞的标准，要求每条染色体上最少有5个基因，实测会过滤不少细胞
  win.size = 25,           # 检测窗口的基因数，默认25个基因，可以在15-150的范围内选择
  KS.cut = 0.05,           # 差异检验的阈值，越大灵敏度越低
  sam.name = "LUAD",       # 样本名称，copykat运行的结果文件会以此为前缀
  distance = "euclidean",  # 聚类的距离，可选c("euclidean", "pearson", "spearman")
  norm.cell.names = "",    # 正常细胞的名称向量
  n.cores = 40)            # 线程数量
save(copykat.res, file = "5.scRNA/copyKat_GSE189357/copyKat_result.RData")

## 3.3 inferCNV----
library(infercnv)

### 筛选分析需要的细胞类型
table(scMerge$Cell_Type)
sc.cnv <- subset(scMerge, idents = c('Epithelial', 'Endothelial', 'T'))
counts <- GetAssayData(sc.cnv, slot = 'counts')
sc.cnv$barcode <- colnames(sc.cnv)

a <- sc.cnv@meta.data[,c('barcode', 'Cell_Type')]
a$Cell_Type=as.character(a$Cell_Type)

sc.epi$barcode <- rownames(sc.epi@meta.data)
b <- sc.epi@meta.data[,c("barcode", "seurat_clusters")]
dd <- dplyr::left_join(a, b, by = "barcode")
dd$seurat_clusters <- as.character(dd$seurat_clusters)
dd[is.na(dd$seurat_clusters),3] <- dd[is.na(dd$seurat_clusters),2]
rownames(dd) = dd$barcode

dd <- dd[,-2]
table(dd$seurat_clusters)
# 0           1          10          11          12          13           2           3           4 
# 1835        1358         254         195         191          51        1332        1179        1032 
# 5           6           7           8           9 Endothelial           T 
# 966         587         576         428         334        4141       26427
identical(sc.cnv@meta.data$barcode, dd$barcode)
dd <- as.data.frame(dd)
sc.cnv$anno <- dd$seurat_clusters
Idents(sc.cnv) <- sc.cnv$anno
anno <- data.frame(Idents(sc.cnv))

# 文件制作方式见https://github.com/broadinstitute/inferCNV/wiki/instructions-create-genome-position-file
# python gtf_to_position_file.py --attribute_name "gene_name" gencode.v43.annotation.gtf human_v43gtf_genes_pos.txt
gene_order <- "human_v43gtf_genes_pos.txt"

infercnv_obj = CreateInfercnvObject(raw_counts_matrix = counts,
                                    annotations_file = anno,
                                    delim = "\t",
                                    gene_order_file = gene_order,
                                    min_max_counts_per_cell = c(100, +Inf),
                                    ref_group_names = c("Endothelial", "T"))
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff = 0.1,
                             out_dir = "5.scRNA/inferCNV_GSE189357",
                             cluster_by_groups = T, plot_steps = F,
                             HMM = FALSE,
                             denoise = TRUE,
                             num_threads = 40,
                             write_expr_matrix = T, output_format = "pdf",
                             analysis_mode = "subclusters")  
saveRDS(infercnv_obj, file = "5.scRNA/inferCNV_GSE189357/infercnv_obj.rds")

### 计算CNV score
infer_CNV_obj <- readRDS('5.scRNA/inferCNV_GSE189357/run.final.infercnv_obj')
expr <- infer_CNV_obj@expr.data
expr[1:4,1:4]
#### 第一种方法：平方和（参考来源：https://blog.csdn.net/qq_38774801/article/details/115841319）
expr2=expr-1
expr2=expr2 ^ 2
CNV_score=as.data.frame(colMeans(expr2))
colnames(CNV_score)="CNV_score"
identical(rownames(sc.cnv@meta.data), rownames(CNV_score))
CNV_score$Cell_clusters <- sc.cnv$anno
table(CNV_score$Cell_clusters)
class(CNV_score$Cell_clusters)
CNV_score$Cell_clusters <- factor(CNV_score$Cell_clusters, levels = c(seq(0,13,1), "T", "Endothelial"))
ggplot(CNV_score, aes(x = Cell_clusters, y = CNV_score, fill = Cell_clusters)) +
  geom_boxplot() +
  scale_fill_manual(values = mycolors) +
  theme_classic(base_size = 14) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(size = 12, angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 12)) +
  ylab("CNV score") + xlab("Cell clusters")

#### 第二种方法：Jimmy计算CNV score（参考来源：https://mp.weixin.qq.com/s/Jb1NhcvNK17MLN3g6HX4SQ; https://mp.weixin.qq.com/s/LzIQMkP6SqhmDIqSvQPyvg）
if(T){
  tmp1 = expr[,infer_CNV_obj@reference_grouped_cell_indices$`T`]
  tmp2 = expr[,infer_CNV_obj@reference_grouped_cell_indices$Endothelial]
  tmp= cbind(tmp1,tmp2)
  down=mean(rowMeans(tmp)) - 2 * mean( apply(tmp, 1, sd))
  up=mean(rowMeans(tmp)) + 2 * mean( apply(tmp, 1, sd))
  oneCopy=up-down
  oneCopy
  a1= down- 2*oneCopy
  a2= down- 1*oneCopy
  down;up
  a3 = up +  1*oneCopy
  a4 = up + 2*oneCopy 
  
  cnv_score_table<-infer_CNV_obj@expr.data
  cnv_score_table[1:4,1:4]
  cnv_score_mat <- as.matrix(cnv_score_table)
  
  # Scoring
  cnv_score_table[cnv_score_mat > 0 & cnv_score_mat < a2] <- "A" #complete loss. 2pts
  cnv_score_table[cnv_score_mat >= a2 & cnv_score_mat < down] <- "B" #loss of one copy. 1pts
  cnv_score_table[cnv_score_mat >= down & cnv_score_mat <  up ] <- "C" #Neutral. 0pts
  cnv_score_table[cnv_score_mat >= up  & cnv_score_mat <= a3] <- "D" #addition of one copy. 1pts
  cnv_score_table[cnv_score_mat > a3  & cnv_score_mat <= a4 ] <- "E" #addition of two copies. 2pts
  cnv_score_table[cnv_score_mat > a4] <- "F" #addition of more than two copies. 2pts
  
  # Check
  table(cnv_score_table[,1])
  # Replace with score 
  cnv_score_table_pts <- cnv_score_mat
  rm(cnv_score_mat)
  # 
  cnv_score_table_pts[cnv_score_table == "A"] <- 2
  cnv_score_table_pts[cnv_score_table == "B"] <- 1
  cnv_score_table_pts[cnv_score_table == "C"] <- 0
  cnv_score_table_pts[cnv_score_table == "D"] <- 1
  cnv_score_table_pts[cnv_score_table == "E"] <- 2
  cnv_score_table_pts[cnv_score_table == "F"] <- 2
  
  cnv_score_table_pts[1:4,1:4]
  str(  as.data.frame(cnv_score_table_pts[1:4,1:4])) 
  cell_scores_CNV <- as.data.frame(colSums(cnv_score_table_pts))
  
  colnames(cell_scores_CNV) <- "cnv_score" 
}

head(cell_scores_CNV)
score=cell_scores_CNV
head(score)
meta <- sc.cnv@meta.data
identical(rownames(meta), rownames(score))
meta$totalCNV = score[match(colnames(sc.cnv),
                            rownames(score)),1]
class(meta$anno)
meta$anno <- factor(meta$anno, levels = c(seq(0,13,1), "T", "Endothelial"))
save(meta, file = "5.scRNA/inferCNV_GSE189357/infercnv_totalCNVscore.RData") #两个方法计算的趋势差不多，采用第二种方法

cnvscore_plot <- ggplot(meta, aes(x = anno, y = totalCNV, fill = anno)) +
  geom_boxplot() +
  scale_fill_manual(values = mycolors) +
  theme_classic(base_size = 14) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(size = 12, angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 12)) +
  ylab("CNV score") + xlab("Cell clusters")
cnvscore_plot
ggsave(filename = "5.scRNA/inferCNV_GSE189357/infercnv_totalCNVscore_plot.pdf",
       cnvscore_plot, height = 4, width = 5, dpi = 300)

## 3.4 把推断结果重新添加到总的meta.data中----
### 3.4.1 copyKAT结果----
pred <- data.table::fread("5.scRNA/copyKat_GSE189357/LUAD_copykat_prediction.txt")
table(pred$copykat.pred)
# aneuploid     diploid not.defined 
# 3294        6642         382 
pred <- as.data.frame(pred)
pred$copykat.pred <- ifelse(pred$copykat.pred == "aneuploid", "Tumor cells",
                            ifelse(pred$copykat.pred == "diploid", "Normal epithelials", "not.defined"))
table(pred$copykat.pred)
metadata_cellorder <- rownames(sc.epi@meta.data)
rownames(pred) <- pred$cell.names
pred <- pred[metadata_cellorder,]
identical(pred$cell.names, colnames(sc.epi))

Idents(sc.epi, cells = colnames(sc.epi)) <- pred$copykat.pred
identical(rownames(pred), rownames(sc.epi@meta.data))
sc.epi$copycat_pred <- pred$copykat.pred
identical(rownames(sc.epi@meta.data), pred$cell.names)
identical(sc.epi@meta.data$copycat_pred, pred$copykat.pred)
table(sc.epi$copycat_pred)

Idents(sc.epi) <- sc.epi$copycat_pred
Idents(scMerge, cells = colnames(sc.epi)) <- Idents(sc.epi)
scMerge$Cells <- Idents(scMerge)

saveRDS(scMerge, file = "5.scRNA/GSE189357_copyKAT_ref.rds")

### 3.4.2 inferCNV结果----
table(sc.epi$RNA_snn_res.0.3)
sc.epi$inferCNV_pred <- ifelse(sc.epi$RNA_snn_res.0.3 %in% c(1,9,10,12), "Normal epithelials", "Tumor cells")
table(sc.epi$inferCNV_pred)

Idents(sc.epi) <- sc.epi$inferCNV_pred
Idents(scMerge, cells = colnames(sc.epi)) <- Idents(sc.epi)
scMerge$Cells2 <- Idents(scMerge)

saveRDS(scMerge, file = "5.scRNA/GSE189357_inferCNV_ref.rds")

### 3.4.3 TCfinder结果----
table(TCfinder_predict_results$cell_type)
TCfinder_predict_results$cell_type <- ifelse(TCfinder_predict_results$cell_type == "tumor", "Tumor cells", "Normal epithelials")
table(TCfinder_predict_results$cell_type)

identical(TCfinder_predict_results$barcode, rownames(sc.epi@meta.data))
sc.epi$TCfinder_pred <- TCfinder_predict_results$cell_type
table(sc.epi$TCfinder_pred)

Idents(sc.epi) <- sc.epi$TCfinder_pred
Idents(scMerge, cells = colnames(sc.epi)) <- Idents(sc.epi)
scMerge$Cells3 <- Idents(scMerge)

saveRDS(scMerge, file = "5.scRNA/GSE189357_TCfinder_ref.rds")

### 3.4.4 判断最终结果----
#### 至少2个算法判断为肿瘤细胞的才是最终的肿瘤细胞
library(dplyr)

tumor_res <- sc.epi@meta.data[,c("copycat_pred", "inferCNV_pred", "TCfinder_pred")]
colnames(tumor_res)[1] <- "copyKat_pred"
table(tumor_res$copyKat_pred)
table(tumor_res$inferCNV_pred)
table(tumor_res$TCfinder_pred)

tumor_res <- tumor_res %>% 
  mutate(final_pred = ifelse(rowSums(select(., copyKat_pred, inferCNV_pred, TCfinder_pred) == "Tumor cells") >= 2, "Tumor cells", "Normal epithelials"))
table(tumor_res$final_pred)
# Normal epithelials        Tumor cells 
# 4601               5717

identical(rownames(tumor_res), rownames(sc.epi@meta.data))
sc.epi$TumorCell_pred <- tumor_res$final_pred

Seurat::Idents(sc.epi) <- sc.epi$TumorCell_pred
Seurat::Idents(scMerge, cells = colnames(sc.epi)) <- Seurat::Idents(sc.epi)
scMerge$Cell_types <- Seurat::Idents(scMerge)
table(scMerge$Cell_types)

saveRDS(scMerge, file = "5.scRNA/GSE189357_final_ref.rds")
saveRDS(sc.epi, file = "5.scRNA/GSE189357_epi.rds")

## 3.5 UMAP可视化----
Idents(sc.epi) <- sc.epi$RNA_snn_res.0.3
p6 <- DimPlot(sc.epi, reduction = "umap", pt.size = 0.6, label = F, label.box = F,
              cols = mycolors) +
  tidydr::theme_dr(xlength = 0.2, 
                   ylength = 0.2,
                   arrow = arrow(length = unit(0.1, "inches"),type = "closed")) + 
  ggtitle("Clusters") +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2, hjust = 0.01),
        legend.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.6, face = "bold", size = 14))
p6
ggsave(filename = '5.scRNA/GSE189357_Epi_Clusters_UMAP.pdf', p6, height = 4.5, width = 5, dpi = 300)

p7 <- DimPlot(sc.epi, reduction = "umap", pt.size = 0.6, label = F, label.box = F,
              cols = c(mycolors[1], mycolors[3], mycolors[2]), group.by = "copycat_pred") +
  tidydr::theme_dr(xlength = 0.2, 
                   ylength = 0.2,
                   arrow = arrow(length = unit(0.1, "inches"),type = "closed")) + 
  ggtitle("CopyKat") +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2, hjust = 0.01),
        legend.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.6, face = "bold", size = 14))
p7
ggsave(filename = '5.scRNA/GSE189357_Epi_CopyKat_UMAP.pdf', p7, height = 4.5, width = 6, dpi = 300)

p8 <- DimPlot(sc.epi, reduction = "umap", pt.size = 0.6, label = F, label.box = F,
              cols = mycolors, group.by = "inferCNV_pred") +
  tidydr::theme_dr(xlength = 0.2, 
                   ylength = 0.2,
                   arrow = arrow(length = unit(0.1, "inches"),type = "closed")) + 
  ggtitle("inferCNV") +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2, hjust = 0.01),
        legend.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.6, face = "bold", size = 14))
p8
ggsave(filename = '5.scRNA/GSE189357_Epi_inferCNV_UMAP.pdf', p8, height = 4.5, width = 6, dpi = 300)

p9 <- DimPlot(sc.epi, reduction = "umap", pt.size = 0.6, label = F, label.box = F,
              cols = mycolors, group.by = "TCfinder_pred") +
  tidydr::theme_dr(xlength = 0.2, 
                   ylength = 0.2,
                   arrow = arrow(length = unit(0.1, "inches"),type = "closed")) + 
  ggtitle("TCfinder") +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2, hjust = 0.01),
        legend.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.6, face = "bold", size = 14))
p9
ggsave(filename = '5.scRNA/GSE189357_Epi_TCfinder_UMAP.pdf', p9, height = 4.5, width = 6, dpi = 300)

p10 <- DimPlot(sc.epi, reduction = "umap", pt.size = 0.6, label = F, label.box = F,
              cols = mycolors, group.by = "TumorCell_pred") +
  tidydr::theme_dr(xlength = 0.2, 
                   ylength = 0.2,
                   arrow = arrow(length = unit(0.1, "inches"),type = "closed")) + 
  ggtitle("Tumor cell prediction") +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2, hjust = 0.01),
        legend.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.6, face = "bold", size = 14))
p10
ggsave(filename = '5.scRNA/GSE189357_Epi_FinalPred_UMAP.pdf', p10, height = 4.5, width = 6, dpi = 300)

save(p6,p7,p8,p9,p10, file = "5.scRNA/GSE189357_Epi_SeuratPlots.RData")
