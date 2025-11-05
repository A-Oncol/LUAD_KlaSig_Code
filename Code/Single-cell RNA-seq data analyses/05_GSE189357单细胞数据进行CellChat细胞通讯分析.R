######目的：GSE189357单细胞数据进行CellChat细胞通讯分析
######作者：申奥
######日期：2024-11-10
######Server: yjx@xiyoucloud.top:34022
######Version: R 4.3.3


# 0.环境设置----
rm(list = ls())
options(stringAsFactors = F)

library(Seurat, lib.loc = "~/bioSoft/seurat_v4/")
library(CellChat) # ‘2.1.2’
library(ggplot2)
library(ggpubr)
library(ggsci)
library(cols4all)


# 1.读取并整理单细胞数据----
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


# 2.CellChat----
## 2.1 构建CellChat对象----
data_input <- sc_obj@assays$RNA@data
meta_data <- sc_obj@meta.data
meta_data = meta_data[!is.na(meta_data$Cell_types),]
data_input = data_input[,row.names(meta_data)]

identical(colnames(data_input), rownames(meta_data))

cellchat <- createCellChat(object = data_input, 
                           meta = meta_data, 
                           group.by = "Cell_types")
levels(cellchat@idents)
unique(cellchat@idents)
table(cellchat@idents)
groupSize <- as.numeric(table(cellchat@idents))

## 2.2 CellChat流程分析----
### 加载CellChatDB数据库
CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
cellchat@DB <- CellChatDB.use
### 提取数据库支持的数据子集
cellchat <- subsetData(cellchat)
### 识别过表达基因
#future::plan("multisession", workers = 4) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
### 识别配体-受体对
cellchat <- identifyOverExpressedInteractions(cellchat)
### 将配体、受体投射到PPI网络
#cellchat <- projectData(cellchat, PPI.human)
### 计算通讯概率，推断细胞通讯网络
cellchat <- computeCommunProb(cellchat, raw.use = T, population.size = F) # 没有进行projectData则raw.use设置为T；population.size = T表示不是分选
### 过滤
cellchat <- filterCommunication(cellchat, min.cells = 10)
### 在信号通路水平推断细胞通讯
cellchat <- computeCommunProbPathway(cellchat)
head(cellchat@net)
head(cellchat@netP)
### 提取预测的细胞通讯网络
df.net <- subsetCommunication(cellchat)
write.table(df.net, file = "5.scRNA/CellChat_GSE189357/GSE189357_CellNet.txt", quote = F, sep = "\t", row.names = F)
### 计算加和的cell-cell通讯网络
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # 计算网络中心性，slot.name = "netP"表示使用聚合后的通信概率网络进行计算，中心性用于评估各细胞在通信网络中的重要性

save(cellchat, file = "5.scRNA/CellChat_GSE189357/cellchat_obj.RData")
saveRDS(cellchat, file = "5.scRNA/CellChat_GSE189357/cellchat_obj.rds")

## 2.3 可视化----
### 2.3.1 circle图----
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, 
                 weight.scale = T, label.edge = F, 
                 title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize,
                 weight.scale = T, label.edge= F,
                 title.name = "Interaction weights/strength") # 高6 宽12

# netVisual_circle(cellchat@net$count, vertex.weight = groupSize, 
#                  weight.scale = T, label.edge= T, sources.use = 'High-upKla_TumorCells',
#                  title.name = "Number of interactions")
# netVisual_circle(cellchat@net$count, vertex.weight = groupSize, 
#                  weight.scale = T, label.edge= T, sources.use = 'Low-upKla_TumorCells',
#                  title.name = "Number of interactions")
# 
# netVisual_circle(cellchat@net$weight, vertex.weight = groupSize,
#                  weight.scale = T, label.edge= F, sources.use = 'High-upKla_TumorCells',
#                  title.name = "Interaction weights/strength")
# netVisual_circle(cellchat@net$weight, vertex.weight = groupSize,
#                  weight.scale = T, label.edge= F, sources.use = 'Low-upKla_TumorCells',
#                  title.name = "Interaction weights/strength")

### 2.3.2 bubble图----
# netVisual_bubble(cellchat, sources.use = c("High-upKla_TumorCells"), angle.x = 45)
# netVisual_bubble(cellchat, sources.use = c("Low-upKla_TumorCells"), angle.x = 45)

bubble_plot <- netVisual_bubble(cellchat, 
                                sources.use = c("High-upKla_TumorCells", "Low-upKla_TumorCells"),
                                angle.x = 45) +
  theme(axis.text = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12))
bubble_plot
ggsave(filename = "5.scRNA/CellChat_GSE189357/GSE189357_CellChat_bubblePlot.pdf", bubble_plot,
       height = 6.5, width = 10, dpi = 300)

### 2.3.3 基于通路的图----
cellchat@netP$pathways
pathways.show <- "MIF"

netVisual_heatmap(cellchat, signaling = pathways.show) # 高4 长5.5
# heatmap_plot
# ggsave(filename = "5.scRNA/CellChat_GSE189357/GSE189357_CellChat_heatmapPlot.pdf", heatmap_plot,
#        height = 4, width = 4, dpi = 300)

#netAnalysis_contribution(cellchat, signaling = pathways.show)
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, color.heatmap = "Blues") # 高4 宽4

#netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
plotGeneExpression(cellchat, signaling = pathways.show, enriched.only = T, type = "violin") # 高4 长5.5

#### 层次结构图
cellchat@netP$pathways
levels(cellchat@idents) 
vertex.receiver = seq(1:2) # a numeric vector
netVisual_aggregate(cellchat, signaling = pathways.show,
                    vertex.receiver = vertex.receiver,layout= "hierarchy")
