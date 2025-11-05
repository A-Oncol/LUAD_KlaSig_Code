######目的：P24的SpaCET分析
######作者：申奥
######日期：2024-09-27
######Server: yjx@xiyoucloud.top:34022
######Version: R 4.3.3


# 0.环境设置----
rm(list = ls())
options(stringAsFactors = F)

library(SpaCET)
library(cols4all)


# 1.读取、质控----
visiumPath <- "1.data/SP/PT24_T1/"
SpaCET_obj <- create.SpaCET.object.10X(visiumPath = visiumPath)
# calculate the QC metrics
SpaCET_obj <- SpaCET.quality.control(SpaCET_obj)
# plot the QC metrics
SpaCET.visualize.spatialFeature(
  SpaCET_obj, 
  spatialType = "QualityControl", 
  spatialFeatures=c("UMI","Gene"),
  imageBg = F
) ## 高度4，宽度8


# 2.deconvolve ST data----
SpaCET_obj <- SpaCET.deconvolution(SpaCET_obj, cancerType = "LUAD", coreNo = 48)
# show the ST deconvolution results
SpaCET_obj@results$deconvolution$propMat[1:13,1:6]
# show the spatial distribution of malignant cells.
c4a_series()
c4a_palettes()

p1 <- SpaCET.visualize.spatialFeature(
  SpaCET_obj, 
  spatialType = "CellFraction", 
  spatialFeatures=c("Malignant"),
  imageBg = F) +
  scale_color_continuous_c4a_seq("hcl.blues2") +
  #labs(color = "Fraction") +
  # theme(legend.title = element_text(size = 12),
  #       legend.text = element_text(size = 11),
  #       title = element_text(size = 12))
  theme(legend.position = "none",
        title = element_text(size = 12))
p1
ggsave(filename = "13.SpaCET/P24_T1/Malignant.pdf", p1, height = 4, width = 4, dpi = 300)

# show the spatial distribution of all cell types.
# SpaCET.visualize.spatialFeature(
#   SpaCET_obj, 
#   spatialType = "CellFraction", 
#   spatialFeatures="All", 
#   pointSize = 0.1, 
#   nrow=5
# ) +
#   scale_color_continuous_c4a_seq("wes.zissou1")

# calculate the cell-cell colocalization.
SpaCET_obj <- SpaCET.CCI.colocalization(SpaCET_obj)
SpaCET.visualize.colocalization(SpaCET_obj)

# # calculate the L-R network score across ST spots.
# SpaCET_obj <- SpaCET.CCI.LRNetworkScore(SpaCET_obj, coreNo = 48)
# SpaCET.visualize.spatialFeature(
#   SpaCET_obj, 
#   spatialType = "LRNetworkScore", 
#   spatialFeatures=c("Network_Score","Network_Score_pv")
# )

# 3.Identify the Tumor-Stroma Interface----
SpaCET_obj <- SpaCET.identify.interface(SpaCET_obj, MalignantCutoff = 0.5)
p2 <- SpaCET.visualize.spatialFeature(SpaCET_obj, spatialType = "Interface", 
                                      spatialFeature = "Interface",
                                      imageBg = F) +
  # theme(legend.title = element_text(size = 12),
  #       legend.text = element_text(size = 11),
  #       title = element_text(size = 12))
  theme(legend.position = "none",
        title = element_text(size = 12))
p2
ggsave(filename = "13.SpaCET/P24_T1/Interface.pdf", p2, height = 4, width = 4, dpi = 300)


# 4.malignant cell states----
SpaCET_obj <- SpaCET.deconvolution.malignant(SpaCET_obj, coreNo = 48, malignantCutoff = 0.7)
# show cancer cell state fraction of the first five spots
SpaCET_obj@results$deconvolution$propMat[c("Malignant cell state A","Malignant cell state B"),1:6]
msa <- SpaCET.visualize.spatialFeature(
  SpaCET_obj, 
  spatialType = "CellFraction", 
  spatialFeatures="Malignant cell state A", 
  imageBg = F
)+
  scale_color_continuous_c4a_seq("hcl.blues2") +
  labs(color = "Fraction") +
  theme(legend.title = element_text(size = 12),
        legend.text = element_text(size = 11),
        title = element_text(size = 12))
msa

msb <- SpaCET.visualize.spatialFeature(
  SpaCET_obj, 
  spatialType = "CellFraction", 
  spatialFeatures="Malignant cell state B", 
  imageBg = F
)+
  scale_color_continuous_c4a_seq("hcl.blues2") +
  labs(color = "Fraction") +
  theme(legend.title = element_text(size = 12),
        legend.text = element_text(size = 11),
        title = element_text(size = 12))
msb

ggsave(filename = "13.SpaCET/P24_T1/Malignant_state_A.pdf", msa, height = 4, width = 4.5, dpi = 300)
ggsave(filename = "13.SpaCET/P24_T1/Malignant_state_B.pdf", msb, height = 4, width = 4.5, dpi = 300)


# 5.计算肿瘤细胞状态评分----
# run gene set calculation
SpaCET_obj <- SpaCET.GeneSetScore(SpaCET_obj, GeneSets="CancerCellState")
# show all gene sets
rownames(SpaCET_obj@results$GeneSetScore)
# visualize two gene sets
SpaCET.visualize.spatialFeature(
  SpaCET_obj, 
  spatialType = "GeneSetScore", 
  spatialFeatures = c("CancerCellState_Cycle","CancerCellState_cEMT"),
  imageBg = F
)


# 6.TLS评分----
# run gene set calculation
SpaCET_obj <- SpaCET.GeneSetScore(SpaCET_obj, GeneSets="TLS")
# visualize TLS
SpaCET.visualize.spatialFeature(
  SpaCET_obj, 
  spatialType = "GeneSetScore", 
  spatialFeatures = c("TLS"),
  imageBg = F
)


# 7.保存----
save(SpaCET_obj, file = "13.SpaCET/P24_T1/P24_T1_SpaCET_obj.RData")
