######目的：对H3K18la ChIP-seq和TCGA-LUAD上调的差异基因进行取交集
######作者：申奥
######日期：2024-10-25
######Server: yjx@xiyoucloud.top:34022
######Version: R 4.3.3


# 0.环境设置----
rm(list = ls())
options(stringAsFactors = F)

library(VennDiagram)


# 1.加载两个基因集并取交集----
load("~/projects/kla/PRJNA857271/6.macs2/H1299_klaGenes.RData")

upDEGs <- read.table("15.diff/upDEG_T_vs_N_0.05_2.txt")
upDEgenes <- rownames(upDEGs)

upKla_genes <- intersect(klaGenes, upDEgenes)
save(upKla_genes, file = "15.diff/upKla_genes.RData")
write.table(upKla_genes, file = "15.diff/upKla_genes.txt", quote = F, sep = "\t", row.names = F)


# 2.韦恩图展示----
klaGenes2 <- klaGenes[!is.na(klaGenes)]
#二元#
venn.diagram(x=list(klaGenes2,upDEgenes),
             scaled = T, # 根据比例显示大小
             alpha= 0.9, #透明度
             lwd=1,lty=1,
             col=c("black", "black"), #圆圈线条粗细、形状、颜色；1 实线, 2 虚线, blank无线条
             label.col ='black' , # 数字颜色abel.col=c('#FFFFCC','#CCFFFF',......)根据不同颜色显示数值颜色
             cex = 2, # 数字大小
             fontface = "bold",  # 字体粗细；加粗bold
             fill=c('#336681','#ba4f4a'), # 填充色 配色https://www.58pic.com/
             category.names = c("H3K18la", "UP") , #标签名
             cat.dist = 0.02, # 标签距离圆圈的远近
             cat.pos = -180, # 标签相对于圆圈的角度cat.pos = c(-10, 10, 135)
             cat.cex = 2, #标签字体大小
             cat.fontface = "bold",  # 标签字体加粗
             cat.col='black' ,   #cat.col=c('#FFFFCC','#CCFFFF',.....)根据相应颜色改变标签颜色
             cat.default.pos = "outer",  # 标签位置, outer内;text 外
             output=TRUE,
             filename='15.diff/venn.tiff',# 文件保存
             imagetype="tiff",  # 类型（tiff png svg）
             resolution = 400,  # 分辨率
             compression = "lzw"# 压缩算法
             
)