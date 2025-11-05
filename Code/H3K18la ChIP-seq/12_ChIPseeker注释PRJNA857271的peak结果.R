######目的：ChIPseeker注释PRJNA857271的peak结果
######作者：申奥
######日期：2024-10-24
######Server: yjx@xiyoucloud.top:34022
######Version: R 4.3.3


# 0.环境设置----
rm(list = ls())
options(stringAsFactors = F)

library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)  # 使用 UCSC 提供的 hg38 基因注释数据库
library(org.Hs.eg.db)  # 基因注释数据库
library(clusterProfiler)  # 可视化工具包
library(ggplot2)
library(openxlsx)


# 1.注释----
peak <- readPeakFile("~/projects/kla/PRJNA857271/6.macs2/H1299_peaks.narrowPeak")

# 载入人类基因注释数据
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# 对峰值进行注释
peakAnno <- annotatePeak(peak, 
                         tssRegion=c(-3000, 3000),  # 设置感兴趣的TSS区域（比如TSS上下游3000bp）
                         TxDb=txdb, 
                         annoDb="org.Hs.eg.db")  # 选择注释数据库
peakAnnoDF <- as.data.frame(peakAnno)


# 2.查看注释结果----
# 显示注释结果的摘要
print(peakAnno)

# 查看头几行的注释信息
head(as.data.frame(peakAnno))

# 画图
pdf("H3K18la_PeakAnnotation_pie_chart.pdf")
plotAnnoPie(peakAnno)
dev.off()


# 3.保存----
klaGenes <- unique(peakAnnoDF$SYMBOL)
save(peakAnno, peakAnnoDF, file = "~/projects/kla/PRJNA857271/6.macs2/H1299_peakAnno.RData")
save(klaGenes, file = "~/projects/kla/PRJNA857271/6.macs2/H1299_klaGenes.RData")
write.table(klaGenes, file = "~/projects/kla/PRJNA857271/6.macs2/H1299_klaGenes.txt", sep = "\t", quote = F)
