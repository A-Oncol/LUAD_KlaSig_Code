######目的：利用50中的RSF预后模型在免疫治疗队列中进行验证
######作者：申奥
######日期：2024-11-28
######Server: yjx@xiyoucloud.top:34022
######Version: R 4.4.1


# 0.环境设置----
rm(list = ls())
options(stringAsFactors = F)

library(readr)
library(stringr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(ggsignif)
library(ggdist)
library(openxlsx)


# 1.TCGA-LUAD队列----
## 1.1 TIDE评分与kla.sig模型评分之间的关系比较----
### 1.1.1 加载risk score评分----
load("12.Prognosis_model/RSF_riskScore.RData")
load("12.Prognosis_model/merge_datasets.RData")

identical(merge_datasets[["TCGA_LUAD"]]$OS, rs[["TCGA_LUAD"]]$OS)
identical(merge_datasets[["TCGA_LUAD"]]$OS.time, rs[["TCGA_LUAD"]]$OS.time)
TCGA_LUAD_riskScore <- data.frame(cbind(Sample_ID = merge_datasets[["TCGA_LUAD"]]$Sample_ID, rs[["TCGA_LUAD"]]))
TCGA_LUAD_riskScore$Risk_group <- ifelse(TCGA_LUAD_riskScore$RS > median(TCGA_LUAD_riskScore$RS), "High-Risk", "Low-Risk")
table(TCGA_LUAD_riskScore$Risk_group)

### 1.1.2 加载28脚本1.3中计算好的TIDE评分----
TIDE_res <- read.csv("11.immunotherapy/TCGA-LUAD_TIDE_output.csv", header = T, row.names = 1, check.names = F, stringsAsFactors = F)
table(TIDE_res$Responder)
# False  True 
# 338   203
samples <- rownames(TIDE_res)
table(str_sub(samples, 16, 16))

TIDE_res <- TIDE_res[which(str_sub(samples, 16, 16) == "A"),]
TIDE_res$Sample_ID <- str_sub(rownames(TIDE_res), 1, 15)

### 1.1.3 合并risk score与TIDE，并作图----
riskScore_TIDE <- inner_join(TCGA_LUAD_riskScore, TIDE_res, by = "Sample_ID")
class(riskScore_TIDE$Risk_group)
riskScore_TIDE$Risk_group <- factor(riskScore_TIDE$Risk_group, levels = c("Low-Risk", "High-Risk"))

save(riskScore_TIDE, file = "11.immunotherapy/riskScore_TIDE.RData")

comparison <- list(c("Low-Risk", "High-Risk"))
#Custom.color <- c("#245892","#877eac","#a9d296")
Custom.color <- c("#79aec1", "#ae3f51")
p1 <- ggplot(riskScore_TIDE, aes(x = Risk_group, y = TIDE, fill = Risk_group)) +
  geom_jitter(mapping = aes(color=Risk_group),width = .05, alpha = 0.5,size=0.9)+ #绘制散点图
  geom_boxplot(position = position_nudge(x = 0.14),width=0.1,outlier.size = 0,outlier.alpha =0)+ #绘制箱线图，并通过position设置偏移
  stat_halfeye(mapping = aes(fill=Risk_group),width = 0.2, .width = 0, justification = -1.2, point_colour = NA,alpha=0.6) + #绘制云雨图，并通过position设置偏移
  scale_fill_manual(values = Custom.color)+   #映射云雨图和箱线图的颜色
  scale_color_manual(values = Custom.color)+  #映射散点的颜色
  #expand_limits(x = c(0.5, 3.8))+ #扩展画板，若显示不全，请根据你的数据范围手动调整或删除此行
  # ylim(11,35)+ #控制y轴显示范围，若显示不全，请根据你的数据范围手动调整或删除此行
  xlab("") +  #设置X轴标题
  ylab("TIDE score") +   #设置Y轴标题
  #ggtitle("NSCLC-GSE135222")+  #设置主标题
  theme(axis.ticks.x = element_line(size = 0,color = "black"),  #自定义主题
        panel.background = element_rect(fill = "white", color = "white"),  #设置画板
        panel.grid.major.x = element_blank(),   #设置网格
        panel.grid.minor.x = element_blank(), #设置网格
        panel.grid.major.y = element_line(color = "gray", size = 0.25), #设置网格
        panel.grid.minor.y = element_blank(), #设置网格
        panel.border = element_rect(color = "black", fill = NA,linewidth = 1), #设置边框
        legend.position = "none", #隐藏图例
        axis.title.x = element_text(size = 13),  #调整X轴标题字体大小
        axis.title.y = element_text(size = 13), #调整Y轴标题字体大小
        axis.text.x = element_text(size = 12,hjust = 0.3), #设置x轴刻度字体偏移，若更换数据，可能需要重新设置
        axis.text.y = element_text(size = 12), #设置Y轴刻度字体大小
        plot.title = element_text(hjust = 0.5)
  )+
  stat_compare_means(comparisons = comparison, method = "wilcox.test", label = "p.format", label.y = 4.5)
p1
ggsave(filename = '11.immunotherapy/TCGA-LUAD_klaSigScore_TIDE_boxplot.pdf', p1, height = 3, width = 2.5, dpi = 300)

p2 <- ggplot(riskScore_TIDE, aes(x = Risk_group, y = Exclusion, fill = Risk_group)) +
  geom_jitter(mapping = aes(color=Risk_group),width = .05, alpha = 0.5,size=0.9)+ #绘制散点图
  geom_boxplot(position = position_nudge(x = 0.14),width=0.1,outlier.size = 0,outlier.alpha =0)+ #绘制箱线图，并通过position设置偏移
  stat_halfeye(mapping = aes(fill=Risk_group),width = 0.2, .width = 0, justification = -1.2, point_colour = NA,alpha=0.6) + #绘制云雨图，并通过position设置偏移
  scale_fill_manual(values = Custom.color)+   #映射云雨图和箱线图的颜色
  scale_color_manual(values = Custom.color)+  #映射散点的颜色
  #expand_limits(x = c(0.5, 3.8))+ #扩展画板，若显示不全，请根据你的数据范围手动调整或删除此行
  # ylim(11,35)+ #控制y轴显示范围，若显示不全，请根据你的数据范围手动调整或删除此行
  xlab("") +  #设置X轴标题
  ylab("T cell Exclusion") +   #设置Y轴标题
  #ggtitle("NSCLC-GSE135222")+  #设置主标题
  theme(axis.ticks.x = element_line(size = 0,color = "black"),  #自定义主题
        panel.background = element_rect(fill = "white", color = "white"),  #设置画板
        panel.grid.major.x = element_blank(),   #设置网格
        panel.grid.minor.x = element_blank(), #设置网格
        panel.grid.major.y = element_line(color = "gray", size = 0.25), #设置网格
        panel.grid.minor.y = element_blank(), #设置网格
        panel.border = element_rect(color = "black", fill = NA,linewidth = 1), #设置边框
        legend.position = "none", #隐藏图例
        axis.title.x = element_text(size = 13),  #调整X轴标题字体大小
        axis.title.y = element_text(size = 13), #调整Y轴标题字体大小
        axis.text.x = element_text(size = 12,hjust = 0.3), #设置x轴刻度字体偏移，若更换数据，可能需要重新设置
        axis.text.y = element_text(size = 12), #设置Y轴刻度字体大小
        plot.title = element_text(hjust = 0.5)
  )+
  stat_compare_means(comparisons = comparison, method = "wilcox.test", label = "p.format", label.y = 4)
p2
ggsave(filename = '11.immunotherapy/TCGA-LUAD_klaSigScore_Exclusion_boxplot.pdf', p2, height = 3, width = 3, dpi = 300)

p3 <- ggplot(riskScore_TIDE, aes(x = Risk_group, y = MDSC, fill = Risk_group)) +
  geom_jitter(mapping = aes(color=Risk_group),width = .05, alpha = 0.5,size=0.9)+ #绘制散点图
  geom_boxplot(position = position_nudge(x = 0.14),width=0.1,outlier.size = 0,outlier.alpha =0)+ #绘制箱线图，并通过position设置偏移
  stat_halfeye(mapping = aes(fill=Risk_group),width = 0.2, .width = 0, justification = -1.2, point_colour = NA,alpha=0.6) + #绘制云雨图，并通过position设置偏移
  scale_fill_manual(values = Custom.color)+   #映射云雨图和箱线图的颜色
  scale_color_manual(values = Custom.color)+  #映射散点的颜色
  #expand_limits(x = c(0.5, 3.8))+ #扩展画板，若显示不全，请根据你的数据范围手动调整或删除此行
  # ylim(11,35)+ #控制y轴显示范围，若显示不全，请根据你的数据范围手动调整或删除此行
  xlab("") +  #设置X轴标题
  ylab("MDSC") +   #设置Y轴标题
  #ggtitle("NSCLC-GSE135222")+  #设置主标题
  theme(axis.ticks.x = element_line(size = 0,color = "black"),  #自定义主题
        panel.background = element_rect(fill = "white", color = "white"),  #设置画板
        panel.grid.major.x = element_blank(),   #设置网格
        panel.grid.minor.x = element_blank(), #设置网格
        panel.grid.major.y = element_line(color = "gray", size = 0.25), #设置网格
        panel.grid.minor.y = element_blank(), #设置网格
        panel.border = element_rect(color = "black", fill = NA,linewidth = 1), #设置边框
        legend.position = "none", #隐藏图例
        axis.title.x = element_text(size = 13),  #调整X轴标题字体大小
        axis.title.y = element_text(size = 13), #调整Y轴标题字体大小
        axis.text.x = element_text(size = 12,hjust = 0.3), #设置x轴刻度字体偏移，若更换数据，可能需要重新设置
        axis.text.y = element_text(size = 12), #设置Y轴刻度字体大小
        plot.title = element_text(hjust = 0.5)
  )+
  stat_compare_means(comparisons = comparison, method = "wilcox.test", label = "p.format", label.y = 0.2)
p3
ggsave(filename = '11.immunotherapy/TCGA-LUAD_klaSigScore_MDSC_boxplot.pdf', p3, height = 3, width = 3, dpi = 300)

## 1.2 risk score分组与TCGA-LUAD队列的TIL浸润的关系----
### 1.2.1 加载TIL评分----
####下载自：https://pmc.ncbi.nlm.nih.gov/articles/PMC5943714/#S34中的附表2
til <- read.xlsx("11.immunotherapy/NIHMS958989-supplement-2.xlsx", sheet = 1)
table(til$Study)

luad_til <- til %>% 
  filter(Study == "LUAD")
luad_til <- luad_til[,c("ParticipantBarcode", "til_percentage", "Slide")]

### 1.2.2 与RiskScore合并----
TCGA_LUAD_riskScore$ParticipantBarcode <- str_sub(TCGA_LUAD_riskScore$Sample_ID, 1, 12)

riskScore_TIL <- inner_join(luad_til, TCGA_LUAD_riskScore, by = "ParticipantBarcode")
class(riskScore_TIL$Risk_group)
riskScore_TIL$Risk_group <- factor(riskScore_TIL$Risk_group, levels = c("Low-Risk", "High-Risk"))

save(riskScore_TIL, file = "11.immunotherapy/riskScore_TIL.RData")

p4 <- ggplot(riskScore_TIL, aes(x = Risk_group, y = til_percentage, fill = Risk_group)) +
  geom_jitter(mapping = aes(color=Risk_group),width = .05, alpha = 0.5,size=0.9)+ #绘制散点图
  geom_boxplot(position = position_nudge(x = 0.14),width=0.1,outlier.size = 0,outlier.alpha =0)+ #绘制箱线图，并通过position设置偏移
  stat_halfeye(mapping = aes(fill=Risk_group),width = 0.2, .width = 0, justification = -1.2, point_colour = NA,alpha=0.6) + #绘制云雨图，并通过position设置偏移
  scale_fill_manual(values = Custom.color)+   #映射云雨图和箱线图的颜色
  scale_color_manual(values = Custom.color)+  #映射散点的颜色
  #expand_limits(x = c(0.5, 3.8))+ #扩展画板，若显示不全，请根据你的数据范围手动调整或删除此行
  # ylim(11,35)+ #控制y轴显示范围，若显示不全，请根据你的数据范围手动调整或删除此行
  xlab("") +  #设置X轴标题
  ylab("TIL percentage") +   #设置Y轴标题
  #ggtitle("NSCLC-GSE135222")+  #设置主标题
  theme(axis.ticks.x = element_line(size = 0,color = "black"),  #自定义主题
        panel.background = element_rect(fill = "white", color = "white"),  #设置画板
        panel.grid.major.x = element_blank(),   #设置网格
        panel.grid.minor.x = element_blank(), #设置网格
        panel.grid.major.y = element_line(color = "gray", size = 0.25), #设置网格
        panel.grid.minor.y = element_blank(), #设置网格
        panel.border = element_rect(color = "black", fill = NA,linewidth = 1), #设置边框
        legend.position = "none", #隐藏图例
        axis.title.x = element_text(size = 13),  #调整X轴标题字体大小
        axis.title.y = element_text(size = 13), #调整Y轴标题字体大小
        axis.text.x = element_text(size = 12,hjust = 0.3), #设置x轴刻度字体偏移，若更换数据，可能需要重新设置
        axis.text.y = element_text(size = 12), #设置Y轴刻度字体大小
        plot.title = element_text(hjust = 0.5)
  )+
  stat_compare_means(comparisons = comparison, method = "wilcox.test", label = "p.format", label.y = 32)
p4
ggsave(filename = '11.immunotherapy/TCGA-LUAD_klaSigScore_TIL_boxplot.pdf', p4, height = 3, width = 2.5, dpi = 300)


# 2.真实免疫治疗队列----
## 2.1 加载RSFfit模型----
load("12.Prognosis_model/RSFfit.RData")
load("12.Prognosis_model/same_genes.RData")
load("15.diff/upKla_genes.RData")
build_model_genes <- intersect(upKla_genes, same_genes) # 318个

## 2.2 GSE126044----
### 2.2.1 读取TPM和注释文件，整理成表达谱后进行scale和risk score计算----
####可以通过https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE126044下载
gse126044_tpm <- read_tsv("11.immunotherapy/GEO/GSE126044_norm_counts_TPM_GRCh38.p13_NCBI.tsv")
annot <- read_tsv('11.immunotherapy/GEO/Human.GRCh38.p13.annot.tsv') %>% 
  select(1,2,5)
table(annot$GeneType)
annot <- annot %>% 
  filter(GeneType %in% c("protein-coding", 'ncRNA')) %>% 
  select(1,2)

gse126044_tpm_annot <- inner_join(annot, gse126044_tpm, by = 'GeneID') %>% 
  select(-1)
which(duplicated(gse126044_tpm_annot$Symbol))
gse126044_tpm_annot <- gse126044_tpm_annot %>% 
  tibble::column_to_rownames("Symbol")

log2_gse126044_tpm <- log2(gse126044_tpm_annot + 1)
log2_gse126044_tpm <- as.data.frame(t(log2_gse126044_tpm))

expr_df <- log2_gse126044_tpm[,build_model_genes]
expr_df <- as.data.frame(scale(expr_df))
riskScore <- predict(fit, newdata = expr_df, type = "risk")$predicted
rs <- data.frame(Sample_ID = rownames(expr_df),
                 RiskScore = riskScore)
rs$RiskGroup <- ifelse(rs$RiskScore > median(rs$RiskScore), "High-Risk", "Low-Risk")
table(rs$RiskGroup)
rs$RiskGroup <- factor(rs$RiskGroup, levels = c("Low-Risk", "High-Risk"))

### 2.2.2 临床信息----
gse126044_clinic <- read.xlsx("11.immunotherapy/GEO/GSE126044_clinicalinfo.xlsx", sheet = 1)
colnames(gse126044_clinic)[4] <- "Response"
table(gse126044_clinic$Response)
gse126044_clinic$Response <- factor(gse126044_clinic$Response, levels = c("NR", "R"))

clinic <- gse126044_clinic

### 2.2.3 合并与画图比较----
risk_clinic <- inner_join(rs, clinic, by = "Sample_ID")

comparison <- list(c("NR", "R"))
#Custom.color <- c("#245892","#877eac","#a9d296")
Custom.color <- c("#79aec1", "#ae3f51")
p5 <- ggplot(risk_clinic, aes(x = Response, y = RiskScore, fill = Response)) +
  geom_jitter(mapping = aes(color=Response),width = .05, alpha = 0.5,size=0.9)+ #绘制散点图
  geom_boxplot(position = position_nudge(x = 0.14),width=0.1,outlier.size = 0,outlier.alpha =0)+ #绘制箱线图，并通过position设置偏移
  stat_halfeye(mapping = aes(fill=Response),width = 0.2, .width = 0, justification = -1.2, point_colour = NA,alpha=0.6) + #绘制云雨图，并通过position设置偏移
  scale_fill_manual(values = Custom.color)+   #映射云雨图和箱线图的颜色
  scale_color_manual(values = Custom.color)+  #映射散点的颜色
  #expand_limits(x = c(0.5, 3.8))+ #扩展画板，若显示不全，请根据你的数据范围手动调整或删除此行
  #ylim(20,110)+ #控制y轴显示范围，若显示不全，请根据你的数据范围手动调整或删除此行
  xlab("") +  #设置X轴标题
  #ylab("Risk Score") +   #设置Y轴标题
  ggtitle("NSCLC-GSE126044")+  #设置主标题
  scale_y_continuous(
    name = "Risk Score",                     # 设置 Y 轴标题
    limits = c(30, 110),                    # 设置 Y 轴显示范围
    breaks = seq(30, 110, by = 20)          # 按步长 10 显示刻度
  ) +
  theme(axis.ticks.x = element_line(size = 0,color = "black"),  #自定义主题
        panel.background = element_rect(fill = "white", color = "white"),  #设置画板
        panel.grid.major.x = element_blank(),   #设置网格
        panel.grid.minor.x = element_blank(), #设置网格
        panel.grid.major.y = element_line(color = "gray", size = 0.25), #设置网格
        panel.grid.minor.y = element_blank(), #设置网格
        panel.border = element_rect(color = "black", fill = NA,linewidth = 1), #设置边框
        legend.position = "none", #隐藏图例
        axis.title.x = element_text(size = 13),  #调整X轴标题字体大小
        axis.title.y = element_text(size = 13), #调整Y轴标题字体大小
        axis.text.x = element_text(size = 12,hjust = 0.3), #设置x轴刻度字体偏移，若更换数据，可能需要重新设置
        axis.text.y = element_text(size = 12), #设置Y轴刻度字体大小
        plot.title = element_text(hjust = 0.5)
  )+
  stat_compare_means(comparisons = comparison, method = "wilcox.test", label = "p.format", label.y = 100)
p5
ggsave(filename = '11.immunotherapy/GSE126044_klaSigScore_response_boxplot.pdf', p5, height = 3, width = 2.5, dpi = 300)

## 2.3 GSE135222----
### 2.3.1 读取TPM和注释文件，整理成表达谱后进行scale和risk score计算----
gse135222_tpm <- read_tsv("11.immunotherapy/GEO/GSE135222_norm_counts_TPM_GRCh38.p13_NCBI.tsv")
annot <- read_tsv('11.immunotherapy/GEO/Human.GRCh38.p13.annot.tsv') %>% 
  select(1,2,5)
table(annot$GeneType)
annot <- annot %>% 
  filter(GeneType %in% c("protein-coding", 'ncRNA')) %>% 
  select(1,2)

gse135222_tpm_annot <- inner_join(annot, gse135222_tpm, by = 'GeneID') %>% 
  select(-1)
which(duplicated(gse135222_tpm_annot$Symbol))
gse135222_tpm_annot <- gse135222_tpm_annot %>% 
  tibble::column_to_rownames("Symbol")

log2_gse135222_tpm <- log2(gse135222_tpm_annot + 1)
log2_gse135222_tpm <- as.data.frame(t(log2_gse135222_tpm))

expr_df <- log2_gse135222_tpm[,build_model_genes]
expr_df <- as.data.frame(scale(expr_df))
riskScore <- predict(fit, newdata = expr_df, type = "risk")$predicted
rs <- data.frame(Sample_ID = rownames(expr_df),
                 RiskScore = riskScore)
rs$RiskGroup <- ifelse(rs$RiskScore > median(rs$RiskScore), "High-Risk", "Low-Risk")
table(rs$RiskGroup)
rs$RiskGroup <- factor(rs$RiskGroup, levels = c("Low-Risk", "High-Risk"))

### 2.3.2 临床信息----
gse135222_clinic <- read.xlsx("11.immunotherapy/GEO/GSE135222_clinical_info.xlsx", sheet = 1)
gse135222_samples <- read.xlsx('11.immunotherapy/GEO/GSE135222_samples.xlsx', sheet = 1)
gse135222_samples$Patient.ID <- str_split_i(gse135222_samples$Sample_Name, ' ', 2)
gse135222_samples$Patient.ID <- as.numeric(gse135222_samples$Patient.ID)

gse135222_response <- inner_join(gse135222_samples, gse135222_clinic, by = 'Patient.ID')
gse135222_response <- gse135222_response %>% 
  select(c('Sample_ID', 'benefit'))
table(gse135222_response$benefit)
colnames(gse135222_response) <- c('Sample_ID', 'Response')
table(gse135222_response$Response)
gse135222_response$Response <- ifelse(gse135222_response$Response == "Y", "R", "NR")
gse135222_response$Response <- factor(gse135222_response$Response, levels = c("NR", "R"))

clinic <- gse135222_response

### 2.3.3 合并与画图比较----
risk_clinic <- inner_join(rs, clinic, by = "Sample_ID")

comparison <- list(c("NR", "R"))
#Custom.color <- c("#245892","#877eac","#a9d296")
Custom.color <- c("#79aec1", "#ae3f51")
p6 <- ggplot(risk_clinic, aes(x = Response, y = RiskScore, fill = Response)) +
  geom_jitter(mapping = aes(color=Response),width = .05, alpha = 0.5,size=0.9)+ #绘制散点图
  geom_boxplot(position = position_nudge(x = 0.14),width=0.1,outlier.size = 0,outlier.alpha =0)+ #绘制箱线图，并通过position设置偏移
  stat_halfeye(mapping = aes(fill=Response),width = 0.2, .width = 0, justification = -1.2, point_colour = NA,alpha=0.6) + #绘制云雨图，并通过position设置偏移
  scale_fill_manual(values = Custom.color)+   #映射云雨图和箱线图的颜色
  scale_color_manual(values = Custom.color)+  #映射散点的颜色
  #expand_limits(x = c(0.5, 3.8))+ #扩展画板，若显示不全，请根据你的数据范围手动调整或删除此行
  #ylim(20,110)+ #控制y轴显示范围，若显示不全，请根据你的数据范围手动调整或删除此行
  xlab("") +  #设置X轴标题
  #ylab("Risk Score") +   #设置Y轴标题
  ggtitle("NSCLC-GSE135222")+  #设置主标题
  scale_y_continuous(
    name = "Risk Score",                     # 设置 Y 轴标题
    limits = c(30, 120),                    # 设置 Y 轴显示范围
    breaks = seq(30, 120, by = 20)          # 按步长 10 显示刻度
  ) +
  theme(axis.ticks.x = element_line(size = 0,color = "black"),  #自定义主题
        panel.background = element_rect(fill = "white", color = "white"),  #设置画板
        panel.grid.major.x = element_blank(),   #设置网格
        panel.grid.minor.x = element_blank(), #设置网格
        panel.grid.major.y = element_line(color = "gray", size = 0.25), #设置网格
        panel.grid.minor.y = element_blank(), #设置网格
        panel.border = element_rect(color = "black", fill = NA,linewidth = 1), #设置边框
        legend.position = "none", #隐藏图例
        axis.title.x = element_text(size = 13),  #调整X轴标题字体大小
        axis.title.y = element_text(size = 13), #调整Y轴标题字体大小
        axis.text.x = element_text(size = 12,hjust = 0.3), #设置x轴刻度字体偏移，若更换数据，可能需要重新设置
        axis.text.y = element_text(size = 12), #设置Y轴刻度字体大小
        plot.title = element_text(hjust = 0.5)
  )+
  stat_compare_means(comparisons = comparison, method = "wilcox.test", label = "p.format", label.y = 110)
p6
ggsave(filename = '11.immunotherapy/GSE135222_klaSigScore_response_boxplot.pdf', p6, height = 3, width = 2.5, dpi = 300)

## 2.4 GSE91061----
### 2.4.1 读取TPM和注释文件，整理成表达谱后进行scale和risk score计算----
gse91061_tpm <- read_tsv("11.immunotherapy/GEO/GSE91061_norm_counts_TPM_GRCh38.p13_NCBI.tsv")
annot <- read_tsv('11.immunotherapy/GEO/Human.GRCh38.p13.annot.tsv') %>% 
  select(1,2,5)
table(annot$GeneType)
annot <- annot %>% 
  filter(GeneType %in% c("protein-coding", 'ncRNA')) %>% 
  select(1,2)

gse91061_tpm_annot <- inner_join(annot, gse91061_tpm, by = 'GeneID') %>% 
  select(-1)
which(duplicated(gse91061_tpm_annot$Symbol))
gse91061_tpm_annot <- gse91061_tpm_annot %>% 
  tibble::column_to_rownames("Symbol")

log2_gse91061_tpm <- log2(gse91061_tpm_annot + 1)
log2_gse91061_tpm <- as.data.frame(t(log2_gse91061_tpm))

expr_df <- log2_gse91061_tpm[,build_model_genes]
expr_df <- as.data.frame(scale(expr_df))
riskScore <- predict(fit, newdata = expr_df, type = "risk")$predicted
rs <- data.frame(Sample_ID = rownames(expr_df),
                 RiskScore = riskScore)
rs$RiskGroup <- ifelse(rs$RiskScore > median(rs$RiskScore), "High-Risk", "Low-Risk")
table(rs$RiskGroup)
rs$RiskGroup <- factor(rs$RiskGroup, levels = c("Low-Risk", "High-Risk"))

### 2.4.2 临床信息----
gse91061_clinic <- read.xlsx("11.immunotherapy/GEO/GSE91061_clinical_info.xlsx", sheet = 1)
table(gse91061_clinic$Response)
table(gse91061_clinic$Treatment)
gse91061_clinic <- gse91061_clinic %>% 
  filter(Response != "UNK")
table(gse91061_clinic$Response)
gse91061_clinic$Response <- ifelse(gse91061_clinic$Response == "PRCR", "R", "NR")
table(gse91061_clinic$Response)
# NR  R 
# 82 23
gse91061_clinic$Response <- factor(gse91061_clinic$Response, levels = c("NR", "R"))

clinic <- gse91061_clinic

### 2.4.3 合并与画图比较----
risk_clinic <- inner_join(rs, clinic, by = "Sample_ID")
table(risk_clinic$RiskGroup, risk_clinic$Response)

comparison <- list(c("NR", "R"))
#Custom.color <- c("#245892","#877eac","#a9d296")
Custom.color <- c("#79aec1", "#ae3f51")
p7 <- ggplot(risk_clinic, aes(x = Response, y = RiskScore, fill = Response)) +
  geom_jitter(mapping = aes(color=Response),width = .05, alpha = 0.5,size=0.9)+ #绘制散点图
  geom_boxplot(position = position_nudge(x = 0.14),width=0.1,outlier.size = 0,outlier.alpha =0)+ #绘制箱线图，并通过position设置偏移
  stat_halfeye(mapping = aes(fill=Response),width = 0.2, .width = 0, justification = -1.2, point_colour = NA,alpha=0.6) + #绘制云雨图，并通过position设置偏移
  scale_fill_manual(values = Custom.color)+   #映射云雨图和箱线图的颜色
  scale_color_manual(values = Custom.color)+  #映射散点的颜色
  #expand_limits(x = c(0.5, 3.8))+ #扩展画板，若显示不全，请根据你的数据范围手动调整或删除此行
  #ylim(20,110)+ #控制y轴显示范围，若显示不全，请根据你的数据范围手动调整或删除此行
  xlab("") +  #设置X轴标题
  #ylab("Risk Score") +   #设置Y轴标题
  ggtitle("Melanoma-GSE91061")+  #设置主标题
  scale_y_continuous(
    name = "Risk Score",                     # 设置 Y 轴标题
    limits = c(30, 110),                    # 设置 Y 轴显示范围
    breaks = seq(30, 110, by = 20)          # 按步长 10 显示刻度
  ) +
  theme(axis.ticks.x = element_line(size = 0,color = "black"),  #自定义主题
        panel.background = element_rect(fill = "white", color = "white"),  #设置画板
        panel.grid.major.x = element_blank(),   #设置网格
        panel.grid.minor.x = element_blank(), #设置网格
        panel.grid.major.y = element_line(color = "gray", size = 0.25), #设置网格
        panel.grid.minor.y = element_blank(), #设置网格
        panel.border = element_rect(color = "black", fill = NA,linewidth = 1), #设置边框
        legend.position = "none", #隐藏图例
        axis.title.x = element_text(size = 13),  #调整X轴标题字体大小
        axis.title.y = element_text(size = 13), #调整Y轴标题字体大小
        axis.text.x = element_text(size = 12,hjust = 0.3), #设置x轴刻度字体偏移，若更换数据，可能需要重新设置
        axis.text.y = element_text(size = 12), #设置Y轴刻度字体大小
        plot.title = element_text(hjust = 0.5)
  )+
  stat_compare_means(comparisons = comparison, method = "wilcox.test", label = "p.format", label.y = 100)
p7
ggsave(filename = '11.immunotherapy/GSE91061_klaSigScore_response_boxplot.pdf', p7, height = 3, width = 2.5, dpi = 300)

## 2.5 GSE106128----
### 2.5.1 整理表达谱并计算评分----
gse106128_exp <- read.table("11.immunotherapy/GEO/GSE106128_series_matrix.txt", header = T, sep = "\t", comment.char = "!")
gpl6947 <- data.table::fread("11.immunotherapy/GEO/GPL6947-13512.txt", skip = 30)
gpl6947 <- gpl6947[,c(1,14)]
colnames(gpl6947)[1] <- "ID_REF"

gse106128_exp_annot <- inner_join(gpl6947, gse106128_exp, by = "ID_REF")
gse106128_exp_annot <- gse106128_exp_annot[,-1]
gse106128_exp_annot <- gse106128_exp_annot %>% 
  filter(Symbol != "")
which(duplicated(gse106128_exp_annot$Symbol))
gse106128_exp_annot <- aggregate(. ~ Symbol, gse106128_exp_annot, mean)
which(duplicated(gse106128_exp_annot$Symbol))

gse106128_exp_annot <- gse106128_exp_annot %>% 
  tibble::column_to_rownames("Symbol")

log2_gse106128_exp <- log2(gse106128_exp_annot + 1)
log2_gse106128_exp <- as.data.frame(t(log2_gse106128_exp))

expr_df <- log2_gse106128_exp[,build_model_genes]
expr_df <- as.data.frame(scale(expr_df))
riskScore <- predict(fit, newdata = expr_df, type = "risk")$predicted
rs <- data.frame(Sample_ID = rownames(expr_df),
                 RiskScore = riskScore)
rs$RiskGroup <- ifelse(rs$RiskScore > median(rs$RiskScore), "High-Risk", "Low-Risk")
table(rs$RiskGroup)
rs$RiskGroup <- factor(rs$RiskGroup, levels = c("Low-Risk", "High-Risk"))

### 2.5.2 临床信息----
gse106128_clinic <- read_tsv("11.immunotherapy/GEO/Melanoma-GSE106128.Response (1).tsv")
gse106128_clinic <- gse106128_clinic[,-1]
gse106128_clinic <- gse106128_clinic[,c(1,5,6,9,10)]
colnames(gse106128_clinic) <- c("Sample_ID", "Treatment", "Response", "OS.time", "OS")

table(gse106128_clinic$Treatment)

gse106128_clinic$OS.time <- gse106128_clinic$OS.time / 365
table(gse106128_clinic$OS)

gse106128_clinic$OS <- ifelse(gse106128_clinic$OS == "Dead", 1, 0)
table(gse106128_clinic$OS)

table(gse106128_clinic$Response)
gse106128_clinic <- gse106128_clinic %>% 
  filter(Response != "UNK")
table(gse106128_clinic$Response)
gse106128_clinic$Response <- factor(gse106128_clinic$Response, levels = c("NR", "R"))

clinic <- gse106128_clinic

### 2.5.3 合并与画图比较----
risk_clinic <- inner_join(rs, clinic, by = "Sample_ID")
table(risk_clinic$RiskGroup, risk_clinic$Response)

comparison <- list(c("NR", "R"))
#Custom.color <- c("#245892","#877eac","#a9d296")
Custom.color <- c("#79aec1", "#ae3f51")
p8 <- ggplot(risk_clinic, aes(x = Response, y = RiskScore, fill = Response)) +
  geom_jitter(mapping = aes(color=Response),width = .05, alpha = 0.5,size=0.9)+ #绘制散点图
  geom_boxplot(position = position_nudge(x = 0.14),width=0.1,outlier.size = 0,outlier.alpha =0)+ #绘制箱线图，并通过position设置偏移
  stat_halfeye(mapping = aes(fill=Response),width = 0.2, .width = 0, justification = -1.2, point_colour = NA,alpha=0.6) + #绘制云雨图，并通过position设置偏移
  scale_fill_manual(values = Custom.color)+   #映射云雨图和箱线图的颜色
  scale_color_manual(values = Custom.color)+  #映射散点的颜色
  #expand_limits(x = c(0.5, 3.8))+ #扩展画板，若显示不全，请根据你的数据范围手动调整或删除此行
  #ylim(20,110)+ #控制y轴显示范围，若显示不全，请根据你的数据范围手动调整或删除此行
  xlab("") +  #设置X轴标题
  #ylab("Risk Score") +   #设置Y轴标题
  ggtitle("Melanoma-GSE106128")+  #设置主标题
  scale_y_continuous(
    name = "Risk Score",                     # 设置 Y 轴标题
    limits = c(30, 110),                    # 设置 Y 轴显示范围
    breaks = seq(30, 110, by = 20)          # 按步长 10 显示刻度
  ) +
  theme(axis.ticks.x = element_line(size = 0,color = "black"),  #自定义主题
        panel.background = element_rect(fill = "white", color = "white"),  #设置画板
        panel.grid.major.x = element_blank(),   #设置网格
        panel.grid.minor.x = element_blank(), #设置网格
        panel.grid.major.y = element_line(color = "gray", size = 0.25), #设置网格
        panel.grid.minor.y = element_blank(), #设置网格
        panel.border = element_rect(color = "black", fill = NA,linewidth = 1), #设置边框
        legend.position = "none", #隐藏图例
        axis.title.x = element_text(size = 13),  #调整X轴标题字体大小
        axis.title.y = element_text(size = 13), #调整Y轴标题字体大小
        axis.text.x = element_text(size = 12,hjust = 0.3), #设置x轴刻度字体偏移，若更换数据，可能需要重新设置
        axis.text.y = element_text(size = 12), #设置Y轴刻度字体大小
        plot.title = element_text(hjust = 0.5)
  )+
  stat_compare_means(comparisons = comparison, method = "wilcox.test", label = "p.format", label.y = 100)
p8
ggsave(filename = '11.immunotherapy/GSE106128_klaSigScore_response_boxplot.pdf', p8, height = 3, width = 2.5, dpi = 300)

### 2.5.4 生存分析----
library(survival)
library(survminer)

surv_object <- Surv(time = risk_clinic$OS.time, event = risk_clinic$OS)
kmfit <- survfit(surv_object ~ RiskGroup, data = risk_clinic)

# 绘制生存曲线
high_color <- pal_npg("nrc")(10)[1]
low_color <- pal_npg("nrc")(10)[2]

km_plot <- ggsurvplot(kmfit, data = risk_clinic,
                      pval = T, # 在图上添加log rank检验的p值
                      pval.method = T,
                      #pval.coord = c(13.5, 0.9), # 调整P值位置，X=10，Y=0.1
                      #pval.method.coord = c(9, 0.9), # 调整P值方法的位置，X=3，Y=0.85
                      conf.int = T,# 添加置信区间
                      risk.table = T, # 在图下方添加风险表
                      legend.labs=c("High","Low"), #表头标签注释分组
                      legend.title="Risk-group",#表头标签
                      #title="Overall survival",#改一下整体名称
                      ylab="Overall Survival (OS)",xlab = " Time (Years)",#修改X轴Y轴名称
                      risk.table.col = "strata", # 风险表加颜色
                      linetype = 1, # 生存曲线的线型
                      surv.median.line = "hv", # 标注出中位生存时间
                      ggtheme = theme_bw(base_size = 12, base_line_size = 1, base_rect_size = 1) + theme(panel.grid = element_blank()), #背景布局
                      palette = c(high_color, low_color))# 图形颜色风格
km_plot

## 2.6 GSE115821----
### 2.6.1 读取TPM和注释文件，整理成表达谱后进行scale和risk score计算----
gse115821_tpm <- read_tsv("11.immunotherapy/GEO/GSE115821_norm_counts_TPM_GRCh38.p13_NCBI.tsv")
annot <- read_tsv('11.immunotherapy/GEO/Human.GRCh38.p13.annot.tsv') %>% 
  select(1,2,5)
table(annot$GeneType)
annot <- annot %>% 
  filter(GeneType %in% c("protein-coding", 'ncRNA')) %>% 
  select(1,2)

gse115821_tpm_annot <- inner_join(annot, gse115821_tpm, by = 'GeneID') %>% 
  select(-1)
which(duplicated(gse115821_tpm_annot$Symbol))
gse115821_tpm_annot <- gse115821_tpm_annot %>% 
  tibble::column_to_rownames("Symbol")

log2_gse115821_tpm <- log2(gse115821_tpm_annot + 1)
log2_gse115821_tpm <- as.data.frame(t(log2_gse115821_tpm))

expr_df <- log2_gse115821_tpm[,build_model_genes]
expr_df <- as.data.frame(scale(expr_df))
riskScore <- predict(fit, newdata = expr_df, type = "risk")$predicted
rs <- data.frame(Sample_ID = rownames(expr_df),
                 RiskScore = riskScore)
rs$RiskGroup <- ifelse(rs$RiskScore > median(rs$RiskScore), "High-Risk", "Low-Risk")
table(rs$RiskGroup)
rs$RiskGroup <- factor(rs$RiskGroup, levels = c("Low-Risk", "High-Risk"))

### 2.6.2 临床信息----
gse115821_clinic <- read.xlsx("11.immunotherapy/GEO/GSE115821_clinical.xlsx", sheet = 1)
gse115821_clinic <- gse115821_clinic[,c(1,2,5)]
colnames(gse115821_clinic) <- c("Sample_ID", "Treatment", "Response")
table(gse115821_clinic$Response)
gse115821_clinic$Response <- factor(gse115821_clinic$Response, levels = c("NR", "R"))

clinic <- gse115821_clinic

### 2.6.3 合并与画图比较----
risk_clinic <- inner_join(rs, clinic, by = "Sample_ID")

comparison <- list(c("NR", "R"))
#Custom.color <- c("#245892","#877eac","#a9d296")
Custom.color <- c("#79aec1", "#ae3f51")
p9 <- ggplot(risk_clinic, aes(x = Response, y = RiskScore, fill = Response)) +
  geom_jitter(mapping = aes(color=Response),width = .05, alpha = 0.5,size=0.9)+ #绘制散点图
  geom_boxplot(position = position_nudge(x = 0.14),width=0.1,outlier.size = 0,outlier.alpha =0)+ #绘制箱线图，并通过position设置偏移
  stat_halfeye(mapping = aes(fill=Response),width = 0.2, .width = 0, justification = -1.2, point_colour = NA,alpha=0.6) + #绘制云雨图，并通过position设置偏移
  scale_fill_manual(values = Custom.color)+   #映射云雨图和箱线图的颜色
  scale_color_manual(values = Custom.color)+  #映射散点的颜色
  #expand_limits(x = c(0.5, 3.8))+ #扩展画板，若显示不全，请根据你的数据范围手动调整或删除此行
  #ylim(20,110)+ #控制y轴显示范围，若显示不全，请根据你的数据范围手动调整或删除此行
  xlab("") +  #设置X轴标题
  #ylab("Risk Score") +   #设置Y轴标题
  ggtitle("Melanoma-GSE115821")+  #设置主标题
  scale_y_continuous(
    name = "Risk Score",                     # 设置 Y 轴标题
    limits = c(30, 110),                    # 设置 Y 轴显示范围
    breaks = seq(30, 110, by = 20)          # 按步长 10 显示刻度
  ) +
  theme(axis.ticks.x = element_line(size = 0,color = "black"),  #自定义主题
        panel.background = element_rect(fill = "white", color = "white"),  #设置画板
        panel.grid.major.x = element_blank(),   #设置网格
        panel.grid.minor.x = element_blank(), #设置网格
        panel.grid.major.y = element_line(color = "gray", size = 0.25), #设置网格
        panel.grid.minor.y = element_blank(), #设置网格
        panel.border = element_rect(color = "black", fill = NA,linewidth = 1), #设置边框
        legend.position = "none", #隐藏图例
        axis.title.x = element_text(size = 13),  #调整X轴标题字体大小
        axis.title.y = element_text(size = 13), #调整Y轴标题字体大小
        axis.text.x = element_text(size = 12,hjust = 0.3), #设置x轴刻度字体偏移，若更换数据，可能需要重新设置
        axis.text.y = element_text(size = 12), #设置Y轴刻度字体大小
        plot.title = element_text(hjust = 0.5)
  )+
  stat_compare_means(comparisons = comparison, method = "wilcox.test", label = "p.format", label.y = 100)
p9
#ggsave(filename = '11.immunotherapy/GSE115821_klaSigScore_response_boxplot.pdf', p9, height = 3, width = 2.5, dpi = 300)

## 2.7 phs000452----
### 2.7.1 整理表达谱并计算评分----
phs000452_expr <- readRDS("11.immunotherapy/Immunotherapy_datasets/Expression/Melanoma-phs000452.Response.Rds")
which(duplicated(phs000452_expr$GENE_SYMBOL))
phs000452_expr <- phs000452_expr %>% 
  tibble::column_to_rownames("GENE_SYMBOL")

log2_phs000452_expr <- log2(phs000452_expr + 1)
log2_phs000452_expr <- as.data.frame(t(log2_phs000452_expr))

expr_df <- log2_phs000452_expr[,build_model_genes]
expr_df <- as.data.frame(scale(expr_df))
riskScore <- predict(fit, newdata = expr_df, type = "risk")$predicted
rs <- data.frame(Sample_ID = rownames(expr_df),
                 RiskScore = riskScore)
rs$RiskGroup <- ifelse(rs$RiskScore > median(rs$RiskScore), "High-Risk", "Low-Risk")
table(rs$RiskGroup)
rs$RiskGroup <- factor(rs$RiskGroup, levels = c("Low-Risk", "High-Risk"))

### 2.7.2 临床信息----
phs000452_clinic <- read_tsv("11.immunotherapy/Immunotherapy_datasets/Clinical/Melanoma-phs000452.Response.tsv")
colnames(phs000452_clinic) <- c("Num", colnames(phs000452_clinic))
phs000452_clinic <- phs000452_clinic[,c(2,7,10,11)]
colnames(phs000452_clinic) <- c("Sample_ID", "Response", "OS.time", "OS")

table(phs000452_clinic$Response)
phs000452_clinic$Response <- ifelse(phs000452_clinic$Response %in% c("PD", "SD"), "NR", "R")
table(phs000452_clinic$Response)
class(phs000452_clinic$Response)
phs000452_clinic$Response <- factor(phs000452_clinic$Response, levels = c("NR", "R"))

phs000452_clinic$OS.time <- phs000452_clinic$OS.time / 365

table(phs000452_clinic$OS)
phs000452_clinic$OS <- ifelse(phs000452_clinic$OS == "Dead", 1, 0)

clinic <- phs000452_clinic

### 2.7.3 合并与画图比较----
risk_clinic <- inner_join(rs, clinic, by = "Sample_ID")

comparison <- list(c("NR", "R"))
#Custom.color <- c("#245892","#877eac","#a9d296")
Custom.color <- c("#79aec1", "#ae3f51")
p10 <- ggplot(risk_clinic, aes(x = Response, y = RiskScore, fill = Response)) +
  geom_jitter(mapping = aes(color=Response),width = .05, alpha = 0.5,size=0.9)+ #绘制散点图
  geom_boxplot(position = position_nudge(x = 0.14),width=0.1,outlier.size = 0,outlier.alpha =0)+ #绘制箱线图，并通过position设置偏移
  stat_halfeye(mapping = aes(fill=Response),width = 0.2, .width = 0, justification = -1.2, point_colour = NA,alpha=0.6) + #绘制云雨图，并通过position设置偏移
  scale_fill_manual(values = Custom.color)+   #映射云雨图和箱线图的颜色
  scale_color_manual(values = Custom.color)+  #映射散点的颜色
  #expand_limits(x = c(0.5, 3.8))+ #扩展画板，若显示不全，请根据你的数据范围手动调整或删除此行
  #ylim(20,110)+ #控制y轴显示范围，若显示不全，请根据你的数据范围手动调整或删除此行
  xlab("") +  #设置X轴标题
  #ylab("Risk Score") +   #设置Y轴标题
  ggtitle("Melanoma-phs000452")+  #设置主标题
  scale_y_continuous(
    name = "Risk Score",                     # 设置 Y 轴标题
    limits = c(30, 110),                    # 设置 Y 轴显示范围
    breaks = seq(30, 110, by = 20)          # 按步长 10 显示刻度
  ) +
  theme(axis.ticks.x = element_line(size = 0,color = "black"),  #自定义主题
        panel.background = element_rect(fill = "white", color = "white"),  #设置画板
        panel.grid.major.x = element_blank(),   #设置网格
        panel.grid.minor.x = element_blank(), #设置网格
        panel.grid.major.y = element_line(color = "gray", size = 0.25), #设置网格
        panel.grid.minor.y = element_blank(), #设置网格
        panel.border = element_rect(color = "black", fill = NA,linewidth = 1), #设置边框
        legend.position = "none", #隐藏图例
        axis.title.x = element_text(size = 13),  #调整X轴标题字体大小
        axis.title.y = element_text(size = 13), #调整Y轴标题字体大小
        axis.text.x = element_text(size = 12,hjust = 0.3), #设置x轴刻度字体偏移，若更换数据，可能需要重新设置
        axis.text.y = element_text(size = 12), #设置Y轴刻度字体大小
        plot.title = element_text(hjust = 0.5)
  )+
  stat_compare_means(comparisons = comparison, method = "wilcox.test", label = "p.format", label.y = 100)
p10
#ggsave(filename = '11.immunotherapy/phs000452_klaSigScore_response_boxplot.pdf', p10, height = 3, width = 2.5, dpi = 300)

### 2.7.4 生存分析----
surv_object <- Surv(time = risk_clinic$OS.time, event = risk_clinic$OS)
kmfit <- survfit(surv_object ~ RiskGroup, data = risk_clinic)

# 绘制生存曲线
high_color <- pal_npg("nrc")(10)[1]
low_color <- pal_npg("nrc")(10)[2]

km_plot <- ggsurvplot(kmfit, data = risk_clinic,
                      pval = T, # 在图上添加log rank检验的p值
                      pval.method = T,
                      #pval.coord = c(13.5, 0.9), # 调整P值位置，X=10，Y=0.1
                      #pval.method.coord = c(9, 0.9), # 调整P值方法的位置，X=3，Y=0.85
                      conf.int = T,# 添加置信区间
                      risk.table = T, # 在图下方添加风险表
                      legend.labs=c("High","Low"), #表头标签注释分组
                      legend.title="Risk-group",#表头标签
                      #title="Overall survival",#改一下整体名称
                      ylab="Overall Survival (OS)",xlab = " Time (Years)",#修改X轴Y轴名称
                      risk.table.col = "strata", # 风险表加颜色
                      linetype = 1, # 生存曲线的线型
                      surv.median.line = "hv", # 标注出中位生存时间
                      ggtheme = theme_bw(base_size = 12, base_line_size = 1, base_rect_size = 1) + theme(panel.grid = element_blank()), #背景布局
                      palette = c(high_color, low_color))# 图形颜色风格
km_plot

## 2.8 PRJEB23709----
### 2.8.1 整理表达谱并计算评分----
PRJEB23709_expr <- readRDS("11.immunotherapy/Immunotherapy_datasets/Expression/Melanoma-PRJEB23709.Response.Rds")
which(duplicated(PRJEB23709_expr$GENE_SYMBOL))
PRJEB23709_expr <- PRJEB23709_expr %>% 
  tibble::column_to_rownames("GENE_SYMBOL")

log2_PRJEB23709_expr <- log2(PRJEB23709_expr + 1)
log2_PRJEB23709_expr <- as.data.frame(t(log2_PRJEB23709_expr))

expr_df <- log2_PRJEB23709_expr[,build_model_genes]
expr_df <- as.data.frame(scale(expr_df))
riskScore <- predict(fit, newdata = expr_df, type = "risk")$predicted
rs <- data.frame(Sample_ID = rownames(expr_df),
                 RiskScore = riskScore)
rs$RiskGroup <- ifelse(rs$RiskScore > median(rs$RiskScore), "High-Risk", "Low-Risk")
table(rs$RiskGroup)
rs$RiskGroup <- factor(rs$RiskGroup, levels = c("Low-Risk", "High-Risk"))

### 2.8.2 临床信息----
PRJEB23709_clinic <- read_tsv("11.immunotherapy/Immunotherapy_datasets/Clinical/Melanoma-PRJEB23709.Response.tsv")
colnames(PRJEB23709_clinic) <- c("Num", colnames(PRJEB23709_clinic))
PRJEB23709_clinic <- PRJEB23709_clinic[,c(2,6,7,10,11)]
colnames(PRJEB23709_clinic) <- c("Sample_ID", "Treatment", "Response", "OS.time", "OS")

table(PRJEB23709_clinic$Response)
PRJEB23709_clinic$Response <- ifelse(PRJEB23709_clinic$Response %in% c("PD", "SD"), "NR", "R")
table(PRJEB23709_clinic$Response)
class(PRJEB23709_clinic$Response)
PRJEB23709_clinic$Response <- factor(PRJEB23709_clinic$Response, levels = c("NR", "R"))

PRJEB23709_clinic$OS.time <- PRJEB23709_clinic$OS.time / 365

table(PRJEB23709_clinic$OS)
PRJEB23709_clinic$OS <- ifelse(PRJEB23709_clinic$OS == "Dead", 1, 0)

clinic <- PRJEB23709_clinic

### 2.8.3 合并与画图比较----
risk_clinic <- inner_join(rs, clinic, by = "Sample_ID")

comparison <- list(c("NR", "R"))
#Custom.color <- c("#245892","#877eac","#a9d296")
Custom.color <- c("#79aec1", "#ae3f51")
p11 <- ggplot(risk_clinic, aes(x = Response, y = RiskScore, fill = Response)) +
  geom_jitter(mapping = aes(color=Response),width = .05, alpha = 0.5,size=0.9)+ #绘制散点图
  geom_boxplot(position = position_nudge(x = 0.14),width=0.1,outlier.size = 0,outlier.alpha =0)+ #绘制箱线图，并通过position设置偏移
  stat_halfeye(mapping = aes(fill=Response),width = 0.2, .width = 0, justification = -1.2, point_colour = NA,alpha=0.6) + #绘制云雨图，并通过position设置偏移
  scale_fill_manual(values = Custom.color)+   #映射云雨图和箱线图的颜色
  scale_color_manual(values = Custom.color)+  #映射散点的颜色
  #expand_limits(x = c(0.5, 3.8))+ #扩展画板，若显示不全，请根据你的数据范围手动调整或删除此行
  #ylim(20,110)+ #控制y轴显示范围，若显示不全，请根据你的数据范围手动调整或删除此行
  xlab("") +  #设置X轴标题
  #ylab("Risk Score") +   #设置Y轴标题
  ggtitle("Melanoma-PRJEB23709")+  #设置主标题
  scale_y_continuous(
    name = "Risk Score",                     # 设置 Y 轴标题
    limits = c(30, 110),                    # 设置 Y 轴显示范围
    breaks = seq(30, 110, by = 20)          # 按步长 10 显示刻度
  ) +
  theme(axis.ticks.x = element_line(size = 0,color = "black"),  #自定义主题
        panel.background = element_rect(fill = "white", color = "white"),  #设置画板
        panel.grid.major.x = element_blank(),   #设置网格
        panel.grid.minor.x = element_blank(), #设置网格
        panel.grid.major.y = element_line(color = "gray", size = 0.25), #设置网格
        panel.grid.minor.y = element_blank(), #设置网格
        panel.border = element_rect(color = "black", fill = NA,linewidth = 1), #设置边框
        legend.position = "none", #隐藏图例
        axis.title.x = element_text(size = 13),  #调整X轴标题字体大小
        axis.title.y = element_text(size = 13), #调整Y轴标题字体大小
        axis.text.x = element_text(size = 12,hjust = 0.3), #设置x轴刻度字体偏移，若更换数据，可能需要重新设置
        axis.text.y = element_text(size = 12), #设置Y轴刻度字体大小
        plot.title = element_text(hjust = 0.5)
  )+
  stat_compare_means(comparisons = comparison, method = "wilcox.test", label = "p.format", label.y = 100)
p11
ggsave(filename = '11.immunotherapy/PRJEB23709_klaSigScore_response_boxplot.pdf', p11, height = 3, width = 2.5, dpi = 300)

### 2.8.4 生存分析----
surv_object <- Surv(time = risk_clinic$OS.time, event = risk_clinic$OS)
kmfit <- survfit(surv_object ~ RiskGroup, data = risk_clinic)

# 绘制生存曲线
high_color <- pal_npg("nrc")(10)[1]
low_color <- pal_npg("nrc")(10)[2]

km_plot <- ggsurvplot(kmfit, data = risk_clinic,
                      pval = T, # 在图上添加log rank检验的p值
                      pval.method = T,
                      #pval.coord = c(13.5, 0.9), # 调整P值位置，X=10，Y=0.1
                      #pval.method.coord = c(9, 0.9), # 调整P值方法的位置，X=3，Y=0.85
                      conf.int = T,# 添加置信区间
                      risk.table = T, # 在图下方添加风险表
                      legend.labs=c("High","Low"), #表头标签注释分组
                      legend.title="Risk-group",#表头标签
                      #title="Overall survival",#改一下整体名称
                      ylab="Overall Survival (OS)",xlab = " Time (Years)",#修改X轴Y轴名称
                      risk.table.col = "strata", # 风险表加颜色
                      linetype = 1, # 生存曲线的线型
                      surv.median.line = "hv", # 标注出中位生存时间
                      ggtheme = theme_bw(base_size = 12, base_line_size = 1, base_rect_size = 1) + theme(panel.grid = element_blank()), #背景布局
                      palette = c(high_color, low_color))# 图形颜色风格
km_plot

## 2.9 Braun----
### 2.9.1 整理表达谱并计算评分----
braun_expr <- readRDS("11.immunotherapy/Immunotherapy_datasets/Expression/RCC-Braun_2020.Response.Rds")
which(duplicated(braun_expr$GENE_SYMBOL))
braun_expr <- braun_expr %>% 
  tibble::column_to_rownames("GENE_SYMBOL")

log2_braun_expr <- log2(braun_expr + 1)
log2_braun_expr <- as.data.frame(t(log2_braun_expr))

expr_df <- log2_braun_expr[,build_model_genes]
expr_df <- as.data.frame(scale(expr_df))
riskScore <- predict(fit, newdata = expr_df, type = "risk")$predicted
rs <- data.frame(Sample_ID = rownames(expr_df),
                 RiskScore = riskScore)
rs$RiskGroup <- ifelse(rs$RiskScore > median(rs$RiskScore), "High-Risk", "Low-Risk")
table(rs$RiskGroup)
rs$RiskGroup <- factor(rs$RiskGroup, levels = c("Low-Risk", "High-Risk"))

### 2.9.2 临床信息----
braun_clinic <- read_tsv("11.immunotherapy/Immunotherapy_datasets/Clinical/RCC-Braun_2020.Response.tsv")
colnames(braun_clinic) <- c("Num", colnames(braun_clinic))
braun_clinic <- braun_clinic[,c(2,6,7,10,11)]
colnames(braun_clinic) <- c("Sample_ID", "Treatment", "Response", "OS.time", "OS")

table(braun_clinic$Response)
braun_clinic$Response <- ifelse(braun_clinic$Response %in% c("PD", "SD"), "NR", 
                                ifelse(braun_clinic$Response == "NE", "UNK", "R"))
table(braun_clinic$Response)
braun_clinic <- braun_clinic %>% 
  filter(Response != "UNK")
table(braun_clinic$Response)

class(braun_clinic$Response)
braun_clinic$Response <- factor(braun_clinic$Response, levels = c("NR", "R"))

braun_clinic$OS.time <- braun_clinic$OS.time / 365

table(braun_clinic$OS)
braun_clinic$OS <- ifelse(braun_clinic$OS == "Dead", 1, 0)

clinic <- braun_clinic

### 2.9.3 合并与画图比较----
risk_clinic <- inner_join(rs, clinic, by = "Sample_ID")

comparison <- list(c("NR", "R"))
#Custom.color <- c("#245892","#877eac","#a9d296")
Custom.color <- c("#79aec1", "#ae3f51")
p12 <- ggplot(risk_clinic, aes(x = Response, y = RiskScore, fill = Response)) +
  geom_jitter(mapping = aes(color=Response),width = .05, alpha = 0.5,size=0.9)+ #绘制散点图
  geom_boxplot(position = position_nudge(x = 0.14),width=0.1,outlier.size = 0,outlier.alpha =0)+ #绘制箱线图，并通过position设置偏移
  stat_halfeye(mapping = aes(fill=Response),width = 0.2, .width = 0, justification = -1.2, point_colour = NA,alpha=0.6) + #绘制云雨图，并通过position设置偏移
  scale_fill_manual(values = Custom.color)+   #映射云雨图和箱线图的颜色
  scale_color_manual(values = Custom.color)+  #映射散点的颜色
  #expand_limits(x = c(0.5, 3.8))+ #扩展画板，若显示不全，请根据你的数据范围手动调整或删除此行
  #ylim(20,110)+ #控制y轴显示范围，若显示不全，请根据你的数据范围手动调整或删除此行
  xlab("") +  #设置X轴标题
  #ylab("Risk Score") +   #设置Y轴标题
  ggtitle("RCC-Braun 2020")+  #设置主标题
  scale_y_continuous(
    name = "Risk Score",                     # 设置 Y 轴标题
    limits = c(30, 110),                    # 设置 Y 轴显示范围
    breaks = seq(30, 110, by = 20)          # 按步长 10 显示刻度
  ) +
  theme(axis.ticks.x = element_line(size = 0,color = "black"),  #自定义主题
        panel.background = element_rect(fill = "white", color = "white"),  #设置画板
        panel.grid.major.x = element_blank(),   #设置网格
        panel.grid.minor.x = element_blank(), #设置网格
        panel.grid.major.y = element_line(color = "gray", size = 0.25), #设置网格
        panel.grid.minor.y = element_blank(), #设置网格
        panel.border = element_rect(color = "black", fill = NA,linewidth = 1), #设置边框
        legend.position = "none", #隐藏图例
        axis.title.x = element_text(size = 13),  #调整X轴标题字体大小
        axis.title.y = element_text(size = 13), #调整Y轴标题字体大小
        axis.text.x = element_text(size = 12,hjust = 0.3), #设置x轴刻度字体偏移，若更换数据，可能需要重新设置
        axis.text.y = element_text(size = 12), #设置Y轴刻度字体大小
        plot.title = element_text(hjust = 0.5)
  )+
  stat_compare_means(comparisons = comparison, method = "wilcox.test", label = "p.format", label.y = 100)
p12
#ggsave(filename = '11.immunotherapy/Braun_klaSigScore_response_boxplot.pdf', p12, height = 3, width = 2.5, dpi = 300)

### 2.9.4 生存分析----
surv_object <- Surv(time = risk_clinic$OS.time, event = risk_clinic$OS)
kmfit <- survfit(surv_object ~ RiskGroup, data = risk_clinic)

# 绘制生存曲线
high_color <- pal_npg("nrc")(10)[1]
low_color <- pal_npg("nrc")(10)[2]

km_plot <- ggsurvplot(kmfit, data = risk_clinic,
                      pval = T, # 在图上添加log rank检验的p值
                      pval.method = T,
                      #pval.coord = c(13.5, 0.9), # 调整P值位置，X=10，Y=0.1
                      #pval.method.coord = c(9, 0.9), # 调整P值方法的位置，X=3，Y=0.85
                      conf.int = T,# 添加置信区间
                      risk.table = T, # 在图下方添加风险表
                      legend.labs=c("High","Low"), #表头标签注释分组
                      legend.title="Risk-group",#表头标签
                      #title="Overall survival",#改一下整体名称
                      ylab="Overall Survival (OS)",xlab = " Time (Years)",#修改X轴Y轴名称
                      risk.table.col = "strata", # 风险表加颜色
                      linetype = 1, # 生存曲线的线型
                      surv.median.line = "hv", # 标注出中位生存时间
                      ggtheme = theme_bw(base_size = 12, base_line_size = 1, base_rect_size = 1) + theme(panel.grid = element_blank()), #背景布局
                      palette = c(high_color, low_color))# 图形颜色风格
km_plot

## 2.10 PRJNA482620----
### 2.10.1 整理表达谱并计算评分----
PRJNA482620_expr <- readRDS("11.immunotherapy/Immunotherapy_datasets/Expression/GBM-PRJNA482620.Response.Rds")
which(duplicated(PRJNA482620_expr$GENE_SYMBOL))
PRJNA482620_expr <- PRJNA482620_expr %>% 
  tibble::column_to_rownames("GENE_SYMBOL")

log2_PRJNA482620_expr <- log2(PRJNA482620_expr + 1)
log2_PRJNA482620_expr <- as.data.frame(t(log2_PRJNA482620_expr))

expr_df <- log2_PRJNA482620_expr[,build_model_genes]
expr_df <- as.data.frame(scale(expr_df))
riskScore <- predict(fit, newdata = expr_df, type = "risk")$predicted
rs <- data.frame(Sample_ID = rownames(expr_df),
                 RiskScore = riskScore)
rs$RiskGroup <- ifelse(rs$RiskScore > median(rs$RiskScore), "High-Risk", "Low-Risk")
table(rs$RiskGroup)
rs$RiskGroup <- factor(rs$RiskGroup, levels = c("Low-Risk", "High-Risk"))

### 2.10.2 临床信息----
PRJNA482620_clinic <- read_tsv("11.immunotherapy/Immunotherapy_datasets/Clinical/GBM-PRJNA482620.Response.tsv")
colnames(PRJNA482620_clinic) <- c("Num", colnames(PRJNA482620_clinic))
PRJNA482620_clinic <- PRJNA482620_clinic[,c(2,6,7,10,11)]
colnames(PRJNA482620_clinic) <- c("Sample_ID", "Treatment", "Response", "OS.time", "OS")

table(PRJNA482620_clinic$Response)
PRJNA482620_clinic$Response <- ifelse(PRJNA482620_clinic$Response == "N", "NR", "R")
table(PRJNA482620_clinic$Response)
class(PRJNA482620_clinic$Response)
PRJNA482620_clinic$Response <- factor(PRJNA482620_clinic$Response, levels = c("NR", "R"))

PRJNA482620_clinic$OS.time <- PRJNA482620_clinic$OS.time / 365

table(PRJNA482620_clinic$OS)
PRJNA482620_clinic$OS <- ifelse(PRJNA482620_clinic$OS == "Dead", 1, 0)

clinic <- PRJNA482620_clinic

### 2.10.3 合并与画图比较----
risk_clinic <- inner_join(rs, clinic, by = "Sample_ID")

comparison <- list(c("NR", "R"))
#Custom.color <- c("#245892","#877eac","#a9d296")
Custom.color <- c("#79aec1", "#ae3f51")
p13 <- ggplot(risk_clinic, aes(x = Response, y = RiskScore, fill = Response)) +
  geom_jitter(mapping = aes(color=Response),width = .05, alpha = 0.5,size=0.9)+ #绘制散点图
  geom_boxplot(position = position_nudge(x = 0.14),width=0.1,outlier.size = 0,outlier.alpha =0)+ #绘制箱线图，并通过position设置偏移
  stat_halfeye(mapping = aes(fill=Response),width = 0.2, .width = 0, justification = -1.2, point_colour = NA,alpha=0.6) + #绘制云雨图，并通过position设置偏移
  scale_fill_manual(values = Custom.color)+   #映射云雨图和箱线图的颜色
  scale_color_manual(values = Custom.color)+  #映射散点的颜色
  #expand_limits(x = c(0.5, 3.8))+ #扩展画板，若显示不全，请根据你的数据范围手动调整或删除此行
  #ylim(20,110)+ #控制y轴显示范围，若显示不全，请根据你的数据范围手动调整或删除此行
  xlab("") +  #设置X轴标题
  #ylab("Risk Score") +   #设置Y轴标题
  ggtitle("GBM-PRJNA482620")+  #设置主标题
  scale_y_continuous(
    name = "Risk Score",                     # 设置 Y 轴标题
    limits = c(30, 110),                    # 设置 Y 轴显示范围
    breaks = seq(30, 110, by = 20)          # 按步长 10 显示刻度
  ) +
  theme(axis.ticks.x = element_line(size = 0,color = "black"),  #自定义主题
        panel.background = element_rect(fill = "white", color = "white"),  #设置画板
        panel.grid.major.x = element_blank(),   #设置网格
        panel.grid.minor.x = element_blank(), #设置网格
        panel.grid.major.y = element_line(color = "gray", size = 0.25), #设置网格
        panel.grid.minor.y = element_blank(), #设置网格
        panel.border = element_rect(color = "black", fill = NA,linewidth = 1), #设置边框
        legend.position = "none", #隐藏图例
        axis.title.x = element_text(size = 13),  #调整X轴标题字体大小
        axis.title.y = element_text(size = 13), #调整Y轴标题字体大小
        axis.text.x = element_text(size = 12,hjust = 0.3), #设置x轴刻度字体偏移，若更换数据，可能需要重新设置
        axis.text.y = element_text(size = 12), #设置Y轴刻度字体大小
        plot.title = element_text(hjust = 0.5)
  )+
  stat_compare_means(comparisons = comparison, method = "wilcox.test", label = "p.format", label.y = 100)
p13
#ggsave(filename = '11.immunotherapy/PRJNA482620_klaSigScore_response_boxplot.pdf', p13, height = 3, width = 2.5, dpi = 300)

### 2.10.4 生存分析----
surv_object <- Surv(time = risk_clinic$OS.time, event = risk_clinic$OS)
kmfit <- survfit(surv_object ~ RiskGroup, data = risk_clinic)

# 绘制生存曲线
high_color <- pal_npg("nrc")(10)[1]
low_color <- pal_npg("nrc")(10)[2]

km_plot <- ggsurvplot(kmfit, data = risk_clinic,
                      pval = T, # 在图上添加log rank检验的p值
                      pval.method = T,
                      #pval.coord = c(13.5, 0.9), # 调整P值位置，X=10，Y=0.1
                      #pval.method.coord = c(9, 0.9), # 调整P值方法的位置，X=3，Y=0.85
                      conf.int = T,# 添加置信区间
                      risk.table = T, # 在图下方添加风险表
                      legend.labs=c("High","Low"), #表头标签注释分组
                      legend.title="Risk-group",#表头标签
                      #title="Overall survival",#改一下整体名称
                      ylab="Overall Survival (OS)",xlab = " Time (Years)",#修改X轴Y轴名称
                      risk.table.col = "strata", # 风险表加颜色
                      linetype = 1, # 生存曲线的线型
                      surv.median.line = "hv", # 标注出中位生存时间
                      ggtheme = theme_bw(base_size = 12, base_line_size = 1, base_rect_size = 1) + theme(panel.grid = element_blank()), #背景布局
                      palette = c(high_color, low_color))# 图形颜色风格
km_plot

## 2.11 PRJEB25780----
### 2.11.1 整理表达谱并计算评分----
PRJEB25780_expr <- readRDS("11.immunotherapy/Immunotherapy_datasets/Expression/STAD-PRJEB25780.Response.Rds")
which(duplicated(PRJEB25780_expr$GENE_SYMBOL))
PRJEB25780_expr <- PRJEB25780_expr %>% 
  tibble::column_to_rownames("GENE_SYMBOL")

log2_PRJEB25780_expr <- log2(PRJEB25780_expr + 1)
log2_PRJEB25780_expr <- as.data.frame(t(log2_PRJEB25780_expr))

expr_df <- log2_PRJEB25780_expr[,build_model_genes]
expr_df <- as.data.frame(scale(expr_df))
riskScore <- predict(fit, newdata = expr_df, type = "risk")$predicted
rs <- data.frame(Sample_ID = rownames(expr_df),
                 RiskScore = riskScore)
rs$RiskGroup <- ifelse(rs$RiskScore > median(rs$RiskScore), "High-Risk", "Low-Risk")
table(rs$RiskGroup)
rs$RiskGroup <- factor(rs$RiskGroup, levels = c("Low-Risk", "High-Risk"))

### 2.11.2 临床信息----
PRJEB25780_clinic <- read_tsv("11.immunotherapy/Immunotherapy_datasets/Clinical/STAD-PRJEB25780.Response.tsv")
colnames(PRJEB25780_clinic) <- c("Num", colnames(PRJEB25780_clinic))
PRJEB25780_clinic <- PRJEB25780_clinic[,c(2,6,7)]
colnames(PRJEB25780_clinic) <- c("Sample_ID", "Treatment", "Response")

table(PRJEB25780_clinic$Treatment)
PRJEB25780_clinic <- PRJEB25780_clinic[1:45,]

table(PRJEB25780_clinic$Response)
PRJEB25780_clinic$Response <- ifelse(PRJEB25780_clinic$Response %in% c("PD", "SD"), "NR", "R")
table(PRJEB25780_clinic$Response)
class(PRJEB25780_clinic$Response)
PRJEB25780_clinic$Response <- factor(PRJEB25780_clinic$Response, levels = c("NR", "R"))

clinic <- PRJEB25780_clinic

### 2.11.3 合并与画图比较----
risk_clinic <- inner_join(rs, clinic, by = "Sample_ID")

comparison <- list(c("NR", "R"))
#Custom.color <- c("#245892","#877eac","#a9d296")
Custom.color <- c("#79aec1", "#ae3f51")
p14 <- ggplot(risk_clinic, aes(x = Response, y = RiskScore, fill = Response)) +
  geom_jitter(mapping = aes(color=Response),width = .05, alpha = 0.5,size=0.9)+ #绘制散点图
  geom_boxplot(position = position_nudge(x = 0.14),width=0.1,outlier.size = 0,outlier.alpha =0)+ #绘制箱线图，并通过position设置偏移
  stat_halfeye(mapping = aes(fill=Response),width = 0.2, .width = 0, justification = -1.2, point_colour = NA,alpha=0.6) + #绘制云雨图，并通过position设置偏移
  scale_fill_manual(values = Custom.color)+   #映射云雨图和箱线图的颜色
  scale_color_manual(values = Custom.color)+  #映射散点的颜色
  #expand_limits(x = c(0.5, 3.8))+ #扩展画板，若显示不全，请根据你的数据范围手动调整或删除此行
  #ylim(20,110)+ #控制y轴显示范围，若显示不全，请根据你的数据范围手动调整或删除此行
  xlab("") +  #设置X轴标题
  #ylab("Risk Score") +   #设置Y轴标题
  ggtitle("STAD-PRJEB25780")+  #设置主标题
  scale_y_continuous(
    name = "Risk Score",                     # 设置 Y 轴标题
    limits = c(30, 110),                    # 设置 Y 轴显示范围
    breaks = seq(30, 110, by = 20)          # 按步长 10 显示刻度
  ) +
  theme(axis.ticks.x = element_line(size = 0,color = "black"),  #自定义主题
        panel.background = element_rect(fill = "white", color = "white"),  #设置画板
        panel.grid.major.x = element_blank(),   #设置网格
        panel.grid.minor.x = element_blank(), #设置网格
        panel.grid.major.y = element_line(color = "gray", size = 0.25), #设置网格
        panel.grid.minor.y = element_blank(), #设置网格
        panel.border = element_rect(color = "black", fill = NA,linewidth = 1), #设置边框
        legend.position = "none", #隐藏图例
        axis.title.x = element_text(size = 13),  #调整X轴标题字体大小
        axis.title.y = element_text(size = 13), #调整Y轴标题字体大小
        axis.text.x = element_text(size = 12,hjust = 0.3), #设置x轴刻度字体偏移，若更换数据，可能需要重新设置
        axis.text.y = element_text(size = 12), #设置Y轴刻度字体大小
        plot.title = element_text(hjust = 0.5)
  )+
  stat_compare_means(comparisons = comparison, method = "wilcox.test", label = "p.format", label.y = 100)
p14
#ggsave(filename = '11.immunotherapy/PRJEB25780_klaSigScore_response_boxplot.pdf', p14, height = 3, width = 2.5, dpi = 300)
