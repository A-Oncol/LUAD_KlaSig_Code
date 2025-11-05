######目的：Kla评分与免疫微环境
######作者：申奥
######日期：2025-02-05
######Server: yjx@xiyoucloud.top:34022
######Version: R 4.4.2


# 0.环境设置----
rm(list = ls())
options(stringAsFactors = F)

library(qs)
library(openxlsx)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(ggsignif)
library(ggdist)
library(IOBR)
packageVersion("IOBR") # ‘0.99.0’


# 1.比较TCGA官方计算好的指标----
## 1.1 加载risk score评分----
load("12.Prognosis_model/RSF_riskScore.RData")
load("12.Prognosis_model/merge_datasets.RData")

identical(merge_datasets[["TCGA_LUAD"]]$OS, rs[["TCGA_LUAD"]]$OS)
identical(merge_datasets[["TCGA_LUAD"]]$OS.time, rs[["TCGA_LUAD"]]$OS.time)
TCGA_LUAD_riskScore <- data.frame(cbind(Sample_ID = merge_datasets[["TCGA_LUAD"]]$Sample_ID, rs[["TCGA_LUAD"]]))
TCGA_LUAD_riskScore$Risk_group <- ifelse(TCGA_LUAD_riskScore$RS > median(TCGA_LUAD_riskScore$RS), "High-Risk", "Low-Risk")
table(TCGA_LUAD_riskScore$Risk_group)

qsave(TCGA_LUAD_riskScore, file = '32.TCGA/TCGA_LUAD_riskScore.qs')

TCGA_LUAD_riskScore$Sample_ID <- str_sub(TCGA_LUAD_riskScore$Sample_ID, 1, 12)

class(TCGA_LUAD_riskScore$Risk_group)
TCGA_LUAD_riskScore$Risk_group <- factor(TCGA_LUAD_riskScore$Risk_group, levels = c("Low-Risk", "High-Risk"))

## 1.2 新抗原----
immunity2018 <- read.xlsx("32.TCGA/NIHMS958212-supplement-2.xlsx", sheet = 1) # Immunity2018的附表2(PMID: 29628290)
table(immunity2018$TCGA.Study)

luad <- immunity2018[which(immunity2018$`TCGA.Study` == "LUAD"),]
colnames(luad)[1] <- "Sample_ID"
luad_SNVneoAg <- na.omit(luad[,c("Sample_ID", "SNV.Neoantigens")])
luad_IndelneoAg <- na.omit(luad[,c("Sample_ID", "Indel.Neoantigens")])

luad_SNVneoAg_riskScore <- inner_join(luad_SNVneoAg, TCGA_LUAD_riskScore, by = "Sample_ID")
luad_IndelneoAg_riskScore <- inner_join(luad_IndelneoAg, TCGA_LUAD_riskScore, by = "Sample_ID")

comparison <- list(c("Low-Risk", "High-Risk"))
#Custom.color <- c("#245892","#877eac","#a9d296")
Custom.color <- c("#79aec1", "#ae3f51")
p1 <- ggplot(luad_SNVneoAg_riskScore, aes(x = Risk_group, y = SNV.Neoantigens, fill = Risk_group)) +
  geom_jitter(mapping = aes(color=Risk_group),width = .05, alpha = 0.5,size=0.9)+ #绘制散点图
  geom_boxplot(position = position_nudge(x = 0.14),width=0.1,outlier.size = 0,outlier.alpha =0)+ #绘制箱线图，并通过position设置偏移
  stat_halfeye(mapping = aes(fill=Risk_group),width = 0.2, .width = 0, justification = -1.2, point_colour = NA,alpha=0.6) + #绘制云雨图，并通过position设置偏移
  scale_fill_manual(values = Custom.color)+   #映射云雨图和箱线图的颜色
  scale_color_manual(values = Custom.color)+  #映射散点的颜色
  #expand_limits(x = c(0.5, 3.8))+ #扩展画板，若显示不全，请根据你的数据范围手动调整或删除此行
  #ylim(20,110)+ #控制y轴显示范围，若显示不全，请根据你的数据范围手动调整或删除此行
  xlab("") +  #设置X轴标题
  #ylab("Risk Score") +   #设置Y轴标题
  #ggtitle("NSCLC-GSE126044")+  #设置主标题
  # scale_y_continuous(
  #   name = "SNV neoantigens",                     # 设置 Y 轴标题
  #   limits = c(30, 110),                    # 设置 Y 轴显示范围
  #   breaks = seq(30, 110, by = 20)          # 按步长 10 显示刻度
  # ) +
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
  stat_compare_means(comparisons = comparison, method = "wilcox.test", label = "p.format")
p1
ggsave(filename = '32.TCGA/TCGA-LUAD_klaSigScore_SNVneoAg_boxplot.pdf', p1, height = 3, width = 2.5, dpi = 300)

p2 <- ggplot(luad_IndelneoAg_riskScore, aes(x = Risk_group, y = Indel.Neoantigens, fill = Risk_group)) +
  geom_jitter(mapping = aes(color=Risk_group),width = .05, alpha = 0.5,size=0.9)+ #绘制散点图
  geom_boxplot(position = position_nudge(x = 0.14),width=0.1,outlier.size = 0,outlier.alpha =0)+ #绘制箱线图，并通过position设置偏移
  stat_halfeye(mapping = aes(fill=Risk_group),width = 0.2, .width = 0, justification = -1.2, point_colour = NA,alpha=0.6) + #绘制云雨图，并通过position设置偏移
  scale_fill_manual(values = Custom.color)+   #映射云雨图和箱线图的颜色
  scale_color_manual(values = Custom.color)+  #映射散点的颜色
  #expand_limits(x = c(0.5, 3.8))+ #扩展画板，若显示不全，请根据你的数据范围手动调整或删除此行
  #ylim(20,110)+ #控制y轴显示范围，若显示不全，请根据你的数据范围手动调整或删除此行
  xlab("") +  #设置X轴标题
  #ylab("Risk Score") +   #设置Y轴标题
  #ggtitle("NSCLC-GSE126044")+  #设置主标题
  # scale_y_continuous(
  #   name = "SNV neoantigens",                     # 设置 Y 轴标题
  #   limits = c(30, 110),                    # 设置 Y 轴显示范围
  #   breaks = seq(30, 110, by = 20)          # 按步长 10 显示刻度
  # ) +
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
  stat_compare_means(comparisons = comparison, method = "wilcox.test", label = "p.format")
p2
ggsave(filename = '32.TCGA/TCGA-LUAD_klaSigScore_IndelneoAg_boxplot.pdf', p2, height = 3, width = 2.5, dpi = 300)

## 1.3 TMB----
rm(luad_SNVneoAg, luad_SNVneoAg_riskScore, luad_IndelneoAg, luad_IndelneoAg_riskScore)

luad_tmb <- na.omit(luad[,c("Sample_ID", "Nonsilent.Mutation.Rate")])
luad_tmb_riskScore <- inner_join(luad_tmb, TCGA_LUAD_riskScore, by = "Sample_ID")

p3 <- ggplot(luad_tmb_riskScore, aes(x = Risk_group, y = Nonsilent.Mutation.Rate, fill = Risk_group)) +
  geom_jitter(mapping = aes(color=Risk_group),width = .05, alpha = 0.5,size=0.9)+ #绘制散点图
  geom_boxplot(position = position_nudge(x = 0.14),width=0.1,outlier.size = 0,outlier.alpha =0)+ #绘制箱线图，并通过position设置偏移
  stat_halfeye(mapping = aes(fill=Risk_group),width = 0.2, .width = 0, justification = -1.2, point_colour = NA,alpha=0.6) + #绘制云雨图，并通过position设置偏移
  scale_fill_manual(values = Custom.color)+   #映射云雨图和箱线图的颜色
  scale_color_manual(values = Custom.color)+  #映射散点的颜色
  #expand_limits(x = c(0.5, 3.8))+ #扩展画板，若显示不全，请根据你的数据范围手动调整或删除此行
  #ylim(20,110)+ #控制y轴显示范围，若显示不全，请根据你的数据范围手动调整或删除此行
  xlab("") +  #设置X轴标题
  #ylab("Risk Score") +   #设置Y轴标题
  #ggtitle("NSCLC-GSE126044")+  #设置主标题
  # scale_y_continuous(
  #   name = "SNV neoantigens",                     # 设置 Y 轴标题
  #   limits = c(30, 110),                    # 设置 Y 轴显示范围
  #   breaks = seq(30, 110, by = 20)          # 按步长 10 显示刻度
  # ) +
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
  stat_compare_means(comparisons = comparison, method = "wilcox.test", label = "p.format")
p3
ggsave(filename = '32.TCGA/TCGA-LUAD_klaSigScore_NonsilentMutationRate_boxplot.pdf', p3, height = 3, width = 2.5, dpi = 300)

## 1.4 ITH----
rm(luad_tmb, luad_tmb_riskScore)

luad_ith <- na.omit(luad[,c("Sample_ID", "Intratumor.Heterogeneity")])
luad_ith_riskScore <- inner_join(luad_ith, TCGA_LUAD_riskScore, by = "Sample_ID")

p4 <- ggplot(luad_ith_riskScore, aes(x = Risk_group, y = Intratumor.Heterogeneity, fill = Risk_group)) +
  geom_jitter(mapping = aes(color=Risk_group),width = .05, alpha = 0.5,size=0.9)+ #绘制散点图
  geom_boxplot(position = position_nudge(x = 0.14),width=0.1,outlier.size = 0,outlier.alpha =0)+ #绘制箱线图，并通过position设置偏移
  stat_halfeye(mapping = aes(fill=Risk_group),width = 0.2, .width = 0, justification = -1.2, point_colour = NA,alpha=0.6) + #绘制云雨图，并通过position设置偏移
  scale_fill_manual(values = Custom.color)+   #映射云雨图和箱线图的颜色
  scale_color_manual(values = Custom.color)+  #映射散点的颜色
  #expand_limits(x = c(0.5, 3.8))+ #扩展画板，若显示不全，请根据你的数据范围手动调整或删除此行
  #ylim(20,110)+ #控制y轴显示范围，若显示不全，请根据你的数据范围手动调整或删除此行
  xlab("") +  #设置X轴标题
  #ylab("Risk Score") +   #设置Y轴标题
  #ggtitle("NSCLC-GSE126044")+  #设置主标题
  # scale_y_continuous(
  #   name = "SNV neoantigens",                     # 设置 Y 轴标题
  #   limits = c(30, 110),                    # 设置 Y 轴显示范围
  #   breaks = seq(30, 110, by = 20)          # 按步长 10 显示刻度
  # ) +
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
  stat_compare_means(comparisons = comparison, method = "wilcox.test", label = "p.format")
p4
ggsave(filename = '32.TCGA/TCGA-LUAD_klaSigScore_ITH_boxplot.pdf', p4, height = 3, width = 2.5, dpi = 300)

## 1.5 倍性----
rm(luad_ith, luad_ith_riskScore)

luad_as <- na.omit(luad[,c("Sample_ID", "Aneuploidy.Score")])
luad_as_riskScore <- inner_join(luad_as, TCGA_LUAD_riskScore, by = "Sample_ID")

p5 <- ggplot(luad_as_riskScore, aes(x = Risk_group, y = Aneuploidy.Score, fill = Risk_group)) +
  geom_jitter(mapping = aes(color=Risk_group),width = .05, alpha = 0.5,size=0.9)+ #绘制散点图
  geom_boxplot(position = position_nudge(x = 0.14),width=0.1,outlier.size = 0,outlier.alpha =0)+ #绘制箱线图，并通过position设置偏移
  stat_halfeye(mapping = aes(fill=Risk_group),width = 0.2, .width = 0, justification = -1.2, point_colour = NA,alpha=0.6) + #绘制云雨图，并通过position设置偏移
  scale_fill_manual(values = Custom.color)+   #映射云雨图和箱线图的颜色
  scale_color_manual(values = Custom.color)+  #映射散点的颜色
  #expand_limits(x = c(0.5, 3.8))+ #扩展画板，若显示不全，请根据你的数据范围手动调整或删除此行
  #ylim(20,110)+ #控制y轴显示范围，若显示不全，请根据你的数据范围手动调整或删除此行
  xlab("") +  #设置X轴标题
  #ylab("Risk Score") +   #设置Y轴标题
  #ggtitle("NSCLC-GSE126044")+  #设置主标题
  scale_y_continuous(
    name = "Aneuploidy.Score",                     # 设置 Y 轴标题
    limits = c(0, 32),                    # 设置 Y 轴显示范围
    breaks = seq(0, 32, by = 10)          # 按步长 10 显示刻度
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
  stat_compare_means(comparisons = comparison, method = "wilcox.test", label = "p.format", label.y = 30)
p5
ggsave(filename = '32.TCGA/TCGA-LUAD_klaSigScore_AneuploidyScore_boxplot.pdf', p5, height = 3, width = 2.5, dpi = 300)

## 1.6 TCR/BCR----
rm(luad_as, luad_as_riskScore)

luad_sub <- na.omit(luad[,c("Sample_ID", "TCR.Shannon")])
luad_sub_riskScore <- inner_join(luad_sub, TCGA_LUAD_riskScore, by = "Sample_ID")

p6 <- ggplot(luad_sub_riskScore, aes(x = Risk_group, y = TCR.Shannon, fill = Risk_group)) +
  geom_jitter(mapping = aes(color=Risk_group),width = .05, alpha = 0.5,size=0.9)+ #绘制散点图
  geom_boxplot(position = position_nudge(x = 0.14),width=0.1,outlier.size = 0,outlier.alpha =0)+ #绘制箱线图，并通过position设置偏移
  stat_halfeye(mapping = aes(fill=Risk_group),width = 0.2, .width = 0, justification = -1.2, point_colour = NA,alpha=0.6) + #绘制云雨图，并通过position设置偏移
  scale_fill_manual(values = Custom.color)+   #映射云雨图和箱线图的颜色
  scale_color_manual(values = Custom.color)+  #映射散点的颜色
  #expand_limits(x = c(0.5, 3.8))+ #扩展画板，若显示不全，请根据你的数据范围手动调整或删除此行
  #ylim(20,110)+ #控制y轴显示范围，若显示不全，请根据你的数据范围手动调整或删除此行
  xlab("") +  #设置X轴标题
  #ylab("Risk Score") +   #设置Y轴标题
  #ggtitle("NSCLC-GSE126044")+  #设置主标题
  scale_y_continuous(
    name = "TCR.Shannon",                     # 设置 Y 轴标题
    limits = c(0, 6),                    # 设置 Y 轴显示范围
    breaks = seq(0, 6, by = 2)          # 按步长 10 显示刻度
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
  stat_compare_means(comparisons = comparison, method = "wilcox.test", label = "p.format", label.y = 5.5)
p6
ggsave(filename = '32.TCGA/TCGA-LUAD_klaSigScore_TCRshannon_boxplot.pdf', p6, height = 3, width = 2.5, dpi = 300)

luad_sub <- na.omit(luad[,c("Sample_ID", "TCR.Richness")])
luad_sub_riskScore <- inner_join(luad_sub, TCGA_LUAD_riskScore, by = "Sample_ID")

p7 <- ggplot(luad_sub_riskScore, aes(x = Risk_group, y = TCR.Richness, fill = Risk_group)) +
  geom_jitter(mapping = aes(color=Risk_group),width = .05, alpha = 0.5,size=0.9)+ #绘制散点图
  geom_boxplot(position = position_nudge(x = 0.14),width=0.1,outlier.size = 0,outlier.alpha =0)+ #绘制箱线图，并通过position设置偏移
  stat_halfeye(mapping = aes(fill=Risk_group),width = 0.2, .width = 0, justification = -1.2, point_colour = NA,alpha=0.6) + #绘制云雨图，并通过position设置偏移
  scale_fill_manual(values = Custom.color)+   #映射云雨图和箱线图的颜色
  scale_color_manual(values = Custom.color)+  #映射散点的颜色
  #expand_limits(x = c(0.5, 3.8))+ #扩展画板，若显示不全，请根据你的数据范围手动调整或删除此行
  #ylim(20,110)+ #控制y轴显示范围，若显示不全，请根据你的数据范围手动调整或删除此行
  xlab("") +  #设置X轴标题
  #ylab("Risk Score") +   #设置Y轴标题
  #ggtitle("NSCLC-GSE126044")+  #设置主标题
  # scale_y_continuous(
  #   name = "TCR.Shannon",                     # 设置 Y 轴标题
  #   limits = c(0, 6),                    # 设置 Y 轴显示范围
  #   breaks = seq(0, 6, by = 2)          # 按步长 10 显示刻度
  # ) +
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
  stat_compare_means(comparisons = comparison, method = "wilcox.test", label = "p.format", label.y = 135)
p7
ggsave(filename = '32.TCGA/TCGA-LUAD_klaSigScore_TCRrichness_boxplot.pdf', p7, height = 3, width = 2.5, dpi = 300)

luad_sub <- na.omit(luad[,c("Sample_ID", "TCR.Evenness")])
luad_sub_riskScore <- inner_join(luad_sub, TCGA_LUAD_riskScore, by = "Sample_ID")

p8 <- ggplot(luad_sub_riskScore, aes(x = Risk_group, y = TCR.Evenness, fill = Risk_group)) +
  geom_jitter(mapping = aes(color=Risk_group),width = .05, alpha = 0.5,size=0.9)+ #绘制散点图
  geom_boxplot(position = position_nudge(x = 0.14),width=0.1,outlier.size = 0,outlier.alpha =0)+ #绘制箱线图，并通过position设置偏移
  stat_halfeye(mapping = aes(fill=Risk_group),width = 0.2, .width = 0, justification = -1.2, point_colour = NA,alpha=0.6) + #绘制云雨图，并通过position设置偏移
  scale_fill_manual(values = Custom.color)+   #映射云雨图和箱线图的颜色
  scale_color_manual(values = Custom.color)+  #映射散点的颜色
  #expand_limits(x = c(0.5, 3.8))+ #扩展画板，若显示不全，请根据你的数据范围手动调整或删除此行
  #ylim(20,110)+ #控制y轴显示范围，若显示不全，请根据你的数据范围手动调整或删除此行
  xlab("") +  #设置X轴标题
  #ylab("Risk Score") +   #设置Y轴标题
  #ggtitle("NSCLC-GSE126044")+  #设置主标题
  scale_y_continuous(
    name = "TCR.Evenness",                     # 设置 Y 轴标题
    limits = c(0.6, 1.1),                    # 设置 Y 轴显示范围
    breaks = seq(0.6, 1.1, by = 0.2)          # 按步长 10 显示刻度
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
  stat_compare_means(comparisons = comparison, method = "wilcox.test", label = "p.format", label.y = 1)
p8
ggsave(filename = '32.TCGA/TCGA-LUAD_klaSigScore_TCRevenness_boxplot.pdf', p8, height = 3, width = 2.5, dpi = 300)

luad_sub <- na.omit(luad[,c("Sample_ID", "BCR.Shannon")])
luad_sub_riskScore <- inner_join(luad_sub, TCGA_LUAD_riskScore, by = "Sample_ID")

p9 <- ggplot(luad_sub_riskScore, aes(x = Risk_group, y = BCR.Shannon, fill = Risk_group)) +
  geom_jitter(mapping = aes(color=Risk_group),width = .05, alpha = 0.5,size=0.9)+ #绘制散点图
  geom_boxplot(position = position_nudge(x = 0.14),width=0.1,outlier.size = 0,outlier.alpha =0)+ #绘制箱线图，并通过position设置偏移
  stat_halfeye(mapping = aes(fill=Risk_group),width = 0.2, .width = 0, justification = -1.2, point_colour = NA,alpha=0.6) + #绘制云雨图，并通过position设置偏移
  scale_fill_manual(values = Custom.color)+   #映射云雨图和箱线图的颜色
  scale_color_manual(values = Custom.color)+  #映射散点的颜色
  #expand_limits(x = c(0.5, 3.8))+ #扩展画板，若显示不全，请根据你的数据范围手动调整或删除此行
  #ylim(20,110)+ #控制y轴显示范围，若显示不全，请根据你的数据范围手动调整或删除此行
  xlab("") +  #设置X轴标题
  #ylab("Risk Score") +   #设置Y轴标题
  #ggtitle("NSCLC-GSE126044")+  #设置主标题
  scale_y_continuous(
    name = "BCR.Shannon",                     # 设置 Y 轴标题
    limits = c(0, 6),                    # 设置 Y 轴显示范围
    breaks = seq(0, 6, by = 2)          # 按步长 10 显示刻度
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
  stat_compare_means(comparisons = comparison, method = "wilcox.test", label = "p.format", label.y = 5.5)
p9
ggsave(filename = '32.TCGA/TCGA-LUAD_klaSigScore_BCRshannon_boxplot.pdf', p9, height = 3, width = 2.5, dpi = 300)

luad_sub <- na.omit(luad[,c("Sample_ID", "BCR.Richness")])
luad_sub_riskScore <- inner_join(luad_sub, TCGA_LUAD_riskScore, by = "Sample_ID")

p10 <- ggplot(luad_sub_riskScore, aes(x = Risk_group, y = BCR.Richness, fill = Risk_group)) +
  geom_jitter(mapping = aes(color=Risk_group),width = .05, alpha = 0.5,size=0.9)+ #绘制散点图
  geom_boxplot(position = position_nudge(x = 0.14),width=0.1,outlier.size = 0,outlier.alpha =0)+ #绘制箱线图，并通过position设置偏移
  stat_halfeye(mapping = aes(fill=Risk_group),width = 0.2, .width = 0, justification = -1.2, point_colour = NA,alpha=0.6) + #绘制云雨图，并通过position设置偏移
  scale_fill_manual(values = Custom.color)+   #映射云雨图和箱线图的颜色
  scale_color_manual(values = Custom.color)+  #映射散点的颜色
  #expand_limits(x = c(0.5, 3.8))+ #扩展画板，若显示不全，请根据你的数据范围手动调整或删除此行
  #ylim(20,110)+ #控制y轴显示范围，若显示不全，请根据你的数据范围手动调整或删除此行
  xlab("") +  #设置X轴标题
  #ylab("Risk Score") +   #设置Y轴标题
  #ggtitle("NSCLC-GSE126044")+  #设置主标题
  # scale_y_continuous(
  #   name = "TCR.Shannon",                     # 设置 Y 轴标题
  #   limits = c(0, 6),                    # 设置 Y 轴显示范围
  #   breaks = seq(0, 6, by = 2)          # 按步长 10 显示刻度
  # ) +
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
  stat_compare_means(comparisons = comparison, method = "wilcox.test", label = "p.format", label.y = 600)
p10
ggsave(filename = '32.TCGA/TCGA-LUAD_klaSigScore_BCRrichness_boxplot.pdf', p10, height = 3, width = 2.5, dpi = 300)

luad_sub <- na.omit(luad[,c("Sample_ID", "BCR.Evenness")])
luad_sub_riskScore <- inner_join(luad_sub, TCGA_LUAD_riskScore, by = "Sample_ID")

p11 <- ggplot(luad_sub_riskScore, aes(x = Risk_group, y = BCR.Evenness, fill = Risk_group)) +
  geom_jitter(mapping = aes(color=Risk_group),width = .05, alpha = 0.5,size=0.9)+ #绘制散点图
  geom_boxplot(position = position_nudge(x = 0.14),width=0.1,outlier.size = 0,outlier.alpha =0)+ #绘制箱线图，并通过position设置偏移
  stat_halfeye(mapping = aes(fill=Risk_group),width = 0.2, .width = 0, justification = -1.2, point_colour = NA,alpha=0.6) + #绘制云雨图，并通过position设置偏移
  scale_fill_manual(values = Custom.color)+   #映射云雨图和箱线图的颜色
  scale_color_manual(values = Custom.color)+  #映射散点的颜色
  #expand_limits(x = c(0.5, 3.8))+ #扩展画板，若显示不全，请根据你的数据范围手动调整或删除此行
  #ylim(20,110)+ #控制y轴显示范围，若显示不全，请根据你的数据范围手动调整或删除此行
  xlab("") +  #设置X轴标题
  #ylab("Risk Score") +   #设置Y轴标题
  #ggtitle("NSCLC-GSE126044")+  #设置主标题
  scale_y_continuous(
    name = "BCR.Evenness",                     # 设置 Y 轴标题
    limits = c(0.2, 1.2),                    # 设置 Y 轴显示范围
    breaks = seq(0.2, 1.2, by = 0.2)          # 按步长 10 显示刻度
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
  stat_compare_means(comparisons = comparison, method = "wilcox.test", label = "p.format", label.y = 1.1)
p11
ggsave(filename = '32.TCGA/TCGA-LUAD_klaSigScore_BCRevenness_boxplot.pdf', p11, height = 3, width = 2.5, dpi = 300)

## 1.7 HRD----
luad_sub <- na.omit(luad[,c("Sample_ID", "Homologous.Recombination.Defects")])
luad_sub_riskScore <- inner_join(luad_sub, TCGA_LUAD_riskScore, by = "Sample_ID")

p12 <- ggplot(luad_sub_riskScore, aes(x = Risk_group, y = Homologous.Recombination.Defects, fill = Risk_group)) +
  geom_jitter(mapping = aes(color=Risk_group),width = .05, alpha = 0.5,size=0.9)+ #绘制散点图
  geom_boxplot(position = position_nudge(x = 0.14),width=0.1,outlier.size = 0,outlier.alpha =0)+ #绘制箱线图，并通过position设置偏移
  stat_halfeye(mapping = aes(fill=Risk_group),width = 0.2, .width = 0, justification = -1.2, point_colour = NA,alpha=0.6) + #绘制云雨图，并通过position设置偏移
  scale_fill_manual(values = Custom.color)+   #映射云雨图和箱线图的颜色
  scale_color_manual(values = Custom.color)+  #映射散点的颜色
  #expand_limits(x = c(0.5, 3.8))+ #扩展画板，若显示不全，请根据你的数据范围手动调整或删除此行
  #ylim(20,110)+ #控制y轴显示范围，若显示不全，请根据你的数据范围手动调整或删除此行
  xlab("") +  #设置X轴标题
  #ylab("Risk Score") +   #设置Y轴标题
  #ggtitle("NSCLC-GSE126044")+  #设置主标题
  scale_y_continuous(
    name = "HRD",                     # 设置 Y 轴标题
    limits = c(0, 80),                    # 设置 Y 轴显示范围
    breaks = seq(0, 80, by = 20)          # 按步长 10 显示刻度
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
  stat_compare_means(comparisons = comparison, method = "wilcox.test", label = "p.format", label.y = 73)
p12
ggsave(filename = '32.TCGA/TCGA-LUAD_klaSigScore_HRD_boxplot.pdf', p12, height = 3, width = 2.5, dpi = 300)


# 2.免疫浸润----
rm(luad_sub, luad_sub_riskScore, luad, immunity2018)

## 2.1 加载表达谱----
load("8.BdyDiffGenes/TCGA-LUAD/TCGA-LUAD_mRNA_tpm.RData") # 28脚本1.2
#load("8.BdyDiffGenes/TCGA-LUAD/TCGA-LUAD_lncRNA_tpm.RData")
which(duplicated(mRNA_tpm_annoted$mRNA_symbol))
luad_tpm <- aggregate(. ~ mRNA_symbol, mRNA_tpm_annoted, mean)
which(duplicated(luad_tpm$mRNA_symbol))
luad_tpm <- luad_tpm %>% 
  column_to_rownames("mRNA_symbol")
colnames(luad_tpm) <- str_sub(colnames(luad_tpm), 1, 16)

a <- str_split_i(colnames(luad_tpm), "-", 4)
table(a)
b <- which(str_sub(colnames(luad_tpm), 14, 14) == 0)

luad_tumor_tpm <- luad_tpm[,b]
log_luad_tumor_tpm <- log2(luad_tumor_tpm + 1)

## 2.2 CIBERSORT----
cibersort <- deconvo_tme(eset = log_luad_tumor_tpm, method = "cibersort", arrays = F, perm = 1000)
qsave(cibersort, file = "32.TCGA/TCGA-LUAD_CIBERSORT_results.qs")

colnames(cibersort)[1] <- "Sample_ID"
cibersort$Sample_ID <- str_sub(cibersort$Sample_ID, 1, 12)
cibersort <- cibersort[,c(1:23)]

dat <- cibersort %>%
  as.data.frame() %>%
  gather(key = Cell_type, value = Proportion, -Sample_ID)
dat$Cell_type <- gsub("_CIBERSORT$", "", dat$Cell_type)

high_risk_samples <- TCGA_LUAD_riskScore[TCGA_LUAD_riskScore$Risk_group == "High-Risk",]$Sample_ID
low_risk_samples <- TCGA_LUAD_riskScore[TCGA_LUAD_riskScore$Risk_group == "Low-Risk",]$Sample_ID

dat$Group <- ifelse(dat$Sample_ID %in% high_risk_samples, 'High-Risk', "Low-Risk")

high_color <- "#ae3f51"
low_color <- "#79aec1"
p <- ggplot(dat, aes(x = Cell_type, y = Proportion, fill = Group)) +
  geom_boxplot(outlier.shape = 21, color = "black", alpha = .6) +
  labs(x = "", y = "Estimated Proportion") +
  theme_bw(base_size = 14, base_line_size = 1, base_rect_size = 1)+
  theme(legend.position = "right") +
  theme(axis.text.x = element_text(angle = 45, hjust= 1, vjust = 1, size = rel(1))) +
  theme(axis.text.y = element_text(size = rel(1.1))) +
  scale_fill_manual(values = c(high_color, low_color)) +
  stat_compare_means(aes(group = Group, label = ..p.signif..), method = "wilcox.test")
p
ggsave(p, filename = "32.TCGA/CIBERSORT.pdf", width = 10, height = 6, dpi = 600)

## 2.3 MCPcounter----
mcp <- deconvo_tme(eset = log_luad_tumor_tpm, method = "mcpcounter")
qsave(mcp, file = "32.TCGA/TCGA-LUAD_MCPcounter_results.qs")

colnames(mcp)[1] <- "Sample_ID"
mcp$Sample_ID <- str_sub(mcp$Sample_ID, 1, 12)

dat <- mcp %>%
  as.data.frame() %>%
  gather(key = Cell_type, value = Proportion, -Sample_ID)
dat$Cell_type <- gsub("_MCPcounter$", "", dat$Cell_type)

high_risk_samples <- TCGA_LUAD_riskScore[TCGA_LUAD_riskScore$Risk_group == "High-Risk",]$Sample_ID
low_risk_samples <- TCGA_LUAD_riskScore[TCGA_LUAD_riskScore$Risk_group == "Low-Risk",]$Sample_ID

dat$Group <- ifelse(dat$Sample_ID %in% high_risk_samples, 'High-Risk', "Low-Risk")

high_color <- "#ae3f51"
low_color <- "#79aec1"
p <- ggplot(dat, aes(x = Cell_type, y = Proportion, fill = Group)) +
  geom_boxplot(outlier.shape = 21, color = "black", alpha = .6) +
  labs(x = "", y = "Estimated Proportion") +
  theme_bw(base_size = 14, base_line_size = 1, base_rect_size = 1)+
  theme(legend.position = "right") +
  theme(axis.text.x = element_text(angle = 45, hjust= 1, vjust = 1, size = rel(1.1))) +
  theme(axis.text.y = element_text(size = rel(1.1))) +
  scale_fill_manual(values = c(high_color, low_color)) +
  stat_compare_means(aes(group = Group, label = ..p.signif..), method = "wilcox.test")
p
ggsave(p, filename = "32.TCGA/MCPcounter.pdf", width = 8, height = 6, dpi = 600)

## 2.4 EPIC----
epic <- deconvo_tme(eset = log_luad_tumor_tpm, method = "epic", arrays = F, tumor = T, scale_mrna = F)
qsave(epic, file = "32.TCGA/TCGA-LUAD_EPIC_results.qs")

colnames(epic)[1] <- "Sample_ID"
epic$Sample_ID <- str_sub(epic$Sample_ID, 1, 12)

dat <- epic %>%
  as.data.frame() %>%
  gather(key = Cell_type, value = Proportion, -Sample_ID)
dat$Cell_type <- gsub("_EPIC$", "", dat$Cell_type)

high_risk_samples <- TCGA_LUAD_riskScore[TCGA_LUAD_riskScore$Risk_group == "High-Risk",]$Sample_ID
low_risk_samples <- TCGA_LUAD_riskScore[TCGA_LUAD_riskScore$Risk_group == "Low-Risk",]$Sample_ID

dat$Group <- ifelse(dat$Sample_ID %in% high_risk_samples, 'High-Risk', "Low-Risk")

high_color <- "#ae3f51"
low_color <- "#79aec1"
p <- ggplot(dat, aes(x = Cell_type, y = Proportion, fill = Group)) +
  geom_boxplot(outlier.shape = 21, color = "black", alpha = .6) +
  labs(x = "", y = "Estimated Proportion") +
  theme_bw(base_size = 14, base_line_size = 1, base_rect_size = 1)+
  theme(legend.position = "right") +
  theme(axis.text.x = element_text(angle = 45, hjust= 1, vjust = 1, size = rel(1.1))) +
  theme(axis.text.y = element_text(size = rel(1.1))) +
  scale_fill_manual(values = c(high_color, low_color)) +
  stat_compare_means(aes(group = Group, label = ..p.signif..), method = "wilcox.test")
p
ggsave(p, filename = "32.TCGA/EPIC.pdf", width = 6, height = 6, dpi = 600)

## 2.5 TIMER----
timer <- deconvo_tme(eset = log_luad_tumor_tpm, method = "timer", group_list = rep("LUAD", dim(log_luad_tumor_tpm)[2]))
qsave(timer, file = "32.TCGA/TCGA-LUAD_TIMER_results.qs")

colnames(timer)[1] <- "Sample_ID"
timer$Sample_ID <- str_sub(timer$Sample_ID, 1, 12)

dat <- timer %>%
  as.data.frame() %>%
  gather(key = Cell_type, value = Proportion, -Sample_ID)
dat$Cell_type <- gsub("_TIMER$", "", dat$Cell_type)

high_risk_samples <- TCGA_LUAD_riskScore[TCGA_LUAD_riskScore$Risk_group == "High-Risk",]$Sample_ID
low_risk_samples <- TCGA_LUAD_riskScore[TCGA_LUAD_riskScore$Risk_group == "Low-Risk",]$Sample_ID

dat$Group <- ifelse(dat$Sample_ID %in% high_risk_samples, 'High-Risk', "Low-Risk")

high_color <- "#ae3f51"
low_color <- "#79aec1"
p <- ggplot(dat, aes(x = Cell_type, y = Proportion, fill = Group)) +
  geom_boxplot(outlier.shape = 21, color = "black", alpha = .6) +
  labs(x = "", y = "Estimated Proportion") +
  theme_bw(base_size = 14, base_line_size = 1, base_rect_size = 1)+
  theme(legend.position = "right") +
  theme(axis.text.x = element_text(angle = 45, hjust= 1, vjust = 1, size = rel(1.1))) +
  theme(axis.text.y = element_text(size = rel(1.1))) +
  scale_fill_manual(values = c(high_color, low_color)) +
  stat_compare_means(aes(group = Group, label = ..p.signif..), method = "wilcox.test")
p
ggsave(p, filename = "32.TCGA/TIMER.pdf", width = 6, height = 4.5, dpi = 600)

## 2.6 quanTIseq----
quantiseq <- deconvo_tme(eset = log_luad_tumor_tpm, tumor = TRUE, arrays = F, scale_mrna = TRUE, method = "quantiseq")
qsave(quantiseq, file = "32.TCGA/TCGA-LUAD_quanTIseq_results.qs")

colnames(quantiseq)[1] <- "Sample_ID"
quantiseq$Sample_ID <- str_sub(quantiseq$Sample_ID, 1, 12)

dat <- quantiseq %>%
  as.data.frame() %>%
  gather(key = Cell_type, value = Proportion, -Sample_ID)
dat$Cell_type <- gsub("_quantiseq$", "", dat$Cell_type)

high_risk_samples <- TCGA_LUAD_riskScore[TCGA_LUAD_riskScore$Risk_group == "High-Risk",]$Sample_ID
low_risk_samples <- TCGA_LUAD_riskScore[TCGA_LUAD_riskScore$Risk_group == "Low-Risk",]$Sample_ID

dat$Group <- ifelse(dat$Sample_ID %in% high_risk_samples, 'High-Risk', "Low-Risk")

high_color <- "#ae3f51"
low_color <- "#79aec1"
p <- ggplot(dat, aes(x = Cell_type, y = Proportion, fill = Group)) +
  geom_boxplot(outlier.shape = 21, color = "black", alpha = .6) +
  labs(x = "", y = "Estimated Proportion") +
  theme_bw(base_size = 14, base_line_size = 1, base_rect_size = 1)+
  theme(legend.position = "right") +
  theme(axis.text.x = element_text(angle = 45, hjust= 1, vjust = 1, size = rel(1.1))) +
  theme(axis.text.y = element_text(size = rel(1.1))) +
  scale_fill_manual(values = c(high_color, low_color)) +
  stat_compare_means(aes(group = Group, label = ..p.signif..), method = "wilcox.test")
p
ggsave(p, filename = "32.TCGA/quanTIseq.pdf", width = 8, height = 6, dpi = 600)
