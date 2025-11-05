######目的：鉴定LUAD高kla.sig评分的药物
######作者：申奥
######日期：2024-11-30
######Server: yjx@xiyoucloud.top:34022
######Version: R 4.4.1


# 0.环境设置----
rm(list = ls())
options(stringAsFactors = F)

library(oncoPredict)
library(readr)
library(stringr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(ggsignif)
library(ggdist)


# 1.加载药物表达谱及响应数据----
### 下载地址：https://osf.io/c6tfx/files/osfstorage
GDSC2_Expr <- readRDS(file="18.drug_targets/oncoPredict/DataFiles/DataFiles/Training Data/GDSC2_Expr (RMA Normalized and Log Transformed).rds")
GDSC2_Res <- readRDS(file = "18.drug_targets/oncoPredict/DataFiles/DataFiles/Training Data/GDSC2_Res.rds")
GDSC2_Res = exp(GDSC2_Res)

dim(GDSC2_Expr) # 行基因，列细胞系
## [1] 17419   805
dim(GDSC2_Res) # 行细胞系，列药物
## [1] 805 198


# 2.加载TCGA-LUAD表达谱数据----
load("12.Prognosis_model/merge_datasets.RData") # 51脚本1.1.1
expr_df <- merge_datasets[["TCGA_LUAD"]] 
expr_df <- expr_df[,-c(2,3)]
which(duplicated(expr_df$Sample_ID))
expr_df <- expr_df %>% 
  tibble::column_to_rownames("Sample_ID")
expr_df <- as.data.frame(t(expr_df))


# 3.预测----
calcPhenotype(trainingExprData = GDSC2_Expr, #训练数据集中的基因表达数据
              trainingPtype = GDSC2_Res, #训练数据集中的药物敏感性数据
              testExprData = as.matrix(expr_df),  #测试数据集中的基因表达数据
              batchCorrect = 'standardize',  # 进行批次效应校正, 根据自己的数据类型不同进行选择。"eb" for array, "standardize" for RNA-seq
              powerTransformPhenotype = TRUE, #这个参数决定是否对药物敏感性数据进行幂变换
              removeLowVaryingGenes = 0.2, #这个参数指定了移除低变异基因的阈值
              minNumSamples = 10,  #这个参数指定了每个组（例如，每种药物敏感性）的最小样本数量
              printOutput = TRUE, #这个参数决定是否在控制台打印输出信息
              removeLowVaringGenesFrom = 'rawData' ) #这个参数指定了从哪个步骤开始移除低变异基因，有的教程用的是“homogenizeData”。
# removeLowVaryingGenesFrom 参数决定从哪个数据处理中阶段开始移除低变异基因。它有以下可选值：
# 'rawData'（默认值）：表示从原始数据（未标准化或处理的数据）中移除低变异基因。
# 'homogenizeData'：表示在数据经过批次校正或标准化（homogenize）之后再移除低变异基因。
#  removeLowVaringGenesFrom参数选择的影响
# 'rawData'：适用于希望在所有后续处理步骤之前就减少低变异基因数量的情况。适合在数据中含有大量低变异基因时使用，可以减少计算复杂度。
# 'homogenizeData'：适合先处理批次效应或其他数据标准化后再进行低变异基因过滤的情况。此设置可能更准确，因为标准化后的数据会更一致，减少噪音对低变异基因判定的影响。
# 使用建议
# 如果数据的批次效应较大或标准化前后差异显著，可以选择 'homogenizeData'，以确保变异计算更符合实际。
# 如果数据本身已经较为一致或对数据处理效率要求较高，可以使用默认的 'rawData'。
# 具体选择取决于数据特性和分析目标。如果不确定，可以先用默认值（'rawData'），观察结果是否合理。


# 4.分析结果----
## 4.1 读取预测结果及risk score分组结果，整合在一起----
pred_res <- read.csv("18.drug_targets/oncoPredict/calcPhenotype_Output/DrugPredictions.csv")
colnames(pred_res)[1] <- "Sample_ID"
which(duplicated(pred_res$Sample_ID))
pred_res2 <- log10(pred_res[,c(2:ncol(pred_res))])
pred_res3 <- cbind(Sample_ID = pred_res$Sample_ID, pred_res2)

load("12.Prognosis_model/RSF_riskScore.RData") # 51脚本1.1.1
identical(merge_datasets[["TCGA_LUAD"]]$OS, rs[["TCGA_LUAD"]]$OS)
identical(merge_datasets[["TCGA_LUAD"]]$OS.time, rs[["TCGA_LUAD"]]$OS.time)
TCGA_LUAD_riskScore <- data.frame(cbind(Sample_ID = merge_datasets[["TCGA_LUAD"]]$Sample_ID, rs[["TCGA_LUAD"]]))
TCGA_LUAD_riskScore$Risk_group <- ifelse(TCGA_LUAD_riskScore$RS > median(TCGA_LUAD_riskScore$RS), "High-Risk", "Low-Risk")
table(TCGA_LUAD_riskScore$Risk_group)
save(TCGA_LUAD_riskScore, file = "12.Prognosis_model/TCGA_LUAD_riskScoreGroup.RData")

merge_df <- inner_join(TCGA_LUAD_riskScore, pred_res3, by = "Sample_ID")
table(merge_df$Risk_group)
merge_df$Risk_group <- factor(merge_df$Risk_group, levels = c("Low-Risk", "High-Risk"))

## 4.2 查看几种常见化疗药物是否有差异----
### 4.2.1 Paclitaxel紫杉醇----
which(str_detect(colnames(merge_df), "Paclitaxel"))
colnames(merge_df)[41]
comparison <- list(c("Low-Risk", "High-Risk"))
#Custom.color <- c("#245892","#877eac","#a9d296")
Custom.color <- c("#79aec1", "#ae3f51")
p1 <- ggplot(merge_df, aes(x = Risk_group, y = Paclitaxel_1080, fill = Risk_group)) +
  geom_jitter(mapping = aes(color=Risk_group),width = .05, alpha = 0.5,size=0.9)+ #绘制散点图
  geom_boxplot(position = position_nudge(x = 0.14),width=0.1,outlier.size = 0,outlier.alpha =0)+ #绘制箱线图，并通过position设置偏移
  stat_halfeye(mapping = aes(fill=Risk_group),width = 0.2, .width = 0, justification = -1.2, point_colour = NA,alpha=0.6) + #绘制云雨图，并通过position设置偏移
  scale_fill_manual(values = Custom.color)+   #映射云雨图和箱线图的颜色
  scale_color_manual(values = Custom.color)+  #映射散点的颜色
  #expand_limits(x = c(0.5, 3.8))+ #扩展画板，若显示不全，请根据你的数据范围手动调整或删除此行
  #ylim(0,0.75)+ #控制y轴显示范围，若显示不全，请根据你的数据范围手动调整或删除此行
  xlab("") +  #设置X轴标题
  ylab(bquote("Log"[10]~"(IC50)")) +   #设置Y轴标题
  ggtitle("Paclitaxel")+  #设置主标题
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
  stat_compare_means(comparisons = comparison, method = "wilcox.test", label = "p.format", label.y = 0)
p1
ggsave(filename = '18.drug_targets/oncoPredict/TCGA-LUAD_Paclitaxel_IC50_boxplot.pdf', p1, height = 3, width = 2.5, dpi = 300)

### 4.2.2 Gefitinib----
which(str_detect(colnames(merge_df), "Gefitinib"))
colnames(merge_df)[11]
comparison <- list(c("Low-Risk", "High-Risk"))
#Custom.color <- c("#245892","#877eac","#a9d296")
Custom.color <- c("#79aec1", "#ae3f51")
p2 <- ggplot(merge_df, aes(x = Risk_group, y = Gefitinib_1010, fill = Risk_group)) +
  geom_jitter(mapping = aes(color=Risk_group),width = .05, alpha = 0.5,size=0.9)+ #绘制散点图
  geom_boxplot(position = position_nudge(x = 0.14),width=0.1,outlier.size = 0,outlier.alpha =0)+ #绘制箱线图，并通过position设置偏移
  stat_halfeye(mapping = aes(fill=Risk_group),width = 0.2, .width = 0, justification = -1.2, point_colour = NA,alpha=0.6) + #绘制云雨图，并通过position设置偏移
  scale_fill_manual(values = Custom.color)+   #映射云雨图和箱线图的颜色
  scale_color_manual(values = Custom.color)+  #映射散点的颜色
  #expand_limits(x = c(0.5, 3.8))+ #扩展画板，若显示不全，请根据你的数据范围手动调整或删除此行
  ylim(0.5,2)+ #控制y轴显示范围，若显示不全，请根据你的数据范围手动调整或删除此行
  xlab("") +  #设置X轴标题
  ylab(bquote("Log"[10]~"(IC50)")) +   #设置Y轴标题
  ggtitle("Gefitinib")+  #设置主标题
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
  stat_compare_means(comparisons = comparison, method = "wilcox.test", label = "p.format", label.y = 1.85)
p2
ggsave(filename = '18.drug_targets/oncoPredict/TCGA-LUAD_Gefitinib_IC50_boxplot.pdf', p2, height = 3, width = 2.5, dpi = 300)

### 4.2.3 Cisplatin----
which(str_detect(colnames(merge_df), "Cisplatin"))
colnames(merge_df)[8]
comparison <- list(c("Low-Risk", "High-Risk"))
#Custom.color <- c("#245892","#877eac","#a9d296")
Custom.color <- c("#79aec1", "#ae3f51")
p3 <- ggplot(merge_df, aes(x = Risk_group, y = Cisplatin_1005, fill = Risk_group)) +
  geom_jitter(mapping = aes(color=Risk_group),width = .05, alpha = 0.5,size=0.9)+ #绘制散点图
  geom_boxplot(position = position_nudge(x = 0.14),width=0.1,outlier.size = 0,outlier.alpha =0)+ #绘制箱线图，并通过position设置偏移
  stat_halfeye(mapping = aes(fill=Risk_group),width = 0.2, .width = 0, justification = -1.2, point_colour = NA,alpha=0.6) + #绘制云雨图，并通过position设置偏移
  scale_fill_manual(values = Custom.color)+   #映射云雨图和箱线图的颜色
  scale_color_manual(values = Custom.color)+  #映射散点的颜色
  #expand_limits(x = c(0.5, 3.8))+ #扩展画板，若显示不全，请根据你的数据范围手动调整或删除此行
  ylim(0,3)+ #控制y轴显示范围，若显示不全，请根据你的数据范围手动调整或删除此行
  xlab("") +  #设置X轴标题
  ylab(bquote("Log"[10]~"(IC50)")) +   #设置Y轴标题
  ggtitle("Cisplatin")+  #设置主标题
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
  stat_compare_means(comparisons = comparison, method = "wilcox.test", label = "p.format", label.y = 2.7)
p3
ggsave(filename = '18.drug_targets/oncoPredict/TCGA-LUAD_Cisplatin_IC50_boxplot.pdf', p3, height = 3, width = 2.5, dpi = 300)

### 4.2.4 Gemcitabine----
which(str_detect(colnames(merge_df), "Gemcitabine"))
colnames(merge_df)[57]
comparison <- list(c("Low-Risk", "High-Risk"))
#Custom.color <- c("#245892","#877eac","#a9d296")
Custom.color <- c("#79aec1", "#ae3f51")
p4 <- ggplot(merge_df, aes(x = Risk_group, y = Gemcitabine_1190, fill = Risk_group)) +
  geom_jitter(mapping = aes(color=Risk_group),width = .05, alpha = 0.5,size=0.9)+ #绘制散点图
  geom_boxplot(position = position_nudge(x = 0.14),width=0.1,outlier.size = 0,outlier.alpha =0)+ #绘制箱线图，并通过position设置偏移
  stat_halfeye(mapping = aes(fill=Risk_group),width = 0.2, .width = 0, justification = -1.2, point_colour = NA,alpha=0.6) + #绘制云雨图，并通过position设置偏移
  scale_fill_manual(values = Custom.color)+   #映射云雨图和箱线图的颜色
  scale_color_manual(values = Custom.color)+  #映射散点的颜色
  #expand_limits(x = c(0.5, 3.8))+ #扩展画板，若显示不全，请根据你的数据范围手动调整或删除此行
  ylim(-2,2)+ #控制y轴显示范围，若显示不全，请根据你的数据范围手动调整或删除此行
  xlab("") +  #设置X轴标题
  ylab(bquote("Log"[10]~"(IC50)")) +   #设置Y轴标题
  ggtitle("Gemcitabine")+  #设置主标题
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
  stat_compare_means(comparisons = comparison, method = "wilcox.test", label = "p.format", label.y = 1.5)
p4
ggsave(filename = '18.drug_targets/oncoPredict/TCGA-LUAD_Gemcitabine_IC50_boxplot.pdf', p4, height = 3, width = 2.5, dpi = 300)

## 4.3 批量分析所有药物----
### 4.3.1 分析药物IC50与risk score的相关性----
cor_results <- data.frame() # 初始化结果数据框

for (drug in colnames(merge_df)[6:ncol(merge_df)]) {
  print(drug)
  test <- cor.test(merge_df$RS, merge_df[[drug]], method = "spearman") # 计算相关性和 p 值
  cor_results <- rbind(cor_results, data.frame(
    Drugs = drug,              # targets_expr 的列名
    Correlation = test$estimate, # 相关性 R
    P_Value = test$p.value     # P 值
  ))
}
cor_results$Adjusted_P <- p.adjust(cor_results$P_Value, method = "BH")
cor_results <- cor_results[order(cor_results$Correlation), ]
save(cor_results, file = "18.drug_targets/TCGA-LUAD_drugIC50_klsSigScore_cor_results.RData")

### 4.3.2 相关性图----
darkblue <- "#0772B9"
lightblue <- "#48C8EF"

sig_cor <- cor_results %>% 
  filter(Correlation < -0.25 & Adjusted_P < 0.05)

p5 <- ggplot(data = sig_cor, aes(Correlation, forcats::fct_reorder(Drugs, Correlation, .desc = T))) +
  geom_segment(aes(xend = 0, yend = Drugs), linetype = 2) +
  geom_point(aes(size = -log10(Adjusted_P)), col = darkblue) +
  scale_size_continuous(range = c(2,8)) +
  scale_x_reverse(breaks = c(0, -0.1, -0.2,-0.3),
                  expand = expansion(mult = c(0.01,.1))) + #左右留空
  theme_classic() +
  labs(x = "Spearman correlation coefficient", y = "", size = bquote("-log"[10]~"("~italic(P)~"-value )")) + 
  theme(legend.position = "bottom", 
        axis.line.y = element_blank(),
        axis.text = element_text(size = 12))
p5
ggsave(filename = '18.drug_targets/oncoPredict/TCGA-LUAD_GDSC2_drugIC50_riskScore_sigCorplot.pdf', p5, height = 4, width = 5, dpi = 300)

### 4.3.3 箱式图----
merge_df2 <- merge_df[,c("Risk_group", sig_cor$Drugs)]

long_df <- merge_df2 %>%
  pivot_longer(-Risk_group, names_to = "Drug", values_to = "IC50")

# 绘制箱式图和散点图，并添加 P 值
custom_colors <- c("High-Risk" = "#ae3f51", 
                   "Low-Risk" = "#79aec1") # 根据需求调整颜色

p6 <- ggplot(long_df, aes(x = Drug, y = IC50, fill = Risk_group)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.8), notch = T, notchwidth = 0.5) + # 箱式图
  #geom_jitter(aes(color = Risk_group), position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), alpha = 0.6) + # 散点图
  stat_compare_means(aes(group = Risk_group), method = "wilcox.test", label = "p.signif", hide.ns = TRUE) + # 添加 P 值
  labs(x = "", y = bquote("Log"[10]~"(IC50)"), fill = "Risk Group", color = "Risk Group") +
  scale_fill_manual(values = custom_colors) + # 自定义填充颜色
  theme_classic2() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 12),
        legend.position = "top", # 图例位置：顶部
        legend.justification = "center") # 居中对齐)
p6
ggsave(filename = '18.drug_targets/oncoPredict/TCGA-LUAD_GDSC2_drugIC50_boxplot.pdf', p6, height = 4, width = 5, dpi = 300)
