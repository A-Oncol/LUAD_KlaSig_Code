######目的：鉴定LUAD高kla.sig评分的可靶向基因
######作者：申奥
######日期：2024-11-30
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
library(biomaRt)
library(openxlsx)


# 1.CPTAC-LUAD队列----
## 1.1 蛋白表达谱整理----
##CPTAC-LUAD数据下载：https://www.linkedomics.org/data_download/CPTAC-pancan-LUAD/
#### 使用wget https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/LUAD/LUAD_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt下载肿瘤样本蛋白组
cptac_pro <- read.table("1.data/CPTAC/LUAD_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt",
                        header = T, sep = "\t")

ensembl_ids <- cptac_pro$idx
ensembl_ids_clean <- sub("\\..*$", "", ensembl_ids)
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
result <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),  # 返回 Ensembl ID 和 Gene Symbol
  filters = "ensembl_gene_id",                       # 输入过滤条件
  values = ensembl_ids_clean,                              # Ensembl ID 列表
  mart = mart                                        # 数据库连接
)
colnames(result)[1] <- "idx"
cptac_pro$idx <- ensembl_ids_clean

pro_expr <- inner_join(result, cptac_pro, by = "idx")
pro_expr <- pro_expr[,-1]
which(duplicated(pro_expr$hgnc_symbol))
pro_expr <- aggregate(. ~ hgnc_symbol, pro_expr, mean)
which(duplicated(pro_expr$hgnc_symbol))
pro_expr <- pro_expr %>% 
  tibble::column_to_rownames("hgnc_symbol")
pro_expr <- as.data.frame(t(pro_expr))

## 1.2 计算kla.sig评分----
load("12.Prognosis_model/RSFfit.RData")
load("12.Prognosis_model/same_genes.RData")
load("15.diff/upKla_genes.RData")
build_model_genes <- intersect(upKla_genes, same_genes) # 318个
write.table(build_model_genes, file = "18.drug_targets/build_model_genes.txt", quote = F, sep = "\t", row.names = F)

expr_df <- pro_expr[,build_model_genes]  # 大多数模型基因不在蛋白谱中，所以不用CPTAC了

## 1.3 看一下Xena下载的CPTAC数据----
# cptac <- read_tsv("18.drug_targets/CPTAC/CPTAC-3.star_tpm.tsv.gz")
# cptac_clinic <- read_tsv("18.drug_targets/CPTAC/CPTAC-3.clinical.tsv.gz")
# table(cptac_clinic$disease_type)
# table(cptac_clinic$primary_diagnosis.diagnoses)
# table(cptac_clinic$primary_site)
# table(cptac_clinic$primary_site.project)
# 
# luad_cptac <- cptac_clinic %>% 
#   filter(primary_site == "Bronchus and lung" & disease_type == "Adenomas and Adenocarcinomas")


# 2.TCGA-LUAD队列----
## 2.1 加载TCGA-LUAD队列表达谱、药物靶点和kla.sig评分----
load("12.Prognosis_model/RSF_riskScore.RData") # 51脚本1.1.1
load("12.Prognosis_model/merge_datasets.RData")
riskScore <- rs[["TCGA_LUAD"]]

expr_df <- merge_datasets[["TCGA_LUAD"]]
identical(riskScore$OS.time, expr_df$OS.time)
# riskScore$Sample_ID <- expr_df$Sample_ID
# write.xlsx(riskScore, file = "18.drug_targets/TCGA-LUAD_riskScore.xlsx") # 只是最后整理文章表格时运行了这两行代码

expr_df <- expr_df[,-c(2,3)]
which(duplicated(expr_df$Sample_ID))
expr_df <- expr_df %>% 
  tibble::column_to_rownames("Sample_ID")

####药物靶点从https://academic.oup.com/bib/article/22/3/bbaa164/5891146?login=false#248046063的Table S10下载
drug_targets <- read.xlsx("18.drug_targets/supplementary_table_bbaa164.xlsx", sheet = 10, startRow = 2)
targets <- unique(drug_targets$Target.genes) # 2249个
drugs <- unique(drug_targets$Agent.names) # 4484个

same_targets <- intersect(colnames(expr_df), targets) # 2196
targets_expr <- expr_df[,same_targets]

## 2.2 计算评分和靶点表达的相关性----
cor_results <- data.frame() # 初始化结果数据框

for (gene in colnames(targets_expr)) {
  print(gene)
  test <- cor.test(riskScore$RS, targets_expr[[gene]], method = "spearman") # 计算相关性和 p 值
  cor_results <- rbind(cor_results, data.frame(
    Target = gene,              # targets_expr 的列名
    Correlation = test$estimate, # 相关性 R
    P_Value = test$p.value     # P 值
  ))
}
cor_results$Adjusted_P <- p.adjust(cor_results$P_Value, method = "BH")
#cor_results$Adjusted_P <- ifelse(cor_results$Adjusted_P == 0, 1e-300, cor_results$Adjusted_P)
save(cor_results, file = "18.drug_targets/TCGA-LUAD_riskScore_drugTargetsExp_cor_results.RData")
write.xlsx(cor_results, file = "18.drug_targets/TCGA-LUAD_riskScore_drugTargetsExp_cor_results.xlsx")

## 2.3 筛选并画图----
cor_results$Significant <- ifelse(cor_results$Correlation > 0.3 & cor_results$Adjusted_P < 0.05, "Significant", "Not Significant")
table(cor_results$Significant)

significant_targets <- cor_results[cor_results$Significant == "Significant",]$Target # 171个
write.table(significant_targets, file = "18.drug_targets/TCGA-LUAD_significant_targets.txt", quote = F, sep = "\t", row.names = F)

# 计算 -log10(Pvalue)
cor_results$logP <- -log10(cor_results$Adjusted_P)
# 绘制图形
p1 <- ggplot(cor_results, aes(x = Correlation, y = logP, color = Significant)) +
  geom_point(size = 1) +  # 散点图
  scale_color_manual(values = c("Significant" = "#ae3f51", "Not Significant" = "grey")) + # 颜色
  geom_vline(xintercept = 0.3, linetype = "dashed", color = "gray") + # 虚线：相关性阈值
  #geom_vline(xintercept = -0.3, linetype = "dashed", color = "black") + # 虚线：负相关性阈值
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray") + # 虚线：P值阈值
  labs(
    #title = "Correlation vs -log10(P-value)",
    x = "Spearman Correlation Coefficient",
    y = "-log10(BH adjusted P-value)",
    color = "Significance"
  ) +
  theme_classic2() +
  theme(axis.text = element_text(size = 12))
p1
ggsave(filename = "18.drug_targets/TCGA-LUAD_riskScore_drugTargetsExp_cor_plot.pdf", p1,
       height = 4, width = 5, dpi = 300)


# # 3.Xena数据库APOLLO-LUAD队列----
# ## 3.1 整理表达谱并计算risk score----
# apollo_tpm <- read_tsv("18.drug_targets/Xena_AppoloLUAD/APOLLO-LUAD.star_tpm.tsv.gz")
# annot <- read.table("18.drug_targets/Xena_AppoloLUAD/gencode.v36.annotation.gtf.gene.probemap", header = T)
# annot <- annot[,c(1,2)]
# colnames(annot) <- c("Ensembl_ID", "Symbol")
# 
# apollo_tpm_annot <- inner_join(annot, apollo_tpm, by = "Ensembl_ID")
# apollo_tpm_annot <- apollo_tpm_annot[,-1]
# which(duplicated(apollo_tpm_annot$Symbol))
# apollo_tpm_annot <- aggregate(. ~ Symbol, apollo_tpm_annot, mean)
# which(duplicated(apollo_tpm_annot$Symbol))
# apollo_tpm_annot <- apollo_tpm_annot %>% 
#   tibble::column_to_rownames("Symbol")
# 
# apollo_tpm_annot <- as.data.frame(t(apollo_tpm_annot))
# 
# expr_df <- apollo_tpm_annot[,build_model_genes]
# expr_df <- as.data.frame(scale(expr_df))
# riskScore <- predict(fit, newdata = expr_df, type = "risk")$predicted
# rs <- data.frame(Sample_ID = rownames(expr_df),
#                  RiskScore = riskScore)
# rs$RiskGroup <- ifelse(rs$RiskScore > median(rs$RiskScore), "High-Risk", "Low-Risk")
# table(rs$RiskGroup)
# rs$RiskGroup <- factor(rs$RiskGroup, levels = c("Low-Risk", "High-Risk"))
# 
# ## 3.2 计算评分和靶点表达的相关性----
# same_targets2 <- intersect(colnames(expr_df), targets) # 84
# targets_expr <- expr_df[,same_targets2]
# 
# cor_results2 <- data.frame() # 初始化结果数据框
# 
# for (gene in colnames(targets_expr)) {
#   print(gene)
#   test <- cor.test(rs$RiskScore, targets_expr[[gene]], method = "spearman") # 计算相关性和 p 值
#   cor_results2 <- rbind(cor_results2, data.frame(
#     Target = gene,              # targets_expr 的列名
#     Correlation = test$estimate, # 相关性 R
#     P_Value = test$p.value     # P 值
#   ))
# }
# cor_results2$Adjusted_P <- p.adjust(cor_results2$P_Value, method = "BH")
# cor_results2$Adjusted_P <- ifelse(cor_results2$Adjusted_P == 0, 1e-300, cor_results$Adjusted_P)
# save(cor_results2, file = "18.drug_targets/Xena-Apollo_riskScore_drugTargetsExp_cor_results.RData")


# 3.CCLE数据----
## 3.1 LUAD细胞系的表达数据----
# cell_expr <- read.csv("18.drug_targets/DepMap/Expression_Public_24Q2_subsetted.csv")
# cell_expr2 <- read.csv("18.drug_targets/DepMap/Batch_corrected_Expression_Public_24Q2_subsetted.csv")
# cell_pro <- read.csv("18.drug_targets/DepMap/Protein_Array_subsetted.csv")
# cell_pro2 <- read.csv("18.drug_targets/DepMap/Proteomics_subsetted_NAsdropped.csv")

## 3.2 CERES数据----
cell_ceres1 <- read.csv("18.drug_targets/DepMap/CRISPR_(DepMap_Public_24Q2+Score,_Chronos)_subsetted.csv")
# cell_ceres2 <- read.csv("18.drug_targets/DepMap/CRISPR_(Project_Score,_CERES)_subsetted.csv")
# cell_ceres3 <- read.csv("18.drug_targets/DepMap/CRISPR_(Project_Score,_Chronos)_subsetted.csv")
# cell_ceres4 <- read.csv("18.drug_targets/DepMap/CRISPR_GeCKO_(DepMap,_CERES)_subsetted.csv")

significant_targets %in% colnames(cell_ceres1)
#cell_ceres_sub <- cell_ceres1[,c("cell_line_display_name",significant_targets)]
cell_ceres_sub <- cell_ceres1[,c("cell_line_display_name", intersect(significant_targets, build_model_genes))]

# 将数据转换为长格式，方便 ggplot2 使用
df_long <- cell_ceres_sub %>%
  pivot_longer(
    cols = -cell_line_display_name,              # 转换第一列以外的列
    names_to = "Gene",         # 新列名为 Gene
    values_to = "CERES_Score"  # CERES 评分
  )

gene_order <- df_long %>%
  group_by(Gene) %>%
  summarise(Mean_CERES = mean(CERES_Score, na.rm = TRUE)) %>%
  arrange(Mean_CERES) %>%
  pull(Gene) # 提取排序后的基因名

# 设置 Gene 列为因子，按均值排序
df_long$Gene <- factor(df_long$Gene, levels = gene_order)

# 绘制图形
p2 <- ggplot(df_long, aes(x = Gene, y = CERES_Score)) +
  geom_boxplot(
    aes(color = Gene), # 箱式图边框按基因着色
    outlier.shape = NA, # 不显示离群点
    alpha = 0.6         # 半透明填充
  ) +
  geom_jitter(
    aes(color = Gene), # 散点按细胞系着色
    width = 0.2,       # 水平方向抖动
    size = 0.8,          # 散点大小
    alpha = 0.7        # 散点透明度
  ) +
  # scale_color_manual(
  #   values = c("Cell1" = "red", "Cell2" = "blue", "Cell3" = "green", 
  #              "Cell4" = "purple", "Cell5" = "orange", # 细胞系颜色
  #              "GeneA" = "darkred", "GeneB" = "darkblue", "GeneC" = "darkgreen") # 基因颜色
  # ) +
  geom_hline(
    yintercept = -1,     # 添加 CERES 评分为 -1 的虚线
    linetype = "dashed", # 虚线样式
    color = "grey"       # 灰色
  ) +
  labs(
    #title = "CERES Scores Sorted by Mean Across Genes",
    x = "Druggable targets", #  (sorted by mean CERES score)
    y = "CERES Score",
    color = "Legend"
  ) +
  theme_pubr() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10), # 横轴基因名倾斜 45°
    legend.position = "none"                          # 图例放在右侧
  )
p2
ggsave(filename = "18.drug_targets/CERES_druggableKlgSigTargets_plot.pdf", p2,
       height = 4, width = 6, dpi = 300)

## 3.3 韦恩图----
library(VennDiagram)
venn.diagram(x=list(significant_targets,build_model_genes),
             scaled = T, # 根据比例显示大小
             alpha= 0.9, #透明度
             lwd=1,lty=1,
             col=c("black", "black"), #圆圈线条粗细、形状、颜色；1 实线, 2 虚线, blank无线条
             label.col ='black' , # 数字颜色abel.col=c('#FFFFCC','#CCFFFF',......)根据不同颜色显示数值颜色
             cex = 2, # 数字大小
             fontface = "bold",  # 字体粗细；加粗bold
             fill=c('#336681','#ba4f4a'), # 填充色 配色https://www.58pic.com/
             category.names = c("Targets", "Kla.Sig") , #标签名
             cat.dist = 0.02, # 标签距离圆圈的远近
             cat.pos = -180, # 标签相对于圆圈的角度cat.pos = c(-10, 10, 135)
             cat.cex = 2, #标签字体大小
             cat.fontface = "bold",  # 标签字体加粗
             cat.col='black' ,   #cat.col=c('#FFFFCC','#CCFFFF',.....)根据相应颜色改变标签颜色
             cat.default.pos = "outer",  # 标签位置, outer内;text 外
             output=TRUE,
             filename='18.drug_targets/venn.tiff',# 文件保存
             imagetype="tiff",  # 类型（tiff png svg）
             resolution = 400,  # 分辨率
             compression = "lzw"# 压缩算法
             
)

final_29genes <- intersect(significant_targets, build_model_genes)
write.table(final_29genes, file = "18.drug_targets/final_29genes.txt", quote = F, sep = "\t", row.names = F)
