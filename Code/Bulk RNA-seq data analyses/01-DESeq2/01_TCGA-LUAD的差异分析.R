######目的：TCGA-LUAD的差异分析
######作者：申奥
######日期：2024-10-25
######Server: yjx@xiyoucloud.top:34022
######Version: R 4.3.3


# 0.环境设置----
rm(list = ls())
options(stringAsFactors = F)

library(DESeq2)
library(tidyverse)


# 1.加载并整理count数据----
load("~/projects/Rprojects/TCGA_data/data/RNA/TCGA-LUAD_mRNA_count.RData")
load("~/projects/Rprojects/TCGA_data/data/RNA/TCGA-LUAD_lncRNA_count.RData")
identical(colnames(lnc_count_annoted), colnames(mRNA_count_annoted))
colnames(lnc_count_annoted)[1] <- "Symbol"
colnames(mRNA_count_annoted)[1] <- "Symbol"
identical(colnames(lnc_count_annoted), colnames(mRNA_count_annoted))

count_exp <- rbind(mRNA_count_annoted, lnc_count_annoted)

which(duplicated(count_exp$Symbol))
count_exp <- aggregate(. ~ Symbol, count_exp, mean)
which(duplicated(count_exp$Symbol))
count_exp <- count_exp %>% 
  column_to_rownames("Symbol")
count_exp <- round(count_exp, digits = 0)
colnames(count_exp) <- str_sub(colnames(count_exp), 1, 16)
save(count_exp, file = "15.diff/TCGA-LUAD_count4diff.RData")


# 2.分组----
dup_sample_cols <- which(duplicated(colnames(count_exp)))
count_exp <- count_exp %>% 
  select(-all_of(dup_sample_cols))
table(str_sub(colnames(count_exp),14,16))
# 01A 01B 01C 02A 11A 11B 
# 513  14   1   2  58   1
group <- ifelse(str_sub(colnames(count_exp),14,14) == "0",
                "tumor", "normal")
table(group)
# normal  tumor 
# 59    530
class(group)
group <- as.factor(group)
class(group)
coldata <- data.frame(row.names = colnames(count_exp),
                      group = group) 


# 3.构建dds对象----
dds <- DESeqDataSetFromMatrix(countData = count_exp,
                              colData = coldata,
                              design= ~ group)


# 4.过滤在每个样本中count均为0的基因----
filt <- rowSums(counts(dds) > 0) > 0
dds_filt <- dds[filt,]


# 5.差异分析----
###设置筛选标准
adj_p <- 0.05
log2FC <- 2

### 差异分析
dds2 <- DESeq(dds_filt) #加参数parallel = T可以发挥计算机多核性能
resultsNames(dds2) #查看是谁比谁: "group_tumor_vs_normal"
res <- results(dds2)
res <- res[order(res$padj),]
head(res)
#summary(res)
res <- na.omit(res)

### 输出差异分析结果
write.table(res, file="15.diff/TCGA-LUAD_All_results.txt", row.names = T, quote = F)
save(res, file="15.diff/TCGA-LUAD_All_results.RData")

table(res$padj < adj_p)


# 6.挑选差异表达基因----
DEGs <- subset(res, ((padj < adj_p) & (abs(log2FoldChange) > log2FC)))
upDEGs <- subset(res, ((padj < adj_p) & (log2FoldChange) > log2FC))
downDEGs <- subset(res, ((padj < adj_p) & (log2FoldChange) < -log2FC))

dim(DEGs)
dim(upDEGs)
dim(downDEGs)
head(DEGs)

write.table(DEGs,file = paste0("15.diff/allDEG_T_vs_N_", adj_p,"_",log2FC,".txt"),
            sep = "\t", quote = F, col.names = NA)
write.table(upDEGs,file = paste0("15.diff/upDEG_T_vs_N_", adj_p,"_",log2FC,".txt"),
            sep = "\t", quote = F, col.names = NA)
write.table(downDEGs,file= paste0("15.diff/downDEG_T_vs_N_", adj_p,"_",log2FC,".txt"),
            sep = "\t", quote = F, col.names = NA)


# 7.绘制火山图----
library(ggplot2)
diff_output <- read.table("15.diff/TCGA-LUAD_All_results.txt")
### 数据分组
up <- diff_output %>% filter(log2FoldChange > log2FC & padj < adj_p)
down <- diff_output %>% filter(log2FoldChange < -log2FC & padj < adj_p)
not <- diff_output %>% filter(abs(log2FoldChange) <= log2FC | padj >= adj_p)

### 增加因子型变量Expression
diff_result <- diff_output %>% 
  mutate(DEGs = case_when(log2FoldChange > log2FC & padj < adj_p ~ "up",
                          log2FoldChange < -log2FC & padj < adj_p ~ "down",
                          abs(log2FoldChange) <= log2FC | padj >= adj_p ~ "not")) %>%
  mutate_at(vars(DEGs), as.factor)

### 绘图
volcano_plot <- ggplot(data = diff_result, aes(log2FoldChange, -log10(padj), color = DEGs))+
  #画点，用DEG变量填充颜色，size调整点的大小，alpha调整透明度
  geom_point(size = 1.5, alpha = 0.6)+
  #手动定义颜色
  scale_color_manual(values = c("#2166AC", "grey", "#B2182B"))+
  ##加竖线，lty = 2 表示虚线
  geom_vline(xintercept = c(-log2FC, log2FC), lty = 2)+
  ##加横线
  geom_hline(yintercept = -log10(adj_p), lty = 2)+
  ##主题样式
  theme_classic(base_size = 15)

volcano_plot

### 保存火山图
ggsave(volcano_plot, filename = "15.diff/volcano.pdf", height = 4, width = 5, dpi = 600)
