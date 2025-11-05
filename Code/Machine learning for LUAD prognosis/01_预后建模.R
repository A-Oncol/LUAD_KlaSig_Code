######目的：预后建模
######作者：申奥
######日期：2024-11-11
######Server: yjx@xiyoucloud.top:34022
######Version: R 4.4.1


# 0.环境设置----
rm(list = ls())
options(stringAsFactors = F)

library(tidyverse)
library(dplyr)
library(data.table)
library(openxlsx)
library(Mime1)
library(survival)
library(survminer)
library(randomForestSRC)
library(glmnet)
library(plsRcox)
library(superpc)
library(gbm)
library(CoxBoost)
library(survivalsvm)
library(tibble)
library(BART)
library(miscTools)
library(compareC)
library(ggplot2)
library(ggsci)
library(tidyr)
library(ggbreak)


# 1.整理各队列的表达谱和生存信息----
## 1.1 TCGA-LUAD----
load("12.Prognosis_model/TCGA-LUAD/TCGA-LUAD_mRNA_tpm.RData")
load("12.Prognosis_model/TCGA-LUAD/TCGA-LUAD_lncRNA_tpm.RData")
head(mRNA_tpm_annoted)[1:4, 1:4]
head(lnc_tpm_annoted)[1:4, 1:4]
colnames(mRNA_tpm_annoted)[1] <- "Symbol"
colnames(lnc_tpm_annoted)[1] <- "Symbol"
identical(colnames(mRNA_tpm_annoted), colnames(lnc_tpm_annoted))

TCGA_LUAD_TPM <- rbind(mRNA_tpm_annoted, lnc_tpm_annoted)
which(duplicated(TCGA_LUAD_TPM$Symbol))
TCGA_LUAD_TPM <- aggregate(. ~ Symbol, TCGA_LUAD_TPM, mean)
which(duplicated(TCGA_LUAD_TPM$Symbol))
TCGA_LUAD_TPM <- TCGA_LUAD_TPM %>% 
  column_to_rownames("Symbol")

colnames(TCGA_LUAD_TPM) <- str_sub(colnames(TCGA_LUAD_TPM), 1, 16)
a <- str_split_i(colnames(TCGA_LUAD_TPM), "-", 4)
table(a)
b <- which(str_sub(colnames(TCGA_LUAD_TPM), 14, 16) == "01A")

luad_tumor_tpm <- TCGA_LUAD_TPM[,b]
log_luad_tumor_tpm <- as.data.frame(t(log2(luad_tumor_tpm + 1)))
which(duplicated(rownames(log_luad_tumor_tpm)))

surv_info <- read.delim("12.Prognosis_model/TCGA-LUAD/Survival_SupplementalTable_S1_20171025_xena_sp", 
                        sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T) 
table(surv_info$`cancer type abbreviation`)
luad_surv_info <- surv_info %>% 
  filter(`cancer type abbreviation` == "LUAD") %>% 
  select(25, 26)
which(duplicated(rownames(luad_surv_info)))
luad_surv_info$Sample_ID <- rownames(luad_surv_info)

log_luad_tumor_tpm$Sample_ID <- str_sub(rownames(log_luad_tumor_tpm), 1, 15)
same_samples <- intersect(luad_surv_info$Sample_ID, log_luad_tumor_tpm$Sample_ID)
TCGA_LUAD <- inner_join(luad_surv_info, log_luad_tumor_tpm, by = "Sample_ID")
dup_samples <- which(duplicated(TCGA_LUAD$Sample_ID))
TCGA_LUAD <- TCGA_LUAD[-dup_samples,]
rownames(TCGA_LUAD) <- NULL
head(TCGA_LUAD)[1:4, 1:4]
TCGA_LUAD <- TCGA_LUAD[,c("Sample_ID", "OS", "OS.time", colnames(TCGA_LUAD)[4:ncol(TCGA_LUAD)])]
head(TCGA_LUAD)[1:4, 1:4]

save(TCGA_LUAD, file = "12.Prognosis_model/TCGA-LUAD/TCGA-LUAD_surv_log2TPM.RData")

## 1.2 GSE14814----
rm(list = ls())

gse14814_exp <- read.table("12.Prognosis_model/GSE14814/GSE14814_series_matrix.txt",
                           header = T, sep = "\t", comment.char = "!")

gpl96_annot <- fread("12.Prognosis_model/GSE14814/GPL96-57554.txt",
                     header = T, sep = "\t", skip = 16)
gpl96_annot <- gpl96_annot[,c(1,11)]
colnames(gpl96_annot) <- c("ID_REF", "Gene_Symbol")

gse14814_exp_annot <- inner_join(gpl96_annot, gse14814_exp, by = "ID_REF")
gse14814_exp_annot <- gse14814_exp_annot[,-1]
gse14814_exp_annot <- gse14814_exp_annot[-which(str_detect(gse14814_exp_annot$Gene_Symbol, "///")),]
which(duplicated(gse14814_exp_annot$Gene_Symbol))
gse14814_exp_annot <- aggregate(. ~ Gene_Symbol, gse14814_exp_annot, mean)
which(duplicated(gse14814_exp_annot$Gene_Symbol))

gse14814_exp_annot <- gse14814_exp_annot %>% 
  filter(Gene_Symbol != "") %>% 
  column_to_rownames("Gene_Symbol")
gse14814_exp_annot <- as.data.frame(t(gse14814_exp_annot))
gse14814_exp_annot$Sample_ID <- rownames(gse14814_exp_annot)

gse14814_clinic <- read.xlsx("12.Prognosis_model/GSE14814/GSE14814_clinical.xlsx", sheet = 1)
table(gse14814_clinic$Histology.type)
# ADC LCUC SQCC 
# 71   10   52 
gse14814_clinic <- gse14814_clinic %>% 
  filter(Histology.type == "ADC") %>% 
  dplyr::select(1,9,10) %>% 
  na.omit() %>% 
  as.data.frame()
gse14814_clinic$OS.time <- as.numeric(gse14814_clinic$OS.time)
table(gse14814_clinic$OS.status)
gse14814_clinic$OS.status <- ifelse(gse14814_clinic$OS.status == "Dead", 1, 0)
table(gse14814_clinic$OS.status)
colnames(gse14814_clinic)[3] <- "OS"
gse14814_clinic <- gse14814_clinic[,c(1,3,2)]

GSE14814 <- inner_join(gse14814_clinic, gse14814_exp_annot, by = "Sample_ID")
save(GSE14814, file = "12.Prognosis_model/GSE14814/GSE14814_surv_exp.RData")

## 1.3 GSE29016----
rm(list = ls())

gse29016_exp <- read.table("12.Prognosis_model/GSE29016/GSE29016_series_matrix.txt",
                           header = T, sep = "\t", comment.char = "!")

gpl6947_annot <- fread("12.Prognosis_model/GSE29016/GPL6947-13512.txt",
                       header = T, sep = "\t", skip = 30)
gpl6947_annot <- gpl6947_annot %>%
  filter(Source == "RefSeq")
gpl6947_annot <- gpl6947_annot[,c(1,7)]
colnames(gpl6947_annot) <- c("ID_REF", "Gene_Symbol")

gse29016_exp_annot <- inner_join(gpl6947_annot, gse29016_exp, by = "ID_REF")
gse29016_exp_annot <- gse29016_exp_annot[,-1]
which(duplicated(gse29016_exp_annot$Gene_Symbol))
gse29016_exp_annot <- aggregate(. ~ Gene_Symbol, gse29016_exp_annot, mean)
which(duplicated(gse29016_exp_annot$Gene_Symbol))

gse29016_exp_annot <- gse29016_exp_annot %>%
  filter(Gene_Symbol != "") %>%
  column_to_rownames("Gene_Symbol")
gse29016_exp_annot <- log2(gse29016_exp_annot + 1)
gse29016_exp_annot <- as.data.frame(t(gse29016_exp_annot))
gse29016_exp_annot$Sample_ID <- rownames(gse29016_exp_annot)

gse29016_clinic <- read.xlsx("12.Prognosis_model/GSE29016/GSE29016_clinical.xlsx", sheet = 1)
table(gse29016_clinic$histology)
# AC LCNEC  SCLC  SqCC
# 38     9     7    12
gse29016_clinic <- gse29016_clinic %>%
  filter(histology == "AC") %>%
  dplyr::select(1,4,5) %>%
  na.omit() %>%
  as.data.frame()
colnames(gse29016_clinic) <- c("Sample_ID", "OS", "OS.time")
gse29016_clinic$OS.time <- as.numeric(gse29016_clinic$OS.time)
table(gse29016_clinic$OS)

GSE29016 <- inner_join(gse29016_clinic, gse29016_exp_annot, by = "Sample_ID")
save(GSE29016, file = "12.Prognosis_model/GSE29016/GSE29016_surv_exp.RData")

## 1.4 GSE30219----
rm(list = ls())

gse30219_exp <- read.table("12.Prognosis_model/GSE30219/GSE30219_series_matrix.txt",
                           header = T, sep = "\t", comment.char = "!")

gpl570 <- fread("12.Prognosis_model/GSE30219/GPL570-55999.txt", 
                header = T, sep = "\t", skip = 16)
gpl570 <- gpl570 %>% 
  select(1,11)
colnames(gpl570) <- c("ID_REF", "Gene_Symbol")

gse30219_exp_annot <- inner_join(gpl570, gse30219_exp, by = "ID_REF")
gse30219_exp_annot <- gse30219_exp_annot[,-1]
gse30219_exp_annot <- gse30219_exp_annot[-which(str_detect(gse30219_exp_annot$Gene_Symbol, "///")),]
which(duplicated(gse30219_exp_annot$Gene_Symbol))
gse30219_exp_annot <- aggregate(. ~ Gene_Symbol, gse30219_exp_annot, mean)
which(duplicated(gse30219_exp_annot$Gene_Symbol))

gse30219_exp_annot <- gse30219_exp_annot %>% 
  filter(Gene_Symbol != "") %>% 
  column_to_rownames("Gene_Symbol")
gse30219_exp_annot <- as.data.frame(t(gse30219_exp_annot))
gse30219_exp_annot$Sample_ID <- rownames(gse30219_exp_annot)

gse30219_clinic <- read.xlsx("12.Prognosis_model/GSE30219/GSE30219_clinical.xlsx", sheet = 1)
table(gse30219_clinic$histology)
# ADC   BAS CARCI   LCC  LCNE   NTL Other   SCC   SQC 
# 85    39    24     3    56    14     4    21    61 
gse30219_clinic <- gse30219_clinic %>% 
  filter(histology == "ADC") %>% 
  dplyr::select(1,10,9) %>% 
  na.omit() %>% 
  as.data.frame()
colnames(gse30219_clinic) <- c("Sample_ID", "OS", "OS.time")
gse30219_clinic$OS.time <- as.numeric(gse30219_clinic$OS.time)
table(gse30219_clinic$OS)
gse30219_clinic$OS <- ifelse(gse30219_clinic$OS == "DEAD", 1, 0)
table(gse30219_clinic$OS)

GSE30219 <- inner_join(gse30219_clinic, gse30219_exp_annot, by = "Sample_ID")
save(GSE30219, file = "12.Prognosis_model/GSE30219/GSE30219_surv_exp.RData")

## 1.5 GSE31210----
rm(list = ls())

gse31210_exp <- read.table("12.Prognosis_model/GSE31210/GSE31210_series_matrix.txt",
                           header = T, sep = "\t", comment.char = "!")

gpl570 <- fread("12.Prognosis_model/GSE31210/GPL570-55999.txt", 
                header = T, sep = "\t", skip = 16)
gpl570 <- gpl570 %>% 
  select(1,11)
colnames(gpl570) <- c("ID_REF", "Gene_Symbol")

gse31210_exp_annot <- inner_join(gpl570, gse31210_exp, by = "ID_REF")
gse31210_exp_annot <- gse31210_exp_annot[,-1]
gse31210_exp_annot <- gse31210_exp_annot[-which(str_detect(gse31210_exp_annot$Gene_Symbol, "///")),]
which(duplicated(gse31210_exp_annot$Gene_Symbol))
gse31210_exp_annot <- aggregate(. ~ Gene_Symbol, gse31210_exp_annot, mean)
which(duplicated(gse31210_exp_annot$Gene_Symbol))

gse31210_exp_annot <- gse31210_exp_annot %>% 
  filter(Gene_Symbol != "") %>% 
  column_to_rownames("Gene_Symbol")
gse31210_exp_annot2 <- log2(gse31210_exp_annot + 1)
gse31210_exp_annot2 <- as.data.frame(t(gse31210_exp_annot2))
gse31210_exp_annot2$Sample_ID <- rownames(gse31210_exp_annot2)

gse31210_clinic <- read.xlsx("12.Prognosis_model/GSE31210/GSE31210_clinical.xlsx", sheet = 1)
gse31210_clinic <- gse31210_clinic %>% 
  select(1,17,18) %>% 
  na.omit() %>% 
  as.data.frame()
colnames(gse31210_clinic) <- c("Sample_ID", "OS", "OS.time")
gse31210_clinic$OS.time <- as.numeric(gse31210_clinic$OS.time)
table(gse31210_clinic$OS)
gse31210_clinic$OS <- ifelse(gse31210_clinic$OS == "dead", 1, 0)
table(gse31210_clinic$OS)

GSE31210 <- inner_join(gse31210_clinic, gse31210_exp_annot2, by = "Sample_ID")
save(GSE31210, file = "12.Prognosis_model/GSE31210/GSE31210_surv_exp.RData")

## 1.6 GSE37745----
rm(list = ls())

gse37745_exp <- read.table("12.Prognosis_model/GSE37745/GSE37745_series_matrix.txt",
                           header = T, sep = "\t", comment.char = "!")

gpl570 <- fread("12.Prognosis_model/GSE37745/GPL570-55999.txt", 
                header = T, sep = "\t", skip = 16)
gpl570 <- gpl570 %>% 
  select(1,11)
colnames(gpl570) <- c("ID_REF", "Gene_Symbol")

gse37745_exp_annot <- inner_join(gpl570, gse37745_exp, by = "ID_REF")
gse37745_exp_annot <- gse37745_exp_annot[,-1]
gse37745_exp_annot <- gse37745_exp_annot[-which(str_detect(gse37745_exp_annot$Gene_Symbol, "///")),]
which(duplicated(gse37745_exp_annot$Gene_Symbol))
gse37745_exp_annot <- aggregate(. ~ Gene_Symbol, gse37745_exp_annot, mean)
which(duplicated(gse37745_exp_annot$Gene_Symbol))

gse37745_exp_annot <- gse37745_exp_annot %>% 
  filter(Gene_Symbol != "") %>% 
  column_to_rownames("Gene_Symbol")
gse37745_exp_annot <- as.data.frame(t(gse37745_exp_annot))
gse37745_exp_annot$Sample_ID <- rownames(gse37745_exp_annot)

gse37745_clinic <- read.xlsx("12.Prognosis_model/GSE37745/GSE37745_clinical.xlsx", sheet = 1)
table(gse37745_clinic$histology)
# adeno    large squamous 
# 106       24       66
gse37745_clinic <- gse37745_clinic %>% 
  filter(histology == "adeno") %>% 
  select(1,3,4) %>% 
  na.omit() %>% 
  as.data.frame()
colnames(gse37745_clinic) <- c("Sample_ID", "OS", "OS.time")
gse37745_clinic$OS.time <- as.numeric(gse37745_clinic$OS.time)
table(gse37745_clinic$OS)
gse37745_clinic$OS <- ifelse(gse37745_clinic$OS == "yes", 1, 0)
table(gse37745_clinic$OS)

GSE37745 <- inner_join(gse37745_clinic, gse37745_exp_annot, by = "Sample_ID")
save(GSE37745, file = "12.Prognosis_model/GSE37745/GSE37745_surv_exp.RData")

## 1.7 GSE42127----
rm(list = ls())

gse42127_exp <- read.table("12.Prognosis_model/GSE42127/GSE42127_series_matrix.txt", 
                           header = T, sep = "\t", comment.char = "!", check.names = F)

gpl6884 <- fread("12.Prognosis_model/GSE42127/GPL6884.annot.gz")
gpl6884 <- gpl6884 %>%
  select(1,3)
colnames(gpl6884) <- c("ID_REF", "Gene_Symbol")

gse42127_exp_annot <- inner_join(gpl6884, gse42127_exp, by = "ID_REF")
gse42127_exp_annot <- gse42127_exp_annot[,-1]
gse42127_exp_annot <- gse42127_exp_annot %>% 
  filter(Gene_Symbol != "")
which(duplicated(gse42127_exp_annot$Gene_Symbol))
gse42127_exp_annot <- aggregate(. ~ Gene_Symbol, gse42127_exp_annot, mean)
which(duplicated(gse42127_exp_annot$Gene_Symbol))
gse42127_exp_annot <- gse42127_exp_annot %>% 
  column_to_rownames("Gene_Symbol")
gse42127_exp_annot <- as.data.frame(t(gse42127_exp_annot))
gse42127_exp_annot$Sample_ID <- rownames(gse42127_exp_annot)

gse42127_clinic <- read.xlsx("12.Prognosis_model/GSE42127/GSE42127_clinical_info.xlsx", sheet = 1)
table(gse42127_clinic$histology)
# Adenocarcionoma        Squamous 
# 133              43
gse42127_clinic <- gse42127_clinic %>% 
  filter(histology == "Adenocarcionoma")
gse42127_clinic <- gse42127_clinic %>% 
  select(1,5,4)
colnames(gse42127_clinic) <- c("Sample_ID", "OS", "OS.time")
gse42127_clinic$OS.time <- as.numeric(gse42127_clinic$OS.time)
table(gse42127_clinic$OS)
gse42127_clinic$OS <- ifelse(gse42127_clinic$OS == "D", 1, 0)
table(gse42127_clinic$OS)

GSE42127 <- inner_join(gse42127_clinic, gse42127_exp_annot, by = "Sample_ID")
save(GSE42127, file = "12.Prognosis_model/GSE42127/GSE42127_surv_exp.RData")

## 1.8 GSE50081----
rm(list = ls())

gse50081_exp <- read.table("12.Prognosis_model/GSE50081/GSE50081_series_matrix.txt", 
                           header = T, sep = "\t", comment.char = "!", check.names = F)

gpl570 <- fread("12.Prognosis_model/GSE50081/GPL570-55999.txt", 
                header = T, sep = "\t", skip = 16)
gpl570 <- gpl570 %>% 
  select(1,11)
colnames(gpl570) <- c("ID_REF", "Gene_Symbol")

gse50081_exp_annot <- inner_join(gpl570, gse50081_exp, by = "ID_REF")
gse50081_exp_annot <- gse50081_exp_annot[,-1]
gse50081_exp_annot <- gse50081_exp_annot[-which(str_detect(gse50081_exp_annot$Gene_Symbol, "///")),]
which(duplicated(gse50081_exp_annot$Gene_Symbol))
gse50081_exp_annot <- aggregate(. ~ Gene_Symbol, gse50081_exp_annot, mean)
which(duplicated(gse50081_exp_annot$Gene_Symbol))

gse50081_exp_annot <- gse50081_exp_annot %>% 
  filter(Gene_Symbol != "") %>% 
  column_to_rownames("Gene_Symbol")
gse50081_exp_annot <- as.data.frame(t(gse50081_exp_annot))
gse50081_exp_annot$Sample_ID <- rownames(gse50081_exp_annot)

gse50081_clinic <- read.xlsx("12.Prognosis_model/GSE50081/GSE50081_clinical.xlsx", sheet = 1)
table(gse50081_clinic$histology)
# adenocarcinoma       adenosquamous carcinoma          large cell carcinoma NSClarge cell carcinoma-mixed    NSCLC-favor adenocarcinoma 
# 127                             2                             7                             1                             1 
# squamous cell carcinoma    squamous cell carcinoma X2 
# 42                             1 
gse50081_clinic <- gse50081_clinic %>% 
  filter(histology == "adenocarcinoma") %>% 
  select(1,12,11) %>% 
  na.omit() %>% 
  as.data.frame()
colnames(gse50081_clinic) <- c("Sample_ID", "OS", "OS.time")
gse50081_clinic$OS.time <- as.numeric(gse50081_clinic$OS.time)
table(gse50081_clinic$OS)
gse50081_clinic$OS <- ifelse(gse50081_clinic$OS == "dead", 1, 0)
table(gse50081_clinic$OS)

GSE50081 <- inner_join(gse50081_clinic, gse50081_exp_annot, by = "Sample_ID")
save(GSE50081, file = "12.Prognosis_model/GSE50081/GSE50081_surv_exp.RData")

## 1.9 GSE68465----
rm(list = ls())

gse68465_exp <- read.table("12.Prognosis_model/GSE68465/GSE68465_series_matrix.txt", 
                           header = T, sep = "\t", comment.char = "!", check.names = F)

gpl96_annot <- fread("12.Prognosis_model/GSE68465/GPL96-57554.txt",
                     header = T, sep = "\t", skip = 16)
gpl96_annot <- gpl96_annot[,c(1,11)]
colnames(gpl96_annot) <- c("ID_REF", "Gene_Symbol")

gse68465_exp_annot <- inner_join(gpl96_annot, gse68465_exp, by = "ID_REF")
gse68465_exp_annot <- gse68465_exp_annot[,-1]
gse68465_exp_annot <- gse68465_exp_annot[-which(str_detect(gse68465_exp_annot$Gene_Symbol, "///")),]
which(duplicated(gse68465_exp_annot$Gene_Symbol))
gse68465_exp_annot <- aggregate(. ~ Gene_Symbol, gse68465_exp_annot, mean)
which(duplicated(gse68465_exp_annot$Gene_Symbol))

gse68465_exp_annot <- gse68465_exp_annot %>% 
  filter(Gene_Symbol != "") %>% 
  column_to_rownames("Gene_Symbol")
gse68465_exp_annot2 <- log2(gse68465_exp_annot + 1)
gse68465_exp_annot2 <- as.data.frame(t(gse68465_exp_annot2))
gse68465_exp_annot2$Sample_ID <- rownames(gse68465_exp_annot2)

gse68465_clinic <- read.xlsx("12.Prognosis_model/GSE68465/GSE68465_clinical.xlsx", sheet = 1)
table(gse68465_clinic$disease_state)
# Lung Adenocarcinoma              Normal 
# 443                  19
gse68465_clinic <- gse68465_clinic %>% 
  filter(disease_state == "Lung Adenocarcinoma") %>% 
  select(1,6,13) %>% 
  na.omit() %>% 
  as.data.frame()
colnames(gse68465_clinic) <- c("Sample_ID", "OS", "OS.time")
gse68465_clinic$OS.time <- as.numeric(gse68465_clinic$OS.time)
table(gse68465_clinic$OS)
gse68465_clinic$OS <- ifelse(gse68465_clinic$OS == "Dead", 1, 0)
table(gse68465_clinic$OS)

GSE68465 <- inner_join(gse68465_clinic, gse68465_exp_annot2, by = "Sample_ID")
save(GSE68465, file = "12.Prognosis_model/GSE68465/GSE68465_surv_exp.RData")

## 1.10 GSE72094----
rm(list = ls())

gse72094_exp <- read.table("12.Prognosis_model/GSE72094/GSE72094_series_matrix.txt", 
                           header = T, sep = "\t", comment.char = "!", check.names = F)

gpl15048 <- fread("12.Prognosis_model/GSE72094/GPL15048.txt", header = T, sep = "\t", skip = 5)
gpl15048 <- gpl15048 %>% 
  select(1,4)
colnames(gpl15048) <- c("ID_REF", "Gene_Symbol")
gpl15048 <- gpl15048 %>% 
  filter(Gene_Symbol != "")

gse72094_exp_annot <- inner_join(gpl15048, gse72094_exp, by = "ID_REF")
gse72094_exp_annot <- gse72094_exp_annot[,-1]
which(duplicated(gse72094_exp_annot$Gene_Symbol))
gse72094_exp_annot <- aggregate(. ~ Gene_Symbol, gse72094_exp_annot, mean)
which(duplicated(gse72094_exp_annot$Gene_Symbol))

gse72094_exp_annot <- gse72094_exp_annot %>% 
  column_to_rownames("Gene_Symbol")
gse72094_exp_annot <- as.data.frame(t(gse72094_exp_annot))
gse72094_exp_annot$Sample_ID <- rownames(gse72094_exp_annot)

gse72094_clinic <- read.xlsx("12.Prognosis_model/GSE72094/GSE72094_clinical.xlsx", sheet = 1)
gse72094_clinic <- gse72094_clinic %>% 
  select(1,11,12) %>% 
  na.omit() %>% 
  as.data.frame()
colnames(gse72094_clinic) <- c("Sample_ID", "OS", "OS.time")
gse72094_clinic$OS.time <- as.numeric(gse72094_clinic$OS.time)
table(gse72094_clinic$OS)
gse72094_clinic$OS <- ifelse(gse72094_clinic$OS == "Dead", 1, 0)
table(gse72094_clinic$OS)

GSE72094 <- inner_join(gse72094_clinic, gse72094_exp_annot, by = "Sample_ID")
save(GSE72094, file = "12.Prognosis_model/GSE72094/GSE72094_surv_exp.RData")

## 1.11 汇总成列表格式----
rm(list = ls())

load("12.Prognosis_model/TCGA-LUAD/TCGA-LUAD_surv_log2TPM.RData")
load("12.Prognosis_model/GSE14814/GSE14814_surv_exp.RData")
load("12.Prognosis_model/GSE29016/GSE29016_surv_exp.RData")
load("12.Prognosis_model/GSE30219/GSE30219_surv_exp.RData")
load("12.Prognosis_model/GSE31210/GSE31210_surv_exp.RData")
load("12.Prognosis_model/GSE37745/GSE37745_surv_exp.RData")
load("12.Prognosis_model/GSE42127/GSE42127_surv_exp.RData")
load("12.Prognosis_model/GSE50081/GSE50081_surv_exp.RData")
load("12.Prognosis_model/GSE68465/GSE68465_surv_exp.RData")
load("12.Prognosis_model/GSE72094/GSE72094_surv_exp.RData")

range(TCGA_LUAD$OS.time)
table(is.na(TCGA_LUAD$OS.time))
TCGA_LUAD <- TCGA_LUAD %>% 
  filter(!is.na(TCGA_LUAD$OS.time))
range(TCGA_LUAD$OS.time)
TCGA_LUAD$OS.time <- TCGA_LUAD$OS.time / 365
range(TCGA_LUAD$OS.time)
table(TCGA_LUAD$OS)

range(GSE14814$OS.time)
table(GSE14814$OS)

range(GSE29016$OS.time)
table(GSE29016$OS)

range(GSE30219$OS.time)
GSE30219$OS.time <- GSE30219$OS.time / 12
range(GSE30219$OS.time)
table(GSE30219$OS)

range(GSE31210$OS.time)
GSE31210$OS.time <- GSE31210$OS.time / 365
range(GSE31210$OS.time)
table(GSE31210$OS)

range(GSE37745$OS.time)
GSE37745$OS.time <- GSE37745$OS.time / 365
range(GSE37745$OS.time)
table(GSE37745$OS)

range(GSE42127$OS.time)
GSE42127$OS.time <- GSE42127$OS.time / 12
range(GSE42127$OS.time)
table(GSE42127$OS)

range(GSE50081$OS.time)
table(GSE50081$OS)

range(GSE68465$OS.time)
table(is.na(GSE68465$OS.time))
GSE68465 <- GSE68465 %>% 
  filter(!is.na(GSE68465$OS.time))
range(GSE68465$OS.time)
GSE68465$OS.time <- GSE68465$OS.time / 12
range(GSE68465$OS.time)
table(GSE68465$OS)

range(GSE72094$OS.time)
GSE72094$OS.time <- GSE72094$OS.time / 365
range(GSE72094$OS.time)
table(GSE72094$OS)

merge_datasets <- list(TCGA_LUAD = TCGA_LUAD,
                       GSE14814 = GSE14814,
                       #GSE26939 = GSE26939,
                       GSE29016 = GSE29016,
                       GSE30219 = GSE30219,
                       GSE31210 = GSE31210,
                       GSE37745 = GSE37745,
                       GSE42127 = GSE42127,
                       GSE50081 = GSE50081,
                       GSE68465 = GSE68465,
                       GSE72094 = GSE72094
)

merge_datasets <- lapply(merge_datasets, function(df) {
  df$OS <- as.integer(df$OS)  # 将目标列因子化
  return(df)  # 返回修改后的数据框
})

same_genes <- Reduce(intersect, lapply(merge_datasets, colnames))

save(merge_datasets, file = "12.Prognosis_model/merge_datasets.RData")
save(same_genes, file = "12.Prognosis_model/same_genes.RData")


# 2.建模前准备工作----
## 2.1 加载用于建模的基因----
load("15.diff/upKla_genes.RData")
load("14.Location_heterogeneity/geneExp/upKla_highUpKlaTumor_sameMal_intersect_markers.RData")
load("14.Location_heterogeneity/geneExp/upKla_highUpKlaTumor_allMal_intersect_markers.RData")

high_upKla_tumorCell_sig <- readRDS("5.scRNA/FindAllMarkers_GSE189357/high_upKla_tumorCell_markers.rds")
high_upKla_tumorCell_sig_genes <- high_upKla_tumorCell_sig$gene
a <- intersect(upKla_genes, high_upKla_tumorCell_sig_genes)

table(upKla_genes %in% same_genes)
# FALSE  TRUE 
# 436   318
table(degs2 %in% same_genes)
# FALSE  TRUE 
# 6    18
table(a %in% same_genes)
# FALSE  TRUE 
# 9    22
build_model_genes <- intersect(upKla_genes, same_genes) # 318个

## 2.2 数据标准化----
mm <- lapply(merge_datasets,function(x){
  x[,-c(1:3)] <- scale(x[,-c(1:3)])
  return(x)})

result <- data.frame()
# TCGA作为训练集
est_data <- mm$TCGA_LUAD
# GEO作为验证集
val_data_list <- mm
pre_var <- build_model_genes
est_dd <- est_data[, c('OS.time', 'OS', pre_var)]
val_dd_list <- lapply(val_data_list, function(x){x[, c('OS.time', 'OS', pre_var)]})


# 3.按照101代码流程建模----
# 设置种子数
seed <- 123

## 3.1 RSF----
rf_nodesize <- 15
set.seed(seed)
fit <- rfsrc(Surv(OS.time,OS)~., data = est_dd, 
             ntree = 1000, nodesize = rf_nodesize,  #该值建议多调整
             splitrule = 'logrank', 
             importance = T, 
             proximity = T, 
             forest = T, 
             seed = seed)
print(fit)
# fit$xvar.names
# importance <- as.data.frame(fit$importance)
rs <- lapply(val_dd_list, function(x){cbind(x[, 1:2], RS  = predict(fit, newdata = x)$predicted)})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- 'RSF'
result <- rbind(result, cc)

## 3.2 RSF + CoxBoost----
set.seed(seed)
fit <- rfsrc(Surv(OS.time, OS)~., data = est_dd, 
             ntree = 1000, nodesize = rf_nodesize,  #该值建议多调整
             splitrule = 'logrank', 
             importance = T, 
             proximity = T, 
             forest = T, 
             seed = seed)
rid <- var.select(object = fit, conservative = "high")
rid <- rid$topvars

est_dd2 <- est_dd[, c('OS.time', 'OS', rid)]
val_dd_list2 <- lapply(val_dd_list, function(x){x[, c('OS.time', 'OS', rid)]})

set.seed(seed)
pen <- optimCoxBoostPenalty(est_dd2[, 'OS.time'], est_dd2[, 'OS'], as.matrix(est_dd2[, -c(1, 2)]), 
                            trace=TRUE, start.penalty = 500, parallel = T)

cv.res <- cv.CoxBoost(est_dd2[, 'OS.time'], est_dd2[, 'OS'], as.matrix(est_dd2[, -c(1, 2)]), 
                      maxstepno = 500, K = 10, type = "verweij",  penalty = pen$penalty)
fit <- CoxBoost(est_dd2[, 'OS.time'], est_dd2[, 'OS'], as.matrix(est_dd2[, -c(1, 2)]), 
                stepno = cv.res$optimal.step, penalty = pen$penalty)
rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, newdata = x[, -c(1, 2)], newtime = x[, 1],  newstatus = x[, 2], type = "lp")))})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('RSF + CoxBoost')
result <- rbind(result, cc)

## 3.3 RSF + Enet----
set.seed(seed)
fit <- rfsrc(Surv(OS.time, OS)~., data = est_dd, 
             ntree = 1000, nodesize = rf_nodesize, #该值建议多调整
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)
rid <- var.select(object = fit, conservative = "high")
rid <- rid$topvars

est_dd2 <- est_dd[, c('OS.time', 'OS', rid)]
val_dd_list2 <- lapply(val_dd_list, function(x){x[, c('OS.time', 'OS', rid)]})

x1 <- as.matrix(est_dd2[, rid])
x2 <- as.matrix(Surv(est_dd2$OS.time, est_dd2$OS))
for (alpha in seq(0.1, 0.9, 0.1)) {
  set.seed(seed)
  fit = cv.glmnet(x1, x2, family = "cox", alpha = alpha, nfolds = 10)
  rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, type = 'link', newx = as.matrix(x[, -c(1, 2)]), s = fit$lambda.min)))})
  cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
    rownames_to_column('ID')
  cc$Model <- paste0('RSF + ', 'Enet', '[α=', alpha, ']')
  result <- rbind(result, cc)
}

## 3.4 RSF + GBM----
set.seed(seed)
fit <- rfsrc(Surv(OS.time, OS)~., data = est_dd,
             ntree = 1000, nodesize = rf_nodesize,  #该值建议多调整
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)
rid <- var.select(object = fit, conservative = "high")
rid <- rid$topvars

est_dd2 <- est_dd[, c('OS.time', 'OS', rid)]
val_dd_list2 <- lapply(val_dd_list, function(x){x[, c('OS.time', 'OS', rid)]})

set.seed(seed)
fit <- gbm(formula = Surv(OS.time, OS)~., data = est_dd2, distribution = 'coxph',
           n.trees = 10000,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10, n.cores = 6)

# find index for number trees with minimum CV error
best <- which.min(fit$cv.error)
set.seed(seed)
fit <- gbm(formula = Surv(OS.time,OS)~., data = est_dd2, distribution = 'coxph',
           n.trees = best,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10,n.cores = 8)
rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, x, n.trees = best, type = 'link')))})
cc <- data.frame(Cindex=sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('RSF + ', 'GBM')
result <- rbind(result, cc)

## 3.5 RSF + Lasso----
set.seed(seed)
fit <- rfsrc(Surv(OS.time, OS)~., data = est_dd,
             ntree = 1000, nodesize = rf_nodesize, #该值建议多调整
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)
rid <- var.select(object = fit, conservative = "high")
rid <- rid$topvars

est_dd2 <- est_dd[, c('OS.time', 'OS', rid)]
val_dd_list2 <- lapply(val_dd_list, function(x){x[, c('OS.time', 'OS', rid)]})

x1 <- as.matrix(est_dd2[, rid])
x2 <- as.matrix(Surv(est_dd2$OS.time, est_dd2$OS))
set.seed(seed)
fit = cv.glmnet(x1, x2,
                nfold = 10, #例文描述：10-fold cross-validation
                family = "binomial", alpha = 1,
                type.measure = "class")
rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, type = 'response', newx = as.matrix(x[, -c(1, 2)]), s = fit$lambda.min)))})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('RSF + ', 'Lasso')
result <- rbind(result, cc)

## 3.6 RSF + plsRcox----
set.seed(seed)
fit <- rfsrc(Surv(OS.time, OS)~., data = est_dd,
             ntree = 1000, nodesize = rf_nodesize, #该值建议多调整
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)
rid <- var.select(object = fit, conservative = "high")
rid <- rid$topvars

est_dd2 <- est_dd[, c('OS.time', 'OS', rid)]
val_dd_list2 <- lapply(val_dd_list, function(x){x[, c('OS.time', 'OS', rid)]})

set.seed(seed)
cv.plsRcox.res = cv.plsRcox(list(x = est_dd2[, rid], time = est_dd2$OS.time, status = est_dd2$OS), nt = 10, verbose = FALSE)
fit <- plsRcox(est_dd2[, rid], time = est_dd2$OS.time, event = est_dd2$OS, nt = as.numeric(cv.plsRcox.res[5]))
rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, type = "lp", newdata = x[, -c(1, 2)])))})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('RSF + ', 'plsRcox')
result <- rbind(result, cc)

## 3.7 RSF + Ridge----
set.seed(seed)
fit <- rfsrc(Surv(OS.time, OS)~., data = est_dd, 
             ntree = 1000, nodesize = rf_nodesize, #该值建议多调整
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)
rid <- var.select(object = fit, conservative = "high")
rid <- rid$topvars

est_dd2 <- est_dd[, c('OS.time', 'OS', rid)]
val_dd_list2 <- lapply(val_dd_list, function(x){x[, c('OS.time', 'OS', rid)]})

x1 <- as.matrix(est_dd2[, rid])
x2 <- as.matrix(Surv(est_dd2$OS.time, est_dd2$OS))
set.seed(seed)
fit = cv.glmnet(x1, x2,
                nfold=10, #例文描述：10-fold cross-validation
                family = "binomial", alpha = 0,
                type.measure = "class")
rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, type = 'response', newx = as.matrix(x[, -c(1, 2)]), s = fit$lambda.min)))})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('RSF + ', 'Ridge')
result <- rbind(result, cc)

## 3.8 RSF + StepCox----
set.seed(seed)
fit <- rfsrc(Surv(OS.time,OS)~., data = est_dd,
             ntree = 1000, nodesize = rf_nodesize, #该值建议多调整
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)
rid <- var.select(object = fit, conservative = "high")
rid <- rid$topvars

est_dd2 <- est_dd[, c('OS.time', 'OS', rid)]
val_dd_list2 <- lapply(val_dd_list, function(x){x[, c('OS.time', 'OS', rid)]})

for (direction in c("both", "backward", "forward")) {
  fit <- step(coxph(Surv(OS.time, OS)~., est_dd2), direction = direction)
  rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS=predict(fit, type = 'risk', newdata = x))})
  cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
    rownames_to_column('ID')
  cc$Model <- paste0('RSF + ', 'StepCox', '[', direction, ']')
  result <- rbind(result, cc)
}

## 3.9 RSF + SuperPC----
set.seed(seed)
fit <- rfsrc(Surv(OS.time,OS)~., data = est_dd,
             ntree = 1000, nodesize = rf_nodesize, ##该值建议多调整
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)
rid <- var.select(object = fit, conservative = "high")
rid <- rid$topvars

est_dd2 <- est_dd[, c('OS.time', 'OS', rid)]
val_dd_list2 <- lapply(val_dd_list, function(x){x[, c('OS.time', 'OS', rid)]})

data <- list(x = t(est_dd2[, -c(1, 2)]), y = est_dd2$OS.time,
             censoring.status = est_dd2$OS, 
             featurenames = colnames(est_dd2)[-c(1, 2)])
set.seed(seed)
fit <- superpc.train(data = data, type = 'survival', s0.perc = 0.5) #default
cv.fit <- superpc.cv(fit, data, n.threshold = 20, #default
                     n.fold = 10, 
                     n.components = 3, 
                     min.features = 5, 
                     max.features = nrow(data$x), 
                     compute.fullcv = TRUE, 
                     compute.preval = TRUE)
rs <- lapply(val_dd_list2, function(w){
  test <- list(x = t(w[, -c(1, 2)]), y = w$OS.time, censoring.status=w$OS, featurenames = colnames(w)[-c(1, 2)])
  ff <- superpc.predict(fit, data, test, threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1, ])], n.components = 1)
  rr <- as.numeric(ff$v.pred)
  rr2 <- cbind(w[, 1:2], RS = rr)
  return(rr2)
})
cc <- data.frame(Cindex=sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('RSF + ', 'SuperPC')
result <- rbind(result, cc)

## 3.10 RSF + survival-SVM----
set.seed(seed)
fit <- rfsrc(Surv(OS.time, OS)~., data = est_dd,
             ntree = 1000, nodesize = rf_nodesize, #该值建议多调整
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)
rid <- var.select(object = fit, conservative = "high")
rid <- rid$topvars

est_dd2 <- est_dd[, c('OS.time', 'OS', rid)]
val_dd_list2 <- lapply(val_dd_list, function(x){x[, c('OS.time', 'OS', rid)]})

fit = survivalsvm(Surv(OS.time, OS)~., data= est_dd2, gamma.mu = 1)
rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS=as.numeric(predict(fit, x)$predicted))})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('RSF + ', 'survival-SVM')
result <- rbind(result,cc)

## 3.11 Enet----
x1 <- as.matrix(est_dd[, pre_var])
x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
for (alpha in seq(0.1, 0.9, 0.1)) {
  set.seed(seed)
  fit = cv.glmnet(x1, x2, family = "cox", alpha = alpha, nfolds = 10)
  rs <- lapply(val_dd_list,function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit,type = 'link', newx = as.matrix(x[,-c(1,2)]), s = fit$lambda.min)))})
  cc <- data.frame(Cindex = sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
    rownames_to_column('ID')
  cc$Model <- paste0('Enet', '[α=', alpha, ']')
  result <- rbind(result, cc)
}

## 3.12 StepCox----
for (direction in c("both", "backward", "forward")) {
  fit <- step(coxph(Surv(OS.time,OS)~., est_dd), direction = direction)
  rs <- lapply(val_dd_list,function(x){cbind(x[, 1:2], RS = predict(fit, type = 'risk', newdata = x))})
  cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox', '[', direction, ']')
  result <- rbind(result, cc)
}

## 3.13 StepCox + 其它----
for (direction in c("both", "backward", "forward")) {
  fit <- step(coxph(Surv(OS.time, OS)~., est_dd), direction = direction)
  rid <- names(coef(fit)) #这里不用卡P值，迭代的结果就是可以纳入的基因
  
  est_dd2 <- est_dd[,c('OS.time', 'OS', rid)]
  val_dd_list2 <- lapply(val_dd_list, function(x){x[, c('OS.time', 'OS', rid)]})
  
  set.seed(seed)
  pen <- optimCoxBoostPenalty(est_dd2[, 'OS.time'], est_dd2[, 'OS'], as.matrix(est_dd2[, -c(1,2)]),
                              trace=TRUE, start.penalty = 500, parallel = T)
  cv.res <- cv.CoxBoost(est_dd2[, 'OS.time'], est_dd2[, 'OS'], as.matrix(est_dd2[, -c(1,2)]),
                        maxstepno = 500, K = 10 , type = "verweij", penalty = pen$penalty)
  fit <- CoxBoost(est_dd2[, 'OS.time'], est_dd2[, 'OS'], as.matrix(est_dd2[, -c(1, 2)]),
                  stepno = cv.res$optimal.step, penalty = pen$penalty)
  rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, newdata = x[, -c(1, 2)], newtime=x[, 1], newstatus=x[,2], type="lp")))})
  cc <- data.frame(Cindex = sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox', '[', direction, ']', ' + CoxBoost')
  result <- rbind(result, cc)
  
  x1 <- as.matrix(est_dd2[, rid])
  x2 <- as.matrix(Surv(est_dd2$OS.time, est_dd2$OS))
  for (alpha in seq(0.1, 0.9, 0.1)) {
    set.seed(seed)
    fit = cv.glmnet(x1, x2, family = "cox",alpha = alpha, nfolds = 10)
    rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, type = 'link', newx = as.matrix(x[, -c(1, 2)]), s = fit$lambda.min)))})
    cc <- data.frame(Cindex = sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
      rownames_to_column('ID')
    cc$Model <- paste0('StepCox', '[', direction, ']', ' + Enet', '[α=', alpha, ']')
    result <- rbind(result, cc)
  }
  
  set.seed(seed)
  fit <- gbm(formula = Surv(OS.time, OS)~., data = est_dd2, distribution = 'coxph',
             n.trees = 10000,
             interaction.depth = 3,
             n.minobsinnode = 10,
             shrinkage = 0.001,
             cv.folds = 10,n.cores = 6)
  # find index for number trees with minimum CV error
  best <- which.min(fit$cv.error)
  set.seed(seed)
  fit <- gbm(formula = Surv(OS.time, OS)~., data = est_dd2, distribution = 'coxph',
             n.trees = best,
             interaction.depth = 3,
             n.minobsinnode = 10,
             shrinkage = 0.001,
             cv.folds = 10,n.cores = 8)
  rs <- lapply(val_dd_list2, function(x){cbind(x[,1:2], RS = as.numeric(predict(fit, x, n.trees = best, type = 'link')))})
  cc <- data.frame(Cindex=sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox', '[', direction, ']', ' + GBM')
  result <- rbind(result, cc)
  
  x1 <- as.matrix(est_dd2[, rid])
  x2 <- as.matrix(Surv(est_dd2$OS.time, est_dd2$OS))
  set.seed(seed)
  fit = cv.glmnet(x1, x2,
                  nfold=10, #例文描述：10-fold cross-validation
                  family = "binomial", alpha = 1,
                  type.measure = "class")
  rs <- lapply(val_dd_list2, function(x){cbind(x[,1:2], RS = as.numeric(predict(fit, type = 'response', newx = as.matrix(x[, -c(1, 2)]), s = fit$lambda.min)))})
  cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox', '[', direction, ']', ' + Lasso')
  result <- rbind(result, cc)
  
  set.seed(seed)
  cv.plsRcox.res = cv.plsRcox(list(x = est_dd2[,rid], time = est_dd2$OS.time, status = est_dd2$OS), nt = 10, verbose = FALSE)
  fit <- plsRcox(est_dd2[, rid], time = est_dd2$OS.time,
                 event = est_dd2$OS, nt = as.numeric(cv.plsRcox.res[5]))
  rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, type = "lp", newdata = x[, -c(1,2)])))})
  cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox', '[', direction, ']', ' + plsRcox')
  result <- rbind(result, cc)
  
  x1 <- as.matrix(est_dd2[, rid])
  x2 <- as.matrix(Surv(est_dd2$OS.time, est_dd2$OS))
  set.seed(seed)
  fit = cv.glmnet(x1, x2,
                  nfold = 10, #例文描述：10-fold cross-validation
                  family = "binomial", alpha = 0,
                  type.measure = "class")
  rs <- lapply(val_dd_list2, function(x){cbind(x[,1:2], RS = as.numeric(predict(fit, type = 'response', newx = as.matrix(x[, -c(1, 2)]), s = fit$lambda.min)))})
  cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox', '[', direction, ']', ' + Ridge')
  result <- rbind(result, cc)
  
  set.seed(seed)
  fit <- rfsrc(Surv(OS.time,OS)~., data = est_dd2,
               ntree = 1000, nodesize = rf_nodesize, #该值建议多调整
               splitrule = 'logrank',
               importance = T,
               proximity = T,
               forest = T,
               seed = seed)
  rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = predict(fit, newdata = x)$predicted)})
  cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox', '[', direction, ']', ' + RSF')
  result <- rbind(result, cc)
  
  data <- list(x = t(est_dd2[, -c(1, 2)]), y = est_dd2$OS.time,
               censoring.status = est_dd2$OS,
               featurenames = colnames(est_dd2)[-c(1,2)])
  set.seed(seed)
  fit <- superpc.train(data = data,type = 'survival', s0.perc = 0.5) #default
  cv.fit <- superpc.cv(fit, data, n.threshold = 20, #default
                       n.fold = 10,
                       n.components = 3,
                       min.features = 5,
                       max.features = nrow(data$x),
                       compute.fullcv = TRUE,
                       compute.preval = TRUE)
  rs <- lapply(val_dd_list2, function(w){
    test <- list(x = t(w[, -c(1,2)]), y = w$OS.time, censoring.status = w$OS, featurenames = colnames(w)[-c(1,2)])
    ff <- superpc.predict(fit, data, test, threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])], n.components = 1)
    rr <- as.numeric(ff$v.pred)
    rr2 <- cbind(w[,1:2], RS = rr)
    return(rr2)
  })
  cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox', '[', direction, ']', ' + SuperPC')
  result <- rbind(result, cc)
  
  fit = survivalsvm(Surv(OS.time,OS)~., data = est_dd2, gamma.mu = 1)
  rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, x)$predicted))})
  cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox', '[', direction, ']', ' + survival-SVM')
  result <- rbind(result, cc)
}

## 3.14 CoxBoost----
set.seed(seed)
pen <- optimCoxBoostPenalty(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]), 
                            trace = TRUE, #start.penalty = 500, 
                            parallel = T)
cv.res <- cv.CoxBoost(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]), maxstepno = 500, K = 10, type = "verweij", penalty = pen$penalty)
fit <- CoxBoost(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                stepno = cv.res$optimal.step, penalty = pen$penalty)
rs <- lapply(val_dd_list, function(x){cbind(x[,1:2], RS = as.numeric(predict(fit, newdata = x[, -c(1,2)], newtime = x[,1], newstatus = x[,2], type = "lp")))})
cc <- data.frame(Cindex = sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('CoxBoost')
result <- rbind(result, cc)

## 3.15 CoxBoost + Enet----
set.seed(seed)
pen <- optimCoxBoostPenalty(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                            trace = TRUE, start.penalty = 500, parallel = T)
cv.res <- cv.CoxBoost(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                      maxstepno = 500, K = 10, type = "verweij", penalty = pen$penalty)
fit <- CoxBoost(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                stepno = cv.res$optimal.step, penalty = pen$penalty)
rid <- as.data.frame(coef(fit))
rid$id <- rownames(rid)
rid <- rid[which(rid$`coef(fit)`!=0), "id"]

est_dd2 <- est_dd[, c('OS.time', 'OS', rid)]
val_dd_list2 <- lapply(val_dd_list, function(x){x[, c('OS.time', 'OS', rid)]})

x1 <- as.matrix(est_dd2[, rid])
x2 <- as.matrix(Surv(est_dd2$OS.time, est_dd2$OS))
for (alpha in seq(0.1, 0.9, 0.1)) {
  set.seed(seed)
  fit = cv.glmnet(x1, x2, family = "cox", alpha = alpha, nfolds = 10)
  rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, type = 'link', newx = as.matrix(x[, -c(1,2)]), s = fit$lambda.min)))})
  cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
    rownames_to_column('ID')
  cc$Model <- paste0('CoxBoost', ' + Enet', '[α=', alpha, ']')
  result <- rbind(result, cc)
}

## 3.16 CoxBoost + GBM----
set.seed(seed)
pen <- optimCoxBoostPenalty(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                            trace = TRUE, start.penalty = 500, parallel = T)
cv.res <- cv.CoxBoost(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                      maxstepno = 500, K= 10, type = "verweij", penalty = pen$penalty)
fit <- CoxBoost(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                stepno = cv.res$optimal.step, penalty = pen$penalty)
rid <- as.data.frame(coef(fit))
rid$id <- rownames(rid)
rid <- rid[which(rid$`coef(fit)`!=0), "id"]

est_dd2 <- est_dd[, c('OS.time', 'OS', rid)]
val_dd_list2 <- lapply(val_dd_list, function(x){x[,c('OS.time', 'OS', rid)]})

set.seed(seed)
fit <- gbm(formula = Surv(OS.time,OS)~., data = est_dd2, distribution = 'coxph',
           n.trees = 10000,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10, n.cores = 6)
# find index for number trees with minimum CV error
best <- which.min(fit$cv.error)
set.seed(seed)
fit <- gbm(formula = Surv(OS.time,OS)~., data = est_dd2, distribution = 'coxph',
           n.trees = best,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10,n.cores = 8)
rs <- lapply(val_dd_list2,function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, x, n.trees = best, type = 'link')))})
cc <- data.frame(Cindex=sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('CoxBoost + ', 'GBM')
result <- rbind(result, cc)

## 3.17 CoxBoost + Lasso----
set.seed(seed)
pen <- optimCoxBoostPenalty(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                            trace = TRUE, start.penalty = 500, parallel = T)
cv.res <- cv.CoxBoost(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                      maxstepno = 500, K = 10, type = "verweij", penalty = pen$penalty)
fit <- CoxBoost(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                stepno=cv.res$optimal.step, penalty=pen$penalty)
rid <- as.data.frame(coef(fit))
rid$id <- rownames(rid)
rid <- rid[which(rid$`coef(fit)` != 0), "id"]

est_dd2 <- est_dd[,c('OS.time', 'OS', rid)]
val_dd_list2 <- lapply(val_dd_list, function(x){x[, c('OS.time', 'OS', rid)]})

x1 <- as.matrix(est_dd2[, rid])
x2 <- as.matrix(Surv(est_dd2$OS.time, est_dd2$OS))
set.seed(seed)
fit = cv.glmnet(x1, x2,
                nfold = 10, #例文描述：10-fold cross-validation
                family = "binomial", alpha = 1,
                type.measure = "class")
rs <- lapply(val_dd_list2, function(x){cbind(x[,1:2], RS = as.numeric(predict(fit, type = 'response', newx = as.matrix(x[, -c(1,2)]), s = fit$lambda.min)))})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('CoxBoost + ', 'Lasso')
result <- rbind(result, cc)

## 3.18 CoxBoost + plsRcox----
set.seed(seed)
pen <- optimCoxBoostPenalty(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                            trace = TRUE, start.penalty = 500, parallel = T)
cv.res <- cv.CoxBoost(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                      maxstepno = 500, K = 10, type = "verweij", penalty = pen$penalty)
fit <- CoxBoost(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                stepno = cv.res$optimal.step, penalty = pen$penalty)
rid <- as.data.frame(coef(fit))
rid$id <- rownames(rid)
rid <- rid[which(rid$`coef(fit)` != 0), "id"]

est_dd2 <- est_dd[,c('OS.time', 'OS', rid)]
val_dd_list2 <- lapply(val_dd_list, function(x){x[,c('OS.time', 'OS', rid)]})

set.seed(seed)
cv.plsRcox.res = cv.plsRcox(list(x = est_dd2[,rid], time = est_dd2$OS.time, status = est_dd2$OS), nt = 10, verbose = FALSE)
fit <- plsRcox(est_dd2[, rid], time = est_dd2$OS.time, event = est_dd2$OS, nt = as.numeric(cv.plsRcox.res[5]))
rs <- lapply(val_dd_list2, function(x){cbind(x[,1:2], RS = as.numeric(predict(fit, type="lp", newdata = x[, -c(1,2)])))})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('CoxBoost + ', 'plsRcox')
result <- rbind(result, cc)

## 3.19 CoxBoost + Ridge----
set.seed(seed)
pen <- optimCoxBoostPenalty(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                            trace = TRUE, start.penalty = 500, parallel = T)
cv.res <- cv.CoxBoost(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                      maxstepno = 500, K=10, type="verweij", penalty = pen$penalty)
fit <- CoxBoost(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                stepno = cv.res$optimal.step, penalty = pen$penalty)
rid <- as.data.frame(coef(fit))
rid$id <- rownames(rid)
rid <- rid[which(rid$`coef(fit)` != 0), "id"]

est_dd2 <- est_dd[,c('OS.time', 'OS', rid)]
val_dd_list2 <- lapply(val_dd_list, function(x){x[,c('OS.time', 'OS', rid)]})

x1 <- as.matrix(est_dd2[, rid])
x2 <- as.matrix(Surv(est_dd2$OS.time, est_dd2$OS))
set.seed(seed)
fit = cv.glmnet(x1, x2,
                nfold=10, #例文描述：10-fold cross-validation
                family = "binomial", alpha = 0,
                type.measure = "class")
rs <- lapply(val_dd_list2, function(x){cbind(x[,1:2], RS = as.numeric(predict(fit, type = 'response', newx = as.matrix(x[, -c(1,2)]), s = fit$lambda.min)))})
cc <- data.frame(Cindex = sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('CoxBoost + ', 'Ridge')
result <- rbind(result, cc)

## 3.20 CoxBoost + StepCox----
set.seed(seed)
pen <- optimCoxBoostPenalty(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                            trace = TRUE, start.penalty = 500, parallel = T)
cv.res <- cv.CoxBoost(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                      maxstepno = 500, K = 10, type = "verweij", penalty = pen$penalty)
fit <- CoxBoost(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                stepno = cv.res$optimal.step, penalty = pen$penalty)
rid <- as.data.frame(coef(fit))
rid$id <- rownames(rid)
rid <- rid[which(rid$`coef(fit)` != 0), "id"]

est_dd2 <- est_dd[,c('OS.time', 'OS', rid)]
val_dd_list2 <- lapply(val_dd_list, function(x){x[,c('OS.time', 'OS', rid)]})

for (direction in c("both", "backward", "forward")) {
  fit <- step(coxph(Surv(OS.time,OS)~., est_dd2), direction = direction)
  rs <- lapply(val_dd_list2, function(x){cbind(x[,1:2], RS = predict(fit, type = 'risk', newdata = x))})
  cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
    rownames_to_column('ID')
  cc$Model <- paste0('CoxBoost + ', 'StepCox', '[', direction, ']')
  result <- rbind(result, cc)
}

## 3.21 CoxBoost + SuperPC----
set.seed(seed)
pen <- optimCoxBoostPenalty(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                            trace = TRUE, start.penalty = 500, parallel = T)
cv.res <- cv.CoxBoost(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[,-c(1,2)]),
                      maxstepno = 500, K= 10, type = "verweij", penalty = pen$penalty)
fit <- CoxBoost(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[,-c(1,2)]),
                stepno = cv.res$optimal.step, penalty = pen$penalty)
rid <- as.data.frame(coef(fit))
rid$id <- rownames(rid)
rid <- rid[which(rid$`coef(fit)` != 0), "id"]

est_dd2 <- est_dd[,c('OS.time', 'OS', rid)]
val_dd_list2 <- lapply(val_dd_list, function(x){x[,c('OS.time', 'OS', rid)]})

data <- list(x = t(est_dd2[, -c(1,2)]), y = est_dd2$OS.time, censoring.status = est_dd2$OS,
             featurenames = colnames(est_dd2)[-c(1,2)])
set.seed(seed)
fit <- superpc.train(data = data, type = 'survival', s0.perc = 0.5) #default
cv.fit <- superpc.cv(fit, data, n.threshold = 20, #default
                     n.fold = 10,
                     n.components = 3,
                     min.features = 5,
                     max.features = nrow(data$x),
                     compute.fullcv = TRUE,
                     compute.preval =TRUE)
rs <- lapply(val_dd_list2, function(w){
  test <- list(x=t(w[, -c(1,2)]), y = w$OS.time, censoring.status = w$OS, featurenames = colnames(w)[-c(1,2)])
  ff <- superpc.predict(fit, data, test, threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])], n.components = 1)
  rr <- as.numeric(ff$v.pred)
  rr2 <- cbind(w[,1:2], RS = rr)
  return(rr2)
})
cc <- data.frame(Cindex = sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('CoxBoost + ', 'SuperPC')
result <- rbind(result, cc)

## 3.22 CoxBoost + survival-SVM----
set.seed(seed)
pen <- optimCoxBoostPenalty(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                            trace = TRUE, start.penalty = 500, parallel = T)
cv.res <- cv.CoxBoost(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                      maxstepno = 500, K = 10, type = "verweij", penalty = pen$penalty)
fit <- CoxBoost(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[,-c(1,2)]),
                stepno = cv.res$optimal.step, penalty = pen$penalty)
rid <- as.data.frame(coef(fit))
rid$id <- rownames(rid)
rid <- rid[which(rid$`coef(fit)` != 0), "id"]

est_dd2 <- est_dd[, c('OS.time', 'OS', rid)]
val_dd_list2 <- lapply(val_dd_list, function(x){x[, c('OS.time', 'OS', rid)]})

fit = survivalsvm(Surv(OS.time, OS)~., data = est_dd2, gamma.mu = 1)
rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, x)$predicted))})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('CoxBoost + ', 'survival-SVM')
result <- rbind(result, cc)

## 3.23 plsRcox----
set.seed(seed)
cv.plsRcox.res = cv.plsRcox(list(x = est_dd[,pre_var], time = est_dd$OS.time, status = est_dd$OS), nt = 10, verbose = FALSE)
fit <- plsRcox(est_dd[,pre_var], time = est_dd$OS.time, event = est_dd$OS, nt = as.numeric(cv.plsRcox.res[5]))
rs <- lapply(val_dd_list, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit,type = "lp", newdata = x[, -c(1, 2)])))})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time,OS) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('plsRcox')
result <- rbind(result, cc)

## 3.24 SuperPC----
data <- list(x = t(est_dd[, -c(1,2)]), y = est_dd$OS.time, censoring.status = est_dd$OS, featurenames = colnames(est_dd)[-c(1, 2)])
set.seed(seed) 
fit <- superpc.train(data = data,type = 'survival', s0.perc = 0.5) #default
cv.fit <- superpc.cv(fit, data, n.threshold = 20, #default
                     n.fold = 10,
                     n.components = 3,
                     min.features = 5,
                     max.features = nrow(data$x),
                     compute.fullcv = TRUE,
                     compute.preval = TRUE)
rs <- lapply(val_dd_list, function(w){
  test <- list(x = t(w[,-c(1,2)]), y = w$OS.time, censoring.status = w$OS, featurenames = colnames(w)[-c(1,2)])
  ff <- superpc.predict(fit, data, test, threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])], n.components = 1)
  rr <- as.numeric(ff$v.pred)
  rr2 <- cbind(w[,1:2], RS = rr)
  return(rr2)
})
cc <- data.frame(Cindex = sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('SuperPC')
result <- rbind(result, cc)

## 3.25 GBM----
set.seed(seed)
fit <- gbm(formula = Surv(OS.time,OS)~., data = est_dd, distribution = 'coxph',
           n.trees = 10000,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10, n.cores = 6)
# find index for number trees with minimum CV error
best <- which.min(fit$cv.error)
set.seed(seed)
fit <- gbm(formula = Surv(OS.time, OS)~., data = est_dd, distribution = 'coxph',
           n.trees = best,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10, n.cores = 8)
rs <- lapply(val_dd_list,function(x){cbind(x[,1:2], RS = as.numeric(predict(fit, x, n.trees = best, type = 'link')))})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('GBM')
result <- rbind(result, cc)

## 3.26 survival-SVM----
fit = survivalsvm(Surv(OS.time,OS)~., data = est_dd, gamma.mu = 0.25, type = "vanbelle1",
                  diff.meth = "makediff1", opt.meth = "quadprog", kernel = "add_kernel")
rs <- lapply(val_dd_list, function(x){cbind(x[,1:2], RS = as.numeric(predict(fit, x)$predicted))})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('survival-SVM')
result <- rbind(result, cc)

## 3.27 Ridge----
x1 <- as.matrix(est_dd[, pre_var])
x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
set.seed(seed)
#fit = glmnet(x1, x2, family = "binomial", alpha = 0, lambda = NULL)
fit = cv.glmnet(x1, x2,
                nfold = 10, #例文描述：10-fold cross-validation
                family = "cox", alpha = 0
)

rs <- lapply(val_dd_list, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, type = 'response', newx = as.matrix(x[, -c(1,2)]), s = fit$lambda.min)))})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time,OS) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('Ridge')
result <- rbind(result, cc)

## 3.28 Lasso----
x1 <- as.matrix(est_dd[, pre_var])
x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
set.seed(seed)
fit = cv.glmnet(x1, x2,
                nfold = 10, #例文描述：10-fold cross-validation
                family = "cox", alpha = 1)
rs <- lapply(val_dd_list, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, type = 'response', newx = as.matrix(x[, -c(1,2)]), s = fit$lambda.min)))})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('Lasso')
result <- rbind(result, cc)

## 3.29 Lasso + CoxBoost----
x1 <- as.matrix(est_dd[, pre_var])
x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
set.seed(seed)
fit = cv.glmnet(x1, x2,
                nfold = 10, #例文描述：10-fold cross-validation
                family = "cox", alpha = 1)
fit$lambda.min
myCoefs <- coef(fit, s = "lambda.min");
rid <- myCoefs@Dimnames[[1]][which(myCoefs != 0 )]
rid <- rid[-1]

est_dd2 <- est_dd[, c('OS.time', 'OS', rid)]
val_dd_list2 <- lapply(val_dd_list, function(x){x[, c('OS.time', 'OS', rid)]})

set.seed(seed)
pen <- optimCoxBoostPenalty(est_dd2[, 'OS.time'], est_dd2[, 'OS'], as.matrix(est_dd2[, -c(1,2)]),
                            trace = TRUE, start.penalty = 500, parallel = T)
cv.res <- cv.CoxBoost(est_dd2[, 'OS.time'], est_dd2[, 'OS'], as.matrix(est_dd2[, -c(1,2)]),
                      maxstepno = 500, K = 10, type = "verweij", penalty = pen$penalty)
fit <- CoxBoost(est_dd2[, 'OS.time'], est_dd2[, 'OS'], as.matrix(est_dd2[, -c(1,2)]),
                stepno = cv.res$optimal.step, penalty = pen$penalty)
rs <- lapply(val_dd_list2, function(x){cbind(x[,1:2], RS = as.numeric(predict(fit, newdata = x[,-c(1,2)], newtime = x[,1], newstatus = x[,2], type = "lp")))})
cc <- data.frame(Cindex=sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('Lasso + CoxBoost')
result <- rbind(result, cc)

## 3.29 Lasso + GBM----
x1 <- as.matrix(est_dd[, pre_var])
x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
set.seed(seed)
fit = cv.glmnet(x1, x2,
                nfold = 10, #例文描述：10-fold cross-validation
                family = "binomial", alpha = 1,
                type.measure = "class")
fit$lambda.min
myCoefs <- coef(fit, s = "lambda.min");
rid <- myCoefs@Dimnames[[1]][which(myCoefs != 0 )]
rid <- rid[-1]

est_dd2 <- est_dd[,c('OS.time', 'OS', rid)]
val_dd_list2 <- lapply(val_dd_list, function(x){x[,c('OS.time', 'OS', rid)]})

set.seed(seed)
fit <- gbm(formula = Surv(OS.time,OS)~., data = est_dd2, distribution = 'coxph',
           n.trees = 10000,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10, n.cores = 6)
# find index for number trees with minimum CV error
best <- which.min(fit$cv.error)
set.seed(seed)
fit <- gbm(formula = Surv(OS.time,OS)~., data = est_dd2, distribution = 'coxph',
           n.trees = best,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10, n.cores = 8)
rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, x, n.trees = best, type = 'link')))})
cc <- data.frame(Cindex = sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('Lasso + ', 'GBM')
result <- rbind(result, cc)

## 3.30 Lasso + plsRcox----
x1 <- as.matrix(est_dd[, pre_var])
x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
set.seed(seed)
fit = cv.glmnet(x1, x2,
                nfold = 10, #例文描述：10-fold cross-validation
                family = "binomial", alpha = 1,
                type.measure = "class")
fit$lambda.min
myCoefs <- coef(fit, s = "lambda.min");
rid <- myCoefs@Dimnames[[1]][which(myCoefs != 0 )]
rid <- rid[-1]

est_dd2 <- est_dd[,c('OS.time', 'OS', rid)]
val_dd_list2 <- lapply(val_dd_list, function(x){x[,c('OS.time', 'OS', rid)]})

set.seed(seed)
cv.plsRcox.res = cv.plsRcox(list(x = est_dd2[, rid], time = est_dd2$OS.time, status = est_dd2$OS), nt = 10, verbose = FALSE)
fit <- plsRcox(est_dd2[, rid], time = est_dd2$OS.time, event = est_dd2$OS, nt = as.numeric(cv.plsRcox.res[5]))
rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, type = "lp", newdata = x[,-c(1,2)])))})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('Lasso + ', 'plsRcox')
result <- rbind(result, cc)

## 3.31 Lasso + RSF----
x1 <- as.matrix(est_dd[, pre_var])
x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
set.seed(seed)
fit = cv.glmnet(x1, x2,
                nfold = 10, #例文描述：10-fold cross-validation
                family = "binomial", alpha = 1,
                type.measure = "class")
fit$lambda.min
myCoefs <- coef(fit, s = "lambda.min");
rid <- myCoefs@Dimnames[[1]][which(myCoefs != 0 )]
rid<-rid[-1]

est_dd2 <- est_dd[,c('OS.time', 'OS', rid)]
val_dd_list2 <- lapply(val_dd_list, function(x){x[, c('OS.time', 'OS', rid)]})

set.seed(seed)
fit <- rfsrc(Surv(OS.time,OS)~., data = est_dd2,
             ntree = 1000, nodesize = rf_nodesize, ##该值建议多调整
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)
rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = predict(fit, newdata = x)$predicted)})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('Lasso', ' + RSF')
result <- rbind(result, cc)

## 3.32 Lasso + stepcox----
x1 <- as.matrix(est_dd[, pre_var])
x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
set.seed(seed)
fit = cv.glmnet(x1, x2,
                nfold = 10, #例文描述：10-fold cross-validation
                family = "binomial", alpha = 1,
                type.measure = "class")
fit$lambda.min
myCoefs <- coef(fit, s = "lambda.min");
rid <- myCoefs@Dimnames[[1]][which(myCoefs != 0 )]
rid <- rid[-1]

est_dd2 <- est_dd[, c('OS.time', 'OS', rid)]
val_dd_list2 <- lapply(val_dd_list, function(x){x[, c('OS.time', 'OS', rid)]})

for (direction in c("both", "backward", "forward")) {
  fit <- step(coxph(Surv(OS.time,OS)~., est_dd2), direction = direction)
  rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = predict(fit, type = 'risk', newdata = x))})
  cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
    rownames_to_column('ID')
  cc$Model <- paste0('Lasso + ', 'StepCox', '[', direction, ']')
  result <- rbind(result, cc)
}

## 3.33 Lasso + SuperPC----
x1 <- as.matrix(est_dd[, pre_var])
x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
set.seed(seed)
fit = cv.glmnet(x1, x2,
                nfold = 10, #例文描述：10-fold cross-validation
                family = "binomial", alpha = 1,
                type.measure = "class")
fit$lambda.min
myCoefs <- coef(fit, s = "lambda.min");
rid <- myCoefs@Dimnames[[1]][which(myCoefs != 0 )]
rid <- rid[-1]

est_dd2 <- est_dd[,c('OS.time', 'OS', rid)]
val_dd_list2 <- lapply(val_dd_list, function(x){x[, c('OS.time', 'OS', rid)]})

data <- list(x = t(est_dd2[,-c(1,2)]), y = est_dd2$OS.time, censoring.status = est_dd2$OS,
             featurenames = colnames(est_dd2)[-c(1,2)])
set.seed(seed)
fit <- superpc.train(data = data,type = 'survival', s0.perc = 0.5) #default
cv.fit <- superpc.cv(fit,data,n.threshold = 20, #default
                     n.fold = 10,
                     n.components = 3,
                     min.features = 5,
                     max.features = nrow(data$x),
                     compute.fullcv = TRUE,
                     compute.preval = TRUE)
rs <- lapply(val_dd_list2, function(w){
  test <- list(x = t(w[,-c(1,2)]), y = w$OS.time, censoring.status = w$OS, featurenames = colnames(w)[-c(1,2)])
  ff <- superpc.predict(fit, data, test, threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])], n.components = 1)
  rr <- as.numeric(ff$v.pred)
  rr2 <- cbind(w[, 1:2], RS = rr)
  return(rr2)
})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('Lasso + ', 'SuperPC')
result <- rbind(result, cc)

## 3.34 Lasso + survival-SVM----
x1 <- as.matrix(est_dd[, pre_var])
x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
set.seed(seed)
fit = cv.glmnet(x1, x2,
                nfold = 10, #例文描述：10-fold cross-validation
                family = "binomial", alpha = 1,
                type.measure = "class")
fit$lambda.min
myCoefs <- coef(fit, s = "lambda.min");
rid <- myCoefs@Dimnames[[1]][which(myCoefs != 0 )]
rid <- rid[-1]

est_dd2 <- est_dd[,c('OS.time', 'OS', rid)]
val_dd_list2 <- lapply(val_dd_list, function(x){x[,c('OS.time', 'OS', rid)]})

fit = survivalsvm(Surv(OS.time,OS)~., data = est_dd2, gamma.mu = 1)
rs <- lapply(val_dd_list2, function(x){cbind(x[,1:2], RS = as.numeric(predict(fit, x)$predicted))})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('Lasso + ', 'survival-SVM')
result <- rbind(result, cc)

## 3.35 汇总结果----
result2 <- result
###将结果的长数据转换为宽数据
dd2 <- pivot_wider(result2, names_from = 'ID', values_from = 'Cindex') %>% as.data.frame()
#将C指数定义为数值型
dd2[,-1] <- apply(dd2[,-1], 2, as.numeric)
#求每个模型的C指数在所有数据集的均值
dd2$All <- apply(dd2[,2:ncol(dd2)], 1, mean)
#求每个模型的C指数在GEO验证集的均值
dd2$GEO <- apply(dd2[,3:11], 1, mean)
###查看每个模型的C指数
head(dd2)

## 3.36 热图绘制----
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(ggsci)
# 根据C指数排序
dd2 <- dd2[order(dd2$All, decreasing = T),]
# 绘制所有数据集的C指数热图
dt <- dd2[, 2:11]
rownames(dt) <- dd2$Model

col_ha <- HeatmapAnnotation(which = "col", Cohort = colnames(dt),
                            #annotation_name_gp = gpar(fontsize = 9),
                            #annotation_name_side = "left",
                            col = list(Cohort=c("TCGA_LUAD" = pal_npg("nrc")(10)[1],
                                                "GSE14814" = pal_npg("nrc")(10)[2],
                                                "GSE29016" = pal_npg("nrc")(10)[3],
                                                "GSE30219" = pal_npg("nrc")(10)[4],
                                                "GSE31210" = pal_npg("nrc")(10)[5],
                                                "GSE37745" = pal_npg("nrc")(10)[6],
                                                "GSE42127" = pal_npg("nrc")(10)[7],
                                                "GSE50081" = pal_npg("nrc")(10)[8],
                                                "GSE68465" = pal_npg("nrc")(10)[9],
                                                "GSE72094" = pal_npg("nrc")(10)[10])),
                            annotation_legend_param = list(Cohort=list(title = "Cohort",
                                                                       title_position = "topleft",
                                                                       title_gp = gpar(fontsize = 12),
                                                                       labels_rot = 0,
                                                                       legend_height = unit(1,"cm"),
                                                                       legend_width = unit(5,"mm"),
                                                                       labels_gp = gpar(fontsize = 8))))

# 行注释
row_ha <- rowAnnotation('Mean Cindex' = anno_barplot(round(rowMeans(dt), 3), bar_width = 1, add_numbers = T,
                                                     labels = c("Mean Cindex"), height = unit(1, "mm"),
                                                     gp = gpar(col = "white", fill = "skyblue1"), numbers_gp = gpar(fontsize = 9),
                                                     axis_param = list(at = c(0, 0.5, 1),
                                                                       labels = c("0", "0.5", "1")),
                                                     width = unit(2.5, "cm")),
                        annotation_name_side = "bottom",
                        annotation_name_gp = gpar(fontsize = 8, fontface = "bold", angle = 90))

# 自定义图形，主要是热图右侧的条形图
cell_fun <- function(j, i, x, y, width, height, fill) {
  grid.text(
    #round(dt[i, j], 3), 
    #format(dt[i, j], nsmall = 3),  # 使用 format 强制显示三位小数
    sprintf("%.3f", dt[i, j]),  # 使用 sprintf 格式化为 3 位小数
    x, y,
    gp = gpar(
      fontsize = 9
    ))
}

# 画出热图
heatmap <- Heatmap(dt, name = " ", #图例标题
                   #参照color_mapping_legend函数设置图例
                   heatmap_legend_param = list(title="",title_position = "topleft", labels_rot = 0,
                                               legend_height = unit(4,"cm"),
                                               legend_width = unit(5,"mm"),
                                               labels_gp = gpar(fontsize = 12)),
                   border = TRUE,
                   column_split = c(colnames(dt)),
                   column_gap = unit(1.5, "mm"),
                   show_column_names = F,
                   show_row_names = T,
                   col = colorRamp2(c(0.5,0.6,0.7), c("#094687", "#DDEAF3", "#79C9C7")), # 选择颜色
                   column_title ="", # 列标题
                   #row_title ="Intersect Gene",
                   #column_title_side = "top", # 列标题位置
                   #column_title_rot = 45, 
                   height = unit(4, "cm"), # 调整热图高度
                   width = unit(12, "cm"),
                   row_title_side = "left",
                   row_title_rot = 90, # 旋转方向
                   column_title_gp = gpar(fontsize = 0), # 颜色，字体，大小
                   #row_title_gp = gpar(fontsize = 15, fontface = "bold",col = "black"),
                   cluster_columns = F,
                   cluster_rows = F,
                   column_order = c(colnames(dt)),
                   show_row_dend = F, # 是否显示聚类树
                   cell_fun = cell_fun,
                   top_annotation = col_ha,
                   right_annotation = row_ha
)
heatmap #宽10 高4

## 3.37 保存结果----
save(result, file = "12.Prognosis_model/model_Cindex_results.RData")
save(dd2, file = "12.Prognosis_model/model_mean-Cindex_results.RData")


# # 4.自己建模
# ## 4.1 RSF
# # 设置随机种子
# set.seed(123)
# 
# # 初始模型构建
# fit <- rfsrc(Surv(OS.time, OS) ~ ., data = est_dd, 
#              ntree = 1000,       # 树的数量
#              nodesize = 15,     # 每个叶节点的最小样本数
#              splitrule = "logrank", 
#              importance = TRUE, 
#              proximity = TRUE)
# print(fit)  # 查看模型结果
# 
# # 对每个 GEO 数据集进行验证，计算 C-index
# cindex_results <- lapply(val_data_list, function(val_data) {
#   # 预测验证集数据
#   pred <- predict(fit, newdata = val_data)
#   
#   # 计算验证集的 C-index
#   cindex <- as.numeric(summary(coxph(Surv(OS.time, OS) ~ pred$predicted, data = val_data))$concordance[1])
#   return(cindex)
# })
# # 打印每个验证集上的 C-index
# names(cindex_results) <- c("TCGA-LUAD", paste0("GEO_", 1:9))
# print(cindex_results)
# 
# # 设置参数搜索范围
# ntree_values <- c(500, 1000, 1500, 2000)    # 树的数量
# nodesize_values <- c(5, 10, 15, 20)   # 每个叶节点的最小样本数
# 
# # 存储结果
# tuning_results <- data.frame(ntree = integer(), nodesize = integer(), avg_cindex = numeric())
# 
# # 网格搜索
# for (ntree in ntree_values) {
#   for (nodesize in nodesize_values) {
#     set.seed(123)
#     fit <- rfsrc(Surv(OS.time, OS) ~ ., data = est_dd, 
#                  ntree = ntree, nodesize = nodesize, 
#                  splitrule = "logrank", 
#                  importance = TRUE, 
#                  proximity = TRUE)
#     
#     # 对每个验证集计算 C-index 并求平均
#     cindex_list <- sapply(val_data_list, function(val_data) {
#       pred <- predict(fit, newdata = val_data)
#       as.numeric(summary(coxph(Surv(OS.time, OS) ~ pred$predicted, data = val_data))$concordance[1])
#     })
#     
#     avg_cindex <- mean(cindex_list)
#     
#     # 记录结果
#     tuning_results <- rbind(tuning_results, data.frame(ntree = ntree, nodesize = nodesize, avg_cindex = avg_cindex))
#   }
# }
# 
# # 查看调参结果
# print(tuning_results)
# 
# # 找到最佳参数组合
# best_params <- tuning_results[which.max(tuning_results$avg_cindex), ]
# print(best_params)
# 
# # 使用最优参数训练模型
# set.seed(123)
# best_fit <- rfsrc(Surv(OS.time, OS) ~ ., data = est_dd, 
#                   ntree = best_params$ntree, 
#                   nodesize = best_params$nodesize, 
#                   splitrule = "logrank", 
#                   importance = TRUE, 
#                   proximity = TRUE)
# 
# # 在每个 GEO 验证集上计算最终 C-index
# final_cindex_results <- sapply(val_data_list, function(val_data) {
#   pred <- predict(best_fit, newdata = val_data)
#   as.numeric(summary(coxph(Surv(OS.time, OS) ~ pred$predicted, data = val_data))$concordance[1])
# })
# 
# # 输出最终 C-index
# names(final_cindex_results) <- c("TCGA-LUAD", paste0("GEO_", 1:9))
# print(final_cindex_results)
# 
# # 提取变量重要性
# var_importance <- data.frame(Variable = names(best_fit$importance), Importance = best_fit$importance)
# 
# # 绘制变量重要性图
# ggplot(var_importance, aes(x = reorder(Variable, Importance), y = Importance)) +
#   geom_bar(stat = "identity") +
#   coord_flip() +
#   labs(title = "Variable Importance in Random Survival Forest",
#        x = "Variable", y = "Importance")


# 4.使用RSF模型计算的risk score进行风险分层、K-M和timeROC作图----
# 设置种子数
seed <- 123

 ## 4.1 重复3.1的RSF得到fit模型----
rf_nodesize <- 15
set.seed(seed)
fit <- rfsrc(Surv(OS.time,OS)~., data = est_dd, 
             ntree = 1000, nodesize = rf_nodesize,  #该值建议多调整
             splitrule = 'logrank', 
             importance = T, 
             proximity = T, 
             forest = T, 
             seed = seed)
print(fit)
# fit$xvar.names
# importance <- as.data.frame(fit$importance)
rs <- lapply(val_dd_list, function(x){cbind(x[, 1:2], RS  = predict(fit, newdata = x)$predicted)})

save(fit, file = "12.Prognosis_model/RSFfit.RData")
save(rs, file = "12.Prognosis_model/RSF_riskScore.RData")

## 4.2 每个队列做K-M图----
### 4.2.1 TCGA-LUAD----
cohort <- rs[["TCGA_LUAD"]]
cohort$Risk_group <- ifelse(cohort$RS > median(cohort$RS), "High-Risk", "Low-Risk")
table(cohort$Risk_group)

# 创建生存对象
surv_object <- Surv(time = cohort$OS.time, event = cohort$OS)

# 根据分组拟合生存曲线
kmfit <- survfit(surv_object ~ Risk_group, data = cohort)

# 绘制生存曲线
high_color <- pal_npg("nrc")(10)[1]
low_color <- pal_npg("nrc")(10)[2]

km_plot <- ggsurvplot(kmfit, data = cohort,
                      pval = T, # 在图上添加log rank检验的p值
                      pval.method = T,
                      pval.coord = c(13.5, 0.9), # 调整P值位置，X=10，Y=0.1
                      pval.method.coord = c(9, 0.9), # 调整P值方法的位置，X=3，Y=0.85
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

pdf(file = "12.Prognosis_model/TCGA-LUAD_RiskScore_KMplot.pdf",height = 5.5, width = 5, onefile = F)
km_plot
dev.off()

### 4.2.2 GSE14814----
cohort <- rs[["GSE14814"]]
cohort$Risk_group <- ifelse(cohort$RS > median(cohort$RS), "High-Risk", "Low-Risk")
table(cohort$Risk_group)

# 创建生存对象
surv_object <- Surv(time = cohort$OS.time, event = cohort$OS)

# 根据分组拟合生存曲线
kmfit <- survfit(surv_object ~ Risk_group, data = cohort)

# 绘制生存曲线
km_plot <- ggsurvplot(kmfit, data = cohort,
                      pval = T, # 在图上添加log rank检验的p值
                      pval.method = T,
                      pval.coord = c(8, 0.9), # 调整P值位置，X=10，Y=0.1
                      pval.method.coord = c(5.8, 0.9), # 调整P值方法的位置，X=3，Y=0.85
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

pdf(file = "12.Prognosis_model/GSE14814_RiskScore_KMplot.pdf",height = 5.5, width = 5, onefile = F)
km_plot
dev.off()

### 4.2.3 GSE29016----
cohort <- rs[["GSE29016"]]
cohort$Risk_group <- ifelse(cohort$RS > median(cohort$RS), "High-Risk", "Low-Risk")
table(cohort$Risk_group)

# 创建生存对象
surv_object <- Surv(time = cohort$OS.time, event = cohort$OS)

# 根据分组拟合生存曲线
kmfit <- survfit(surv_object ~ Risk_group, data = cohort)

# 绘制生存曲线
km_plot <- ggsurvplot(kmfit, data = cohort,
                      pval = T, # 在图上添加log rank检验的p值
                      pval.method = T,
                      pval.coord = c(11, 0.9), # 调整P值位置，X=10，Y=0.1
                      pval.method.coord = c(7.5, 0.9), # 调整P值方法的位置，X=3，Y=0.85
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

pdf(file = "12.Prognosis_model/GSE29016_RiskScore_KMplot.pdf",height = 5.5, width = 5, onefile = F)
km_plot
dev.off()

### 4.2.4 GSE30219----
cohort <- rs[["GSE30219"]]
cohort$Risk_group <- ifelse(cohort$RS > median(cohort$RS), "High-Risk", "Low-Risk")
table(cohort$Risk_group)

# 创建生存对象
surv_object <- Surv(time = cohort$OS.time, event = cohort$OS)

# 根据分组拟合生存曲线
kmfit <- survfit(surv_object ~ Risk_group, data = cohort)

# 绘制生存曲线
km_plot <- ggsurvplot(kmfit, data = cohort,
                      pval = T, # 在图上添加log rank检验的p值
                      pval.method = T,
                      pval.coord = c(15, 0.9), # 调整P值位置，X=10，Y=0.1
                      pval.method.coord = c(10, 0.9), # 调整P值方法的位置，X=3，Y=0.85
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

pdf(file = "12.Prognosis_model/GSE30219_RiskScore_KMplot.pdf",height = 5.5, width = 5, onefile = F)
km_plot
dev.off()

### 4.2.5 GSE31210----
cohort <- rs[["GSE31210"]]
cohort$Risk_group <- ifelse(cohort$RS > median(cohort$RS), "High-Risk", "Low-Risk")
table(cohort$Risk_group)

# 创建生存对象
surv_object <- Surv(time = cohort$OS.time, event = cohort$OS)

# 根据分组拟合生存曲线
kmfit <- survfit(surv_object ~ Risk_group, data = cohort)

# 绘制生存曲线
km_plot <- ggsurvplot(kmfit, data = cohort,
                      pval = T, # 在图上添加log rank检验的p值
                      pval.method = T,
                      pval.coord = c(7.5, 0.4), # 调整P值位置，X=10，Y=0.1
                      pval.method.coord = c(5, 0.4), # 调整P值方法的位置，X=3，Y=0.85
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

pdf(file = "12.Prognosis_model/GSE31210_RiskScore_KMplot.pdf",height = 5.5, width = 5, onefile = F)
km_plot
dev.off()

### 4.2.6 GSE37745----
cohort <- rs[["GSE37745"]]
cohort$Risk_group <- ifelse(cohort$RS > median(cohort$RS), "High-Risk", "Low-Risk")
table(cohort$Risk_group)

# 创建生存对象
surv_object <- Surv(time = cohort$OS.time, event = cohort$OS)

# 根据分组拟合生存曲线
kmfit <- survfit(surv_object ~ Risk_group, data = cohort)

# 绘制生存曲线
km_plot <- ggsurvplot(kmfit, data = cohort,
                      pval = T, # 在图上添加log rank检验的p值
                      pval.method = T,
                      pval.coord = c(12, 0.9), # 调整P值位置，X=10，Y=0.1
                      pval.method.coord = c(8, 0.9), # 调整P值方法的位置，X=3，Y=0.85
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

pdf(file = "12.Prognosis_model/GSE37745_RiskScore_KMplot.pdf",height = 5.5, width = 5, onefile = F)
km_plot
dev.off()

### 4.2.7 GSE42127----
cohort <- rs[["GSE42127"]]
cohort$Risk_group <- ifelse(cohort$RS > median(cohort$RS), "High-Risk", "Low-Risk")
table(cohort$Risk_group)

# 创建生存对象
surv_object <- Surv(time = cohort$OS.time, event = cohort$OS)

# 根据分组拟合生存曲线
kmfit <- survfit(surv_object ~ Risk_group, data = cohort)

# 绘制生存曲线
km_plot <- ggsurvplot(kmfit, data = cohort,
                      pval = T, # 在图上添加log rank检验的p值
                      pval.method = T,
                      pval.coord = c(3, 0.3), # 调整P值位置，X=10，Y=0.1
                      pval.method.coord = c(0.5, 0.3), # 调整P值方法的位置，X=3，Y=0.85
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

pdf(file = "12.Prognosis_model/GSE42127_RiskScore_KMplot.pdf",height = 5.5, width = 5, onefile = F)
km_plot
dev.off()

### 4.2.8 GSE50081----
cohort <- rs[["GSE50081"]]
cohort$Risk_group <- ifelse(cohort$RS > median(cohort$RS), "High-Risk", "Low-Risk")
table(cohort$Risk_group)

# 创建生存对象
surv_object <- Surv(time = cohort$OS.time, event = cohort$OS)

# 根据分组拟合生存曲线
kmfit <- survfit(surv_object ~ Risk_group, data = cohort)

# 绘制生存曲线
km_plot <- ggsurvplot(kmfit, data = cohort,
                      pval = T, # 在图上添加log rank检验的p值
                      pval.method = T,
                      pval.coord = c(8, 0.9), # 调整P值位置，X=10，Y=0.1
                      pval.method.coord = c(5.5, 0.9), # 调整P值方法的位置，X=3，Y=0.85
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

pdf(file = "12.Prognosis_model/GSE50081_RiskScore_KMplot.pdf",height = 5.5, width = 5, onefile = F)
km_plot
dev.off()

### 4.2.9 GSE68465----
cohort <- rs[["GSE68465"]]
cohort$Risk_group <- ifelse(cohort$RS > median(cohort$RS), "High-Risk", "Low-Risk")
table(cohort$Risk_group)

# 创建生存对象
surv_object <- Surv(time = cohort$OS.time, event = cohort$OS)

# 根据分组拟合生存曲线
kmfit <- survfit(surv_object ~ Risk_group, data = cohort)

# 绘制生存曲线
km_plot <- ggsurvplot(kmfit, data = cohort,
                      pval = T, # 在图上添加log rank检验的p值
                      pval.method = T,
                      pval.coord = c(12, 0.9), # 调整P值位置，X=10，Y=0.1
                      pval.method.coord = c(8, 0.9), # 调整P值方法的位置，X=3，Y=0.85
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

pdf(file = "12.Prognosis_model/GSE68465_RiskScore_KMplot.pdf",height = 5.5, width = 5, onefile = F)
km_plot
dev.off()

### 4.2.10 GSE72094----
cohort <- rs[["GSE72094"]]
cohort$Risk_group <- ifelse(cohort$RS > median(cohort$RS), "High-Risk", "Low-Risk")
table(cohort$Risk_group)

# 创建生存对象
surv_object <- Surv(time = cohort$OS.time, event = cohort$OS)

# 根据分组拟合生存曲线
kmfit <- survfit(surv_object ~ Risk_group, data = cohort)

# 绘制生存曲线
km_plot <- ggsurvplot(kmfit, data = cohort,
                      pval = T, # 在图上添加log rank检验的p值
                      pval.method = T,
                      pval.coord = c(4.5, 0.9), # 调整P值位置，X=10，Y=0.1
                      pval.method.coord = c(3, 0.9), # 调整P值方法的位置，X=3，Y=0.85
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

pdf(file = "12.Prognosis_model/GSE72094_RiskScore_KMplot.pdf",height = 5.5, width = 5, onefile = F)
km_plot
dev.off()

## 4.3 每个队列画时间依赖ROC----
library(timeROC)

### 4.3.1 TCGA-LUAD----
cohort <- rs[["TCGA_LUAD"]]

roc_1year <- timeROC(
  T = cohort$OS.time,              # 生存时间
  delta = cohort$OS,        # 生存状态
  marker = cohort$RS, # 预测分数
  cause = 1,                # 指定事件
  times = 1,                # 指定时间点 (1 年)
  iid = TRUE
)

roc_3year <- timeROC(
  T = cohort$OS.time,              # 生存时间
  delta = cohort$OS,        # 生存状态
  marker = cohort$RS, # 预测分数
  cause = 1,                # 指定事件
  times = 3,                # 指定时间点 (1 年)
  iid = TRUE
)

roc_5year <- timeROC(
  T = cohort$OS.time,              # 生存时间
  delta = cohort$OS,        # 生存状态
  marker = cohort$RS, # 预测分数
  cause = 1,                # 指定事件
  times = 5,                # 指定时间点 (1 年)
  iid = TRUE
)

print(roc_1year)
print(roc_3year)
print(roc_5year)

# 查看 AUC 值
auc_1year <- roc_1year$AUC[2]
auc_3year <- roc_3year$AUC[2]
auc_5year <- roc_5year$AUC[2]

# 绘制 ROC 曲线
plot(roc_1year, time = 1, col = pal_npg("nrc")(10)[1], title = "Time-dependent ROC", lwd = 2)
plot(roc_3year, time = 3, add = TRUE, col = pal_npg("nrc")(10)[2], lwd = 2)
plot(roc_5year, time = 5, add = TRUE, col = pal_npg("nrc")(10)[3], lwd = 2)

legend("bottomright", 
       legend = c(paste0("1 Year AUC: ", sprintf("%.3f", auc_1year)),
                  paste0("3 Year AUC: ", sprintf("%.3f", auc_3year)),
                  paste0("5 Year AUC: ", sprintf("%.3f", auc_5year))),
       col = c(pal_npg("nrc")(10)[1], pal_npg("nrc")(10)[2], pal_npg("nrc")(10)[3]), 
       lty = 1, lwd = 2) # 宽 5.5 高 5

### 4.3.2 GSE14814----
cohort <- rs[["GSE14814"]]

roc_1year <- timeROC(
  T = cohort$OS.time,              # 生存时间
  delta = cohort$OS,        # 生存状态
  marker = cohort$RS, # 预测分数
  cause = 1,                # 指定事件
  times = 1,                # 指定时间点 (1 年)
  iid = TRUE
)

roc_3year <- timeROC(
  T = cohort$OS.time,              # 生存时间
  delta = cohort$OS,        # 生存状态
  marker = cohort$RS, # 预测分数
  cause = 1,                # 指定事件
  times = 3,                # 指定时间点 (1 年)
  iid = TRUE
)

roc_5year <- timeROC(
  T = cohort$OS.time,              # 生存时间
  delta = cohort$OS,        # 生存状态
  marker = cohort$RS, # 预测分数
  cause = 1,                # 指定事件
  times = 5,                # 指定时间点 (1 年)
  iid = TRUE
)

print(roc_1year)
print(roc_3year)
print(roc_5year)

# 查看 AUC 值
auc_1year <- roc_1year$AUC[2]
auc_3year <- roc_3year$AUC[2]
auc_5year <- roc_5year$AUC[2]

# 绘制 ROC 曲线
plot(roc_1year, time = 1, col = pal_npg("nrc")(10)[1], title = "Time-dependent ROC", lwd = 2)
plot(roc_3year, time = 3, add = TRUE, col = pal_npg("nrc")(10)[2], lwd = 2)
plot(roc_5year, time = 5, add = TRUE, col = pal_npg("nrc")(10)[3], lwd = 2)

legend("bottomright", 
       legend = c(paste0("1 Year AUC: ", sprintf("%.3f", auc_1year)),
                  paste0("3 Year AUC: ", sprintf("%.3f", auc_3year)),
                  paste0("5 Year AUC: ", sprintf("%.3f", auc_5year))),
       col = c(pal_npg("nrc")(10)[1], pal_npg("nrc")(10)[2], pal_npg("nrc")(10)[3]), 
       lty = 1, lwd = 2) # 宽 5.5 高 5

### 4.3.3 GSE29016----
cohort <- rs[["GSE29016"]]

roc_1year <- timeROC(
  T = cohort$OS.time,              # 生存时间
  delta = cohort$OS,        # 生存状态
  marker = cohort$RS, # 预测分数
  cause = 1,                # 指定事件
  times = 1,                # 指定时间点 (1 年)
  iid = TRUE
)

roc_3year <- timeROC(
  T = cohort$OS.time,              # 生存时间
  delta = cohort$OS,        # 生存状态
  marker = cohort$RS, # 预测分数
  cause = 1,                # 指定事件
  times = 3,                # 指定时间点 (1 年)
  iid = TRUE
)

roc_5year <- timeROC(
  T = cohort$OS.time,              # 生存时间
  delta = cohort$OS,        # 生存状态
  marker = cohort$RS, # 预测分数
  cause = 1,                # 指定事件
  times = 5,                # 指定时间点 (1 年)
  iid = TRUE
)

print(roc_1year)
print(roc_3year)
print(roc_5year)

# 查看 AUC 值
auc_1year <- roc_1year$AUC[2]
auc_3year <- roc_3year$AUC[2]
auc_5year <- roc_5year$AUC[2]

# 绘制 ROC 曲线
plot(roc_1year, time = 1, col = pal_npg("nrc")(10)[1], title = "Time-dependent ROC", lwd = 2)
plot(roc_3year, time = 3, add = TRUE, col = pal_npg("nrc")(10)[2], lwd = 2)
plot(roc_5year, time = 5, add = TRUE, col = pal_npg("nrc")(10)[3], lwd = 2)

legend("bottomright", 
       legend = c(paste0("1 Year AUC: ", sprintf("%.3f", auc_1year)),
                  paste0("3 Year AUC: ", sprintf("%.3f", auc_3year)),
                  paste0("5 Year AUC: ", sprintf("%.3f", auc_5year))),
       col = c(pal_npg("nrc")(10)[1], pal_npg("nrc")(10)[2], pal_npg("nrc")(10)[3]), 
       lty = 1, lwd = 2) # 宽 5.5 高 5

### 4.3.4 GSE30219----
cohort <- rs[["GSE30219"]]

roc_1year <- timeROC(
  T = cohort$OS.time,              # 生存时间
  delta = cohort$OS,        # 生存状态
  marker = cohort$RS, # 预测分数
  cause = 1,                # 指定事件
  times = 1,                # 指定时间点 (1 年)
  iid = TRUE
)

roc_3year <- timeROC(
  T = cohort$OS.time,              # 生存时间
  delta = cohort$OS,        # 生存状态
  marker = cohort$RS, # 预测分数
  cause = 1,                # 指定事件
  times = 3,                # 指定时间点 (1 年)
  iid = TRUE
)

roc_5year <- timeROC(
  T = cohort$OS.time,              # 生存时间
  delta = cohort$OS,        # 生存状态
  marker = cohort$RS, # 预测分数
  cause = 1,                # 指定事件
  times = 5,                # 指定时间点 (1 年)
  iid = TRUE
)

print(roc_1year)
print(roc_3year)
print(roc_5year)

# 查看 AUC 值
auc_1year <- roc_1year$AUC[2]
auc_3year <- roc_3year$AUC[2]
auc_5year <- roc_5year$AUC[2]

# 绘制 ROC 曲线
plot(roc_1year, time = 1, col = pal_npg("nrc")(10)[1], title = "Time-dependent ROC", lwd = 2)
plot(roc_3year, time = 3, add = TRUE, col = pal_npg("nrc")(10)[2], lwd = 2)
plot(roc_5year, time = 5, add = TRUE, col = pal_npg("nrc")(10)[3], lwd = 2)

legend("bottomright", 
       legend = c(paste0("1 Year AUC: ", sprintf("%.3f", auc_1year)),
                  paste0("3 Year AUC: ", sprintf("%.3f", auc_3year)),
                  paste0("5 Year AUC: ", sprintf("%.3f", auc_5year))),
       col = c(pal_npg("nrc")(10)[1], pal_npg("nrc")(10)[2], pal_npg("nrc")(10)[3]), 
       lty = 1, lwd = 2) # 宽 5.5 高 5

### 4.3.5 GSE31210----
cohort <- rs[["GSE31210"]]

roc_1year <- timeROC(
  T = cohort$OS.time,              # 生存时间
  delta = cohort$OS,        # 生存状态
  marker = cohort$RS, # 预测分数
  cause = 1,                # 指定事件
  times = 1,                # 指定时间点 (1 年)
  iid = TRUE
)

roc_3year <- timeROC(
  T = cohort$OS.time,              # 生存时间
  delta = cohort$OS,        # 生存状态
  marker = cohort$RS, # 预测分数
  cause = 1,                # 指定事件
  times = 3,                # 指定时间点 (1 年)
  iid = TRUE
)

roc_5year <- timeROC(
  T = cohort$OS.time,              # 生存时间
  delta = cohort$OS,        # 生存状态
  marker = cohort$RS, # 预测分数
  cause = 1,                # 指定事件
  times = 5,                # 指定时间点 (1 年)
  iid = TRUE
)

print(roc_1year)
print(roc_3year)
print(roc_5year)

# 查看 AUC 值
auc_1year <- roc_1year$AUC[2]
auc_3year <- roc_3year$AUC[2]
auc_5year <- roc_5year$AUC[2]

# 绘制 ROC 曲线
plot(roc_1year, time = 1, col = pal_npg("nrc")(10)[1], title = "Time-dependent ROC", lwd = 2)
plot(roc_3year, time = 3, add = TRUE, col = pal_npg("nrc")(10)[2], lwd = 2)
plot(roc_5year, time = 5, add = TRUE, col = pal_npg("nrc")(10)[3], lwd = 2)

legend("bottomright", 
       legend = c(paste0("1 Year AUC: ", sprintf("%.3f", auc_1year)),
                  paste0("3 Year AUC: ", sprintf("%.3f", auc_3year)),
                  paste0("5 Year AUC: ", sprintf("%.3f", auc_5year))),
       col = c(pal_npg("nrc")(10)[1], pal_npg("nrc")(10)[2], pal_npg("nrc")(10)[3]), 
       lty = 1, lwd = 2) # 宽 5.5 高 5

### 4.3.6 GSE37745----
cohort <- rs[["GSE37745"]]

roc_1year <- timeROC(
  T = cohort$OS.time,              # 生存时间
  delta = cohort$OS,        # 生存状态
  marker = cohort$RS, # 预测分数
  cause = 1,                # 指定事件
  times = 1,                # 指定时间点 (1 年)
  iid = TRUE
)

roc_3year <- timeROC(
  T = cohort$OS.time,              # 生存时间
  delta = cohort$OS,        # 生存状态
  marker = cohort$RS, # 预测分数
  cause = 1,                # 指定事件
  times = 3,                # 指定时间点 (1 年)
  iid = TRUE
)

roc_5year <- timeROC(
  T = cohort$OS.time,              # 生存时间
  delta = cohort$OS,        # 生存状态
  marker = cohort$RS, # 预测分数
  cause = 1,                # 指定事件
  times = 5,                # 指定时间点 (1 年)
  iid = TRUE
)

print(roc_1year)
print(roc_3year)
print(roc_5year)

# 查看 AUC 值
auc_1year <- roc_1year$AUC[2]
auc_3year <- roc_3year$AUC[2]
auc_5year <- roc_5year$AUC[2]

# 绘制 ROC 曲线
plot(roc_1year, time = 1, col = pal_npg("nrc")(10)[1], title = "Time-dependent ROC", lwd = 2)
plot(roc_3year, time = 3, add = TRUE, col = pal_npg("nrc")(10)[2], lwd = 2)
plot(roc_5year, time = 5, add = TRUE, col = pal_npg("nrc")(10)[3], lwd = 2)

legend("bottomright", 
       legend = c(paste0("1 Year AUC: ", sprintf("%.3f", auc_1year)),
                  paste0("3 Year AUC: ", sprintf("%.3f", auc_3year)),
                  paste0("5 Year AUC: ", sprintf("%.3f", auc_5year))),
       col = c(pal_npg("nrc")(10)[1], pal_npg("nrc")(10)[2], pal_npg("nrc")(10)[3]), 
       lty = 1, lwd = 2) # 宽 5.5 高 5

### 4.3.7 GSE42127----
cohort <- rs[["GSE42127"]]

roc_1year <- timeROC(
  T = cohort$OS.time,              # 生存时间
  delta = cohort$OS,        # 生存状态
  marker = cohort$RS, # 预测分数
  cause = 1,                # 指定事件
  times = 1,                # 指定时间点 (1 年)
  iid = TRUE
)

roc_3year <- timeROC(
  T = cohort$OS.time,              # 生存时间
  delta = cohort$OS,        # 生存状态
  marker = cohort$RS, # 预测分数
  cause = 1,                # 指定事件
  times = 3,                # 指定时间点 (1 年)
  iid = TRUE
)

roc_5year <- timeROC(
  T = cohort$OS.time,              # 生存时间
  delta = cohort$OS,        # 生存状态
  marker = cohort$RS, # 预测分数
  cause = 1,                # 指定事件
  times = 5,                # 指定时间点 (1 年)
  iid = TRUE
)

print(roc_1year)
print(roc_3year)
print(roc_5year)

# 查看 AUC 值
auc_1year <- roc_1year$AUC[2]
auc_3year <- roc_3year$AUC[2]
auc_5year <- roc_5year$AUC[2]

# 绘制 ROC 曲线
plot(roc_1year, time = 1, col = pal_npg("nrc")(10)[1], title = "Time-dependent ROC", lwd = 2)
plot(roc_3year, time = 3, add = TRUE, col = pal_npg("nrc")(10)[2], lwd = 2)
plot(roc_5year, time = 5, add = TRUE, col = pal_npg("nrc")(10)[3], lwd = 2)

legend("bottomright", 
       legend = c(paste0("1 Year AUC: ", sprintf("%.3f", auc_1year)),
                  paste0("3 Year AUC: ", sprintf("%.3f", auc_3year)),
                  paste0("5 Year AUC: ", sprintf("%.3f", auc_5year))),
       col = c(pal_npg("nrc")(10)[1], pal_npg("nrc")(10)[2], pal_npg("nrc")(10)[3]), 
       lty = 1, lwd = 2) # 宽 5.5 高 5

### 4.3.8 GSE50081----
cohort <- rs[["GSE50081"]]

roc_1year <- timeROC(
  T = cohort$OS.time,              # 生存时间
  delta = cohort$OS,        # 生存状态
  marker = cohort$RS, # 预测分数
  cause = 1,                # 指定事件
  times = 1,                # 指定时间点 (1 年)
  iid = TRUE
)

roc_3year <- timeROC(
  T = cohort$OS.time,              # 生存时间
  delta = cohort$OS,        # 生存状态
  marker = cohort$RS, # 预测分数
  cause = 1,                # 指定事件
  times = 3,                # 指定时间点 (1 年)
  iid = TRUE
)

roc_5year <- timeROC(
  T = cohort$OS.time,              # 生存时间
  delta = cohort$OS,        # 生存状态
  marker = cohort$RS, # 预测分数
  cause = 1,                # 指定事件
  times = 5,                # 指定时间点 (1 年)
  iid = TRUE
)

print(roc_1year)
print(roc_3year)
print(roc_5year)

# 查看 AUC 值
auc_1year <- roc_1year$AUC[2]
auc_3year <- roc_3year$AUC[2]
auc_5year <- roc_5year$AUC[2]

# 绘制 ROC 曲线
plot(roc_1year, time = 1, col = pal_npg("nrc")(10)[1], title = "Time-dependent ROC", lwd = 2)
plot(roc_3year, time = 3, add = TRUE, col = pal_npg("nrc")(10)[2], lwd = 2)
plot(roc_5year, time = 5, add = TRUE, col = pal_npg("nrc")(10)[3], lwd = 2)

legend("bottomright", 
       legend = c(paste0("1 Year AUC: ", sprintf("%.3f", auc_1year)),
                  paste0("3 Year AUC: ", sprintf("%.3f", auc_3year)),
                  paste0("5 Year AUC: ", sprintf("%.3f", auc_5year))),
       col = c(pal_npg("nrc")(10)[1], pal_npg("nrc")(10)[2], pal_npg("nrc")(10)[3]), 
       lty = 1, lwd = 2) # 宽 5.5 高 5

### 4.3.9 GSE68465----
cohort <- rs[["GSE68465"]]

roc_1year <- timeROC(
  T = cohort$OS.time,              # 生存时间
  delta = cohort$OS,        # 生存状态
  marker = cohort$RS, # 预测分数
  cause = 1,                # 指定事件
  times = 1,                # 指定时间点 (1 年)
  iid = TRUE
)

roc_3year <- timeROC(
  T = cohort$OS.time,              # 生存时间
  delta = cohort$OS,        # 生存状态
  marker = cohort$RS, # 预测分数
  cause = 1,                # 指定事件
  times = 3,                # 指定时间点 (1 年)
  iid = TRUE
)

roc_5year <- timeROC(
  T = cohort$OS.time,              # 生存时间
  delta = cohort$OS,        # 生存状态
  marker = cohort$RS, # 预测分数
  cause = 1,                # 指定事件
  times = 5,                # 指定时间点 (1 年)
  iid = TRUE
)

print(roc_1year)
print(roc_3year)
print(roc_5year)

# 查看 AUC 值
auc_1year <- roc_1year$AUC[2]
auc_3year <- roc_3year$AUC[2]
auc_5year <- roc_5year$AUC[2]

# 绘制 ROC 曲线
plot(roc_1year, time = 1, col = pal_npg("nrc")(10)[1], title = "Time-dependent ROC", lwd = 2)
plot(roc_3year, time = 3, add = TRUE, col = pal_npg("nrc")(10)[2], lwd = 2)
plot(roc_5year, time = 5, add = TRUE, col = pal_npg("nrc")(10)[3], lwd = 2)

legend("bottomright", 
       legend = c(paste0("1 Year AUC: ", sprintf("%.3f", auc_1year)),
                  paste0("3 Year AUC: ", sprintf("%.3f", auc_3year)),
                  paste0("5 Year AUC: ", sprintf("%.3f", auc_5year))),
       col = c(pal_npg("nrc")(10)[1], pal_npg("nrc")(10)[2], pal_npg("nrc")(10)[3]), 
       lty = 1, lwd = 2) # 宽 5.5 高 5

### 4.3.10 GSE72094----
cohort <- rs[["GSE72094"]]

roc_1year <- timeROC(
  T = cohort$OS.time,              # 生存时间
  delta = cohort$OS,        # 生存状态
  marker = cohort$RS, # 预测分数
  cause = 1,                # 指定事件
  times = 1,                # 指定时间点 (1 年)
  iid = TRUE
)

roc_3year <- timeROC(
  T = cohort$OS.time,              # 生存时间
  delta = cohort$OS,        # 生存状态
  marker = cohort$RS, # 预测分数
  cause = 1,                # 指定事件
  times = 3,                # 指定时间点 (1 年)
  iid = TRUE
)

roc_5year <- timeROC(
  T = cohort$OS.time,              # 生存时间
  delta = cohort$OS,        # 生存状态
  marker = cohort$RS, # 预测分数
  cause = 1,                # 指定事件
  times = 5,                # 指定时间点 (1 年)
  iid = TRUE
)

print(roc_1year)
print(roc_3year)
print(roc_5year)

# 查看 AUC 值
auc_1year <- roc_1year$AUC[2]
auc_3year <- roc_3year$AUC[2]
auc_5year <- roc_5year$AUC[2]

# 绘制 ROC 曲线
plot(roc_1year, time = 1, col = pal_npg("nrc")(10)[1], title = "Time-dependent ROC", lwd = 2)
plot(roc_3year, time = 3, add = TRUE, col = pal_npg("nrc")(10)[2], lwd = 2)
plot(roc_5year, time = 5, add = TRUE, col = pal_npg("nrc")(10)[3], lwd = 2)

legend("bottomright", 
       legend = c(paste0("1 Year AUC: ", sprintf("%.3f", auc_1year)),
                  paste0("3 Year AUC: ", sprintf("%.3f", auc_3year)),
                  paste0("5 Year AUC: ", sprintf("%.3f", auc_5year))),
       col = c(pal_npg("nrc")(10)[1], pal_npg("nrc")(10)[2], pal_npg("nrc")(10)[3]), 
       lty = 1, lwd = 2) # 宽 5.5 高 5