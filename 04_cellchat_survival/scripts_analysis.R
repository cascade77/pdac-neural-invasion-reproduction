
# Person 4 Analysis Script
# PDAC Neural Invasion Paper Reproduction

# === SURVIVAL ANALYSIS ===
library(survival)
library(survminer)
library(dplyr)

# Load data
expr <- read.table("TCGA.PAAD.sampleMap_HiSeqV2.gz", header=TRUE, sep="\t", row.names=1)
surv_clean <- read.table("https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/survival%2FPAAD_survival.txt", header=TRUE, sep="\t")

# TGFBI expression
colnames(expr) <- gsub("\\.", "-", colnames(expr))
merged <- merge(surv_clean, data.frame(sample=colnames(expr), tgfbi=as.numeric(expr["TGFBI",])), by="sample")
merged$group <- ifelse(merged$tgfbi > median(merged$tgfbi), "High", "Low")

# Kaplan-Meier plot
fit <- survfit(Surv(OS.time, OS) ~ group, data=merged)
ggsurvplot(fit, data=merged, pval=TRUE, palette=c("#E24B4A","#185FA5"),
           title="Overall Survival by TGFBI Expression (TCGA-PAAD)")

# === CELLCHAT ANALYSIS ===
library(CellChat)
library(Seurat)
# ... (CellChat pipeline run on seurat_annotated.rds)

# === MONOCLE3 TRAJECTORY ===
library(monocle3)
library(SeuratWrappers)
# ... (Trajectory analysis on Schwann cells)
