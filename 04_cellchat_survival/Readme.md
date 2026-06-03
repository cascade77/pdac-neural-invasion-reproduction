# 04 CellChat, Trajectory, and Survival Analysis

---

## Objectives

Trying to investigate how a certain gene (TGFBI) and certain nerve cells (Schwann cells) behave and impact pancreatic cancer:

- To download TCGA pancreatic cancer patient data and determine whether patients with higher or lower TGFBI gene expression survive longer
- To run CellChat and find out which cells are sending TGF-beta signals to Schwann cells
- To run Monocle3 to find out how Schwann cells change their state over time during pancreatic cancer

---

## Background

Pancreatic cancer is one of the most deadly cancers mainly because it is usually detected very late and does not respond well to most treatments. One of the main reasons is neural invasion, basically when cancer cells start to move and invade the nerves running through the pancreas. Patients whose tumors show high neural invasion show worse results. But before recently, no one knew exactly which cell types were responsible for enabling or resisting this nerve invasion.

**Schwann cells** wrap around and protect nerves and are found along nerves throughout the body. When cancerous cells invade nerves, Schwann cells are right at the site of invasion. The paper we are reproducing found that Schwann cells are not just passive bystanders but they actively interact with cancer cells and other immune cells.

**TGFBI** stands for Transforming Growth Factor Beta Induced. It is a gene that gets switched on in response to TGF-beta signaling. The survival analysis showed that pancreatic cancer patients with high TGFBI expression in their tumors survive for shorter periods of time.

---

## Dataset

- TCGA-PAAD gene expression and survival data (183 patients) from UCSC Xena and GDC portal
- Annotated Seurat object from `02_clustering_annotation` (8344 cells, 5 cell types)


---

## Platform

- **OS:** Ubuntu (Linux)
- **Language:** R
- **Environment:** Run via `Rscript` from terminal

---

## Dependencies

```bash
Rscript -e "install.packages('BiocManager', repos='https://cloud.r-project.org')"
Rscript -e "BiocManager::install(c('CellChat', 'monocle3'))"
Rscript -e "install.packages(c('survival', 'survminer', 'ggplot2', 'dplyr'), repos='https://cloud.r-project.org')"
```

---

## Workflow

### Survival Analysis

```r
library(survival)
library(survminer)

tgfbi_expr <- as.numeric(expr["TGFBI", ])
clinical$tgfbi_group <- ifelse(tgfbi_expr > median(tgfbi_expr), "High", "Low")
fit <- survfit(Surv(overall_survival_days, vital_status) ~ tgfbi_group, data = clinical)
ggsurvplot(fit, data = clinical, pval = TRUE)
```

TGFBI expression for all 183 patients was retrieved and split into high and low categories. A Kaplan-Meier curve was drawn to compare survival between the two groups using R.

### CellChat Analysis

```r
library(CellChat)

cellchat <- createCellChat(object = seurat_obj, group.by = "cell_type")
cellchat@DB <- CellChatDB.human
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- computeCommunProb(cellchat)
```

CellChat was run on the single-cell data from `02_clustering_annotation` to find out which cells are sending TGF-beta signals to Schwann cells.

### Trajectory Analysis

```r
library(monocle3)

schwann_cds <- as.cell_data_set(subset(seurat_obj, idents = "Schwann cells"))
schwann_cds <- cluster_cells(schwann_cds)
schwann_cds <- learn_graph(schwann_cds)
schwann_cds <- order_cells(schwann_cds)
```

Monocle3 pseudotime analysis was run by clustering cells, learning the trajectory graph via R, and ordering cells by pseudotime. This shows how Schwann cells change their state over time along a trajectory.

---

## Results

### Survival Analysis

![TGFBI survival plot](results/tgfbi_survival_plot.jpeg)

The Kaplan-Meier plot shows overall survival of 183 pancreatic cancer patients split by TGFBI expression. Patients with high TGFBI (red line) show slightly worse survival compared to low TGFBI (blue line). The p-value is 0.083 which means the trend exists but is not strongly significant, this could be because the sample size is small (183 patients). This still supports the paper's finding that TGFBI is linked to worse outcomes.

### CellChat TGF-beta Signaling

![CellChat TGF-beta plot](results/cellchat_tgfb_plot.jpeg)

The circle plot shows TGF-beta signaling between 5 cell types. Tumor cells and macrophages send the strongest TGF-beta signals (thick lines) to other cells including Schwann cells. Schwann cells also signal back to themselves (the loop). This confirms the paper's finding that TGF-beta signaling from the tumor microenvironment reaches Schwann cells and likely drives the TGFBI+ Schwann cell subtype.

### Monocle3 Trajectory

![Monocle3 trajectory](results/monocle3_trajectory.jpeg)

The UMAP plot shows Schwann cells colored by pseudotime (dark blue = early, yellow = late). The trajectory line shows Schwann cells start from one state (dark blue, top left) and gradually transition to a different state (yellow, bottom). This suggests Schwann cells go through a developmental change inside the tumor, possibly transitioning toward the TGFBI+ subtype as pseudotime increases.

---

## Reference

Chen MM, Gao Q, Ning H, et al. Integrated single-cell and spatial transcriptomics uncover distinct cellular subtypes involved in neural invasion in pancreatic cancer. *Cancer Cell.* 2025. [https://doi.org/10.1016/j.ccell.2025.04.014](https://doi.org/10.1016/j.ccell.2025.04.014)
