# PDAC Neural Invasion Reproduction

## Project Overview

This project reproduces key analyses from the paper:

**"Integrated Single-Cell and Spatial Transcriptomics Uncover Distinct Cellular Subtypes Involved in Neural Invasion in Pancreatic Cancer"**

The study investigates how different cell populations within the pancreatic ductal adenocarcinoma (PDAC) tumor microenvironment contribute to neural invasion, a process where cancer cells infiltrate surrounding nerves and promote disease progression.

Using single-cell RNA sequencing (scRNA-seq) analysis, we reproduced major clustering and cell annotation steps described in the publication and identified cell populations associated with neural invasion.

---

## Objectives

- Reproduce the single-cell transcriptomic workflow from the original paper.
- Identify major cell populations present in PDAC samples.
- Perform dimensionality reduction and clustering.
- Annotate clusters using known marker genes.
- Investigate Schwann-cell heterogeneity.
- Explore expression of **TGFBI**, a key marker associated with neural invasion.

---

## Dataset

The analysis was performed using pancreatic cancer single-cell transcriptomic data described in the original publication.

The original study combined:

- Single-cell RNA sequencing (scRNA-seq)
- Single-nucleus RNA sequencing (snRNA-seq)
- Spatial transcriptomics

to investigate cellular populations involved in neural invasion.

---

# Analysis Workflow

## 1. Data Preprocessing

- Quality control filtering
- Normalization
- Identification of highly variable genes
- Scaling of expression values

## 2. Dimensionality Reduction

Principal Component Analysis (PCA) was used to reduce data dimensionality while preserving biological variation.

## 3. Clustering

Graph-based clustering was performed to identify transcriptionally distinct cell populations.

## 4. Cell Type Annotation

Clusters were annotated using canonical marker genes reported in the literature and the original paper.

## 5. Schwann Cell Subclustering

Schwann cells were isolated and re-clustered to investigate their internal heterogeneity and identify subpopulations linked to neural invasion.

## 6. Marker Validation

Feature plots and dot plots were used to validate cluster annotations and marker gene expression patterns.

---

# My Contribution (P2)

### Clustering and Cell-Type Annotation

My primary responsibilities included:

- Performing clustering analysis on processed single-cell data.
- Generating UMAP visualizations.
- Identifying cluster-specific marker genes.
- Annotating clusters using known biological markers.
- Performing Schwann-cell subclustering.
- Investigating TGFBI expression patterns.
- Interpreting cluster identities and biological relevance.

---

# Results

## 1. UMAP of Major Cell Types

The UMAP visualization revealed distinct cellular populations within the PDAC microenvironment.

Identified cell populations included:

- T cells
- B cells
- Macrophages
- Fibroblasts
- Schwann cells
- Tumor cells

The clear separation between clusters indicates successful clustering and biologically meaningful cell identities.

### Biological Significance

The tumor microenvironment consists of multiple interacting immune, stromal, neural, and malignant cell populations that collectively influence tumor progression and neural invasion.

---

## 2. Schwann Cell Subclustering

A second clustering analysis was performed specifically on Schwann cells.

Multiple Schwann-cell subtypes were identified, demonstrating that Schwann cells are not a homogeneous population.

### Biological Significance

The original paper showed that different Schwann-cell states play different roles in neural invasion, with specific subtypes contributing to tumor progression and nerve infiltration.

---

## 3. TGFBI Feature Plot

Expression of **TGFBI** was visualized across Schwann-cell populations.

TGFBI was not uniformly expressed, indicating enrichment within specific cellular subpopulations.

### Biological Significance

TGFBI is one of the most important markers reported in the original study.

The paper identified a **TGFBI-positive Schwann-cell subtype** associated with:

- Cancer-cell migration
- Neural invasion
- Poor patient prognosis

The observed expression pattern supports the existence of biologically distinct TGFBI-expressing Schwann-cell populations.

---

## 4. Marker Gene Dot Plot

A dot plot was generated to validate cluster annotations using canonical marker genes.

| Cell Type | Marker Gene |
|------------|------------|
| T Cells | CD3D |
| B Cells | CD79A |
| Macrophages | C1QC |
| Fibroblasts | COL1A1 |
| Schwann Cells | SOX10 |
| Tumor Cells | EPCAM |

Dot size represents the percentage of cells expressing a gene, while color intensity represents average expression level.

### Biological Significance

The marker distribution confirms that clusters correspond to biologically meaningful cell types and validates the annotation strategy.

---

# Key Findings

- Multiple cell populations exist within the PDAC tumor microenvironment.
- Schwann cells exhibit significant heterogeneity.
- Distinct Schwann-cell subtypes can be identified through subclustering.
- TGFBI expression is enriched in specific Schwann-cell populations.
- Findings are consistent with the original study's conclusion that TGFBI-positive Schwann cells contribute to neural invasion.

---

# Technologies Used

- R
- Seurat
- ggplot2
- dplyr
- Single-cell RNA Sequencing (scRNA-seq)

---

# Repository Structure

```text
project/
│
├── data/
├── scripts/
├── results/
│   ├── UMAPs/
│   ├── Marker_Analysis/
│   ├── Schwann_Subclustering/
│   └── Feature_Plots/
│
├── figures/
├── notebook/
└── README.md
```

---

# Reference

Chen MM, Gao Q, Ning H, et al.

**Integrated Single-Cell and Spatial Transcriptomics Uncover Distinct Cellular Subtypes Involved in Neural Invasion in Pancreatic Cancer.**

*Cancer Cell*, 2025.

---

## Authors

This repository was developed as part of a reproduction study of the PDAC neural invasion paper.

**P2 Role:** Clustering Analysis, Cell-Type Annotation, Schwann Cell Subclustering, Marker Gene Interpretation, and Biological Result Analysis.
