# PDAC Neural Invasion Reproduction



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
<img width="1600" height="1200" alt="umap_celltypes-1" src="https://github.com/user-attachments/assets/cc22c2ef-c431-475e-a05d-fa7b1575c797" />

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
<img width="1200" height="1000" alt="umap_schwann_subtypes-1" src="https://github.com/user-attachments/assets/6b4767ca-a7f5-4448-9843-43b788dfbb2c" />

A second clustering analysis was performed specifically on Schwann cells.

Multiple Schwann-cell subtypes were identified, demonstrating that Schwann cells are not a homogeneous population.

### Biological Significance

The original paper showed that different Schwann-cell states play different roles in neural invasion, with specific subtypes contributing to tumor progression and nerve infiltration.

---

## 3. TGFBI Feature Plot
<img width="1000" height="1000" alt="featureplot_TGFBI-1" src="https://github.com/user-attachments/assets/c2f97245-1e77-4fef-b854-1a247fc1a44a" />

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
<img width="1600" height="1200" alt="umap_celltypes-1" src="https://github.com/user-attachments/assets/dbec4bf9-6252-4881-b64f-7769bac1fd62" />

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
