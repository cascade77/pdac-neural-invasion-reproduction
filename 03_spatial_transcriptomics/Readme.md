# Spatial Transcriptomics Analysis of Pancreatic Cancer Neural Invasion

Data: https://drive.google.com/drive/folders/1eALXvYj_y0KnL6aL_TokZauuRtQzlywi?usp=sharing

Notebook: https://drive.google.com/file/d/1Aa0cUpYYrjbE5ruqFfO6JwKwECdjkSfH/view?usp=sharing

---

## Background

Pancreatic Cancer is one of the deadliest cancers because tumor cells physically invade the nerves around them, this is called Neural Invasion (NI). Some patients have high neural invasion and some have low neural invasion. The powerful technology used here is Spatial Transcriptomics.

---

## Dataset

The data is from the Gene Expression Omnibus database under accession number GSE278694. Four tissue slides were selected for analysis:

| Sample ID | Patient | Neural Invasion |
|-----------|---------|-----------------|
| GSM8552946 | PA05 | High NI |
| GSM8552950 | PA11 | High NI |
| GSM8552951 | PA12 | Low NI |
| GSM8552953 | PA22 | Low NI |

Two High-NI and two Low-NI slides were selected to enable direct biological comparison.

---

## Platform

- **OS:** Ubuntu (Linux)
- **Language:** Python
- **Environment:** Jupyter Notebook

---

## Dependencies

```bash
pip install scanpy anndata pandas numpy matplotlib seaborn
```

---

## Workflow

### 1. Data Acquisition

Downloaded data from GEO dataset and extracted 4 Visium slides.

### 2. Data Loading

```python
import scanpy as sc

adata = sc.read_visium(path)
adata.obs["sample_id"] = sample_id
adata.obs["condition"] = "High_NI" if "High" in sample_id else "Low_NI"
```

Loaded all 4 slides into Python via Scanpy. Each slide is represented as a matrix where rows are tissue spots (roughly 3,000 to 5,000 per slide) and columns are genes (roughly 33,000 total). Each slide was labeled with its NI status (High/Low) for downstream comparison.

### 3. Quality Control

```python
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata = adata[adata.obs["pct_counts_mt"] < 25]
```

Low quality spots and noisy genes are removed using the following thresholds: spots with fewer than 200 detected genes are removed, spots with more than 25% mitochondrial reads are removed as they are likely dead or damaged cells, and genes detected in fewer than 3 spots are removed.

### 4. Normalization

```python
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
```

Normalized all spots to 10,000 total counts in order to remove bias of sequencing depth.

### 5. Feature Selection

```python
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
```

Selected the top 2000 most variable genes.

### 6. Dimensionality Reduction

```python
sc.tl.pca(adata, n_comps=20)
sc.tl.umap(adata)
```

PCA reduced 2000 genes to 20 principal components capturing dominant variation. UMAP projected 20 PCs to a 2D visualization.

UMAP Tumor Low NI:

![UMAP Tumor Low NI PA12](results/umapmarkers_Tumor_LowNI_PA12.png)

UMAP Tumor High NI:

![UMAP Tumor High NI PA11](results/umapmarkers_Tumor_HighNI_PA11.png)

### 7. Clustering

```python
sc.pp.neighbors(adata, n_pcs=20)
sc.tl.leiden(adata, resolution=0.5)
```

Used Leiden Clustering to group spots with similar expression profiles. Each cluster represents a distinct cell type or tissue region.

![Clusters Low NI PA12](results/umapclusters_LowNI_PA12.png)

### 8. Cell Type Annotation

```python
marker_genes = {
    "Schwann cells": ["TGFBI", "S100B", "MPZ"],
    "Macrophages":   ["CD68", "MRC1"],
    "Fibroblasts":   ["COL1A1", "ACTA2"],
    "Tumor cells":   ["EPCAM", "KRT19"],
    "T cells":       ["CD3D", "CD8A"]
}
```

Clusters were annotated using known marker genes:

| Cell Type | Marker Genes |
|-----------|-------------|
| Schwann cells | TGFBI, S100B, MPZ |
| Macrophages | CD68, MRC1 |
| Fibroblasts | COL1A1, ACTA2 |
| Tumor cells | EPCAM, KRT19 |
| T cells | CD3D, CD8A |

![Schwann TGFBI Low NI PA12](results/umapmarkers_Schwann_TGFBI_LowNI_PA12.png)

![Macrophage High NI PA11](results/umapmarkers_Macrophage_HighNI_PA11.png)

### 9. High-NI vs Low-NI Comparison

The 4 slides were merged into 1 dataset. UMAP visualization colored by NI status revealed distinct separation between High and Low NI samples which reveals meaningful biological variations.

![Combined UMAP](results/umapcombined_umap.png)

### 10. Differential Expression Analysis

```python
sc.tl.rank_genes_groups(adata_combined, groupby="condition", method="wilcoxon")
```

Run Wilcoxon rank-sum test between High-NI and Low-NI spots in order to identify genes specifically upregulated or downregulated during neural invasion. Results are visualized as a volcano plot where the x-axis is log2 fold change (High-NI vs Low-NI), the y-axis is negative log10 adjusted p-value, red dots are significantly upregulated in High-NI, and blue dots are significantly downregulated in High-NI.

![Volcano plot](results/volcano_plot.png)

![Dotplot DE High NI vs Low NI](results/dotplot_DE_highNI_vs_lowNI.png)

### 11. Cluster Proportion Analysis

```python
cluster_counts = adata.obs.groupby(["sample_id", "leiden"]).size().reset_index()
```

Compared the relative abundance of each cell cluster across High NI and Low NI samples.

![Cluster counts High NI PA05](results/cluster_counts_HighNI_PA05.png)

![Cluster counts Low NI PA22](results/cluster_counts_LowNI_PA22.png)

![Cluster proportions comparison](results/cluster_proportions_comparison.png)

---

## All Result Figures

### Cluster UMAPs

![Clusters High NI PA05](results/umapclusters_HighNI_PA05.png)

![Clusters High NI PA11](results/umapclusters_HighNI_PA11.png)

![Clusters Low NI PA22](results/umapclusters_LowNI_PA22.png)

### Tumor Markers

![Tumor High NI PA05](results/umapmarkers_Tumor_HighNI_PA05.png)

![Tumor Low NI PA22](results/umapmarkers_Tumor_LowNI_PA22.png)

### Schwann TGFBI Markers

![Schwann TGFBI High NI PA05](results/umapmarkers_Schwann_TGFBI_HighNI_PA05.png)

![Schwann TGFBI High NI PA11](results/umapmarkers_Schwann_TGFBI_HighNI_PA11.png)

![Schwann TGFBI Low NI PA22](results/umapmarkers_Schwann_TGFBI_LowNI_PA22.png)

### Macrophage Markers

![Macrophage High NI PA05](results/umapmarkers_Macrophage_HighNI_PA05.png)

![Macrophage Low NI PA12](results/umapmarkers_Macrophage_LowNI_PA12.png)

![Macrophage Low NI PA22](results/umapmarkers_Macrophage_LowNI_PA22.png)

### Fibroblast Markers

![Fibroblast High NI PA05](results/umapmarkers_Fibroblast_HighNI_PA05.png)

![Fibroblast High NI PA11](results/umapmarkers_Fibroblast_HighNI_PA11.png)

![Fibroblast Low NI PA12](results/umapmarkers_Fibroblast_LowNI_PA12.png)

![Fibroblast Low NI PA22](results/umapmarkers_Fibroblast_LowNI_PA22.png)

### Cluster Counts

![Cluster counts High NI PA11](results/cluster_counts_HighNI_PA11.png)

![Cluster counts Low NI PA12](results/cluster_counts_LowNI_PA12.png)

---

## Key Findings

High-NI and Low-NI tumors show distinct cellular compositions on UMAP. Schwann cells and macrophages are spatially concentrated in specific tumor regions. Differential expression analysis identifies candidate genes driving neural invasion. Cell type abundance differs between High-NI and Low-NI samples, consistent with the study's central hypothesis.

---

## Reference

Chen MM, Gao Q, Ning H, et al. Integrated single-cell and spatial transcriptomics uncover distinct cellular subtypes involved in neural invasion in pancreatic cancer. *Cancer Cell.* 2025. [https://doi.org/10.1016/j.ccell.2025.04.014](https://doi.org/10.1016/j.ccell.2025.04.014)
