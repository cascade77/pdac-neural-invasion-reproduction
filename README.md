# PDAC Neural Invasion Reproduction


# 01 Preprocessing

This is the starting point of the whole reproduction pipeline. Before any group member can do clustering, cell type annotation, or spatial analysis (any step of the pipeline), the raw single-cell RNA sequencing data needs to be cleaned up and made comparable across samples. That is what this part does.

---

## Dataset

The data comes from four pancreatic ductal adenocarcinoma (PDAC) tumor samples submitted on GEO under accession [GSE278694](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE278694). The original paper profiled patients with varying degrees of neural invasion (NI), which is basically how much cancer has spread into surrounding nerves. Two samples come from high neural invasion patients and two from low neural invasion patients:

| Sample ID | GEO Accession | Condition |
|-----------|---------------|-----------|
| PA01 | GSM8552940 | High NI |
| PA02 | GSM8552941 | High NI |
| PA21 | GSM8552952 | Low NI |
| PA22 | GSM8552953 | Low NI |

Each file is an HDF5 (.h5) filtered feature-barcode matrix produced by Cell Ranger, which is the standard 10x Genomics output format. These contain raw UMI counts per gene per cell barcode.

---

## Platform

- **OS:** Ubuntu 24 (Linux)
- **Language:** R
- **Environment:** Run via `Rscript` from terminal (not RStudio, not Colab)

---

## Dependencies

```bash
Rscript -e "install.packages('BiocManager', repos='https://cloud.r-project.org'); BiocManager::install('Seurat')"
Rscript -e "install.packages(c('harmony', 'hdf5r', 'ggplot2'), repos='https://cloud.r-project.org')"
```

Seurat is the main framework for single-cell analysis. hdf5r is needed to read the .h5 files. ggplot2 handles the QC plots. harmony is used for batch correction across the four samples.

---

## Workflow

### Step 1: Load the raw data

```r
library(Seurat)
library(harmony)
library(ggplot2)

data_dir <- "~/Documents/Repos/pdac-neural-invasion-reproduction/data/raw"
out_dir  <- "~/Documents/Repos/pdac-neural-invasion-reproduction/01_preprocessing/results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

sample_files <- list(
    high_NI_1 = file.path(data_dir, "GSM8552940_PA01_filtered_feature_bc_matrix.h5"),
    high_NI_2 = file.path(data_dir, "GSM8552941_PA02_filtered_feature_bc_matrix.h5"),
    low_NI_1  = file.path(data_dir, "GSM8552952_PA21_filtered_feature_bc_matrix.h5"),
    low_NI_2  = file.path(data_dir, "GSM8552953_PA22_filtered_feature_bc_matrix.h5")
)
```

Each .h5 file is a sparse matrix of gene counts across thousands of cell barcodes. We  start with labeling each sample with its ID and condition so that information is preserved throughout the analysis. `Read10X_h5` reads the compressed matrix and `CreateSeuratObject` covers it into a Seurat object that all downstream analysis functions work with.

### Step 2: Compute per-cell QC metrics

```r
load_sample <- function(path, sample_id)
{
    counts <- Read10X_h5(path)
    obj    <- CreateSeuratObject(counts, project = sample_id, min.cells = 3, min.features = 200)
    obj$sample_id <- sample_id
    obj$condition <- ifelse(grepl("high", sample_id), "high_NI", "low_NI")
    obj[["pct_mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
    obj
}

seurat_list <- mapply(load_sample, sample_files, names(sample_files), SIMPLIFY = FALSE)
```

The most really important metric here is the mitochondrial gene percentage (`pct_mt`). Mitochondrial genes start with `MT-` in human data. Cells with very high mitochondrial content are usually damaged or dying, because when a cell lyses, cytoplasmic mRNA leaks out but mitochondria stay intact, so you end up sequencing a lot of mitochondrial reads relative to the real transcriptome content. We also look at `nFeature_RNA` (number of detected genes) and `nCount_RNA` (total UMIs). Low gene counts often mean empty droplets or dead cells. Very high gene counts can mean two cells were captured in one droplet (doublets).

### Step 3: Merge samples and visualize QC before filtering

```r
merged <- merge(seurat_list[[1]], y = seurat_list[-1], add.cell.ids = names(seurat_list))

qc_plot <- VlnPlot(merged, features = c("nFeature_RNA", "nCount_RNA", "pct_mt"), ncol = 3)
ggsave(file.path(out_dir, "qc_before_filter.jpeg"), qc_plot, width = 12, height = 5)
```

We merge all 4 samples into a single object before filtering so the violin plots show the distribution across all samples together. This lets you see the overall quality landscape and pick the sensible cutoffs.

### Step 4: Filter low-quality cells

```r
merged <- subset(merged, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & pct_mt < 25)
```

Cells with fewer than 200 detected genes are almost certainly empty droplets. Cells with more than 6000 detected genes are likely doublets. Cells with more than 25% mitochondrial reads are likely dead or damaged. These thresholds are consistent with what the original paper used and are standard in the field.

### Step 5: Normalize and find variable features

```r
merged <- NormalizeData(merged)
merged <- FindVariableFeatures(merged, nfeatures = 2000)
merged <- ScaleData(merged)
merged <- RunPCA(merged, npcs = 30)
```

Raw counts are not directly comparable across cells because each cell has a different total number of reads captured. Normalization divides each cell's counts by its total, then multiplies by 10,000 (log normalization), making values comparable. We then identify the 2000 most variably expressed genes, which are the informative ones that differ meaningfully between cell types. Scaling centers and scales expression so that highly expressed genes do not dominate the PCA. PCA then reduces the data from thousands of gene dimensions down to 30 principal components that capture the main sources of variation.

### Step 6: Harmony batch correction and save output

```r
merged <- RunHarmony(merged, group.by.vars = "sample_id", dims.use = 1:20)

saveRDS(merged, file.path(data_dir, "../processed_subset.rds"))
```

Even after normalization, cells from different patients will cluster by patient rather than by cell type, because there are patient-specific technical and biological differences. Harmony corrects for this by iteratively adjusting the PCA embedding until cells group by biology rather than by which patient they came from. We use the first 20 PCs as input to Harmony. The final processed Seurat object is saved as an .rds file, which is what the downstream analysis (Person 2, clustering) loads directly.

---

## Results

### QC Distributions Before Filtering

![QC violin plot before filtering](results/qc_before_filter.jpeg)

Violin plot showing the distribution of number of detected genes (`nFeature_RNA`), total UMI counts (`nCount_RNA`), and mitochondrial percentage (`pct_mt`) across all four merged samples before any filtering was applied. The distributions look reasonable, with no extreme outliers that would indicate a completely failed sample.

---

## Output

The main output gets passed to the next stage is `data/processed_subset.rds`, a 216MB Seurat object containing all four samples merged, filtered, normalized, and Harmony-integrated. For clustering step, the file is loaded with:

```r
merged <- readRDS("path/to/processed_subset.rds")
```

This file is too large for GitHub (above the 100MB limit) so it is shared separately via [Google Drive](https://drive.google.com/file/d/160K5EYsv9vKkPDBecVerEKo63Ih-o8BS/view?usp=sharing).

---


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



**Integrated Single-Cell and Spatial Transcriptomics Uncover Distinct Cellular Subtypes Involved in Neural Invasion in Pancreatic Cancer.**

*Cancer Cell*, 2025.

---

## Authors

This repository was developed as part of a reproduction study of the PDAC neural invasion paper.

**P2 Role:** Clustering Analysis, Cell-Type Annotation, Schwann Cell Subclustering, Marker Gene Interpretation, and Biological Result Analysis.


# Spatial Transcriptomics Analysis of Pancreatic Cancer Neural Invasion
Data:

https://drive.google.com/drive/folders/1eALXvYj_y0KnL6aL_TokZauuRtQzlywi?usp=sharing

Notebook:

https://drive.google.com/file/d/1Aa0cUpYYrjbE5ruqFfO6JwKwECdjkSfH/view?usp=sharing


 **Pancreatic Cancer:**
Pancreatic Cancer is one of the deadliest cancers beacuse tumor cellls physically invade the nerves around them this is called Neural Invasaion(NI).Some patients have high neural invasion and some have low neural invasion. 
The powerful technologies used is:
Spatial Transcriptomics

**Data Used:**
The data is from database Gene Expression Omnibus.
Accession number: GSE278694
Selected 4 tissue slides for analysis 
 Sample ID | Patient | Neural Invasion |
|-----------|---------|-----------------|
| GSM8552946 | PA#05 | High NI |
| GSM8552950 | PA#11 | High NI |
| GSM8552951 | PA#12 | Low NI |
| GSM8552953 | PA#22 | Low NI |

Two High-NI and two Low-NI slides were selected to enable direct biological comparison.
<img width="700" height="900" alt="image" src="https://github.com/user-attachments/assets/d1e9be6e-9253-4d64-967c-f9048e0f3abd" />

## Analysis Pipeline

### 1. Data Acquisition
Dwnloaded data from GEO daatset and extracted 4 visium slides.

### 2. Data Loading
Loaded all 4 slides into Python via Scanpy. Each slide is represented as a matrix:
- **Rows** tissue spots (~3,000–5,000 per slide)
- **Columns** genes (~33,000 total)

Each slide was labeled with its NI status (High/Low) for downstream comparison.

### 3. Quality Control
The low quality spots are removed and noisy genes using following thresholds 
1. Spots with **<200 detected genes** are removed
2. spots with **<25% mitochondrial reads** are removed which are likely dead are damaged cells
3. Genes detected in **<3 spots** are removed.

### 4. Normalization
Normalized all spots to **10,000 total counts** in order to remove bias of sequencing depth

### 5. Feature Selection
Select top 2000 most variable genes

### 6. Dimensionality Reduction
1.PCA reduced 2000 genes to 20 principal components capturing dominant variation
2.UMAP projected 20 PCs to 2D viusalizations 

UMAP Tumor Low NI:
<img width="700" height="400" alt="umapmarkers_Tumor_LowNI_PA12" src="https://github.com/user-attachments/assets/48e9480b-c5f0-4903-9e27-f73dd574a000" />

UMAP Tumor High NI:
<img width="700" height="400" alt="umapmarkers_Fibroblast_HighNI_PA11" src="https://github.com/user-attachments/assets/a7077315-ef7c-442f-b820-4fb88148dfc3" />

### 7. Clustering
Used **Leiden Clustering** to group spots with similar expression profiles.Each cluster represents a distinct cell type or tissue region.
<img width="700" height="400" alt="umapclusters_LowNI_PA12" src="https://github.com/user-attachments/assets/841fa153-491c-4907-a8a5-165d8ad1080b" />

### 8. Cell Type Annotation
Clusters were annotated using known marker genes:

| Cell Type | Marker Genes |
|-----------|-------------|
| Schwann cells | TGFBI, S100B, MPZ |
| Macrophages | CD68, MRC1 |
| Fibroblasts | COL1A1, ACTA2 |
| Tumor cells | EPCAM, KRT19 |
| T cells | CD3D, CD8A |

<img width="700" height="400" alt="umapmarkers_Schwann_TGFBI_LowNI_PA12" src="https://github.com/user-attachments/assets/5a40f35e-cbd5-4004-9a31-d2dc0c3e3aa2" />

<img width="700" height="400" alt="umapmarkers_Macrophage_HighNI_PA11" src="https://github.com/user-attachments/assets/1fb00778-3465-4ae8-a456-55f4889230a5" />

### 9. High-NI vs Low-NI Comparison
The 4 slides were merged into 1 dataset. UMAP visualization colored by NI tstaue revealed distict seperation between Hogh and low NI samples which reveals meaningful biological variations.

### 10. Differential Expression Analysis
Run **Wilcoxon rank-sum test** between High-NI and Low-NI spots in order to identify genes specifically upregulated or downregulated during neural invasion.

Results are visualized as a **volcano plot**:
- X-axis: log2 fold change (High-NI vs Low-NI)
- Y-axis: −log10 adjusted p-value
- Red dots: significantly upregulated in High-NI
- Blue dots: significantly downregulated in High-NI

<img width="1956" height="1555" alt="image" src="https://github.com/user-attachments/assets/86ecfef6-9877-46c5-ba22-ace9f8c806f5" />

### 11. Cluster Proportion Analysis
Compared the reletive abundance of each celll cluster across High NI and low NI samples.  
Cluster Counts High NI:
<img width="700" height="400" alt="cluster_counts_HighNI_PA05" src="https://github.com/user-attachments/assets/d0da3100-9610-42cb-96fd-4e00b0fdc24f" />
Cluster Counts Low NI:
<img width="700" height="400" alt="cluster_counts_LowNI_PA22" src="https://github.com/user-attachments/assets/e19d22e7-c050-495a-9cac-cc32f37224d3" />

## Output Figures

| Figure | Description |
|--------|-------------|
| Cluster UMAP | Cell populations across all tumor samples |
| Marker gene plots | Spatial distribution of key cell type markers |
| Combined UMAP | High-NI vs Low-NI sample separation |
| Differential expression dotplot | Genes activated in nerve-invaded tissue |
| Volcano plot | Significantly upregulated genes in High-NI tumors |
| Cluster proportions | Cell type balance between High and Low NI |
| Marker gene score maps | Spatial cell type abundance across tissue |


## Key Findings

- High-NI and Low-NI tumors show **distinct cellular compositions** on UMAP
- **Schwann cells and macrophages** are spatially concentrated in specific tumor regions
- Differential expression analysis identifies candidate genes driving neural invasion
- Cell type abundance differs between High-NI and Low-NI samples, consistent with the study's central hypothesis


## Objectives:
Trying to investigate how a certain gene(TGFBI) certain nerve cells ( Schwann cells) behave /imapcat in pancreatic cancer

-To download TCGA pancreatic cancer patient data and derive wether patients with higher or lower TGFBI gene expression survve longer<br>
-To run cell chat and find out which cells are sending TGF-beta signals to schwannn cells<br>
-To run Monocle3 to find out how schwann cells change their state ver time during  pancreatihow oc cancer<br>

## Background

Pancreatic cancer is one of the most deadly cancers  mainly because it is usually
detected very late and does not respond well to most treatments.

One of the main reasons is neural invasion. This is basically
when cancer cells start to move and invade the nerves running through the pancreas.
Patients whose tumors show high neural invasion show worse results. But befre recently, noone knew exactly which cell types
were responsible for enabling or resisting this nerve invasion.

### Schwann cells:
Wrap around and protect nerves. These are
found along nerves throughout the body. When
cancerous cells invade nerves, Schwann cells are right at the site of invasion.
The paper we are trying to reproduce found that Schwann cells are not just passive bystanders but they
actively interact with cancer cells and other immune cells.

### TGFBI stands for Transforming Growth Factor Beta Induced.
It is a gene that gets switched on in response to **TGF-beta** signaling. My survival analysis showed that pancreatic
cancer patients with high TGFBI expression in their tumors survive for
shorter periods of time.

## Methodology

### Data Used
- TCGA-PAAD gene expression and survival data (183 patients) from UCSC Xena
- Annotated Seurat object from Person 2 (8344 cells, 5 cell types)

### Survival Analysis
I got TGFBI expression for all 183 patients and split them into  certain categories (high and
low). Then I drew a Kaplan-Meier curve to compare
survival between the two groups using R langauge.

### CellChat Analysis
I ran CellChat on the single cell data from Person 2 to find out which cells
are sending TGF-beta signals to Schwann cells.

### Trajectory Analysis
I took ran Monocle3 pseudotime analysis so basically  I ran clustering, learned the
trajectory graph via R, and ordered cells by pseudotime. This shows how Schwann cells change their state over time along a trajectory

## Results

### Survival Analysis
The Kaplan-Meier plot shows overall survival of 183 pancreatic cancer
patients split by TGFBI expression. Patients with high TGFBI (red line)
show slightly worse survival compared to low TGFBI (blue line). The
p-value is 0.083 which means the trend exists but is not strongly
significant — this could be because the sample size is small (183 patients).
This still supports the paper's finding that TGFBI is linked to worse outcomes.

### CellChat TGF-beta Signaling
The circle plot shows TGF-beta signaling between 5 cell types. Tumor cells
and Macrophages send the strongest TGF-beta signals (thick lines) to other
cells including Schwann cells. Schwann cells also signal back to themselves
(the loop). This confirms the paper's finding that TGF-beta signaling
from the tumor microenvironment reaches Schwann cells and likely drives
the TGFBI+ Schwann cell subtype.

### Monocle3 Trajectory
The UMAP plot shows Schwann cells colored by pseudotime (dark blue = early,
yellow = late). The trajectory line shows Schwann cells start from one
state (dark blue, top left) and gradually transition to a different state
(yellow, bottom). This suggests Schwann cells go through a developmental
change inside the tumor — possibly transitioning toward the TGFBI+ subtype
as pseudotime increases.


## Reference

Chen MM, Gao Q, Ning H, et al. Integrated single-cell and spatial transcriptomics uncover distinct cellular subtypes involved in neural invasion in pancreatic cancer. *Cancer Cell.* 2025. [https://doi.org/10.1016/j.ccell.2025.04.014](https://doi.org/10.1016/j.ccell.2025.04.01

Chen MM, Gao Q, Ning H, et al.
