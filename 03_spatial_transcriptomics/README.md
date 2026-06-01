# Spatial Transcriptomics Analysis of Pancreatic Cancer Neural Invasion

https://drive.google.com/drive/folders/1eALXvYj_y0KnL6aL_TokZauuRtQzlywi?usp=drive_link

 **Pancreatic Cancer:**
Pancreatic Cancer is one of the deadliest cancers beacuse tumor cellls physically invade the nerves around them this is called Neural Invasaion(NI).Some patients have high neural invasion and some have low neural invasion. 
Tow powerful technologies used are:
1. Single cell RNA sequencing
2. Spatial Transcriptomics

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

##Analysis Pipeline

###1. Data Acquisistion:
Dwnloaded data from GEO daatset and extracted 4 visium slides.

###2.Data Loading:
Loaded all 4 slides into Python via Scanpy. Each slide is represented as a matrix:
- **Rows** — tissue spots (~3,000–5,000 per slide)
- **Columns** — genes (~33,000 total)

Each slide was labeled with its NI status (High/Low) for downstream comparison.

###3.Quality Control:
The low quality spots are removed and noisy genes using following thresholds 
1. Spots with **<200 detected genes** are removed
2. spots with **<25% mitochondrial reads** are removed which are likely dead are damaged cells
3. Genes detected in **<3 spots** are removed.

###4.Normalization:
Normalized all spots to **10,000 total counts** in order to remove bias of sequencing depth

##5.Feature Selection:
Select top 2000 most variable genes

###6.Dimensionality Reduction:
1.PCA reduced 2000 genes to 20 principal components capturing dominant variation
2.UMAP projected 20 PCs to 2D viusalizations 

UMAP Tumor Low NI:
<img width="1225" height="562" alt="umapmarkers_Tumor_LowNI_PA12" src="https://github.com/user-attachments/assets/48e9480b-c5f0-4903-9e27-f73dd574a000" />

UMAP Tumor High NI:
<img width="1249" height="562" alt="umapmarkers_Fibroblast_HighNI_PA11" src="https://github.com/user-attachments/assets/a7077315-ef7c-442f-b820-4fb88148dfc3" />

###7.Clustering:
Used **Leiden Clustering** to group spots with similar expression profiles.Each cluster represents a distinct cell type or tissue region.
<img width="627" height="556" alt="umapclusters_LowNI_PA12" src="https://github.com/user-attachments/assets/841fa153-491c-4907-a8a5-165d8ad1080b" />

###8.Cell Type Annotation:
Clusters were annotated using known marker genes:

| Cell Type | Marker Genes |
|-----------|-------------|
| Schwann cells | TGFBI, S100B, MPZ |
| Macrophages | CD68, MRC1 |
| Fibroblasts | COL1A1, ACTA2 |
| Tumor cells | EPCAM, KRT19 |
| T cells | CD3D, CD8A |

<img width="1805" height="562" alt="umapmarkers_Schwann_TGFBI_LowNI_PA12" src="https://github.com/user-attachments/assets/5a40f35e-cbd5-4004-9a31-d2dc0c3e3aa2" />

<img width="1249" height="562" alt="umapmarkers_Macrophage_HighNI_PA11" src="https://github.com/user-attachments/assets/1fb00778-3465-4ae8-a456-55f4889230a5" />

###9.High NI vs Low NI comaprison:
The 4 slides were merged into 1 dataset. UMAP visualization colored by NI tstaue revealed distict seperation between Hogh and low NI samples which reveals meaningful biological variations.

###10.Differential Expression Anlaysis:
Run **Wilcoxon rank-sum test** between High-NI and Low-NI spots in order to identify genes specifically upregulated or downregulated during neural invasion.

Results are visualized as a **volcano plot**:
- X-axis: log2 fold change (High-NI vs Low-NI)
- Y-axis: −log10 adjusted p-value
- Red dots: significantly upregulated in High-NI
- Blue dots: significantly downregulated in High-NI

<img width="1956" height="1555" alt="image" src="https://github.com/user-attachments/assets/86ecfef6-9877-46c5-ba22-ace9f8c806f5" />

###11.Cluster Proportion Analysis:
Compared the reletive abundance of each celll cluster across High NI and low NI samples.  
Cluster Counts High NI:
<img width="1200" height="750" alt="cluster_counts_HighNI_PA05" src="https://github.com/user-attachments/assets/d0da3100-9610-42cb-96fd-4e00b0fdc24f" />
Cluster Counts Low NI:
<img width="1200" height="750" alt="cluster_counts_LowNI_PA22" src="https://github.com/user-attachments/assets/e19d22e7-c050-495a-9cac-cc32f37224d3" />


