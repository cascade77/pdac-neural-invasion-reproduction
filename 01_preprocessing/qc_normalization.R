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

load_sample <- function(path, sample_id)
{
    counts <- Read10X_h5(path)
    obj    <- CreateSeuratObject(counts, project = sample_id, min.cells = 3, min.features = 200)
    obj$sample_id <- sample_id
    obj$condition <- ifelse(grepl("high", sample_id), "high_NI", "low_NI")
    obj[["pct_mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
    obj
}

cat("loading samples...\n")
seurat_list <- mapply(load_sample, sample_files, names(sample_files), SIMPLIFY = FALSE)

merged <- merge(seurat_list[[1]], y = seurat_list[-1], add.cell.ids = names(seurat_list))

qc_plot <- VlnPlot(merged, features = c("nFeature_RNA", "nCount_RNA", "pct_mt"), ncol = 3)
ggsave(file.path(out_dir, "qc_before_filter.jpeg"), qc_plot, width = 12, height = 5)
cat("QC plot saved\n")

merged <- subset(merged, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & pct_mt < 25)

cat("normalizing...\n")
merged <- NormalizeData(merged)
merged <- FindVariableFeatures(merged, nfeatures = 2000)
merged <- ScaleData(merged)
merged <- RunPCA(merged, npcs = 30)

cat("running harmony integration...\n")
merged <- RunHarmony(merged, group.by.vars = "sample_id", dims.use = 1:20)

saveRDS(merged, file.path(data_dir, "../processed_subset.rds"))
cat("done! processed object saved.\n")
