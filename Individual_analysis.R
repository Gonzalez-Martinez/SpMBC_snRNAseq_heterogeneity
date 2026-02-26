suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(future)
  library(dplyr)
  library(ggplot2)
  library(openxlsx)
  library(DoubletFinder)
  library(infercnv)
})

set.seed(1234)
options(future.globals.maxSize = 6 * 1024^3)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript MBC1.R <path_to_filtered_feature_bc_matrix.h5>", call. = FALSE)
in_h5 <- args[1]
if (!file.exists(in_h5)) stop("Input file not found: ", in_h5, call. = FALSE)

outdir <- "results_MBC1"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

qc_min_features <- 200
qc_min_log10gumi <- 0.8
qc_max_percent_mt <- 5
min_cells_per_gene <- 10

pcs_df <- 1:15
pcs_umap <- 1:12
res_main <- 0.6

doublet_rate <- 0.075
pN <- 0.25

run_infercnv <- FALSE
gene_order_file <- "gene_order_file.txt"   # placeholder (add to .gitignore or provide separately)
infercnv_ref <- "T cells"
infercnv_out <- file.path(outdir, "infercnv_MBC1")

hdf5_obj <- Read10X_h5(in_h5, use.names = TRUE, unique.features = TRUE)
MBC1 <- CreateSeuratObject(hdf5_obj, min.cells = 0, min.features = 0)
MBC1$orig.ident <- "MBC1"

MBC1[["percent_mt"]] <- PercentageFeatureSet(MBC1, pattern = "^MT-")
MBC1$log10GenesPerUMI <- log10(MBC1$nFeature_RNA) / log10(MBC1$nCount_RNA)

qc_vln <- VlnPlot(MBC1, features = c("nFeature_RNA","nCount_RNA","percent_mt"), pt.size = 0)
ggsave(file.path(outdir, "QC_violin_raw.pdf"), qc_vln, width = 8, height = 4)

MBC1 <- subset(
  MBC1,
  subset = nFeature_RNA > qc_min_features &
    log10GenesPerUMI > qc_min_log10gumi &
    percent_mt < qc_max_percent_mt
)

counts <- GetAssayData(MBC1, slot = "counts")
keep_genes <- Matrix::rowSums(counts > 0) >= min_cells_per_gene
MBC1 <- MBC1[keep_genes, ]

MBC1 <- SCTransform(MBC1, verbose = FALSE)
MBC1 <- RunPCA(MBC1, npcs = 30, verbose = FALSE)
MBC1 <- RunUMAP(MBC1, dims = pcs_umap, verbose = FALSE)
MBC1 <- FindNeighbors(MBC1, dims = pcs_umap, verbose = FALSE)
MBC1 <- FindClusters(MBC1, resolution = res_main, verbose = FALSE)

save(MBC1, file = file.path(outdir, "MBC1_preDF.RData"))

sweep.res.list <- paramSweep(MBC1, PCs = pcs_df, sct = TRUE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)

pK <- bcmvn %>%
  filter(BCmetric == max(BCmetric)) %>%
  pull(pK) %>% as.character() %>% as.numeric()

annotations <- MBC1@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(doublet_rate * nrow(MBC1@meta.data))
nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))

MBC1_D <- doubletFinder(MBC1, PCs = pcs_umap, pN = pN, pK = pK, nExp = nExp_poi, reuse.pANN = NULL, sct = TRUE)
pANN_col <- colnames(MBC1_D@meta.data)[grep("^pANN_", colnames(MBC1_D@meta.data))][1]
MBC1_D <- doubletFinder(MBC1_D, PCs = pcs_umap, pN = pN, pK = pK, nExp = nExp_poi.adj, reuse.pANN = pANN_col, sct = TRUE)

df_cols <- colnames(MBC1_D@meta.data)[grep("^DF.classifications_", colnames(MBC1_D@meta.data))]
df_class_col <- df_cols[length(df_cols)]
colnames(MBC1_D@meta.data)[colnames(MBC1_D@meta.data) == df_class_col] <- "DF_class"

MBC1_singlets <- subset(MBC1_D, subset = DF_class == "Singlet")
save(MBC1_singlets, file = file.path(outdir, "MBC1_singlets.RData"))

MBC1_singlets <- RunPCA(MBC1_singlets, npcs = 30, verbose = FALSE)
MBC1_singlets <- RunUMAP(MBC1_singlets, dims = 1:20, verbose = FALSE)
MBC1_singlets <- FindNeighbors(MBC1_singlets, dims = 1:20, verbose = FALSE)
MBC1_singlets <- FindClusters(MBC1_singlets, resolution = res_main, verbose = FALSE)

p_umap <- DimPlot(MBC1_singlets, reduction = "umap", label = TRUE) + NoLegend()
ggsave(file.path(outdir, "UMAP_MBC1_singlets.pdf"), p_umap, width = 6.5, height = 6)

DefaultAssay(MBC1_singlets) <- "SCT"
markers <- FindAllMarkers(MBC1_singlets, only.pos = TRUE)
save(markers, file = file.path(outdir, "markers_MBC1_singlets.RData"))

top100 <- markers %>% group_by(cluster) %>% slice_max(order_by = avg_log2FC, n = 100)
write.xlsx(top100, file = file.path(outdir, "top100_markers_by_cluster.xlsx"), rowNames = FALSE)

save(MBC1_singlets, file = file.path(outdir, "MBC1_singlets_final.RData"))

if (run_infercnv) {
  if (!file.exists(gene_order_file)) stop("Missing gene_order_file: ", gene_order_file, call. = FALSE)
  
  Object_to_CNVs <- MBC1_singlets
  DefaultAssay(Object_to_CNVs) <- "RNA"
  counts_matrix <- GetAssayData(Object_to_CNVs, layer = "counts")
  sample_annotation <- Object_to_CNVs@meta.data[, "seurat_clusters", drop = FALSE]
  
  infercnv_obj <- CreateInfercnvObject(
    raw_counts_matrix = counts_matrix,
    annotations_file = sample_annotation,
    delim = "\t",
    gene_order_file = gene_order_file,
    ref_group_names = infercnv_ref
  )
  
  infercnv::run(
    infercnv_obj,
    cutoff = 0.1,
    out_dir = infercnv_out,
    denoise = TRUE,
    HMM = TRUE
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(future)
  library(dplyr)
  library(ggplot2)
  library(openxlsx)
  library(DoubletFinder)
  library(infercnv)
})

set.seed(1234)
options(future.globals.maxSize = 6 * 1024^3)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript MBC2.R <path_to_filtered_feature_bc_matrix.h5>", call. = FALSE)
in_h5 <- args[1]
if (!file.exists(in_h5)) stop("Input file not found: ", in_h5, call. = FALSE)

outdir <- "results_MBC2"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

qc_min_features <- 200
qc_min_log10gumi <- 0.8
qc_max_percent_mt <- 5
min_cells_per_gene <- 10

pcs_df <- 1:15
pcs_umap <- 1:12
res_main <- 0.6

doublet_rate <- 0.075
pN <- 0.25

run_infercnv <- FALSE
gene_order_file <- "gene_order_file.txt"
infercnv_ref <- "T cells"
infercnv_out <- file.path(outdir, "infercnv_MBC2")

hdf5_obj <- Read10X_h5(in_h5, use.names = TRUE, unique.features = TRUE)
MBC2 <- CreateSeuratObject(hdf5_obj, min.cells = 0, min.features = 0)
MBC2$orig.ident <- "MBC2"

MBC2[["percent_mt"]] <- PercentageFeatureSet(MBC2, pattern = "^MT-")
MBC2$log10GenesPerUMI <- log10(MBC2$nFeature_RNA) / log10(MBC2$nCount_RNA)

qc_vln <- VlnPlot(MBC2, features = c("nFeature_RNA","nCount_RNA","percent_mt"), pt.size = 0)
ggsave(file.path(outdir, "QC_violin_raw.pdf"), qc_vln, width = 8, height = 4)

MBC2 <- subset(
  MBC2,
  subset = nFeature_RNA > qc_min_features &
    log10GenesPerUMI > qc_min_log10gumi &
    percent_mt < qc_max_percent_mt
)

counts <- GetAssayData(MBC2, slot = "counts")
keep_genes <- Matrix::rowSums(counts > 0) >= min_cells_per_gene
MBC2 <- MBC2[keep_genes, ]

MBC2 <- SCTransform(MBC2, verbose = FALSE)
MBC2 <- RunPCA(MBC2, npcs = 30, verbose = FALSE)
MBC2 <- RunUMAP(MBC2, dims = pcs_umap, verbose = FALSE)
MBC2 <- FindNeighbors(MBC2, dims = pcs_umap, verbose = FALSE)
MBC2 <- FindClusters(MBC2, resolution = res_main, verbose = FALSE)

save(MBC2, file = file.path(outdir, "MBC2_preDF.RData"))

sweep.res.list <- paramSweep(MBC2, PCs = pcs_df, sct = TRUE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)

pK <- bcmvn %>%
  filter(BCmetric == max(BCmetric)) %>%
  pull(pK) %>% as.character() %>% as.numeric()

annotations <- MBC2@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(doublet_rate * nrow(MBC2@meta.data))
nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))

MBC2_D <- doubletFinder(MBC2, PCs = pcs_umap, pN = pN, pK = pK, nExp = nExp_poi, reuse.pANN = NULL, sct = TRUE)
pANN_col <- colnames(MBC2_D@meta.data)[grep("^pANN_", colnames(MBC2_D@meta.data))][1]
MBC2_D <- doubletFinder(MBC2_D, PCs = pcs_umap, pN = pN, pK = pK, nExp = nExp_poi.adj, reuse.pANN = pANN_col, sct = TRUE)

df_cols <- colnames(MBC2_D@meta.data)[grep("^DF.classifications_", colnames(MBC2_D@meta.data))]
df_class_col <- df_cols[length(df_cols)]
colnames(MBC2_D@meta.data)[colnames(MBC2_D@meta.data) == df_class_col] <- "DF_class"

MBC2_singlets <- subset(MBC2_D, subset = DF_class == "Singlet")
save(MBC2_singlets, file = file.path(outdir, "MBC2_singlets.RData"))

MBC2_singlets <- RunPCA(MBC2_singlets, npcs = 30, verbose = FALSE)
MBC2_singlets <- RunUMAP(MBC2_singlets, dims = 1:20, verbose = FALSE)
MBC2_singlets <- FindNeighbors(MBC2_singlets, dims = 1:20, verbose = FALSE)
MBC2_singlets <- FindClusters(MBC2_singlets, resolution = res_main, verbose = FALSE)

p_umap <- DimPlot(MBC2_singlets, reduction = "umap", label = TRUE) + NoLegend()
ggsave(file.path(outdir, "UMAP_MBC2_singlets.pdf"), p_umap, width = 6.5, height = 6)

DefaultAssay(MBC2_singlets) <- "SCT"
markers <- FindAllMarkers(MBC2_singlets, only.pos = TRUE)
save(markers, file = file.path(outdir, "markers_MBC2_singlets.RData"))

top100 <- markers %>% group_by(cluster) %>% slice_max(order_by = avg_log2FC, n = 100)
write.xlsx(top100, file = file.path(outdir, "top100_markers_by_cluster.xlsx"), rowNames = FALSE)

save(MBC2_singlets, file = file.path(outdir, "MBC2_singlets_final.RData"))

if (run_infercnv) {
  if (!file.exists(gene_order_file)) stop("Missing gene_order_file: ", gene_order_file, call. = FALSE)
  
  Object_to_CNVs <- MBC2_singlets
  DefaultAssay(Object_to_CNVs) <- "RNA"
  counts_matrix <- GetAssayData(Object_to_CNVs, layer = "counts")
  sample_annotation <- Object_to_CNVs@meta.data[, "seurat_clusters", drop = FALSE]
  
  infercnv_obj <- CreateInfercnvObject(
    raw_counts_matrix = counts_matrix,
    annotations_file = sample_annotation,
    delim = "\t",
    gene_order_file = gene_order_file,
    ref_group_names = infercnv_ref
  )
  
  infercnv::run(
    infercnv_obj,
    cutoff = 0.1,
    out_dir = infercnv_out,
    denoise = TRUE,
    HMM = TRUE
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(future)
  library(dplyr)
  library(ggplot2)
  library(openxlsx)
  library(DoubletFinder)
  library(infercnv)
})

set.seed(1234)
options(future.globals.maxSize = 6 * 1024^3)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript MBC3.R <path_to_filtered_feature_bc_matrix.h5>", call. = FALSE)
in_h5 <- args[1]
if (!file.exists(in_h5)) stop("Input file not found: ", in_h5, call. = FALSE)

outdir <- "results_MBC3"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

qc_min_features <- 200
qc_min_log10gumi <- 0.8
qc_max_percent_mt <- 5
min_cells_per_gene <- 10

pcs_df <- 1:15
pcs_umap <- 1:12
res_main <- 0.6

doublet_rate <- 0.075
pN <- 0.25

run_infercnv <- FALSE
gene_order_file <- "gene_order_file.txt"
infercnv_ref <- "T cells"
infercnv_out <- file.path(outdir, "infercnv_MBC3")

hdf5_obj <- Read10X_h5(in_h5, use.names = TRUE, unique.features = TRUE)
MBC3 <- CreateSeuratObject(hdf5_obj, min.cells = 0, min.features = 0)
MBC3$orig.ident <- "MBC3"

MBC3[["percent_mt"]] <- PercentageFeatureSet(MBC3, pattern = "^MT-")
MBC3$log10GenesPerUMI <- log10(MBC3$nFeature_RNA) / log10(MBC3$nCount_RNA)

qc_vln <- VlnPlot(MBC3, features = c("nFeature_RNA","nCount_RNA","percent_mt"), pt.size = 0)
ggsave(file.path(outdir, "QC_violin_raw.pdf"), qc_vln, width = 8, height = 4)

MBC3 <- subset(
  MBC3,
  subset = nFeature_RNA > qc_min_features &
    log10GenesPerUMI > qc_min_log10gumi &
    percent_mt < qc_max_percent_mt
)

counts <- GetAssayData(MBC3, slot = "counts")
keep_genes <- Matrix::rowSums(counts > 0) >= min_cells_per_gene
MBC3 <- MBC3[keep_genes, ]

MBC3 <- SCTransform(MBC3, verbose = FALSE)
MBC3 <- RunPCA(MBC3, npcs = 30, verbose = FALSE)
MBC3 <- RunUMAP(MBC3, dims = pcs_umap, verbose = FALSE)
MBC3 <- FindNeighbors(MBC3, dims = pcs_umap, verbose = FALSE)
MBC3 <- FindClusters(MBC3, resolution = res_main, verbose = FALSE)

save(MBC3, file = file.path(outdir, "MBC3_preDF.RData"))

sweep.res.list <- paramSweep(MBC3, PCs = pcs_df, sct = TRUE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)

pK <- bcmvn %>%
  filter(BCmetric == max(BCmetric)) %>%
  pull(pK) %>% as.character() %>% as.numeric()

annotations <- MBC3@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(doublet_rate * nrow(MBC3@meta.data))
nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))

MBC3_D <- doubletFinder(MBC3, PCs = pcs_umap, pN = pN, pK = pK, nExp = nExp_poi, reuse.pANN = NULL, sct = TRUE)
pANN_col <- colnames(MBC3_D@meta.data)[grep("^pANN_", colnames(MBC3_D@meta.data))][1]
MBC3_D <- doubletFinder(MBC3_D, PCs = pcs_umap, pN = pN, pK = pK, nExp = nExp_poi.adj, reuse.pANN = pANN_col, sct = TRUE)

df_cols <- colnames(MBC3_D@meta.data)[grep("^DF.classifications_", colnames(MBC3_D@meta.data))]
df_class_col <- df_cols[length(df_cols)]
colnames(MBC3_D@meta.data)[colnames(MBC3_D@meta.data) == df_class_col] <- "DF_class"

MBC3_singlets <- subset(MBC3_D, subset = DF_class == "Singlet")
save(MBC3_singlets, file = file.path(outdir, "MBC3_singlets.RData"))

MBC3_singlets <- RunPCA(MBC3_singlets, npcs = 30, verbose = FALSE)
MBC3_singlets <- RunUMAP(MBC3_singlets, dims = 1:20, verbose = FALSE)
MBC3_singlets <- FindNeighbors(MBC3_singlets, dims = 1:20, verbose = FALSE)
MBC3_singlets <- FindClusters(MBC3_singlets, resolution = res_main, verbose = FALSE)

p_umap <- DimPlot(MBC3_singlets, reduction = "umap", label = TRUE) + NoLegend()
ggsave(file.path(outdir, "UMAP_MBC3_singlets.pdf"), p_umap, width = 6.5, height = 6)

DefaultAssay(MBC3_singlets) <- "SCT"
markers <- FindAllMarkers(MBC3_singlets, only.pos = TRUE)
save(markers, file = file.path(outdir, "markers_MBC3_singlets.RData"))

top100 <- markers %>% group_by(cluster) %>% slice_max(order_by = avg_log2FC, n = 100)
write.xlsx(top100, file = file.path(outdir, "top100_markers_by_cluster.xlsx"), rowNames = FALSE)

save(MBC3_singlets, file = file.path(outdir, "MBC3_singlets_final.RData"))

if (run_infercnv) {
  if (!file.exists(gene_order_file)) stop("Missing gene_order_file: ", gene_order_file, call. = FALSE)
  
  Object_to_CNVs <- MBC3_singlets
  DefaultAssay(Object_to_CNVs) <- "RNA"
  counts_matrix <- GetAssayData(Object_to_CNVs, layer = "counts")
  sample_annotation <- Object_to_CNVs@meta.data[, "seurat_clusters", drop = FALSE]
  
  infercnv_obj <- CreateInfercnvObject(
    raw_counts_matrix = counts_matrix,
    annotations_file = sample_annotation,
    delim = "\t",
    gene_order_file = gene_order_file,
    ref_group_names = infercnv_ref
  )
  
  infercnv::run(
    infercnv_obj,
    cutoff = 0.1,
    out_dir = infercnv_out,
    denoise = TRUE,
    HMM = TRUE
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(future)
  library(dplyr)
  library(ggplot2)
  library(openxlsx)
  library(DoubletFinder)
  library(infercnv)
})

set.seed(1234)
options(future.globals.maxSize = 6 * 1024^3)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript MBC4.R <path_to_filtered_feature_bc_matrix.h5>", call. = FALSE)
in_h5 <- args[1]
if (!file.exists(in_h5)) stop("Input file not found: ", in_h5, call. = FALSE)

outdir <- "results_MBC4"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

qc_min_features <- 200
qc_min_log10gumi <- 0.8
qc_max_percent_mt <- 5
min_cells_per_gene <- 10

pcs_df <- 1:15
pcs_umap <- 1:12
res_main <- 0.6

doublet_rate <- 0.075
pN <- 0.25

run_infercnv <- FALSE
gene_order_file <- "gene_order_file.txt"
infercnv_ref <- "T cells"
infercnv_out <- file.path(outdir, "infercnv_MBC4")

hdf5_obj <- Read10X_h5(in_h5, use.names = TRUE, unique.features = TRUE)
MBC4 <- CreateSeuratObject(hdf5_obj, min.cells = 0, min.features = 0)
MBC4$orig.ident <- "MBC4"

MBC4[["percent_mt"]] <- PercentageFeatureSet(MBC4, pattern = "^MT-")
MBC4$log10GenesPerUMI <- log10(MBC4$nFeature_RNA) / log10(MBC4$nCount_RNA)

qc_vln <- VlnPlot(MBC4, features = c("nFeature_RNA","nCount_RNA","percent_mt"), pt.size = 0)
ggsave(file.path(outdir, "QC_violin_raw.pdf"), qc_vln, width = 8, height = 4)

MBC4 <- subset(
  MBC4,
  subset = nFeature_RNA > qc_min_features &
    log10GenesPerUMI > qc_min_log10gumi &
    percent_mt < qc_max_percent_mt
)

counts <- GetAssayData(MBC4, slot = "counts")
keep_genes <- Matrix::rowSums(counts > 0) >= min_cells_per_gene
MBC4 <- MBC4[keep_genes, ]

MBC4 <- SCTransform(MBC4, verbose = FALSE)
MBC4 <- RunPCA(MBC4, npcs = 30, verbose = FALSE)
MBC4 <- RunUMAP(MBC4, dims = pcs_umap, verbose = FALSE)
MBC4 <- FindNeighbors(MBC4, dims = pcs_umap, verbose = FALSE)
MBC4 <- FindClusters(MBC4, resolution = res_main, verbose = FALSE)

save(MBC4, file = file.path(outdir, "MBC4_preDF.RData"))

sweep.res.list <- paramSweep(MBC4, PCs = pcs_df, sct = TRUE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)

pK <- bcmvn %>%
  filter(BCmetric == max(BCmetric)) %>%
  pull(pK) %>% as.character() %>% as.numeric()

annotations <- MBC4@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(doublet_rate * nrow(MBC4@meta.data))
nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))

MBC4_D <- doubletFinder(MBC4, PCs = pcs_umap, pN = pN, pK = pK, nExp = nExp_poi, reuse.pANN = NULL, sct = TRUE)
pANN_col <- colnames(MBC4_D@meta.data)[grep("^pANN_", colnames(MBC4_D@meta.data))][1]
MBC4_D <- doubletFinder(MBC4_D, PCs = pcs_umap, pN = pN, pK = pK, nExp = nExp_poi.adj, reuse.pANN = pANN_col, sct = TRUE)

df_cols <- colnames(MBC4_D@meta.data)[grep("^DF.classifications_", colnames(MBC4_D@meta.data))]
df_class_col <- df_cols[length(df_cols)]
colnames(MBC4_D@meta.data)[colnames(MBC4_D@meta.data) == df_class_col] <- "DF_class"

MBC4_singlets <- subset(MBC4_D, subset = DF_class == "Singlet")
save(MBC4_singlets, file = file.path(outdir, "MBC4_singlets.RData"))

MBC4_singlets <- RunPCA(MBC4_singlets, npcs = 30, verbose = FALSE)
MBC4_singlets <- RunUMAP(MBC4_singlets, dims = 1:20, verbose = FALSE)
MBC4_singlets <- FindNeighbors(MBC4_singlets, dims = 1:20, verbose = FALSE)
MBC4_singlets <- FindClusters(MBC4_singlets, resolution = res_main, verbose = FALSE)

p_umap <- DimPlot(MBC4_singlets, reduction = "umap", label = TRUE) + NoLegend()
ggsave(file.path(outdir, "UMAP_MBC4_singlets.pdf"), p_umap, width = 6.5, height = 6)

DefaultAssay(MBC4_singlets) <- "SCT"
markers <- FindAllMarkers(MBC4_singlets, only.pos = TRUE)
save(markers, file = file.path(outdir, "markers_MBC4_singlets.RData"))

top100 <- markers %>% group_by(cluster) %>% slice_max(order_by = avg_log2FC, n = 100)
write.xlsx(top100, file = file.path(outdir, "top100_markers_by_cluster.xlsx"), rowNames = FALSE)

save(MBC4_singlets, file = file.path(outdir, "MBC4_singlets_final.RData"))

if (run_infercnv) {
  if (!file.exists(gene_order_file)) stop("Missing gene_order_file: ", gene_order_file, call. = FALSE)
  
  Object_to_CNVs <- MBC4_singlets
  DefaultAssay(Object_to_CNVs) <- "RNA"
  counts_matrix <- GetAssayData(Object_to_CNVs, layer = "counts")
  sample_annotation <- Object_to_CNVs@meta.data[, "seurat_clusters", drop = FALSE]
  
  infercnv_obj <- CreateInfercnvObject(
    raw_counts_matrix = counts_matrix,
    annotations_file = sample_annotation,
    delim = "\t",
    gene_order_file = gene_order_file,
    ref_group_names = infercnv_ref
  )
  
  infercnv::run(
    infercnv_obj,
    cutoff = 0.1,
    out_dir = infercnv_out,
    denoise = TRUE,
    HMM = TRUE
  )
}
