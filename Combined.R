suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(stringr)
  library(future)
  library(ggplot2)
  library(openxlsx)
  library(here)
})

set.seed(1234)
options(future.globals.maxSize = 8 * 1024^3)

obj_dir <- here("data", "objects")

dir.create(here("results", "objects"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("results", "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("results", "qc"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("figures", "figure1"), recursive = TRUE, showWarnings = FALSE)

sample_files <- c(
  "NST1_singlets.RData",
  "NST2_singlets.RData",
  "NST3_singlets.RData",
  "NST4_singlets.RData",
  "NST5_singlets.RData",
  "NST6_singlets.RData",
  "MBC1_singlets.RData",
  "MBC2_singlets.RData",
  "MBC3_singlets.RData",
  "MBC4_singlets.RData"
)

sample_ids <- c(
  "NST1","NST2","NST3","NST4","NST5","NST6",
  "MBC1","MBC2","MBC3","MBC4"
)

histology <- ifelse(grepl("^NST", sample_ids), "NST", "MBC")

prep_for_merge <- function(o, sid, htype) {

  md <- o@meta.data
  for (cn in colnames(md)) {
    x <- md[[cn]]
    if (is.matrix(x) && ncol(x) == 1) md[[cn]] <- as.vector(x)
  }
  o@meta.data <- md

  if (!"RNA" %in% Assays(o)) {
    old <- DefaultAssay(o)
    o <- RenameAssay(o, old, "RNA")
  }
  DefaultAssay(o) <- "RNA"

  if ("SCT" %in% Assays(o)) o[["SCT"]] <- NULL

  o <- DietSeurat(
    o,
    assays = "RNA",
    layers = c("counts", "data"),
    dimreducs = NULL,
    graphs = NULL,
    misc = FALSE
  )

  o$Sample_ID <- sid
  o$Histological_type <- htype
  o$orig.ident <- sid

  return(o)
}

objs <- list()

for (i in seq_along(sample_files)) {

  fp <- file.path(obj_dir, sample_files[i])
  if (!file.exists(fp)) stop("Missing file: ", fp)

  load(fp)

  obj_name <- ls()[sapply(ls(), function(x) inherits(get(x), "Seurat"))][1]
  o <- get(obj_name)

  o <- prep_for_merge(o, sample_ids[i], histology[i])
  objs[[sample_ids[i]]] <- o

  rm(list = obj_name)
}

obj_all <- merge(
  x = objs[[1]],
  y = objs[-1],
  add.cell.ids = names(objs),
  project = "MBC_NST_combined"
)

lev_nst <- str_sort(unique(obj_all$Sample_ID[grepl("^NST", obj_all$Sample_ID)]), numeric = TRUE)
lev_mbc <- str_sort(unique(obj_all$Sample_ID[grepl("^MBC", obj_all$Sample_ID)]), numeric = TRUE)

obj_all$Sample_ID <- factor(obj_all$Sample_ID, levels = c(lev_nst, lev_mbc))
obj_all$Histological_type <- factor(obj_all$Histological_type, levels = c("NST","MBC"))

saveRDS(obj_all, here("results","objects","obj_all_merged.rds"))

split_list <- SplitObject(obj_all, split.by = "Sample_ID")

split_list <- lapply(split_list, function(x) {
  DefaultAssay(x) <- "RNA"
  SCTransform(x, vst.flavor = "v2", verbose = FALSE)
})

obj_sct <- merge(split_list[[1]], y = split_list[-1])

features <- SelectIntegrationFeatures(split_list, nfeatures = 3000)
DefaultAssay(obj_sct) <- "SCT"
VariableFeatures(obj_sct) <- features

obj_sct <- RunPCA(obj_sct, assay = "SCT", verbose = FALSE)
obj_sct <- RunUMAP(obj_sct, dims = 1:16)
obj_sct <- FindNeighbors(obj_sct, dims = 1:16)
obj_sct <- FindClusters(obj_sct, resolution = 0.3)

Idents(obj_sct) <- "SCT_snn_res.0.3"
obj_sct$seurat_clusters <- obj_sct$SCT_snn_res.0.3

saveRDS(obj_sct, here("results","objects","obj_sct_clusters.rds"))

message("Combined object ready.")

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(openxlsx)
  library(here)
})

dir.create(here("results", "objects"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("results", "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("figures", "figure1"), recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(here("results", "objects", "obj_sct_clusters.rds"))

Idents(obj) <- "SCT_snn_res.0.3"

cluster_to_celltype_raw <- c(
  "0"  = "CAFs",
  "1"  = "TAMs",
  "2"  = "T-Epi 1",
  "3"  = "T-Epi 2",
  "4"  = "T cells",
  "5"  = "Mes",
  "6"  = "Endo",
  "7"  = "Basal",
  "8"  = "TAMs (SPP1+)",
  "9"  = "T-Myo",
  "10" = "T-LCycling",
  "11" = "Pericytes",
  "12" = "Mast cells",
  "13" = "B cells",
  "14" = "Neutrophils",
  "15" = "N-L"
)

names(cluster_to_celltype_raw) <- levels(obj)
obj <- RenameIdents(obj, cluster_to_celltype_raw)
obj$Cell_type <- as.character(Idents(obj))

celltype_final_map <- c(
  "T-Epi 1" = "T-L1",
  "T-Epi 2" = "T-L2",
  "Mes" = "T-Mes",
  "TAMs" = "TAMs 1",
  "TAMs (SPP1+)" = "TAMs 2"
)

obj$Cell_type <- dplyr::recode(obj$Cell_type, !!!celltype_final_map)

ordered_levels <- c(
  "N-L",
  "T-L2",
  "T-L1",
  "T-LCycling",
  "T-Ker",
  "T-Myo",
  "T-Mes",
  "Basal",
  "CAFs",
  "Endo",
  "Pericytes",
  "TAMs 1",
  "TAMs 2",
  "T cells",
  "B cells",
  "Mast cells",
  "Neutrophils"
)

obj$Cell_type <- factor(obj$Cell_type, levels = ordered_levels)
Idents(obj) <- obj$Cell_type

saveRDS(obj, here("results", "objects", "obj_sct_annotated_final.rds"))

DimPlot2(obj, label=TRUE, repel = T, label.size = 7, theme = theme_umap_arrows())
ClusterDistrBar(origin = obj$Sample_ID, cluster = obj$Cell_type, flip = FALSE, reverse_order = FALSE)
                          
grouped_features_fig1 <- list(
  "Lum"      = c("ESR1","PGR","GATA3","FOXA1","ANKRD30A","CDH1","EPCAM"),
  "Cycl"     = c("MKI67","TOP2A"),
  "Basal"    = c("KRT5","KRT14","KRT17","CDH3","TP63","ITGA6"),
  "Kert"     = c("IVL","SPRR1B"),
  "Mes"      = c("COL1A1","COL1A2","COL5A1","COL6A1","LOX","SPARC","CTHRC1","SNAI2","TGFBI"),
  "Endo"     = c("CDH5","VWF"),
  "Peri"     = c("PDGFRB","RGS5"),
  "TAMs"     = c("CSF1R","CD163","SPP1","APOE"),
  "T cells"  = c("CD3D","CD8A"),
  "B cells"  = c("MS4A1","CD79A"),
  "Mast"     = c("KIT","MS4A2"),
  "Neutr"    = c("S100A8","CSF3R")
)

DefaultAssay(obj) <- "SCT"
DotPlot2(obj, features = grouped_features_fig1) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(face = "italic", size = 10)
  )
