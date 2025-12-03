############################################################
## 2.CellType_Basic.R
## 对指定 BigGroup 细胞类型做：子集预处理 + 聚类 + markers + 基础图
## 适用 celltype: "EPITHELIAL", "IMMUNE", "VASCULAR", "MESENCHYMAL"
############################################################

source(file.path("analysis", "0.BASIC.R"))

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

# 映射 celltype 到对象文件名（在 1.0.All_Preprocessing.R 已经保存）
.celltype_to_objfile <- function(celltype) {
  switch(
    toupper(celltype),
    "EPITHELIAL"  = "sc_Obj_EPITHELIAL.rds",
    "IMMUNE"      = "sc_Obj_IMMUNE.rds",
    "VASCULAR"    = "sc_Obj_VASCULAR.rds",
    "MESENCHYMAL" = "sc_Obj_MESENCHYMAL.rds",
    stop("Unknown celltype: ", celltype)
  )
}

# 输出目录: results/Basic/<celltype>/
.basic_output_dir <- function(celltype) {
  out_dir <- file.path(project_dir, "results", "Basic", toupper(celltype))
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  out_dir
}

# 主函数 ---------------------------------------------------------------
run_celltype_basic <- function(celltype,
                               resolution = 0.8,
                               harmony_var = "GSM") {
  celltype <- toupper(celltype)
  msg("[Basic] Celltype = ", celltype)

  # 1. 载入对象 ---------------------------------------------------------
  obj_file <- .celltype_to_objfile(celltype)
  so <- load_rds(obj_file)
  msg("[Basic] Loaded: ", obj_file,
      " (cells = ", ncol(so), ", features = ", nrow(so), ")")

  out_dir <- .basic_output_dir(celltype)

  # 2. 预处理：Normalize / HVG / Scale / PCA / Harmony / UMAP / TSNE / Cluster
  msg("[Basic] Preprocessing & clustering ...")

  so <- NormalizeData(so, normalization.method = "LogNormalize", scale.factor = 10000) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
    ScaleData() %>%
    RunPCA(npcs = 50) %>%
    RunHarmony(harmony_var, plot_convergence = FALSE) %>%
    RunUMAP(reduction = "harmony", dims = 1:30) %>%
    RunTSNE(reduction = "harmony", dims = 1:30) %>%
    FindNeighbors(reduction = "harmony", dims = 1:30) %>%
    FindClusters(resolution = resolution) %>%
    identity()

  # 保存预处理后的对象
  save_rds(so, paste0("sc_Obj_", celltype, "_BasicProcessed.rds"))

  # 3. Cluster markers ---------------------------------------------------
  msg("[Basic] Finding cluster markers ...")

  Idents(so) <- so$seurat_clusters

  markers <- FindAllMarkers(
    so,
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = 0.25
  )

  save_rds(markers, paste0("DE_markers_", celltype, "_clusters.rds"))
  write.csv(
    markers,
    file.path(project_dir, "results", "tables",
              paste0("DE_markers_", celltype, "_clusters.csv")),
    row.names = FALSE
  )

  # 4. SampleType 内的 DEG (Normal vs Primary vs Metastasis) ------------
  if ("SampleType" %in% colnames(so@meta.data)) {
    msg("[Basic] SampleType DEG inside ", celltype, " ...")
    Idents(so) <- so$SampleType

    if (all(c("Normal", "Primary") %in% levels(so$SampleType))) {
      DE_N_vs_P <- FindMarkers(so, ident.1 = "Normal", ident.2 = "Primary",
                               logfc.threshold = 0.25, min.pct = 0.2)
      save_rds(DE_N_vs_P, paste0("DE_", celltype, "_Normal_vs_Primary.rds"))
    }

    if (all(c("Primary", "Metastasis") %in% levels(so$SampleType))) {
      DE_P_vs_M <- FindMarkers(so, ident.1 = "Primary", ident.2 = "Metastasis",
                               logfc.threshold = 0.25, min.pct = 0.2)
      save_rds(DE_P_vs_M, paste0("DE_", celltype, "_Primary_vs_Metastasis.rds"))
    }

    if (all(c("Normal", "Metastasis") %in% levels(so$SampleType))) {
      DE_N_vs_M <- FindMarkers(so, ident.1 = "Normal", ident.2 = "Metastasis",
                               logfc.threshold = 0.25, min.pct = 0.2)
      save_rds(DE_N_vs_M, paste0("DE_", celltype, "_Normal_vs_Metastasis.rds"))
    }
  } else {
    msg("[Basic] No SampleType column found, skip within-celltype DEG.")
  }

  # 5. 基础图：TSNE 按 cluster / SampleType ------------------------------
  msg("[Basic] Plotting basic TSNE ...")

  color_sampletype <- c(
    Normal     = "#1B6CA8",
    Primary    = "#E89A1B",
    Metastasis = "#D64550"
  )

  # tSNE by cluster
  p_cluster <- DimPlot(
    so,
    reduction = "tsne",
    group.by  = "seurat_clusters",
    label     = TRUE,
    repel     = TRUE,
    raster    = FALSE,
    pt.size   = 0.2
  ) + theme(aspect.ratio = 1)

  ggsave(
    filename = file.path(out_dir, paste0("TSNE_", celltype, "_cluster.pdf")),
    plot     = p_cluster,
    width    = 6,
    height   = 6
  )

  # tSNE by SampleType
  if ("SampleType" %in% colnames(so@meta.data)) {
    p_sampletype <- DimPlot(
      so,
      reduction = "tsne",
      group.by  = "SampleType",
      cols      = color_sampletype,
      raster    = FALSE,
      pt.size   = 0.2
    ) + theme(aspect.ratio = 1)

    ggsave(
      filename = file.path(out_dir, paste0("TSNE_", celltype, "_SampleType.pdf")),
      plot     = p_sampletype,
      width    = 6,
      height   = 6
    )
  }

  msg("[Basic] Done for celltype = ", celltype)
  invisible(so)
}
