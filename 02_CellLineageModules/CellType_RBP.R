############################################################
## 3.CellType_RBP.R
## 对指定细胞大类做 RBP 相关基因集准备 & 交集统计
############################################################

source(file.path("analysis", "0.BASIC.R"))

library(Seurat)
library(dplyr)
library(readr)
library(VennDiagram)

# RBP 基因清单路径（人类）
rbp_gene_file <- file.path(project_dir, "data", "raw", "Homo_sapiens_RBPs.csv")

.rbp_output_dir <- function(celltype) {
  out_dir <- file.path(project_dir, "results", "RBP", toupper(celltype))
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  out_dir
}

.celltype_to_basic_objfile <- function(celltype) {
  paste0("sc_Obj_", toupper(celltype), "_BasicProcessed.rds")
}

run_celltype_rbp <- function(celltype) {
  celltype <- toupper(celltype)
  msg("[RBP] Celltype = ", celltype)

  out_dir <- .rbp_output_dir(celltype)

  # 1. 载入预处理后的对象 ---------------------------------------------
  obj_file <- .celltype_to_basic_objfile(celltype)
  so <- load_rds(obj_file)

  # 2. 载入 RBP 列表 ---------------------------------------------------
  if (!file.exists(rbp_gene_file)) {
    stop("[RBP] RBP gene file not found: ", rbp_gene_file)
  }
  rbp_tbl <- read_csv(rbp_gene_file, show_col_types = FALSE)
  rbp_genes <- unique(rbp_tbl$`Gene symbol`)

  features_all  <- rownames(so)
  features_hvg  <- VariableFeatures(so)

  genes_rbp_all <- intersect(features_all, rbp_genes)
  genes_rbp_hvg <- intersect(features_hvg, rbp_genes)

  msg("[RBP] Total features: ", length(features_all),
      "; HVG: ", length(features_hvg),
      "; RBP all-intersect: ", length(genes_rbp_all),
      "; RBP HVG-intersect: ", length(genes_rbp_hvg))

  # 保存基因集
  save_rds(genes_rbp_all, paste0("RBP_genes_", celltype, "_allFeatures.rds"))
  save_rds(genes_rbp_hvg, paste0("RBP_genes_", celltype, "_HVG.rds"))

  write.csv(
    data.frame(Gene = genes_rbp_all),
    file.path(project_dir, "results", "tables",
              paste0("RBP_genes_", celltype, "_allFeatures.csv")),
    row.names = FALSE
  )

  write.csv(
    data.frame(Gene = genes_rbp_hvg),
    file.path(project_dir, "results", "tables",
              paste0("RBP_genes_", celltype, "_HVG.csv")),
    row.names = FALSE
  )

  # 3. 画一个简单的 Venn（All / HVG / RBP） ---------------------------
  venn_input <- list(
    AllFeatures = features_all,
    HVG         = features_hvg,
    RBP         = rbp_genes
  )

  venn_plot <- venn.diagram(
    x = venn_input,
    filename = NULL,
    fill = c("#B0C4DE", "#F5DEB3", "#FFC0CB"),
    alpha = 0.5,
    cex   = 1.4,
    cat.cex = 1.2,
    margin = 0.1
  )

  pdf(file.path(out_dir, paste0("Venn_RBP_", celltype, ".pdf")),
      width = 6, height = 6)
  grid::grid.draw(venn_plot)
  dev.off()

  msg("[RBP] RBP gene sets + Venn saved for ", celltype)
  invisible(list(
    genes_rbp_all  = genes_rbp_all,
    genes_rbp_hvg  = genes_rbp_hvg
  ))
}
