############################################################
## SupplyFig2_CustomPlots.R
## 补充图（可按需要扩展）—— dotplot / violin / heatmap
############################################################

source(file.path("analysis", "0.BASIC.R"))

library(Seurat)
library(ggplot2)
library(pheatmap)

# dotplot --------------------------------------------------------------
plot_dot_custom <- function(
  sc_obj_file,
  genes,
  group_var = "seurat_clusters",
  output_prefix = "SupplyFig2_dot"
) {
  so <- load_rds(sc_obj_file)

  p <- DotPlot(so, features = genes, group.by = group_var) +
    RotatedAxis() +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      aspect.ratio = 1
    )

  ggsave(
    file.path(project_dir, "results", "figures",
              paste0(output_prefix, ".pdf")),
    p, width = 6, height = 4
  )
}

# violin ---------------------------------------------------------------
plot_violin_custom <- function(
  sc_obj_file,
  genes,
  group_var = "SampleType",
  output_prefix = "SupplyFig2_violin"
) {
  so <- load_rds(sc_obj_file)

  p <- VlnPlot(so, features = genes, group.by = group_var, pt.size = 0) +
    theme_classic() +
    theme(aspect.ratio = 1)

  ggsave(
    file.path(project_dir, "results", "figures",
              paste0(output_prefix, ".pdf")),
    p, width = 6, height = 4
  )
}

# heatmap ---------------------------------------------------------------
plot_heatmap_custom <- function(
  mat,                     # 基因 × 样本矩阵
  output_prefix = "SupplyFig2_heatmap"
) {
  pdf(
    file.path(project_dir, "results", "figures",
              paste0(output_prefix, ".pdf")),
    width = 4, height = 6
  )
  pheatmap(
    mat,
    scale = "row",
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_rownames = TRUE,
    show_colnames = TRUE,
    color = colorRampPalette(c("navy", "white", "firebrick3"))(100)
  )
  dev.off()
}
