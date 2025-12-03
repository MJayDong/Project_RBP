############################################################
## Fig1A_Tsne_UMAP.R
## 绘制 tSNE / UMAP（适合作为论文 Figure）
############################################################

source(file.path("analysis", "0.BASIC.R"))

library(Seurat)
library(ggplot2)

plot_tsne_umap <- function(
  sc_obj_file,
  color_by = "BigGroup",
  plot_tsne = TRUE,
  plot_umap = TRUE,
  output_prefix = "Fig1A"
) {
  so <- load_rds(sc_obj_file)

  # 主题（适合出版）
  theme_pub <- theme_classic() +
    theme(
      axis.text  = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      legend.title = element_blank(),
      legend.position = "right",
      aspect.ratio = 1
    )

  if (plot_tsne) {
    p_tsne <- DimPlot(
      so,
      reduction = "tsne",
      group.by = color_by,
      pt.size = 0.2,
      raster = FALSE
    ) + theme_pub

    ggsave(
      filename = file.path(project_dir, "results", "figures",
                           paste0(output_prefix, "_TSNE_", color_by, ".pdf")),
      plot = p_tsne,
      width = 4.5, height = 4.5
    )
  }

  if (plot_umap) {
    p_umap <- DimPlot(
      so,
      reduction = "umap",
      group.by = color_by,
      pt.size = 0.2,
      raster = FALSE
    ) + theme_pub

    ggsave(
      filename = file.path(project_dir, "results", "figures",
                           paste0(output_prefix, "_UMAP_", color_by, ".pdf")),
      plot = p_umap,
      width = 4.5, height = 4.5
    )
  }

  msg("[Fig1A] Done.")
}
