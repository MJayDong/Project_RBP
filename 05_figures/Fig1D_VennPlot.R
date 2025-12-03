############################################################
## Fig1D_VennPlot.R
## 绘制高质量 Venn 图（适合论文）
############################################################

source(file.path("analysis", "0.BASIC.R"))

library(VennDiagram)
library(ggplot2)

plot_venn_sets <- function(
  list_input,
  output_prefix = "Fig1D_Venn"
) {
  venn_plot <- venn.diagram(
    x = list_input,
    filename = NULL,
    fill = c("#1B6CA8", "#E89A1B", "#D64550"),
    alpha = 0.55,
    lty = 0,
    cex = 1.5,
    cat.cex = 1.4,
    cat.fontface = "bold"
  )

  pdf(
    file.path(project_dir, "results", "figures",
              paste0(output_prefix, ".pdf")),
    width = 4, height = 4
  )
  grid::grid.draw(venn_plot)
  dev.off()

  msg("[Fig1D] Done.")
}
