############################################################
## 1.2.All_Visualization.R
## All cells — 所有可视化图（聚类图 / 细胞组成 / marker / top5 / venn）
## Author: Mingjie
############################################################

source(file.path("analysis", "0.BASIC.R"))

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(VennDiagram)

msg("[All_Vis] Start visualization module ...")

## ============================
## 0. 创建输出目录
## ============================
fig_dir <- file.path(project_dir, "results", "figures", "All")
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)


## ============================
## 1. 载入对象和表格
## ============================
sc_all  <- load_rds("sc_Obj_Allcells_withBigGroup.rds")
sample_big_count <- load_rds("All_sample_biggroup_counts.rds")
sample_big_pct   <- load_rds("All_sample_biggroup_percentages.rds")

top5_all <- load_rds("DE_top5_allclusters.rds")

DE_N_vs_P <- load_rds("DE_Normal_vs_Primary.rds")
DE_P_vs_M <- load_rds("DE_Primary_vs_Metastasis.rds")
DE_N_vs_M <- load_rds("DE_Normal_vs_Metastasis.rds")

msg("[All_Vis] Data loaded.")


## ============================
## 2. 自定义配色方案（可统一修改）
## ============================
color_SampleType <- c(
  Normal     = "#1B6CA8",
  Primary    = "#E89A1B",
  Metastasis = "#D64550"
)

color_BigGroup <- c(
  EPITHELIAL       = "#498EA4",
  IMMUNE           = "#E54924",
  MESENCHYMAL      = "#8F7E4F",
  VASCULAR-RELATED = "#6A5ACD",
  OTHER            = "#999999"
)


## ============================
## 3. 聚类图（DimPlot）
## ============================
msg("[All_Vis] Plotting DimPlot ...")

fig_tsne_sampletype <- DimPlot(
  sc_all,
  reduction = "tsne",
  group.by  = "SampleType",
  cols      = color_SampleType,
  raster    = FALSE,
  pt.size   = 0.5
) + theme(aspect.ratio = 1)

save_fig(fig_tsne_sampletype, "TSNE_by_SampleType.pdf", 6, 6)


fig_tsne_biggroup <- DimPlot(
  sc_all,
  reduction = "tsne",
  group.by  = "BigGroup",
  cols      = color_BigGroup,
  raster    = FALSE,
  pt.size   = 0.5
) + theme(aspect.ratio = 1)

save_fig(fig_tsne_biggroup, "TSNE_by_BigGroup.pdf", 6, 6)

msg("[All_Vis] DimPlot saved.")


## ============================
## 4. 细胞组成图（SampleType × BigGroup）
## ============================

## --- 堆叠百分比条形图（按 SampleType） ---
p_pct <- ggplot(sample_big_pct, aes(
  x = SampleType,
  y = Percentage,
  fill = BigGroup
)) +
  geom_bar(stat = "identity", position = "fill", color = "grey20") +
  scale_fill_manual(values = color_BigGroup) +
  theme_minimal() +
  ylab("Percentage") +
  theme(
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12)
  )

save_fig(p_pct, "Composition_Percentage.pdf", 6, 6)


## --- 计数条形图（按 SampleType） ---
p_count <- ggplot(sample_big_count, aes(
  x = SampleType,
  y = CellNumber,
  fill = BigGroup
)) +
  geom_bar(stat = "identity", position = "dodge", color = "grey10") +
  scale_fill_manual(values = color_BigGroup) +
  theme_minimal() +
  ylab("Cell Number") +
  theme(
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12)
  )

save_fig(p_count, "Composition_CellNumber.pdf", 6, 6)

msg("[All_Vis] Composition plots saved.")


## ============================
## 5. Marker / Top5 的 DotPlot
## ============================
msg("[All_Vis] Plotting Top5 marker DotPlot ...")

features_top5 <- unique(top5_all$gene)

fig_top5 <- DotPlot(
  sc_all,
  features = features_top5,
  cols = c("#498EA4", "#E54924")  # 蓝-红
) +
  RotatedAxis() +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 7, angle = 90),
    axis.text.y = element_text(size = 8)
  )

save_fig(fig_top5, "Top5_Markers_DotPlot.pdf", 14, 7)

msg("[All_Vis] Top5 marker dotplot saved.")


## ============================
## 6. 三组差异基因的 Venn diagram
## ============================

msg("[All_Vis] Plotting DEG Venn diagram ...")

gene_NP <- rownames(DE_N_vs_P)
gene_PM <- rownames(DE_P_vs_M)
gene_NM <- rownames(DE_N_vs_M)

venn_deg <- venn.diagram(
  x = list(
    Normal_vs_Primary     = gene_NP,
    Primary_vs_Metastasis = gene_PM,
    Normal_vs_Metastasis  = gene_NM
  ),
  filename = NULL,
  fill = c("#1B6CA8", "#E89A1B", "#D64550"),
  alpha = 0.6,
  cex = 1.4,
  cat.fontface = "bold",
  margin = 0.1
)

pdf(file.path(fig_dir, "DEG_Venn.pdf"), width = 6, height = 6)
grid::grid.draw(venn_deg)
dev.off()

msg("[All_Vis] Venn diagram saved.")


## ============================
## END
## ============================
msg("[All_Vis] All visualization complete.")
