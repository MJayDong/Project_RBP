############################################################
## 8.Cluster_GVis.R
## 使用 ClusterGVis 对 cluster/group 富集结果做统一分析（不出图）
############################################################

source(file.path("analysis", "0.BASIC.R"))

library(Seurat)
library(ClusterGVis)
library(org.Hs.eg.db)
library(dplyr)

# 输入：sc_Obj_Allcells_withBigGroup.rds + markers 表
# 你可以根据需要换成其他 markers 表

run_cluster_gvis <- function(markers_file = "DE_markers_allclusters.csv") {
  msg("[ClusterGVis] Start cluster/group enrichment ...")

  sc_all <- load_rds("sc_Obj_Allcells_withBigGroup.rds")

  # 构建 Group（例如 BigGroup_SampleType）
  sc_all$Group <- paste0(sc_all$BigGroup, "_", sc_all$SampleType)

  Idents(sc_all) <- sc_all$Group

  # 读入 marker 表
  mf <- file.path(project_dir, "results", "tables", markers_file)
  if (!file.exists(mf)) stop("[ClusterGVis] Markers file not found: ", mf)

  markers <- read.csv(mf)

  # markers 表需要至少包含：gene, cluster/group 等信息
  # 为简单起见，这里直接把 cluster 当作 group 名
  if (!all(c("gene", "cluster") %in% colnames(markers))) {
    stop("[ClusterGVis] markers must contain 'gene' and 'cluster' columns.")
  }

  # 使用 ClusterGVis 内置的 enrichment 函数
  # 文档中常用 groupGSEA / clusterGroupEnrichment 等函数
  # 这里假定使用 runGSEA 函数接口
  enrich <- ClusterGVis::runGSEA(
    gene     = markers$gene,
    cluster  = markers$cluster,
    type     = "BP",
    organism = "hsa",
    pvalueCutoff = 0.05,
    topn     = 20,
    seed     = 5201314
  )

  # enrich 是一个长表，包含 group/cluster, pathway, NES 等
  save_rds(enrich, "ClusterGVis_enrichment_table.rds")
  write.csv(
    enrich,
    file.path(project_dir, "results", "tables", "ClusterGVis_enrichment_table.csv"),
    row.names = FALSE
  )

  # 构建 cluster × pathway 的 NES 矩阵
  nes_mat <- enrich %>%
    select(group, term, NES) %>%
    tidyr::pivot_wider(names_from = term, values_from = NES)

  save_rds(nes_mat, "ClusterGVis_NES_matrix.rds")
  write.csv(
    nes_mat,
    file.path(project_dir, "results", "tables", "ClusterGVis_NES_matrix.csv"),
    row.names = FALSE
  )

  msg("[ClusterGVis] Enrichment completed.")
  invisible(list(enrich = enrich, NES_matrix = nes_mat))
}
