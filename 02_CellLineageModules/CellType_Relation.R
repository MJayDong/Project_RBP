############################################################
## 4.CellType_Relation.R
## TF–Target–RBP 关系的交集分析（scMLnet + SCENIC + RBP）
############################################################

source(file.path("analysis", "0.BASIC.R"))

library(readr)
library(dplyr)

rbp_gene_file <- file.path(project_dir, "data", "raw", "Homo_sapiens_RBPs.csv")

.relation_output_dir <- function(celltype) {
  out_dir <- file.path(project_dir, "results", "Relation", toupper(celltype))
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  out_dir
}

# 这个模块主要是“读文件 → overlap → 导出表”
# 由于不同细胞类型的 scMLnet 文件名不一样，所以用参数传路径

run_celltype_relation <- function(celltype,
                                  scmlnet_file,
                                  scenic_file) {
  celltype <- toupper(celltype)
  out_dir  <- .relation_output_dir(celltype)

  msg("[Relation] Celltype = ", celltype)
  msg("[Relation] scMLnet file: ", scmlnet_file)
  msg("[Relation] SCENIC file: ", scenic_file)

  if (!file.exists(scmlnet_file)) stop("scMLnet file not found: ", scmlnet_file)
  if (!file.exists(scenic_file))  stop("SCENIC file not found: ", scenic_file)
  if (!file.exists(rbp_gene_file)) stop("RBP gene file not found: ", rbp_gene_file)

  # 1. 读取 RBP 列表
  rbp_tbl   <- read_csv(rbp_gene_file, show_col_types = FALSE)
  rbp_genes <- unique(rbp_tbl$`Gene symbol`)

  # 2. 读取 scMLnet (TF_Target) ---------------------------------------
  scMLnet <- read.table(scmlnet_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  colnames(scMLnet) <- "Pair"

  scMLnet$TF         <- gsub("_.*", "", scMLnet$Pair)
  scMLnet$TargetGene <- gsub(".*_", "", scMLnet$Pair)

  # 3. 读取 SCENIC TF regulon 信息 -------------------------------------
  scenic_ctx <- read_csv(scenic_file, skip = 2, col_names = FALSE, show_col_types = FALSE)
  colnames(scenic_ctx) <- c(
    "TF", "MotifID", "AUC", "NES", "MotifSimilarityQ", "OrthologousIdentity",
    "Annotation", "Context", "TargetGenes", "RankAtMax"
  )

  # 4. Overlap: TF & Target & RBP --------------------------------------
  overlap_TF <- intersect(scMLnet$TF, scenic_ctx$TF)

  msg("[Relation] Overlap TF count: ", length(overlap_TF))

  overlap_regulon <- scMLnet[scMLnet$TF %in% overlap_TF, ]
  unique_TF <- unique(overlap_regulon$TF)

  msg("[Relation] Regulon TFs: ", paste(unique_TF, collapse = ", "))

  rbp_targets <- intersect(overlap_regulon$TargetGene, rbp_genes)

  msg("[Relation] Overlap TF–Target–RBP genes: ", length(rbp_targets))

  res_tbl <- overlap_regulon %>%
    filter(TargetGene %in% rbp_targets) %>%
    select(TF, TargetGene)

  # 保存结果表 --------------------------------------------------------
  save_rds(res_tbl, paste0("Relation_TF_RBP_", celltype, ".rds"))
  write.csv(
    res_tbl,
    file.path(project_dir, "results", "tables",
              paste0("Relation_TF_RBP_", celltype, ".csv")),
    row.names = FALSE
  )

  msg("[Relation] Relation table saved for ", celltype)
  invisible(res_tbl)
}
