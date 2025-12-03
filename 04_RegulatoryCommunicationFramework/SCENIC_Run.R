############################################################
## 11.SCENIC_Run.R
## SCENIC 全流程（无图）—— 输出 Regulon、AUC、ctx
############################################################

source(file.path("analysis", "0.BASIC.R"))

library(SCENIC)
library(AUCell)
library(RcisTarget)
library(Seurat)
library(dplyr)

# --------------------------------------------------------------------
# 输入：
#   sc_obj_file: Seurat 对象（建议为某个 celltype 子集）
#   org: "mgi" / "hgnc"
#   db_dir: SCENIC 数据库路径
# 输出：Regulons, AUC, ctx文件
# --------------------------------------------------------------------

run_SCENIC <- function(
  sc_obj_file,
  org       = "hgnc",
  db_dir    = file.path(project_dir, "SCENIC", "databases"),
  output_prefix = "SCENIC"
) {
  so <- load_rds(sc_obj_file)

  exprMat <- as.matrix(so@assays$RNA@counts)
  genesKept <- geneFiltering(exprMat, minCountsPerGene = 3, minSamples = 0.03 * ncol(exprMat))
  exprMat <- exprMat[genesKept, ]

  scenicOptions <- initializeScenic(
    org = org,
    dbDir = db_dir,
    nCores = 8
  )
  saveRDS(scenicOptions, paste0(output_prefix, "_scenicOptions_step1.rds"))

  # Step1: 共表达 + GENIE3
  runCorrelation(exprMat, scenicOptions)
  exprMat_filtered <- exprMat
  runGenie3(exprMat_filtered, scenicOptions)

  scenicOptions <- readRDS(paste0(output_prefix, "_scenicOptions_step1.rds"))

  # Step2: RcisTarget 获得 regulons
  runSCENIC_1_coexNetwork2modules(scenicOptions)
  runSCENIC_2_createRegulons(scenicOptions)
  runSCENIC_3_scoreCells(scenicOptions, exprMat)

  aucell_regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
  regulons           <- loadInt(scenicOptions, "regulons")
  ctx                <- loadInt(scenicOptions, "regulonTargetsInfo")

  save_rds(regulons,        paste0(output_prefix, "_regulons.rds"))
  save_rds(aucell_regulonAUC, paste0(output_prefix, "_regulonAUC.rds"))
  save_rds(ctx,             paste0(output_prefix, "_ctx_table.rds"))

  write.csv(ctx,
    file.path(project_dir, "results", "tables",
              paste0(output_prefix, "_ctx_table.csv")),
    row.names = FALSE
  )

  msg("[SCENIC] Done.")
  invisible(list(regulons = regulons, AUC = aucell_regulonAUC, ctx = ctx))
}
