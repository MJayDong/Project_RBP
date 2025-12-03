############################################################
## 13.MultiNicheNet_Run.R
## MultiNicheNet 分析（无图）—— 输出 ligand 活性 + 贡献矩阵
############################################################

source(file.path("analysis", "0.BASIC.R"))

library(Seurat)
library(dplyr)
library(nichenetr)
library(tibble)

run_MultiNicheNet <- function(
  sc_obj_file,
  receiver_celltype,   # 目标细胞类型，如 “EPITHELIAL”
  sender_celltypes,    # 提供信号的多个细胞类型
  output_prefix = "MultiNicheNet"
) {
  so <- load_rds(sc_obj_file)

  # Receiver subset
  so_r <- subset(so, subset = BigGroup == receiver_celltype)

  # Sender subsets
  so_s_list <- lapply(sender_celltypes, function(ct) {
    subset(so, subset = BigGroup == ct)
  })
  names(so_s_list) <- sender_celltypes

  # 获取表达矩阵
  rec_exp <- as.matrix(so_r@assays$RNA@data)

  sender_exp_list <- lapply(so_s_list, function(obj) {
    as.matrix(obj@assays$RNA@data)
  })

  # MultiNicheNet 主流程：Ligand activity + target prediction
  ligand_activities <- lapply(sender_exp_list, function(exp) {
    predict_ligand_activities(
      receiver_expression = rec_exp,
      sender_expression   = exp
    )
  })

  # 合并 ligand 活性
  ligand_activity_table <- bind_rows(
    lapply(names(ligand_activities), function(sender) {
      df <- ligand_activities[[sender]]
      df$Sender <- sender
      df
    })
  )

  save_rds(ligand_activity_table, paste0(output_prefix, "_ligand_activity.rds"))
  write.csv(
    ligand_activity_table,
    file.path(project_dir, "results", "tables",
              paste0(output_prefix, "_ligand_activity.csv")),
    row.names = FALSE
  )

  # 多 celltype 的 ligand–target 推断
  LT_list <- lapply(sender_exp_list, function(exp) {
    nichenetr::infer_ligand_target_links(
      receiver_expression = rec_exp,
      sender_expression   = exp
    )
  })
  names(LT_list) <- sender_celltypes

  # 合并 ligand-target 贡献表
  LT_all <- bind_rows(
    lapply(names(LT_list), function(sender) {
      d <- LT_list[[sender]]
      d$Sender <- sender
      d
    })
  )

  save_rds(LT_all, paste0(output_prefix, "_ligand_target_links.rds"))
  write.csv(
    LT_all,
    file.path(project_dir, "results", "tables",
              paste0(output_prefix, "_ligand_target_links.csv")),
    row.names = FALSE
  )

  msg("[MultiNicheNet] Done.")
  invisible(list(
    ligand_activity = ligand_activity_table,
    ligand_target   = LT_all
  ))
}
