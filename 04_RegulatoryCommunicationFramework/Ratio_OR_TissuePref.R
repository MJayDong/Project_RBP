############################################################
## 10.Ratio_OR_TissuePref.R
## 细胞比例 & OR 富集分析（无图）
############################################################

source(file.path("analysis", "0.BASIC.R"))

library(Seurat)
library(dplyr)

# 1) 细胞比例计算 ---------------------------------------------------

# 比如：state_var = "StateGroup", split_var = "SampleType" 或 "PatientID"
compute_cell_ratio <- function(so, state_var, split_var) {
  if (!state_var %in% colnames(so@meta.data))
    stop("state_var not in meta.data: ", state_var)
  if (!split_var %in% colnames(so@meta.data))
    stop("split_var not in meta.data: ", split_var)

  tab <- table(so@meta.data[[state_var]], so@meta.data[[split_var]])
  ratio <- prop.table(tab, margin = 2)  # 每列标准化

  ratio_df <- as.data.frame(ratio)
  colnames(ratio_df) <- c("State", "Split", "Ratio")
  ratio_df
}

# 2) 简单 OR 计算（2×2）---------------------------------------------

compute_or_2x2 <- function(a, b, c, d) {
  # 2×2：
  #         condition+
  # group+     a      b
  # group-     c      d
  # OR = (a*d)/(b*c)
  or <- (a * d) / (b * c)
  # Fisher 检验 p 值
  p  <- fisher.test(matrix(c(a, b, c, d), nrow = 2))$p.value
  list(OR = or, p_value = p)
}

# 3) 基于 meta.data 的 OR 批量计算（按状态 vs 其他）------------------

# 举例：
#  - group_var: "SampleType" (比如 Normal/Primary/Metastasis)
#  - target_group: "Metastasis"
#  - state_var: "StateGroup" 或 "CellName"
run_ratio_or_analysis <- function(
  sc_obj_file = "sc_Obj_Allcells_withBigGroup.rds",
  state_var   = "StateGroup",
  group_var   = "SampleType",
  target_group = "Metastasis"
) {
  so <- load_rds(sc_obj_file)

  if (!all(c(state_var, group_var) %in% colnames(so@meta.data))) {
    stop("state_var or group_var not found in meta.data.")
  }

  ratio_df <- compute_cell_ratio(so, state_var, group_var)
  save_rds(ratio_df, paste0("Ratio_", state_var, "_by_", group_var, ".rds"))
  write.csv(
    ratio_df,
    file.path(project_dir, "results", "tables",
              paste0("Ratio_", state_var, "_by_", group_var, ".csv")),
    row.names = FALSE
  )

  # 对每个状态计算：目标组 vs 非目标组的 OR
  states <- unique(ratio_df$State)
  or_list <- lapply(states, function(st) {
    # 构建 2×2
    # a: 目标组 & 该状态
    # b: 目标组 & 非该状态
    # c: 非目标组 & 该状态
    # d: 非目标组 & 非该状态
    md <- so@meta.data
    cond_state <- md[[state_var]] == st
    cond_group <- md[[group_var]] == target_group

    a <- sum(cond_state & cond_group)
    b <- sum(!cond_state & cond_group)
    c <- sum(cond_state & !cond_group)
    d <- sum(!cond_state & !cond_group)

    if (min(a, b, c, d) == 0) {
      # 避免 0，简单加 0.5 校正
      a <- a + 0.5; b <- b + 0.5; c <- c + 0.5; d <- d + 0.5
    }

    res <- compute_or_2x2(a, b, c, d)
    data.frame(
      State      = st,
      GroupVar   = group_var,
      Target     = target_group,
      a = a, b = b, c = c, d = d,
      OR         = res$OR,
      p_value    = res$p_value,
      stringsAsFactors = FALSE
    )
  }) %>% bind_rows()

  or_list$FDR <- p.adjust(or_list$p_value, method = "BH")

  save_rds(or_list, paste0("OR_", state_var, "_vs_", target_group, ".rds"))
  write.csv(
    or_list,
    file.path(project_dir, "results", "tables",
              paste0("OR_", state_var, "_vs_", target_group, ".csv")),
    row.names = FALSE
  )

  msg("[Ratio/OR] Done.")
  invisible(list(ratio = ratio_df, OR = or_list))
}
