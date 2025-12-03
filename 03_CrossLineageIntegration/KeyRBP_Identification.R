############################################################
## 7.KeyRBP_Identification.R
## 识别全局关键 RBP：汇总 DEG 特征 + 网络特征 → 生成 RBP 排名表
############################################################

source(file.path("analysis", "0.BASIC.R"))

library(dplyr)
library(readr)
library(stringr)
library(randomForest)

rbp_gene_file <- file.path(project_dir, "data", "raw", "Homo_sapiens_RBPs.csv")

# 读取所有与 RBP / DEG / Relation 相关的表 ------------------------------

.load_rbp_genes <- function() {
  rbp_tbl <- read_csv(rbp_gene_file, show_col_types = FALSE)
  unique(rbp_tbl$`Gene symbol`)
}

# 加载前面模块生成的 Relation_TF_RBP_* 表
.load_relation_tables <- function() {
  rel_dir <- file.path(project_dir, "results", "tables")
  files   <- list.files(rel_dir, pattern = "^Relation_TF_RBP_.*\\.csv$", full.names = TRUE)
  if (length(files) == 0) return(NULL)
  rel_list <- lapply(files, read.csv)
  names(rel_list) <- gsub("^Relation_TF_RBP_|\\.csv$", "", basename(files))
  bind_rows(
    lapply(names(rel_list), function(ct) {
      df <- rel_list[[ct]]
      df$CellType <- ct
      df
    })
  )
}

# 加载所有与 RBP 相关的 DEG 表（例如 DE_markers_*_clusters.csv 等）
.load_deg_tables <- function() {
  tab_dir <- file.path(project_dir, "results", "tables")
  files   <- list.files(tab_dir, pattern = "^DE_.*\\.csv$", full.names = TRUE)
  if (length(files) == 0) return(NULL)
  deg_list <- lapply(files, read.csv)
  names(deg_list) <- basename(files)
  bind_rows(
    lapply(names(deg_list), function(fn) {
      df <- deg_list[[fn]]
      df$Source <- fn
      df
    })
  )
}

# 主函数：计算 RBP 特征 & 排名 -------------------------------------

run_keyRBP_identification <- function() {
  msg("[KeyRBP] Start identification ...")

  rbp_genes <- .load_rbp_genes()
  msg("[KeyRBP] RBP gene number: ", length(rbp_genes))

  deg_all <- .load_deg_tables()
  if (is.null(deg_all)) stop("[KeyRBP] No DEG tables found in results/tables.")

  rel_all <- .load_relation_tables()

  # DEG 特征：以 gene 为单位汇总
  # 假定 DEG 表至少有列：gene, avg_log2FC, p_val / p_val_adj
  if (!"gene" %in% colnames(deg_all)) {
    stop("[KeyRBP] DEG tables must contain 'gene' column.")
  }

  deg_all <- deg_all %>%
    mutate(is_RBP = gene %in% rbp_genes)

  deg_rbp <- deg_all %>%
    filter(is_RBP)

  msg("[KeyRBP] DEG-RBP rows: ", nrow(deg_rbp))

  deg_summary <- deg_rbp %>%
    group_by(gene) %>%
    summarise(
      n_DE_context    = n_distinct(Source),
      mean_log2FC     = mean(avg_log2FC, na.rm = TRUE),
      max_log2FC      = max(avg_log2FC, na.rm = TRUE),
      min_p_val       = suppressWarnings(min(p_val, na.rm = TRUE)),
      min_p_val_adj   = if ("p_val_adj" %in% colnames(deg_rbp))
        suppressWarnings(min(p_val_adj, na.rm = TRUE)) else NA_real_,
      .groups = "drop"
    )

  # 网络特征：TF–RBP 关系的出现次数
  if (!is.null(rel_all)) {
    net_summary <- rel_all %>%
      group_by(TargetGene) %>%
      summarise(
        n_TF_link    = n(),
        n_TF_unique  = n_distinct(TF),
        .groups = "drop"
      ) %>%
      rename(gene = TargetGene)
  } else {
    net_summary <- data.frame(gene = character(0), n_TF_link = integer(0), n_TF_unique = integer(0))
  }

  # 合并特征
  features <- full_join(
    deg_summary,
    net_summary,
    by = "gene"
  ) %>%
    mutate(
      n_TF_link   = ifelse(is.na(n_TF_link), 0, n_TF_link),
      n_TF_unique = ifelse(is.na(n_TF_unique), 0, n_TF_unique)
    )

  # 只保留真实 RBP
  features <- features %>%
    filter(gene %in% rbp_genes)

  # 构建一个简单的综合分数（可按需要改）
  # 分数 = 标准化(mean_log2FC, n_DE_context, n_TF_link, n_TF_unique)之和
  scale_safe <- function(x) {
    if (all(is.na(x))) return(rep(0, length(x)))
    as.numeric(scale(x))
  }

  features <- features %>%
    mutate(
      score_logFC   = scale_safe(mean_log2FC),
      score_nDE     = scale_safe(n_DE_context),
      score_nTF     = scale_safe(n_TF_link),
      score_nTFuniq = scale_safe(n_TF_unique),
      KeyRBP_Score  = score_logFC + score_nDE + score_nTF + score_nTFuniq
    ) %>%
    arrange(desc(KeyRBP_Score))

  # 保存结果
  save_rds(features, "KeyRBP_features_allcelltypes.rds")
  write.csv(
    features,
    file.path(project_dir, "results", "tables", "KeyRBP_features_allcelltypes.csv"),
    row.names = FALSE
  )

  msg("[KeyRBP] Features table saved. Top 10 RBP:")
  print(head(features %>% select(gene, KeyRBP_Score), 10))

  invisible(features)
}
