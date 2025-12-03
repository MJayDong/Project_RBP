############################################################
## 0.BASIC.R
## 初始化项目环境、目录结构、依赖包、基础函数
## Author: Mingjie
## Project: Project_3 (ER+ Breast Cancer RBP Project)
############################################################

## -------- 1. 项目根目录（参数化路径） -------- ##
project_dir <- "path/to/Project_3"
setwd(project_dir)

## 创建标准化目录结构（如不存在则自动创建）
dir_list <- c(
  "data/raw",
  "data/processed",
  "results",
  "results/figures",
  "results/tables",
  "scripts",
  "scripts/functions",
  "logs"
)
invisible(lapply(dir_list, function(x) {
  if (!dir.exists(x)) dir.create(x, recursive = TRUE)
}))

## -------- 2. 加载依赖包 -------- ##
pkg_list <- c(
  "Seurat", "dplyr", "ggplot2", "Matrix", "data.table",
  "patchwork", "stringr", "magrittr", "tidyr",
  "readr", "purrr"
)

for (pkg in pkg_list) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

message("[BASIC] All packages loaded.")

## -------- 3. 常用功能函数 -------- ##

## 3.1 统一的消息打印格式
msg <- function(...) {
  cat("[", format(Sys.time(), "%H:%M:%S"), "] ", ..., "\n", sep = "")
}

## 3.2 标准化保存函数
save_rds <- function(obj, filename) {
  full_path <- file.path(project_dir, "data/processed", filename)
  saveRDS(obj, full_path)
  msg("Saved RDS:", full_path)
}

load_rds <- function(filename) {
  full_path <- file.path(project_dir, "data/processed", filename)
  msg("Loaded RDS:", full_path)
  readRDS(full_path)
}

## 3.3 标准化绘图导出
save_fig <- function(fig, filename, width = 6, height = 5) {
  ggsave(
    filename = file.path(project_dir, "results/figures", filename),
    plot = fig,
    width = width,
    height = height
  )
  msg("Saved figure:", filename)
}

## -------- 4. metadata 统一处理函数（如需要） -------- ##

merge_metadata <- function(meta_df) {
  meta_df$SampleID <- paste0(meta_df$PatientID, "_", meta_df$Region)
  meta_df$CellState <- paste0(meta_df$CellName, "_", meta_df$SampleType)
  return(meta_df)
}

msg("[BASIC] Initialization complete.")
