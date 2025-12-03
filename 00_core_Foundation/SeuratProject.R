############################################################
## 0.SeuratProject.R
## 从 GSE161529 原始矩阵构建 & 整合 Seurat 对象 (All cells)
## Author: Mingjie
############################################################

## 依赖：先运行 0.BASIC.R，完成 project_dir 等基础设置
source(file.path("analysis", "0.BASIC.R"))

## -------- 1. 加载额外包 -------- ##
library(Seurat)
library(Matrix)
library(harmony)
library(stringr)
library(dplyr)
library(ggplot2)

msg("[SeuratProject] Start building Seurat project for GSE161529 ...")

## -------- 2. 定义原始数据目录 & 读入 10X 矩阵 -------- ##
## 假设你把 GSE161529 的 matrix/barcode/features 放在：
##   data/raw/GSE161529/
raw_dir <- file.path(project_dir, "data", "raw", "GSE161529")
if (!dir.exists(raw_dir)) {
  stop("Raw data directory not found: ", raw_dir)
}

setwd(raw_dir)

## 读取所有以 GSM 开头的文件（barcodes / matrix）
ls_file <- list.files("./", "^GSM")

## 这里沿用你原来的设计：每两个文件是一对 barcodes / mtx，
## genes 使用同一个 GSE161529_features.tsv.gz
scList <- vector("list", length = length(ls_file) / 2)
for (i in seq_len(length(ls_file) / 2)) {
  scList[[i]] <- read10X(
    mtx     = ls_file[2 * i],
    genes   = "GSE161529_features.tsv.gz",
    barcodes= ls_file[2 * i - 1],
    DGEList = TRUE
  )
}
names(scList) <- unique(str_split(ls_file, "_", simplify = TRUE)[, 1])

## -------- 3. 定义样本分组信息 -------- ##
## 这里直接用你原脚本里的 GSM 分组向量（请从原脚本 copy 过来）
## ↓↓↓ 这一块你用原始代码替换掉即可 ↓↓↓
Normal        <- c(
  # "GSM4909253", "GSM4909254", ...
)
Primary       <- c(
  # "GSM4909296", "GSM4909297", ...
)
Metastasis    <- c(
  # "GSM4909307", "GSM4909308", ...
)
InvolvedLymph <- c(
  # "GSM4909308","GSM4909310","GSM4909312","GSM4909314","GSM4909316","GSM4909318","GSM4909321"
)
Patient       <- c(
  # "0092", "0019", ...  # 跟 names(scList) 的顺序对应
)
## ↑↑↑ 这一块你用原始代码替换掉即可 ↑↑↑

if (length(Patient) != length(scList)) {
  warning("Length of Patient vector does not match number of samples.")
}

## -------- 4. 处理重复基因符号 & 重命名细胞列名 -------- ##
## 去除 duplicated Symbol，保留 counts 总和更大的那一行
for (i in seq_along(scList)) {
  judge <- identical(rownames(scList[[i]]$counts), rownames(scList[[i]]$genes))
  if (!judge) {
    stop("Row names of counts and genes not identical for sample: ", names(scList)[i])
  }
  
  idx_dup   <- which(duplicated(scList[[i]]$genes$Symbol))
  v_drop    <- numeric(length(idx_dup))
  
  for (j in seq_along(idx_dup)) {
    idx_cur   <- idx_dup[j]
    idx_first <- match(scList[[i]]$genes$Symbol[idx_cur], scList[[i]]$genes$Symbol)
    if (sum(scList[[i]]$counts[idx_cur, ]) >= sum(scList[[i]]$counts[idx_first, ])) {
      v_drop[j] <- idx_first
    } else {
      v_drop[j] <- idx_cur
    }
  }
  
  scList[[i]]$counts <- scList[[i]]$counts[-v_drop, ]
  scList[[i]]$genes  <- scList[[i]]$genes[-v_drop, ]
  rownames(scList[[i]]$counts) <- scList[[i]]$genes$Symbol
  
  msg("Sample ", names(scList)[i], " duplicated genes cleaned.")
}

## 为每个 sample 赋 SampleType & CellSource，并重命名列：
## SampleType_GSM_CellSource_00001
for (i in seq_along(scList)) {
  gsm_id <- names(scList)[i]
  
  if (gsm_id %in% Normal) {
    sample_type <- "Normal"
  } else if (gsm_id %in% Primary) {
    sample_type <- "Primary"
  } else if (gsm_id %in% Metastasis) {
    sample_type <- "Metastasis"
  } else {
    sample_type <- "Undetermined"
  }
  
  if (gsm_id %in% InvolvedLymph) {
    cell_source <- "InvolvedLymph"
  } else if (gsm_id %in% Metastasis) {
    cell_source <- "Metastasis"
  } else if (gsm_id %in% Primary) {
    cell_source <- "Primary"
  } else if (gsm_id %in% Normal) {
    cell_source <- "Normal"
  } else {
    cell_source <- "Undetermined"
  }
  
  n_cells <- ncol(scList[[i]]$counts)
  new_names <- sprintf(
    "%s_%s_%s_%05d",
    sample_type, gsm_id, cell_source, seq_len(n_cells)
  )
  colnames(scList[[i]]$counts) <- new_names
  
  msg("Renamed cells for sample ", gsm_id, " as ", sample_type, "_", cell_source)
}

## -------- 5. 构建每个样本的 Seurat 对象 + QC + 归一化 -------- ##
sc_Obj <- lapply(
  X = seq_along(scList),
  FUN = function(i) {
    counts_mat <- scList[[i]]$counts
    obj <- CreateSeuratObject(
      counts    = counts_mat,
      min.cells = 3,
      project   = names(scList)[i]
    )
    
    obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
    
    ## 你原脚本的 QC 标准：
    ## percent.mt < 20
    ## 200 < nFeature_RNA < 10000
    ## nCount_RNA > 500
    obj <- subset(
      obj,
      subset = percent.mt < 20 &
               nFeature_RNA > 200 & nFeature_RNA < 10000 &
               nCount_RNA > 500
    )
    
    obj <- NormalizeData(obj)
    obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
    
    obj
  }
)

names(sc_Obj) <- names(scList)

## 保存中间对象（可选）
save_rds(sc_Obj, "sc_Obj_list.rds")

## -------- 6. 合并所有样本 & 批次校正 & 聚类 -------- ##
msg("[SeuratProject] Merging all samples ...")

sc_ObjCombined <- merge(
  x = sc_Obj[[1]],
  y = sc_Obj[2:length(sc_Obj)]
)

## 加入 meta 信息：SampleType / CellSource / GSM / Patient
meta <- sc_ObjCombined@meta.data
meta$CellID <- rownames(meta)

meta$GSM <- sapply(strsplit(meta$CellID, "_"), `[`, 2)
meta$SampleType  <- sapply(strsplit(meta$CellID, "_"), `[`, 1)
meta$CellSource  <- sapply(strsplit(meta$CellID, "_"), `[`, 3)

## 映射 PatientID（假定 Patient 向量顺序与 names(scList) 对应）
sample_patient_df <- data.frame(
  GSM       = names(scList),
  PatientID = Patient,
  stringsAsFactors = FALSE
)

meta <- meta %>%
  left_join(sample_patient_df, by = c("GSM" = "GSM"))

sc_ObjCombined@meta.data <- meta

## PCA
sc_ObjCombined <- ScaleData(sc_ObjCombined, verbose = FALSE)
sc_ObjCombined <- RunPCA(sc_ObjCombined, npcs = 30, verbose = FALSE)

## Harmony 按 SampleType 去批次
msg("[SeuratProject] Running Harmony ...")
sc_ObjCombined <- RunHarmony(sc_ObjCombined, "SampleType")

## t-SNE + 邻居 + 聚类
msg("[SeuratProject] Running t-SNE & clustering ...")
sc_ObjCombined <- RunTSNE(sc_ObjCombined, reduction = "harmony", dims = 1:20, seed.use = 2022)
sc_ObjCombined <- FindNeighbors(sc_ObjCombined, reduction = "harmony", dims = 1:20)
sc_ObjCombined <- FindClusters(sc_ObjCombined, resolution = 0.5)

## -------- 7. 基础 QC 图：检查批次效应 -------- ##
fig_pca <- DimPlot(
  sc_ObjCombined,
  reduction = "pca",
  group.by  = "SampleType"
)
save_fig(fig_pca, "batch_effect_pca.pdf", width = 6, height = 5)

fig_harmony <- DimPlot(
  sc_ObjCombined,
  reduction = "harmony",
  group.by  = "SampleType"
)
save_fig(fig_harmony, "batch_effect_harmony.pdf", width = 6, height = 5)

## -------- 8. 保存整合后的 Seurat 对象 -------- ##
save_rds(sc_ObjCombined, "sc_ObjCombined_allcells.rds")
msg("[SeuratProject] Done. Saved sc_ObjCombined_allcells.rds.")
