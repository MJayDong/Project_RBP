############################################################
## 1.0.All_Preprocessing.R
## All cells：元数据规范 + BigGroup 构建 + 组成统计
## Author: Mingjie
############################################################

## 依赖基础环境
source(file.path("analysis", "0.BASIC.R"))

library(Seurat)
library(dplyr)
library(stringr)

msg("[All_Preprocessing] Load integrated Seurat object ...")

## 1. 读取整合后的对象 -------------------------------------------
sc_ObjCombined <- load_rds("sc_ObjCombined_allcells.rds")

## 2. 规范 SampleType 因子顺序 -----------------------------------
## 如果之前有拼写错误（如 Matastasis），这里一并修正
if ("SampleType" %in% colnames(sc_ObjCombined@meta.data)) {
  sc_ObjCombined$SampleType <- gsub(
    pattern = "Matastasis",
    replacement = "Metastasis",
    x = sc_ObjCombined$SampleType
  )
  sc_ObjCombined$SampleType <- factor(
    sc_ObjCombined$SampleType,
    levels = c("Normal", "Primary", "Metastasis")
  )
} else {
  stop("[All_Preprocessing] meta.data 中缺少 SampleType 列。")
}

msg("[All_Preprocessing] SampleType levels: ",
    paste(levels(sc_ObjCombined$SampleType), collapse = ", "))

## 3. 从 CellName 构建 BigGroup -----------------------------------
## 假设 CellName 已经是你在后续步骤用来表示细胞类型的那一列
if (!"CellName" %in% colnames(sc_ObjCombined@meta.data)) {
  stop("[All_Preprocessing] meta.data 中缺少 CellName，无法构建 BigGroup。")
}

big_group <- sc_ObjCombined$CellName

## EPITHELIAL 相关
big_group <- gsub("Breast cancer cells",   "EPITHELIAL", big_group)
big_group <- gsub("Breast luminal cells",  "EPITHELIAL", big_group)
big_group <- gsub("Myoepithelial cells",   "EPITHELIAL", big_group)
big_group <- gsub("Epithelial cells",      "EPITHELIAL", big_group)
big_group <- gsub("Epidermal cells",       "EPITHELIAL", big_group)

## MESENCHYMAL 相关
big_group <- gsub("Fibroblasts",           "MESENCHYMAL", big_group)
big_group <- gsub("Pericytes",             "MESENCHYMAL", big_group)

## VASCULAR-RELATED 相关
big_group <- gsub("Endothelial cells",     "VASCULAR-RELATED", big_group)

## IMMUNE 相关
big_group <- gsub("T cells",               "IMMUNE", big_group)
big_group <- gsub("B cells",               "IMMUNE", big_group)
big_group <- gsub("Plasma cells",          "IMMUNE", big_group)
big_group <- gsub("Macrophages",           "IMMUNE", big_group)
big_group <- gsub("Mast cells",            "IMMUNE", big_group)

## 挂回对象
sc_ObjCombined$BigGroup <- big_group

## 只保留我们关心的四大类，其他统称为 OTHER（如果存在）
valid_groups <- c("EPITHELIAL", "IMMUNE", "MESENCHYMAL", "VASCULAR-RELATED")
sc_ObjCombined$BigGroup[!(sc_ObjCombined$BigGroup %in% valid_groups)] <- "OTHER"

sc_ObjCombined$BigGroup <- factor(
  sc_ObjCombined$BigGroup,
  levels = c("EPITHELIAL", "IMMUNE", "MESENCHYMAL", "VASCULAR-RELATED", "OTHER")
)

msg("[All_Preprocessing] BigGroup table:")
print(table(sc_ObjCombined$BigGroup))

## 4. 细胞组成统计表（SampleType × BigGroup） ---------------------

## 4.1 计数表
sample_big_table_count <- as.data.frame(
  table(sc_ObjCombined$SampleType, sc_ObjCombined$BigGroup),
  stringsAsFactors = FALSE
)
colnames(sample_big_table_count) <- c("SampleType", "BigGroup", "CellNumber")

## 4.2 每个 SampleType 内百分比
sample_big_table_pct <- sample_big_table_count %>%
  group_by(SampleType) %>%
  mutate(Percentage = 100 * CellNumber / sum(CellNumber)) %>%
  ungroup()

## 保存结果，方便 1.2.All_Visualization.R 直接使用
save_rds(sample_big_table_count, "All_sample_biggroup_counts.rds")
save_rds(sample_big_table_pct,   "All_sample_biggroup_percentages.rds")

msg("[All_Preprocessing] Saved composition tables: All_sample_biggroup_*")

## 5. 按 BigGroup 切分对象并保存（供后续 RBP/DE 模块用） ------------

Idents(sc_ObjCombined) <- sc_ObjCombined$BigGroup

sc_Obj_EPITHELIAL   <- subset(sc_ObjCombined, idents = "EPITHELIAL")
sc_Obj_IMMUNE       <- subset(sc_ObjCombined, idents = "IMMUNE")
sc_Obj_MESENCHYMAL  <- subset(sc_ObjCombined, idents = "MESENCHYMAL")
sc_Obj_VASCULAR     <- subset(sc_ObjCombined, idents = "VASCULAR-RELATED")

save_rds(sc_ObjCombined,  "sc_Obj_Allcells_withBigGroup.rds")
save_rds(sc_Obj_EPITHELIAL,  "sc_Obj_EPITHELIAL.rds")
save_rds(sc_Obj_IMMUNE,      "sc_Obj_IMMUNE.rds")
save_rds(sc_Obj_MESENCHYMAL, "sc_Obj_MESENCHYMAL.rds")
save_rds(sc_Obj_VASCULAR,    "sc_Obj_VASCULAR.rds")

msg("[All_Preprocessing] Saved subset Seurat objects by BigGroup.")
msg("[All_Preprocessing] Done.")
