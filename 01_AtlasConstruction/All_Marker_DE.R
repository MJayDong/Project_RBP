############################################################
## 1.1.All_Marker_DE.R
## All cells：Marker / DEG 计算（无可视化）
## Author: Mingjie
############################################################

## 基础环境
source(file.path("analysis", "0.BASIC.R"))

library(Seurat)
library(dplyr)
library(stringr)

msg("[All_DE] Start DEG and marker analysis ...")

## ---------------- 1. 读取对象 ---------------------------
sc_all  <- load_rds("sc_Obj_Allcells_withBigGroup.rds")

sc_epi  <- load_rds("sc_Obj_EPITHELIAL.rds")
sc_imm  <- load_rds("sc_Obj_IMMUNE.rds")
sc_mes  <- load_rds("sc_Obj_MESENCHYMAL.rds")
sc_vas  <- load_rds("sc_Obj_VASCULAR.rds")

msg("[All_DE] Objects loaded.")

## ---------------- 2. All cluster markers -----------------
msg("[All_DE] Finding all-cluster markers ...")

Idents(sc_all) <- sc_all$seurat_clusters

markers_allclusters <- FindAllMarkers(
  sc_all,
  only.pos     = TRUE,
  min.pct      = 0.25,
  logfc.threshold = 0.25
)

save_rds(markers_allclusters, "DE_markers_allclusters.rds")
write.csv(markers_allclusters,
          file.path(project_dir,"results","tables","DE_markers_allclusters.csv"),
          row.names = FALSE)

msg("[All_DE] Cluster markers saved: DE_markers_allclusters.*")


## ---------------- 3. 三组比较：Normal / Primary / Metastasis -----------------

msg("[All_DE] Normal vs Primary vs Metastasis ...")

Idents(sc_all) <- sc_all$SampleType

DE_N_vs_P <- FindMarkers(
  sc_all,
  ident.1 = "Normal",
  ident.2 = "Primary",
  logfc.threshold = 0.25,
  min.pct = 0.2
)

DE_P_vs_M <- FindMarkers(
  sc_all,
  ident.1 = "Primary",
  ident.2 = "Metastasis",
  logfc.threshold = 0.25,
  min.pct = 0.2
)

DE_N_vs_M <- FindMarkers(
  sc_all,
  ident.1 = "Normal",
  ident.2 = "Metastasis",
  logfc.threshold = 0.25,
  min.pct = 0.2
)

## 保存
save_rds(DE_N_vs_P, "DE_Normal_vs_Primary.rds")
save_rds(DE_P_vs_M, "DE_Primary_vs_Metastasis.rds")
save_rds(DE_N_vs_M, "DE_Normal_vs_Metastasis.rds")

write.csv(DE_N_vs_P, file.path(project_dir,"results/tables","DE_Normal_vs_Primary.csv"))
write.csv(DE_P_vs_M, file.path(project_dir,"results/tables","DE_Primary_vs_Metastasis.csv"))
write.csv(DE_N_vs_M, file.path(project_dir,"results/tables","DE_Normal_vs_Metastasis.csv"))

msg("[All_DE] DEG saved for N-P, P-M, N-M.")


## ---------------- 4. 各大类（EPITHELIAL / IMMUNE / MESENCHYMAL / VASCULAR）内部 marker -----------------

## 一个辅助函数
run_DE_for_object <- function(obj, group_label, prefix) {
  Idents(obj) <- obj$seurat_clusters
  
  msg("[All_DE] Running FindAllMarkers for ", group_label)
  
  mk <- FindAllMarkers(
    obj,
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = 0.25
  )
  
  outfile_rds <- paste0("DE_markers_", prefix, ".rds")
  outfile_csv <- paste0("DE_markers_", prefix, ".csv")
  
  save_rds(mk, outfile_rds)
  write.csv(mk,
            file.path(project_dir,"results","tables",outfile_csv),
            row.names = FALSE)
  
  msg("[All_DE] Saved: ", outfile_rds)
  return(mk)
}

mk_epi <- run_DE_for_object(sc_epi, "EPITHELIAL", "EPITHELIAL")
mk_imm <- run_DE_for_object(sc_imm, "IMMUNE", "IMMUNE")
mk_mes <- run_DE_for_object(sc_mes, "MESENCHYMAL", "MESENCHYMAL")
mk_vas <- run_DE_for_object(sc_vas, "VASCULAR", "VASCULAR")

msg("[All_DE] Subgroup marker detection complete.")


## ---------------- 5. 导出 Top5 marker（为可视化准备） --------------
msg("[All_DE] Export Top5 per cluster ...")

extract_top5 <- function(marker_tbl) {
  marker_tbl %>%
    group_by(cluster) %>%
    top_n(5, avg_log2FC) %>%
    arrange(cluster, desc(avg_log2FC))
}

top5_all <- extract_top5(markers_allclusters)

save_rds(top5_all, "DE_top5_allclusters.rds")
write.csv(top5_all,
          file.path(project_dir,"results","tables","DE_top5_allclusters.csv"),
          row.names = FALSE)

msg("[All_DE] Top5 markers saved.")

## ---------------- END ----------------
msg("[All_DE] All DEG and marker analysis done.")
