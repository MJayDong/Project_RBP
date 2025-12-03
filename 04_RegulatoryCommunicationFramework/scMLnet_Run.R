############################################################
## 12.scMLnet_Run.R
## 运行 scMLnet（无图）—— 输出 TF–Target–L–R 网络
############################################################

source(file.path("analysis", "0.BASIC.R"))

library(scMLnet)
library(Seurat)
library(dplyr)

run_scMLnet <- function(
  sc_obj_file,
  sender_group,   # 如 "EPITHELIAL"
  receiver_group, # 如 "IMMUNE"
  output_prefix = "scMLnet"
) {
  so <- load_rds(sc_obj_file)

  # sender & receiver 子集
  so_sender   <- subset(so, subset = BigGroup == sender_group)
  so_receiver <- subset(so, subset = BigGroup == receiver_group)

  sender_exp   <- as.matrix(so_sender@assays$RNA@counts)
  receiver_exp <- as.matrix(so_receiver@assays$RNA@counts)

  # scMLnet 主分析
  mlnet_res <- scMLnet(
    sender_exp,
    receiver_exp,
    species = "Human",
    L_R_pairs = "Human",
    TF_target = "Human",
    ncores = 8
  )

  # 输出四层网络
  L_R  <- mlnet_res$Ligand_Receptor
  R_TF <- mlnet_res$Receptor_TF
  TF_T <- mlnet_res$TF_Target

  save_rds(L_R,   paste0(output_prefix, "_LigandReceptor.rds"))
  save_rds(R_TF,  paste0(output_prefix, "_ReceptorTF.rds"))
  save_rds(TF_T,  paste0(output_prefix, "_TFtarget.rds"))

  write.csv(
    TF_T,
    file.path(project_dir, "results", "tables",
              paste0(output_prefix, "_TFtarget.csv")),
    row.names = FALSE
  )

  msg("[scMLnet] Done.")
  invisible(list(L_R = L_R, R_TF = R_TF, TF_T = TF_T))
}
