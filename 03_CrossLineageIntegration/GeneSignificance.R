############################################################
## 9.GeneSignificance.R
## 基因在不同 group 间的显著性分析（无图）
############################################################

source(file.path("analysis", "0.BASIC.R"))

library(Seurat)
library(dplyr)

# group_var: 比如 "SampleType", "BigGroup", "seurat_clusters", "group"
# genes: 一个向量，指定关心的基因

run_gene_significance <- function(
  sc_obj_file = "sc_Obj_Allcells_withBigGroup.rds",
  genes,
  group_var = "SampleType"
) {
  msg("[GeneSig] File: ", sc_obj_file, "  group_var: ", group_var)

  so <- load_rds(sc_obj_file)

  if (!group_var %in% colnames(so@meta.data)) {
    stop("[GeneSig] group_var not found in meta.data: ", group_var)
  }

  so <- SetIdent(so, value = group_var)

  # 平均表达矩阵
  avg_exp <- AverageExpression(so, features = genes, group.by = group_var)$RNA
  avg_exp <- as.data.frame(avg_exp)
  avg_exp$Gene <- rownames(avg_exp)

  save_rds(avg_exp, paste0("GeneSig_avgExp_", group_var, ".rds"))
  write.csv(
    avg_exp,
    file.path(project_dir, "results", "tables",
              paste0("GeneSig_avgExp_", group_var, ".csv")),
    row.names = FALSE
  )

  # 显著性检验：多组用 Kruskal-Wallis，两组用 Wilcoxon
  group_levels <- levels(as.factor(so@meta.data[[group_var]]))
  n_group <- length(group_levels)

  gene_sig <- lapply(genes, function(g) {
    expr_vec <- FetchData(so, vars = g)[, 1]
    grp_vec  <- as.factor(so@meta.data[[group_var]])

    if (n_group == 2) {
      test_res <- wilcox.test(expr_vec ~ grp_vec)
      data.frame(
        Gene = g,
        Test = "Wilcoxon",
        p_value = test_res$p.value,
        stringsAsFactors = FALSE
      )
    } else {
      test_res <- kruskal.test(expr_vec ~ grp_vec)
      data.frame(
        Gene = g,
        Test = "Kruskal",
        p_value = test_res$p.value,
        stringsAsFactors = FALSE
      )
    }
  }) %>%
    bind_rows()

  gene_sig$FDR <- p.adjust(gene_sig$p_value, method = "BH")

  save_rds(gene_sig, paste0("GeneSig_stats_", group_var, ".rds"))
  write.csv(
    gene_sig,
    file.path(project_dir, "results", "tables",
              paste0("GeneSig_stats_", group_var, ".csv")),
    row.names = FALSE
  )

  msg("[GeneSig] Done.")
  invisible(list(avgExp = avg_exp, stats = gene_sig))
}
