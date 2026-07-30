suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(stringr); library(Matrix)
  library(SingleCellExperiment); library(slingshot); library(irlba)
  library(ggplot2); library(svglite); library(matrixStats)
})
set.seed(1234)
# ================== Colors / cluster levels ==================
cluster_levels <- as.character(0:11)
cluster_colors <- c(
  "0"="#1f77b4","1"="#ff7f0e","2"="#aec7e8","3"="#ffbb78",
  "4"="#2ca02c","5"="#ff9896","6"="#98df8a","7"="#d62728",
  "8"="#8c564b","9"="#c49c94","10"="#e377c2","11"="#bcbd22"
)
# ================== 1) Read and merge (according to your file format) ==================
data_dir <- "C:/Users/巨噬细胞重新聚类/upp5/表达矩阵/gene"  # ← change to your path
files <- list.files(data_dir, pattern="\\.xlsx$", full.names=TRUE)
stopifnot("No xlsx files found" = length(files) > 0)
read_expr_xlsx <- function(path){
  df <- as.data.frame(readxl::read_excel(path), check.names = FALSE)
  stopifnot(ncol(df) >= 5)
  gene_id <- as.character(df[[1]]); gene_symbol <- as.character(df[[2]])
  gene_symbol[is.na(gene_symbol) | gene_symbol==""] <- gene_id
  genes <- make.unique(ifelse(is.na(gene_symbol)|gene_symbol=="","NA_gene",gene_symbol))
  expr_part <- df[,5:ncol(df), drop=FALSE]
  for(j in seq_len(ncol(expr_part))) expr_part[[j]] <- suppressWarnings(as.numeric(expr_part[[j]]))
  m <- as.matrix(expr_part); rownames(m) <- genes; colnames(m) <- make.unique(colnames(m))
  m[is.na(m)] <- 0; storage.mode(m) <- "double"; m
}
expr_list <- list(); meta_list <- list()
for(f in files){
  sample_id <- str_trim(str_remove(basename(f), "\\.xlsx$"))     # e.g., Ctrl 1 cluster 0
  clx <- stringr::str_match(sample_id, "cluster\\s*(\\d+)")[,2]
  group_id <- str_trim(str_replace(sample_id, "cluster\\s*\\d+", ""))
  mat <- read_expr_xlsx(f); if(ncol(mat)==0) next
  new_cols <- paste0(colnames(mat), "_", sample_id); colnames(mat) <- new_cols
  expr_list[[f]] <- mat
  meta <- data.frame(cell_id=new_cols, sample_file=sample_id, group=group_id,
                     cluster_original=ifelse(is.na(clx), NA, clx), stringsAsFactors=FALSE)
  rownames(meta) <- meta$cell_id; meta_list[[f]] <- meta
}
stopifnot(length(expr_list) > 0)
# Union alignment + merge
all_genes <- unique(unlist(lapply(expr_list, rownames)))
expr_list_reindexed <- lapply(expr_list, function(m){
  miss <- setdiff(all_genes, rownames(m))
  if(length(miss)>0) m <- rbind(m, matrix(0, nrow=length(miss), ncol=ncol(m),
                                          dimnames=list(miss, colnames(m))))
  m <- m[all_genes,,drop=FALSE]; m[is.na(m)] <- 0; storage.mode(m) <- "double"; m
})
expr_all <- do.call(cbind, expr_list_reindexed)
expr_all <- Matrix(as.matrix(expr_all), sparse=TRUE)
cell_meta <- do.call(rbind, meta_list)
stopifnot(ncol(expr_all) == nrow(cell_meta))
cell_meta$cluster_original <- factor(as.character(cell_meta$cluster_original), levels=cluster_levels)
rownames(cell_meta) <- cell_meta$cell_id
# Remove all-zero genes/cells (optional)
#keep_gene <- Matrix::rowSums(expr_all) > 0
#keep_cell <- Matrix::colSums(expr_all) > 0
#expr_all  <- expr_all[keep_gene, keep_cell, drop=FALSE]
#cell_meta <- cell_meta[keep_cell,, drop=FALSE]
# ================== 2) Read the updated M1/M2 gene list (two columns: Gene + Type) ==================
gene_list_file <- "C:/Users/monica/Desktop/M1-M2 score genelist.csv"  # ← change to actual path
gene_tbl <- read.csv(gene_list_file, stringsAsFactors = FALSE)
gene_tbl <- gene_tbl %>%
  dplyr::mutate(
    Type = toupper(trimws(Type)),   # unify to uppercase and remove spaces
    Gene = trimws(Gene)
  ) %>%
  dplyr::filter(Type %in% c("M1", "M2"))
m1_genes <- gene_tbl %>% dplyr::filter(Type == "M1") %>% dplyr::pull(Gene) %>% unique()
m2_genes <- gene_tbl %>% dplyr::filter(Type == "M2") %>% dplyr::pull(Gene) %>% unique()
message(sprintf("Read M1 genes: %d, M2 genes: %d", length(m1_genes), length(m2_genes)))
suppressPackageStartupMessages({
  library(Matrix); library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
})

# ===== A) Standardize gene names: remove leading/trailing quotes, spaces, uppercase =====
rownames(expr_all) <- gsub("^'|'$", "", rownames(expr_all))
rownames(expr_all) <- trimws(rownames(expr_all))
rownames(expr_all) <- toupper(rownames(expr_all))

m1_genes <- toupper(trimws(m1_genes))
m2_genes <- toupper(trimws(m2_genes))
m1_genes <- m1_genes[!is.na(m1_genes) & m1_genes != ""]
m2_genes <- m2_genes[!is.na(m2_genes) & m2_genes != ""]

m1_use <- intersect(m1_genes, rownames(expr_all))
m2_use <- intersect(m2_genes, rownames(expr_all))
message(sprintf("Matched: M1=%d, M2=%d", length(m1_use), length(m2_use)))
stopifnot(length(m1_use) > 0, length(m2_use) > 0)

# ===== B) Simple normalization then compute per-cell scores =====
# Use library size normalization + CPM(1e4) + log1p to mitigate sequencing depth effects
libsize <- Matrix::colSums(expr_all)
expr_cpm <- t(t(expr_all) / pmax(libsize, 1)) * 1e4
expr_log <- log1p(expr_cpm)           # log1p works directly on sparse matrix

M1_score <- Matrix::colMeans(expr_log[m1_use, , drop = FALSE])
M2_score <- Matrix::colMeans(expr_log[m2_use, , drop = FALSE])

# ===== C) Attach Ctrl / LPS / ACD labels + generate violin plots =====
# Extract "condition" label from cell_meta$group (robust to case and variations)
score_df <- data.frame(
  cell_id  = colnames(expr_all),
  M1_score = as.numeric(M1_score),
  M2_score = as.numeric(M2_score)
) %>% 
  left_join(cell_meta, by = "cell_id") %>%
  mutate(
    condition = case_when(
      str_detect(toupper(group), "\\bCTRL?\\b|CONTROL") ~ "Ctrl",
      str_detect(toupper(group), "\\bLPS\\b") ~ "LPS",
      str_detect(toupper(group), "A[-_ ]?CDS?|ACD") ~ "ACD",
      TRUE ~ "Other"
    ),
    condition = factor(condition, levels = c("Ctrl","LPS","ACD","Other"))
  )

# Optional: remove "Other" if none exist
score_df <- score_df %>% filter(condition %in% c("Ctrl","LPS","ACD"))

# Colors (for three groups)
cond_colors <- c(Ctrl="#aec7e8", LPS="#ff7f0e", ACD="#17becf")

# —— Violin plots (overall, not faceted by cluster) ——
p_m1_vln <- ggplot(score_df, aes(x = condition, y = M1_score, fill = condition)) +
  geom_violin(trim = TRUE, scale = "width") +
  geom_boxplot(width = 0.15, outlier.size = 0.3, color = "black", alpha = 0.8) +
  scale_fill_manual(values = cond_colors) +
  labs(title = "M1 score by condition", x = NULL, y = "M1 score (log1p CPM mean)") +
  theme_minimal(base_size = 12) + theme(legend.position = "none")

p_m2_vln <- ggplot(score_df, aes(x = condition, y = M2_score, fill = condition)) +
  geom_violin(trim = TRUE, scale = "width") +
  geom_boxplot(width = 0.15, outlier.size = 0.3, color = "black", alpha = 0.8) +
  scale_fill_manual(values = cond_colors) +
  labs(title = "M2 score by condition", x = NULL, y = "M2 score (log1p CPM mean)") +
  theme_minimal(base_size = 12) + theme(legend.position = "none")

p_violin <- p_m1_vln / p_m2_vln + patchwork::plot_layout(heights = c(1,1))
print(p_violin)

# —— Optional: faceted violin plots by cluster (to visualize group differences per cluster) ——
p_m1_vln_fac <- ggplot(score_df, aes(x = condition, y = M1_score, fill = condition)) +
  geom_violin(trim = TRUE, scale = "width") +
  geom_boxplot(width = 0.15, outlier.size = 0.2, color = "black", alpha = 0.8) +
  scale_fill_manual(values = cond_colors) +
  facet_wrap(~ cluster_original, ncol = 4, drop = FALSE) +
  labs(title = "M1 score by condition (faceted by cluster)", x = NULL, y = "M1 score") +
  theme_minimal(base_size = 11) + theme(legend.position = "none")

p_m2_vln_fac <- ggplot(score_df, aes(x = condition, y = M2_score, fill = condition)) +
  geom_violin(trim = TRUE, scale = "width") +
  geom_boxplot(width = 0.15, outlier.size = 0.2, color = "black", alpha = 0.8) +
  scale_fill_manual(values = cond_colors) +
  facet_wrap(~ cluster_original, ncol = 4, drop = FALSE) +
  labs(title = "M2 score by condition (faceted by cluster)", x = NULL, y = "M2 score") +
  theme_minimal(base_size = 11) + theme(legend.position = "none")

# ===== D) Export: PDF/SVG + summary tables =====
out_dir <- file.path(data_dir, "M1_M2_scores"); dir.create(out_dir, showWarnings = FALSE)

# Overall violin plots (three groups)
ggsave(file.path(out_dir, "M1_M2_scores_violin_by_condition.pdf"), p_violin,
       device = cairo_pdf, width = 180, height = 140, units = "mm")
svglite::svglite(file.path(out_dir, "M1_M2_scores_violin_by_condition.svg"),
                 width = 180/25.4, height = 140/25.4); print(p_violin); dev.off()

# Faceted (optional)
ggsave(file.path(out_dir, "M1_score_violin_by_condition_facet_cluster.pdf"), p_m1_vln_fac,
       device = cairo_pdf, width = 210, height = 220, units = "mm")
ggsave(file.path(out_dir, "M2_score_violin_by_condition_facet_cluster.pdf"), p_m2_vln_fac,
       device = cairo_pdf, width = 210, height = 220, units = "mm")

# Export per‑cell scores (with condition)
write.csv(score_df[, c("cell_id","condition","cluster_original","M1_score","M2_score")],
          file.path(out_dir, "M1_M2_scores_per_cell_with_condition.csv"), row.names = FALSE)

# Summary by condition (median/IQR)
cond_summary <- score_df %>%
  group_by(condition) %>%
  summarise(
    n_cells   = n(),
    M1_median = median(M1_score, na.rm = TRUE),
    M1_IQR    = IQR(M1_score, na.rm = TRUE),
    M2_median = median(M2_score, na.rm = TRUE),
    M2_IQR    = IQR(M2_score, na.rm = TRUE),
    .groups = "drop"
  )
write.csv(cond_summary, file.path(out_dir, "M1_M2_scores_by_condition.csv"), row.names = FALSE)

# Detailed summary by condition and cluster (useful for supplementary tables)
cond_cluster_summary <- score_df %>%
  group_by(condition, cluster_original) %>%
  summarise(
    n_cells   = n(),
    M1_median = median(M1_score, na.rm = TRUE),
    M1_IQR    = IQR(M1_score, na.rm = TRUE),
    M2_median = median(M2_score, na.rm = TRUE),
    M2_IQR    = IQR(M2_score, na.rm = TRUE),
    .groups = "drop"
  )
write.csv(cond_cluster_summary,
          file.path(out_dir, "M1_M2_scores_by_condition_and_cluster.csv"),
          row.names = FALSE)
suppressPackageStartupMessages({
  library(ggpubr)   # ← added for significance testing and annotation
  library(rstatix)  # ← added for statistical result tables
})

# ===== E) Statistical tests among three groups (Kruskal + pairwise comparisons) =====
## 1) Kruskal–Wallis test for overall differences
kw_M1 <- kruskal_test(score_df, M1_score ~ condition)
kw_M2 <- kruskal_test(score_df, M2_score ~ condition)
kw_M1; kw_M2   # output to console

## 2) If overall difference exists, perform pairwise comparisons (Dunn's test with Bonferroni correction)
pairwise_M1 <- dunn_test(score_df, M1_score ~ condition, p.adjust.method = "bonferroni")
pairwise_M2 <- dunn_test(score_df, M2_score ~ condition, p.adjust.method = "bonferroni")

## 3) Save result tables
write.csv(kw_M1,  file.path(out_dir, "Kruskal_M1_score.csv"), row.names = FALSE)
write.csv(kw_M2,  file.path(out_dir, "Kruskal_M2_score.csv"), row.names = FALSE)
write.csv(pairwise_M1, file.path(out_dir, "Dunn_pairwise_M1_score.csv"), row.names = FALSE)
write.csv(pairwise_M2, file.path(out_dir, "Dunn_pairwise_M2_score.csv"), row.names = FALSE)
# ==== Add significance annotations to plots ====
# M1
p_m1_vln_sig <- ggviolin(score_df, x = "condition", y = "M1_score",
                         fill = "condition", palette = cond_colors,
                         add = "boxplot", width = 0.9, trim = TRUE) +
  stat_compare_means(method = "kruskal.test", label.y = max(score_df$M1_score)*1.05) +  # overall p‑value
  stat_compare_means(comparisons = list(c("Ctrl","LPS"), c("LPS","ACD"), c("Ctrl","ACD")),
                     method = "wilcox.test", label = "p.signif", hide.ns = TRUE) +
  labs(title = "M1 score by condition", y = "M1 score (log1p CPM mean)") +
  theme_minimal(base_size = 12) + theme(legend.position = "none")

# M2
p_m2_vln_sig <- ggviolin(score_df, x = "condition", y = "M2_score",
                         fill = "condition", palette = cond_colors,
                         add = "boxplot", width = 0.9, trim = TRUE) +
  stat_compare_means(method = "kruskal.test", label.y = max(score_df$M2_score)*1.05) +
  stat_compare_means(comparisons = list(c("Ctrl","LPS"), c("LPS","ACD"), c("Ctrl","ACD")),
                     method = "wilcox.test", label = "p.signif", hide.ns = TRUE) +
  labs(title = "M2 score by condition", y = "M2 score (log1p CPM mean)") +
  theme_minimal(base_size = 12) + theme(legend.position = "none")

p_violin_sig <- p_m1_vln_sig / p_m2_vln_sig + patchwork::plot_layout(heights = c(1,1))
print(p_violin_sig)

# Save violin plots with significance
ggsave(file.path(out_dir, "M1_M2_scores_violin_by_condition_with_significance.pdf"),
       p_violin_sig, device = cairo_pdf, width = 180, height = 150, units = "mm")
svglite::svglite(file.path(out_dir, "M1_M2_scores_violin_by_condition_with_significance.svg"),
                 width = 180/25.4, height = 150/25.4); print(p_violin_sig); dev.off()
