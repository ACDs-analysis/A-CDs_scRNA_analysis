suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(stringr); library(Matrix)
  library(SingleCellExperiment); library(slingshot); library(irlba)
  library(ggplot2); library(svglite); library(matrixStats)
})
set.seed(1234)
# ================== 颜色/cluster级别 ==================
cluster_levels <- as.character(0:11)
cluster_colors <- c(
  "0"="#1f77b4","1"="#ff7f0e","2"="#aec7e8","3"="#ffbb78",
  "4"="#2ca02c","5"="#ff9896","6"="#98df8a","7"="#d62728",
  "8"="#8c564b","9"="#c49c94","10"="#e377c2","11"="#bcbd22"
)
# ================== 1) 读取并合并（按你文件格式） ==================
data_dir <- "E:/科研/单细胞测序/gene"   # ←改成你的路径
files <- list.files(data_dir, pattern="\\.xlsx$", full.names=TRUE)
stopifnot("未找到xlsx文件" = length(files) > 0)
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
  sample_id <- str_trim(str_remove(basename(f), "\\.xlsx$"))     # 例：Ctrl 1 cluster 0
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
# 并集对齐 + 合并
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
# 去掉全零基因/细胞
#keep_gene <- Matrix::rowSums(expr_all) > 0
#keep_cell <- Matrix::colSums(expr_all) > 0
#expr_all  <- expr_all[keep_gene, keep_cell, drop=FALSE]
#cell_meta <- cell_meta[keep_cell,, drop=FALSE]
# ================== 2) 选ECM关键基因 & 预处理 ==================
# 统一基因名为大写，便于稳健匹配
rownames(expr_all) <- gsub("^'|'$", "", rownames(expr_all))
rownames(expr_all) <- trimws(rownames(expr_all))
rownames(expr_all) <- toupper(rownames(expr_all))

ecm_core_genes <- c(
  "IRF1","IRF5","IRF7","KLF2","KLF4","STAT1","NFKB1","RELA","NR1H3","PPARG","ARG1"
)
ecm_core_genes <- unique(toupper(ecm_core_genes))

# 从 sample_file 中提取条件为 Ctrl/LPS/ACD（保留原样本标签）
# 你的 sample_file 类似 "Ctrl 1 cluster 0"、"ACD 3 cluster 4"
cell_meta$condition <- dplyr::case_when(
  grepl("^\\s*CTRL\\b",  toupper(cell_meta$group)) ~ "Ctrl",
  grepl("^\\s*LPS\\b",   toupper(cell_meta$group)) ~ "LPS",
  grepl("^\\s*A\\s*CDS?\\b|^\\s*ACD\\b", toupper(cell_meta$group)) ~ "ACD",
  TRUE ~ NA_character_
)
cell_meta <- cell_meta[!is.na(cell_meta$condition), , drop = FALSE]
cell_meta$condition <- factor(cell_meta$condition, levels = c("Ctrl","LPS","ACD"))

# 同步表达矩阵到保留的细胞
expr_all <- expr_all[, rownames(cell_meta), drop = FALSE]

# 计算 log1p-CPM（每细胞归一化到100万）
calc_log1p_cpm <- function(sparse_counts) {
  libsize <- Matrix::colSums(sparse_counts)
  libsize[libsize == 0] <- 1  # 防止除零
  # CPM = counts / libsize * 1e6
  cpm <- t(t(sparse_counts) / libsize) * 1e6
  log1p(cpm)
}
expr_log1p_cpm <- calc_log1p_cpm(expr_all)

# 检查丢失基因
present <- ecm_core_genes[ecm_core_genes %in% rownames(expr_log1p_cpm)]
missing <- setdiff(ecm_core_genes, present)
if (length(missing) > 0) {
  message("下列ECM基因在矩阵中未找到（请核对别名/物种符号）：\n",
          paste(missing, collapse = ", "))
}
stopifnot(length(present) > 0)

# 抽取这些基因表达，整理成长表
suppressPackageStartupMessages({ library(tidyr) })
expr_ecm <- as.matrix(expr_log1p_cpm[present, , drop = FALSE])
df_long <- as.data.frame(t(expr_ecm), check.names = FALSE)
df_long$cell_id <- rownames(df_long)
df_long <- tidyr::pivot_longer(df_long, cols = all_of(present),
                               names_to = "gene", values_to = "log1pCPM")
plot_df <- dplyr::left_join(df_long, 
                            cell_meta[, c("cell_id","condition")], 
                            by = "cell_id") %>%
  dplyr::filter(!is.na(condition))

# ========= 小提琴图（分面：gene） =========
cond_colors <- c("Ctrl"="#aec7e8","LPS"="#ff7f0e","ACD"="#17becf")

p_ecm10 <- ggplot(plot_df, aes(x = condition, y = log1pCPM, fill = condition)) +
  geom_violin(width = 0.9, trim = TRUE, scale = "width",
              color = "black", linewidth = 0.2) +
  geom_boxplot(width = 0.15, outlier.size = 0.4, outlier.stroke = 0.2,
               fill = "white", color = "black", linewidth = 0.2) +
  facet_wrap(~ gene, scales = "free_y", ncol = 5) +
  scale_fill_manual(values = cond_colors, drop = FALSE) +
  labs(x = NULL, y = "Expression (log1p CPM)",
       title = "ECM (heatmap-selected) genes in macrophages") +
  theme_minimal(base_size = 10) +
  
theme(
  legend.position = "none",
  strip.text = element_text(face = "bold"),
  panel.grid.major = element_blank(),  # 添加这行：移除主要网格线
  panel.grid.minor = element_blank(),   # 移除次要网格线
  panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4)
)
# ========= 统计学标注：Kruskal + 两两Wilcoxon (BH) =========
suppressPackageStartupMessages({ library(ggpubr); library(rstatix) })

valid_genes10 <- plot_df %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise(n_groups = dplyr::n_distinct(condition[!is.na(log1pCPM)]), .groups="drop") %>%
  dplyr::filter(n_groups >= 2) %>% dplyr::pull(gene)

plot_df_stat <- plot_df %>% dplyr::filter(gene %in% valid_genes10)

# 整体（可选）
kruskal_ecm <- plot_df_stat %>%
  group_by(gene) %>%
  rstatix::kruskal_test(log1pCPM ~ condition) %>%
  ungroup()


pairwise_comparisons <- list(c("Ctrl","LPS"), c("LPS","ACD"), c("Ctrl","ACD"))

pw_ecm <- plot_df_stat %>%
  dplyr::group_by(gene) %>%
  rstatix::pairwise_wilcox_test(
    log1pCPM ~ condition,
    comparisons = pairwise_comparisons,
    p.adjust.method = "BH",
    paired = FALSE
  ) %>%
  dplyr::ungroup()

# 给每个基因分配不重叠的 y.position
y_top_tbl <- plot_df_stat %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise(y_top = max(log1pCPM, na.rm = TRUE), .groups = "drop")

offset <- 0.08
pw_ecm <- pw_ecm %>%
  dplyr::left_join(y_top_tbl, by = "gene") %>%
  dplyr::group_by(gene) %>%
  dplyr::mutate(rank_in_gene = dplyr::row_number(),
                y.position = y_top + rank_in_gene * (offset * pmax(y_top, 1) + 0.05)) %>%
  dplyr::ungroup()

# 仅保留真实存在的比较（避免某些基因缺组别）
present_levels_tbl <- plot_df_stat %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise(levels_present = list(unique(as.character(condition))), .groups = "drop")

pw_ecm <- pw_ecm %>%
  dplyr::left_join(present_levels_tbl, by = "gene") %>%
  dplyr::rowwise() %>%
  dplyr::filter(all(c(group1, group2) %in% levels_present)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(p.label = dplyr::if_else(!is.na(p.adj.signif), p.adj.signif, rstatix::p_format(p.adj)))

# 在 p_ecm10 上叠加显著性
p_ecm10_sig <- p_ecm10 +
  ggpubr::stat_pvalue_manual(
    pw_ecm,
    label = "p.label",
    y.position = "y.position",
    xmin = "group1", xmax = "group2",
    tip.length = 0.01,
    step.increase = 0,
    hide.ns = TRUE
  )

print(p_ecm10_sig)


# ========= 保存 =========
dir.create("figs_ecm", showWarnings = FALSE)
ggsave(
  filename = file.path("figs_ecm", "TF.pdf"),
  plot = p_ecm10_sig,
  width = 240, height = 160, units = "mm",
  device = grDevices::cairo_pdf # 或者省略 device, 让 ggsave 按扩展名选 pdf()
)
