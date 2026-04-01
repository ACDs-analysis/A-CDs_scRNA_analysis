suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(stringr); library(Matrix)
  library(SingleCellExperiment); library(slingshot); library(irlba)
  library(ggplot2); library(svglite); library(matrixStats)
  library(tidyr); library(ggpubr); library(rstatix)
})
set.seed(1234)

# ================== 颜色/cluster级别 ==================
cluster_levels <- as.character(0:11)
cluster_colors <- c(
  "0"="#1f77b4","1"="#ff7f0e","2"="#aec7e8","3"="#ffbb78",
  "4"="#2ca02c","5"="#ff9896","6"="#98df8a","7"="#d62728",
  "8"="#8c564b","9"="#c49c94","10"="#e377c2","11"="#bcbd22"
)

# ================== 1) 内皮细胞：读取并合并（CSV，每个文件一个样本/cluster） ==================
data_dir <- "E:/科研/单细胞测序/EC gene"   # ← 内皮细胞路径
files <- list.files(data_dir, pattern="\\.csv$", full.names=TRUE)
stopifnot("未找到csv文件" = length(files) > 0)

read_expr_csv <- function(path){
  df <- tryCatch(
    as.data.frame(utils::read.csv(path, check.names = FALSE), stringsAsFactors = FALSE),
    error = function(e) as.data.frame(utils::read.csv(path, check.names = FALSE, fileEncoding="UTF-8-BOM"), stringsAsFactors = FALSE)
  )
  stopifnot(ncol(df) >= 2)
  
  # 假设第1列是gene_id，第2列是gene_symbol（若为空用gene_id代替）
  gene_id <- as.character(df[[1]])
  gene_symbol <- as.character(df[[2]])
  gene_symbol[is.na(gene_symbol) | gene_symbol==""] <- gene_id
  genes <- make.unique(ifelse(is.na(gene_symbol)|gene_symbol=="","NA_gene",gene_symbol))
  
  # 尽量兼容不同导出格式：优先从第5列开始，否则第3列，否则第2列
  expr_start <- if (ncol(df) >= 5) 5 else if (ncol(df) >= 3) 3 else 2
  expr_part <- df[, expr_start:ncol(df), drop = FALSE]
  
  # 强制数值化（非数值→NA→0）
  for(j in seq_len(ncol(expr_part))){
    expr_part[[j]] <- suppressWarnings(as.numeric(expr_part[[j]]))
  }
  m <- as.matrix(expr_part)
  m[is.na(m)] <- 0
  rownames(m) <- genes
  colnames(m) <- make.unique(colnames(m))
  storage.mode(m) <- "double"
  m
}

expr_list <- list(); meta_list <- list()
for(f in files){
  sample_id <- str_trim(str_remove(basename(f), "\\.csv$"))   # 例：Ctrl 1 Cluster 6
  # 兼容大小写 "Cluster 6" / "cluster 6"
  clx <- stringr::str_match(sample_id, "(?i)cluster\\s*(\\d+)")[,2]
  # 从sample_id去掉"cluster x"得到 group（Ctrl 1 / LPS 2 / ACD 3 …）
  group_id <- str_trim(stringr::str_replace(sample_id, "(?i)cluster\\s*\\d+", ""))
  
  mat <- read_expr_csv(f); if(ncol(mat)==0) next
  new_cols <- paste0(colnames(mat), "_", sample_id)
  colnames(mat) <- new_cols
  expr_list[[f]] <- mat
  
  meta <- data.frame(
    cell_id = new_cols,
    sample_file = sample_id,
    group = group_id,
    cluster_original = ifelse(is.na(clx), NA, clx),
    stringsAsFactors = FALSE
  )
  rownames(meta) <- meta$cell_id
  meta_list[[f]] <- meta
}
stopifnot(length(expr_list) > 0)

# —— 并集对齐 + 合并
all_genes <- unique(unlist(lapply(expr_list, rownames)))
expr_list_reindexed <- lapply(expr_list, function(m){
  miss <- setdiff(all_genes, rownames(m))
  if(length(miss)>0){
    m <- rbind(m, matrix(0, nrow=length(miss), ncol=ncol(m),
                         dimnames=list(miss, colnames(m))))
  }
  m <- m[all_genes,,drop=FALSE]
  m[is.na(m)] <- 0
  storage.mode(m) <- "double"
  m
})
expr_all <- do.call(cbind, expr_list_reindexed)
expr_all <- Matrix(as.matrix(expr_all), sparse=TRUE)

cell_meta <- do.call(rbind, meta_list)
stopifnot(ncol(expr_all) == nrow(cell_meta))
cell_meta$cluster_original <- factor(as.character(cell_meta$cluster_original), levels=cluster_levels)
rownames(cell_meta) <- cell_meta$cell_id

# ================== 2) 选ECM关键基因 & 预处理（内皮） ==================
# 统一基因名为大写，便于稳健匹配
rownames(expr_all) <- gsub("^'|'$", "", rownames(expr_all))
rownames(expr_all) <- trimws(rownames(expr_all))
rownames(expr_all) <- toupper(rownames(expr_all))

ecm_core_genes <- c(
  "ACKR2","TNFRSF1A","ITGA9","ITGA2","ITGA1","ITGB1",
  "KDR","FLT1","ICAM2",
  "PECAM1","OSMR","IL6ST"
)
ecm_core_genes <- unique(toupper(ecm_core_genes))   # 统一转大写

# 从 group 中提取条件为 Ctrl/LPS/ACD（保留原样本标签）
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
  libsize[libsize == 0] <- 1
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

# 抽取表达，整理成长表
expr_ecm <- as.matrix(expr_log1p_cpm[present, , drop = FALSE])
df_long <- as.data.frame(t(expr_ecm), check.names = FALSE)
df_long$cell_id <- rownames(df_long)
df_long <- tidyr::pivot_longer(df_long, cols = all_of(present),
                               names_to = "gene", values_to = "log1pCPM")
plot_df <- dplyr::left_join(df_long, 
                            cell_meta[, c("cell_id","condition")], 
                            by = "cell_id") %>%
  dplyr::filter(!is.na(condition))

# ========= 小提琴图（分面：gene，5列） =========
cond_colors <- c("Ctrl"="#aec7e8","LPS"="#ff7f0e","ACD"="#17becf")

p_ecm10 <- ggplot(plot_df, aes(x = condition, y = log1pCPM, fill = condition)) +
  geom_violin(width = 0.9, trim = TRUE, scale = "width",
              color = "black", linewidth = 0.2) +
  geom_boxplot(width = 0.15, outlier.size = 0.4, outlier.stroke = 0.2,
               fill = "white", color = "black", linewidth = 0.2) +
  facet_wrap(~ gene, scales = "free_y", ncol = 5) +
  scale_fill_manual(values = cond_colors, drop = FALSE) +
  labs(x = NULL, y = "Expression (log1p CPM)",
       title = "Endothelial cells: ECM-related genes") +
  theme_minimal(base_size = 10) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold"),
    panel.grid.major = element_blank(),  # 添加这行：移除主要网格线
    panel.grid.minor = element_blank(),   # 移除次要网格线
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4)
  )

# ========= 统计学标注：Kruskal + 两两Wilcoxon (BH) =========
valid_genes10 <- plot_df %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise(n_groups = dplyr::n_distinct(condition[!is.na(log1pCPM)]), .groups="drop") %>%
  dplyr::filter(n_groups >= 2) %>% dplyr::pull(gene)

plot_df_stat <- plot_df %>% dplyr::filter(gene %in% valid_genes10)

# 整体（可选）
kruskal_ecm <- plot_df_stat %>%
  dplyr::group_by(gene) %>%
  rstatix::kruskal_test(log1pCPM ~ condition) %>%
  dplyr::ungroup()

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

# 叠加显著性
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

# ========= 保存（5列布局；可以按需要调整整图尺寸）=========
dir.create("figs_ecm", showWarnings = FALSE)
ggsave(
  filename = file.path("figs_ecm", "EC_Endothelial_ECM8_violin_withSig.pdf"),
  plot = p_ecm10_sig,
  width = 240, height = 160, units = "mm",
  device = grDevices::cairo_pdf
)
