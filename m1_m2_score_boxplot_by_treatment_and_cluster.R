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
data_dir <- "C:/Users/monica/Documents/我的文件/科研/博后/课题选择/附子碳点与脓毒症/单细胞测序/华大/巨噬细胞重新聚类/upp5/表达矩阵/gene"  # ←改成你的路径
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
# ================== 2) 读取新版 M1/M2 基因列表（Gene + Type 两列） ==================
gene_list_file <- "C:/Users/monica/Desktop/M1-M2 score genelist.csv"  # ←改成实际路径
gene_tbl <- read.csv(gene_list_file, stringsAsFactors = FALSE)
gene_tbl <- gene_tbl %>%
  dplyr::mutate(
    Type = toupper(trimws(Type)),   # 统一大写、去空格
    Gene = trimws(Gene)
  ) %>%
  dplyr::filter(Type %in% c("M1", "M2"))
m1_genes <- gene_tbl %>% dplyr::filter(Type == "M1") %>% dplyr::pull(Gene) %>% unique()
m2_genes <- gene_tbl %>% dplyr::filter(Type == "M2") %>% dplyr::pull(Gene) %>% unique()
message(sprintf("读取 M1 基因数: %d, M2 基因数: %d", length(m1_genes), length(m2_genes)))
suppressPackageStartupMessages({
  library(Matrix); library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
})

# ===== A) 规范基因名：去首尾引号、去空格、统一大写 =====
rownames(expr_all) <- gsub("^'|'$", "", rownames(expr_all))
rownames(expr_all) <- trimws(rownames(expr_all))
rownames(expr_all) <- toupper(rownames(expr_all))

m1_genes <- toupper(trimws(m1_genes))
m2_genes <- toupper(trimws(m2_genes))
m1_genes <- m1_genes[!is.na(m1_genes) & m1_genes != ""]
m2_genes <- m2_genes[!is.na(m2_genes) & m2_genes != ""]

m1_use <- intersect(m1_genes, rownames(expr_all))
m2_use <- intersect(m2_genes, rownames(expr_all))
message(sprintf("匹配到: M1=%d, M2=%d", length(m1_use), length(m2_use)))
stopifnot(length(m1_use) > 0, length(m2_use) > 0)

# ===== B) 简单归一化后计算每个细胞的 score =====
# 用库大小归一化 + CPM(1e4) + log1p，避免细胞测序深度差异干扰
libsize <- Matrix::colSums(expr_all)
expr_cpm <- t(t(expr_all) / pmax(libsize, 1)) * 1e4
expr_log <- log1p(expr_cpm)           # 稀疏矩阵也可直接 log1p

M1_score <- Matrix::colMeans(expr_log[m1_use, , drop = FALSE])
M2_score <- Matrix::colMeans(expr_log[m2_use, , drop = FALSE])

# ===== C) 贴上 Ctrl / LPS / ACD 标签 + 生成小提琴图 =====
# 从 cell_meta$group 提取“实验组”标签（鲁棒匹配大小写与可能写法）
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

# 可选：如果没有“Other”，可以去掉
score_df <- score_df %>% filter(condition %in% c("Ctrl","LPS","ACD"))

# 颜色（按三组）
cond_colors <- c(Ctrl="#aec7e8", LPS="#ff7f0e", ACD="#17becf")

# —— 小提琴图（整体，不分 cluster）——
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

# —— 可选：按 cluster 分面的小提琴（看组间在各 cluster 的差异）——
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

# ===== D) 导出：PDF/SVG + 统计汇总 =====
out_dir <- file.path(data_dir, "M1_M2_scores"); dir.create(out_dir, showWarnings = FALSE)

# 整体三组小提琴
ggsave(file.path(out_dir, "M1_M2_scores_violin_by_condition.pdf"), p_violin,
       device = cairo_pdf, width = 180, height = 140, units = "mm")
svglite::svglite(file.path(out_dir, "M1_M2_scores_violin_by_condition.svg"),
                 width = 180/25.4, height = 140/25.4); print(p_violin); dev.off()

# 分面（可选）
ggsave(file.path(out_dir, "M1_score_violin_by_condition_facet_cluster.pdf"), p_m1_vln_fac,
       device = cairo_pdf, width = 210, height = 220, units = "mm")
ggsave(file.path(out_dir, "M2_score_violin_by_condition_facet_cluster.pdf"), p_m2_vln_fac,
       device = cairo_pdf, width = 210, height = 220, units = "mm")

# 导出每细胞分数（含 condition）
write.csv(score_df[, c("cell_id","condition","cluster_original","M1_score","M2_score")],
          file.path(out_dir, "M1_M2_scores_per_cell_with_condition.csv"), row.names = FALSE)

# 三组汇总（中位数/IQR）
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

# 三组 × cluster 细分汇总（用于论文补充表）
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
  library(ggpubr)   # ← 新增，用于显著性检验与标注
  library(rstatix)  # ← 新增，用于统计结果表格
})

# ===== E) 统计三组间差异（Kruskal + 两两比较） =====
## 1) Kruskal–Wallis 检验整体差异
kw_M1 <- kruskal_test(score_df, M1_score ~ condition)
kw_M2 <- kruskal_test(score_df, M2_score ~ condition)
kw_M1; kw_M2   # 控制台输出

## 2) 若整体有差异，则做两两比较（Dunn's test 带 Bonferroni 校正）
pairwise_M1 <- dunn_test(score_df, M1_score ~ condition, p.adjust.method = "bonferroni")
pairwise_M2 <- dunn_test(score_df, M2_score ~ condition, p.adjust.method = "bonferroni")

## 3) 保存结果表
write.csv(kw_M1,  file.path(out_dir, "Kruskal_M1_score.csv"), row.names = FALSE)
write.csv(kw_M2,  file.path(out_dir, "Kruskal_M2_score.csv"), row.names = FALSE)
write.csv(pairwise_M1, file.path(out_dir, "Dunn_pairwise_M1_score.csv"), row.names = FALSE)
write.csv(pairwise_M2, file.path(out_dir, "Dunn_pairwise_M2_score.csv"), row.names = FALSE)
# ==== 在图上加显著性标注 ====
# M1
p_m1_vln_sig <- ggviolin(score_df, x = "condition", y = "M1_score",
                         fill = "condition", palette = cond_colors,
                         add = "boxplot", width = 0.9, trim = TRUE) +
  stat_compare_means(method = "kruskal.test", label.y = max(score_df$M1_score)*1.05) +  # 整体p值
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

# 保存带显著性的小提琴图
ggsave(file.path(out_dir, "M1_M2_scores_violin_by_condition_with_significance.pdf"),
       p_violin_sig, device = cairo_pdf, width = 180, height = 150, units = "mm")
svglite::svglite(file.path(out_dir, "M1_M2_scores_violin_by_condition_with_significance.svg"),
                 width = 180/25.4, height = 150/25.4); print(p_violin_sig); dev.off()

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
  library(rstatix)
  library(svglite)
})

# ============ 0) 准备：标准化 condition、配色、cluster 因子水平 ============
if(!"condition" %in% names(score_df)) {
  score_df <- score_df %>%
    mutate(
      condition = case_when(
        stringr::str_detect(toupper(group), "\\bCTRL?\\b|CONTROL") ~ "Ctrl",
        stringr::str_detect(toupper(group), "\\bLPS\\b") ~ "LPS",
        stringr::str_detect(toupper(group), "A[-_ ]?CDS?|\\bACD\\b") ~ "ACD",
        TRUE ~ "Other"
      )
    )
}

score_df <- score_df %>%
  filter(condition %in% c("Ctrl","LPS","ACD")) %>%
  mutate(
    condition = factor(condition, levels = c("Ctrl","LPS","ACD")),
    cluster_original = factor(as.character(cluster_original), levels = as.character(0:11))
  ) %>% droplevels()

# 配色：避免 ggpubr 与 scale_* 冲突，这里只用 ggpubr 的 palette，传“无名向量”
cond_colors_named <- c(Ctrl="#aec7e8", LPS="#ff7f0e", ACD="#17becf")
cond_colors <- unname(cond_colors_named)

# 你想标注的三对比较：
comparisons_list <- list(c("Ctrl","LPS"), c("LPS","ACD"), c("Ctrl","ACD"))

# 一个小工具：为分组后的两两比较结果生成每个 cluster 的 y.position
.make_ypos <- function(df_scores, y_col, by_col = "cluster_original", bump = 0.02) {
  # 每个 cluster 的最大 y
  ymax_tbl <- df_scores %>%
    group_by(.data[[by_col]]) %>%
    summarise(ymax = max(.data[[y_col]], na.rm = TRUE), .groups = "drop")
  # 为每个 cluster 的每条比较递增抬高位置（避免文字重叠）
  df_out <- df_scores %>%
    group_by(.data[[by_col]]) %>%
    mutate(.rank = row_number()) %>%
    ungroup() %>%
    left_join(ymax_tbl, by = by_col) %>%
    mutate(y.position = ymax + .rank * bump * pmax(ymax, 1)) %>% # bump 按范围比例抬升
    select(-.rank)
  df_out
}

# ============ 1) 各 cluster 内统计：Kruskal（整体）+ Dunn（两两） ============
# —— M1：Kruskal（稳健抓取列名）
kw_M1_by_cluster <- score_df %>%
  dplyr::group_by(cluster_original) %>%
  rstatix::kruskal_test(M1_score ~ condition) %>%
  tibble::as_tibble() %>%                    # ← 先转 tibble
  dplyr::rename(
    statistic_M1 = .data$statistic,
    p_M1         = .data$p
  ) %>%
  dplyr::select(cluster_original, statistic_M1, p_M1) %>%
  dplyr::ungroup()

kw_M2_by_cluster <- score_df %>%
  dplyr::group_by(cluster_original) %>%
  rstatix::kruskal_test(M2_score ~ condition) %>%
  tibble::as_tibble() %>%
  dplyr::rename(
    statistic_M2 = .data$statistic,
    p_M2         = .data$p
  ) %>%
  dplyr::select(cluster_original, statistic_M2, p_M2) %>%
  dplyr::ungroup()
# 用 score_df 计算每个 cluster 的 ymax，再给 pairwise 表添加 y.position
.make_ypos <- function(pw_tbl, score_df, y_col, by_col = "cluster_original", bump = 0.03) {
  ymax_tbl <- score_df %>%
    dplyr::group_by(.data[[by_col]]) %>%
    dplyr::summarise(ymax = max(.data[[y_col]], na.rm = TRUE), .groups = "drop")
  
  pw_tbl %>%
    dplyr::group_by(.data[[by_col]]) %>%
    dplyr::mutate(.rank = dplyr::row_number()) %>%
    dplyr::ungroup() %>%
    dplyr::left_join(ymax_tbl, by = by_col) %>%
    dplyr::mutate(y.position = ymax + .rank * bump * pmax(ymax, 1)) %>%
    dplyr::select(-.rank)
}
# M1
pw_M1_by_cluster <- score_df %>%
  dplyr::group_by(cluster_original) %>%
  rstatix::dunn_test(M1_score ~ condition, p.adjust.method = "bonferroni") %>%
  tibble::as_tibble() %>%
  dplyr::ungroup() %>%
  dplyr::filter(paste(group1, group2) %in% c("Ctrl LPS","LPS ACD","Ctrl ACD")) %>%
  rstatix::add_significance(p.col = "p.adj") %>%
  dplyr::group_by(cluster_original) %>%
  dplyr::mutate(comp_order = dplyr::row_number()) %>%
  dplyr::ungroup()

pw_M1_by_cluster <- .make_ypos(
  pw_tbl   = pw_M1_by_cluster,
  score_df = score_df,
  y_col    = "M1_score",
  by_col   = "cluster_original",
  bump     = 0.03
)

# M2
pw_M2_by_cluster <- score_df %>%
  dplyr::group_by(cluster_original) %>%
  rstatix::dunn_test(M2_score ~ condition, p.adjust.method = "bonferroni") %>%
  tibble::as_tibble() %>%
  dplyr::ungroup() %>%
  dplyr::filter(paste(group1, group2) %in% c("Ctrl LPS","LPS ACD","Ctrl ACD")) %>%
  rstatix::add_significance(p.col = "p.adj") %>%
  dplyr::group_by(cluster_original) %>%
  dplyr::mutate(comp_order = dplyr::row_number()) %>%
  dplyr::ungroup()

pw_M2_by_cluster <- .make_ypos(
  pw_tbl   = pw_M2_by_cluster,
  score_df = score_df,
  y_col    = "M2_score",
  by_col   = "cluster_original",
  bump     = 0.03
)

# —— 统一生成 p.signif（M1）——
pw_M1_by_cluster <- pw_M1_by_cluster %>%
  dplyr::mutate(p.adj = suppressWarnings(as.numeric(p.adj)))

if(!"p.signif" %in% names(pw_M1_by_cluster)) {
  if("p.adj.signif" %in% names(pw_M1_by_cluster)) {
    pw_M1_by_cluster <- dplyr::rename(pw_M1_by_cluster, p.signif = p.adj.signif)
  } else {
    pw_M1_by_cluster <- pw_M1_by_cluster %>%
      dplyr::mutate(
        p.signif = dplyr::case_when(
          is.na(p.adj)        ~ NA_character_,
          p.adj <= 1e-4       ~ "****",
          p.adj <= 1e-3       ~ "***",
          p.adj <= 1e-2       ~ "**",
          p.adj <= 0.05       ~ "*",
          TRUE                ~ "ns"
        )
      )
  }
}

# —— 统一生成 p.signif（M2）——
pw_M2_by_cluster <- pw_M2_by_cluster %>%
  dplyr::mutate(p.adj = suppressWarnings(as.numeric(p.adj)))

if(!"p.signif" %in% names(pw_M2_by_cluster)) {
  if("p.adj.signif" %in% names(pw_M2_by_cluster)) {
    pw_M2_by_cluster <- dplyr::rename(pw_M2_by_cluster, p.signif = p.adj.signif)
  } else {
    pw_M2_by_cluster <- pw_M2_by_cluster %>%
      dplyr::mutate(
        p.signif = dplyr::case_when(
          is.na(p.adj)        ~ NA_character_,
          p.adj <= 1e-4       ~ "****",
          p.adj <= 1e-3       ~ "***",
          p.adj <= 1e-2       ~ "**",
          p.adj <= 0.05       ~ "*",
          TRUE                ~ "ns"
        )
      )
  }
}

# ============ 2) 分面小提琴 + 显著性标注（每个 cluster 单独面板） ============
# 注意：不要再叠加 scale_fill_manual，否则会与 ggpubr 的 palette 冲突
p_M1_clusters_sig <- ggviolin(
  score_df, x = "condition", y = "M1_score",
  fill = "condition", palette = cond_colors,
  add = "boxplot", width = 0.9, trim = TRUE
) +
  facet_wrap(~ cluster_original, ncol = 4, drop = FALSE) +
  ggpubr::stat_pvalue_manual(
    pw_M1_by_cluster,
    label = "p.signif",
    y.position = "y.position",
    tip.length = 0.01,
    hide.ns = TRUE
  ) +
  labs(title = "M1 score by condition (per cluster)", x = NULL, y = "M1 score (log1p CPM mean)") +
  theme_minimal(base_size = 10) +
  theme(legend.position = "none")

p_M2_clusters_sig <- ggviolin(
  score_df, x = "condition", y = "M2_score",
  fill = "condition", palette = cond_colors,
  add = "boxplot", width = 0.9, trim = TRUE
) +
  facet_wrap(~ cluster_original, ncol = 4, drop = FALSE) +
  ggpubr::stat_pvalue_manual(
    pw_M2_by_cluster,
    label = "p.signif",
    y.position = "y.position",
    tip.length = 0.01,
    hide.ns = TRUE
  ) +
  labs(title = "M2 score by condition (per cluster)", x = NULL, y = "M2 score (log1p CPM mean)") +
  theme_minimal(base_size = 10) +
  theme(legend.position = "none")

# ============ 3) 保存结果 ============
out_dir <- file.path(data_dir, "M1_M2_scores"); dir.create(out_dir, showWarnings = FALSE)

# 表格
write.csv(kw_M1_by_cluster, file.path(out_dir, "Kruskal_M1_by_cluster.csv"), row.names = FALSE)
write.csv(kw_M2_by_cluster, file.path(out_dir, "Kruskal_M2_by_cluster.csv"), row.names = FALSE)
write.csv(pw_M1_by_cluster, file.path(out_dir, "Dunn_pairwise_M1_by_cluster.csv"), row.names = FALSE)
write.csv(pw_M2_by_cluster, file.path(out_dir, "Dunn_pairwise_M2_by_cluster.csv"), row.names = FALSE)

# 图形（PDF+SVG）
ggsave(file.path(out_dir, "M1_score_violin_by_condition_facet_cluster_with_sig.pdf"),
       p_M1_clusters_sig, device = cairo_pdf, width = 210, height = 220, units = "mm")
svglite::svglite(file.path(out_dir, "M1_score_violin_by_condition_facet_cluster_with_sig.svg"),
                 width = 210/25.4, height = 220/25.4); print(p_M1_clusters_sig); dev.off()

ggsave(file.path(out_dir, "M2_score_violin_by_condition_facet_cluster_with_sig.pdf"),
       p_M2_clusters_sig, device = cairo_pdf, width = 210, height = 220, units = "mm")
svglite::svglite(file.path(out_dir, "M2_score_violin_by_condition_facet_cluster_with_sig.svg"),
                 width = 210/25.4, height = 220/25.4); print(p_M2_clusters_sig); dev.off()
