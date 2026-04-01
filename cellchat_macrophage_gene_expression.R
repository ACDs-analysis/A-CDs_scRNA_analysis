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
data_dir <- "E:/科研/单细胞测序/gene"   # ← change to your path
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
# ================== 2) Select key ECM genes & preprocessing ==================
# Convert gene names to uppercase for robust matching
rownames(expr_all) <- gsub("^'|'$", "", rownames(expr_all))
rownames(expr_all) <- trimws(rownames(expr_all))
rownames(expr_all) <- toupper(rownames(expr_all))

ecm_core_genes <- c(
  "CCL5","TNF","SPP1","FN1","VEGFA","ITGAL","ITGB2","PECAM1","OSM","COL1A2","COL4A2","COL1A1" 
)
ecm_core_genes <- unique(toupper(ecm_core_genes))

# Extract conditions from sample_file as Ctrl/LPS/ACD (preserve original sample labels)
# Your sample_file looks like "Ctrl 1 cluster 0", "ACD 3 cluster 4", etc.
cell_meta$condition <- dplyr::case_when(
  grepl("^\\s*CTRL\\b",  toupper(cell_meta$group)) ~ "Ctrl",
  grepl("^\\s*LPS\\b",   toupper(cell_meta$group)) ~ "LPS",
  grepl("^\\s*A\\s*CDS?\\b|^\\s*ACD\\b", toupper(cell_meta$group)) ~ "ACD",
  TRUE ~ NA_character_
)
cell_meta <- cell_meta[!is.na(cell_meta$condition), , drop = FALSE]
cell_meta$condition <- factor(cell_meta$condition, levels = c("Ctrl","LPS","ACD"))

# Synchronize expression matrix to retained cells
expr_all <- expr_all[, rownames(cell_meta), drop = FALSE]

# Compute log1p-CPM (per cell normalization to 1e6)
calc_log1p_cpm <- function(sparse_counts) {
  libsize <- Matrix::colSums(sparse_counts)
  libsize[libsize == 0] <- 1  # prevent division by zero
  # CPM = counts / libsize * 1e6
  cpm <- t(t(sparse_counts) / libsize) * 1e6
  log1p(cpm)
}
expr_log1p_cpm <- calc_log1p_cpm(expr_all)

# Check for missing genes
present <- ecm_core_genes[ecm_core_genes %in% rownames(expr_log1p_cpm)]
missing <- setdiff(ecm_core_genes, present)
if (length(missing) > 0) {
  message("The following ECM genes were not found in the matrix (please check aliases/species symbols):\n",
          paste(missing, collapse = ", "))
}
stopifnot(length(present) > 0)

# Extract expression for these genes and reshape into long format
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

# ========= Violin plot (faceted by gene) =========
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
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),   # Remove minor grid lines
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4)
  )


# ========= Statistical annotation: Kruskal + pairwise Wilcoxon (BH) =========
suppressPackageStartupMessages({ library(ggpubr); library(rstatix) })

valid_genes10 <- plot_df %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise(n_groups = dplyr::n_distinct(condition[!is.na(log1pCPM)]), .groups="drop") %>%
  dplyr::filter(n_groups >= 2) %>% dplyr::pull(gene)

plot_df_stat <- plot_df %>% dplyr::filter(gene %in% valid_genes10)

# Overall (optional)
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

# Assign non-overlapping y.positions for each gene
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

# Keep only comparisons that actually exist (some genes may lack certain groups)
present_levels_tbl <- plot_df_stat %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise(levels_present = list(unique(as.character(condition))), .groups = "drop")

pw_ecm <- pw_ecm %>%
  dplyr::left_join(present_levels_tbl, by = "gene") %>%
  dplyr::rowwise() %>%
  dplyr::filter(all(c(group1, group2) %in% levels_present)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(p.label = dplyr::if_else(!is.na(p.adj.signif), p.adj.signif, rstatix::p_format(p.adj)))

# Add significance annotations to p_ecm10
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


# ========= Save =========
dir.create("figs_ecm", showWarnings = FALSE)
ggsave(
  filename = file.path("figs_ecm", "ECM_10genes_violin_withSig.pdf"),
  plot = p_ecm10_sig,
  width = 240, height = 160, units = "mm",
  device = grDevices::cairo_pdf # or omit device, let ggsave choose based on extension
)
