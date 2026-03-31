suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(stringr); library(Matrix)
  library(SingleCellExperiment); library(slingshot); library(irlba)
  library(ggplot2); library(svglite); library(matrixStats)
  library(tidyr); library(ggpubr); library(rstatix)
})
set.seed(1234)

# ========= Robust cluster number parsing function (strongly recommended) =========
parse_cluster_number <- function(x){
  x <- as.character(x)
  
  # Remove invisible spaces and collapse multiple spaces
  x <- gsub("\u00A0", " ", x, fixed = TRUE)
  x <- stringr::str_squish(x)
  
  # Strategy 1: explicit keyword + number (cluster / clus / cl / ec / endo)
  m1 <- stringr::str_match(
    x,
    "(?i)(cluster|clus|cl|endo|ec)[^0-9]*([0-9]+)"
  )[,3]
  
  # Strategy 2: standalone number at the end of the string
  m2 <- stringr::str_match(
    x,
    "(?<![0-9])([0-9]{1,3})(?![0-9])$"
  )[,2]
  
  cl <- ifelse(!is.na(m1), m1, m2)
  cl <- sub("^0+", "", cl)   # Remove leading zeros
  cl[cl == ""] <- NA
  cl
}


# ================== 1) Endothelial cells: read and merge (CSV, one file per sample/cluster) ==================
data_dir <- "E:/科研/单细胞测序/EC gene"   # ← Endothelial cell path
files <- list.files(data_dir, pattern="\\.csv$", full.names=TRUE)
stopifnot("No CSV files found" = length(files) > 0)

read_expr_csv <- function(path){
  df <- tryCatch(
    as.data.frame(utils::read.csv(path, check.names = FALSE), stringsAsFactors = FALSE),
    error = function(e) as.data.frame(utils::read.csv(path, check.names = FALSE, fileEncoding="UTF-8-BOM"), stringsAsFactors = FALSE)
  )
  stopifnot(ncol(df) >= 2)
  
  # Assume column 1 is gene_id, column 2 is gene_symbol (use gene_id if empty)
  gene_id <- as.character(df[[1]])
  gene_symbol <- as.character(df[[2]])
  idx <- which(is.na(gene_symbol) | gene_symbol == "")
  gene_symbol[idx] <- gene_id[idx]
  
  genes <- make.unique(ifelse(is.na(gene_symbol)|gene_symbol=="","NA_gene",gene_symbol))
  
  # Accommodate different export formats: prefer starting from column 5, else column 3, else column 2
  expr_start <- if (ncol(df) >= 5) 5 else if (ncol(df) >= 3) 3 else 2
  expr_part <- df[, expr_start:ncol(df), drop = FALSE]
  
  # Force numeric conversion (non-numeric -> NA -> 0)
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
  sample_id <- str_trim(str_remove(basename(f), "\\.csv$"))   # e.g., Ctrl 1 Cluster 6
  # Accommodate case "Cluster 6" / "cluster 6"
  sample_id <- sub("\\.[cC][sS][vV]$", "", basename(f))
  sample_id <- gsub("\u00A0", " ", sample_id, fixed = TRUE)
  sample_id <- stringr::str_squish(sample_id)
  
  clx <- parse_cluster_number(sample_id)
  
  # Remove "cluster x" from sample_id to get group (Ctrl 1 / LPS 2 / ACD 3 …)
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

# —— Union alignment and merging
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
cluster_levels <- sort(unique(na.omit(as.character(cell_meta$cluster_original))))
cell_meta$cluster_original <- factor(as.character(cell_meta$cluster_original), levels = cluster_levels)

rownames(cell_meta) <- cell_meta$cell_id

# ================== 2) Select key leak genes & preprocessing (endothelial) ==================
# Convert gene names to uppercase for robust matching
rownames(expr_all) <- gsub("^'|'$", "", rownames(expr_all))
rownames(expr_all) <- trimws(rownames(expr_all))
rownames(expr_all) <- toupper(rownames(expr_all))

leak_genes <- c(
  "ANGPT2","MYLK","SELE","VCAM1","SELE","SELP","MMP9","IL1B","IL6","TNF","CXCL2"
) 
# Alternative candidates: "ICAM1","VCAM1","SELE","SELP","ANGPT2","EDN1","PTGS2","IL6ST","PLVAP"
leak_genes <- unique(toupper(leak_genes))



# Extract conditions from group: Ctrl/LPS/ACD (preserve original sample labels)
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
  libsize[libsize == 0] <- 1
  cpm <- t(t(sparse_counts) / libsize) * 1e6
  log1p(cpm)
}
expr_log1p_cpm <- calc_log1p_cpm(expr_all)

# Check for missing genes
present <- leak_genes[leak_genes %in% rownames(expr_log1p_cpm)]
missing <- setdiff(leak_genes, present)
if (length(missing) > 0) {
  message("The following leak genes were not found in the matrix:\n",
          paste(missing, collapse = ", "))
}
stopifnot(length(present) > 0)

expr_leak <- as.matrix(expr_log1p_cpm[present, , drop = FALSE])

df_long <- as.data.frame(t(expr_leak), check.names = FALSE)
df_long$cell_id <- rownames(df_long)

df_long <- tidyr::pivot_longer(
  df_long,
  cols = all_of(present),
  names_to = "gene",
  values_to = "log1pCPM"
)

plot_df <- dplyr::left_join(
  df_long,
  cell_meta[, c("cell_id","condition")],
  by = "cell_id"
) %>%
  dplyr::filter(!is.na(condition))



# ========= Violin plot (faceted by gene) data preparation =========
df_long <- as.data.frame(t(expr_leak), check.names = FALSE)
df_long$cell_id <- rownames(df_long)

df_long <- tidyr::pivot_longer(
  df_long,
  cols = all_of(present),
  names_to = "gene",
  values_to = "log1pCPM",
  values_drop_na = TRUE
)

plot_df <- dplyr::left_join(
  df_long,
  cell_meta[, c("cell_id", "condition", "cluster_original")],
  by = "cell_id"
) %>%
  dplyr::filter(!is.na(condition))


cond_colors <- c("Ctrl"="#aec7e8","LPS"="#ff7f0e","ACD"="#17becf")

p1 <- ggplot(plot_df, aes(condition, log1pCPM, fill = condition)) +
  geom_violin(trim=TRUE, scale="width") +
  geom_boxplot(width=0.15, fill="white", outlier.size=0.3) +
  facet_wrap(~gene, scales="free_y", ncol=4) +
  scale_fill_manual(values = cond_colors) +
  theme_minimal() +
  labs(title = "Endothelial Leak genes",
       y = "log1p CPM", x = NULL)

print(p1)


# ================== Leak score ==================
leak_score <- Matrix::colMeans(expr_leak)

leak_df <- data.frame(
  cell_id   = colnames(expr_leak),
  LeakScore = as.numeric(leak_score),
  stringsAsFactors = FALSE
) %>%
  dplyr::left_join(
    cell_meta[, c("cell_id","condition","cluster_original")],
    by = "cell_id"
  ) %>%
  dplyr::filter(!is.na(condition), !is.na(cluster_original)) %>%
  dplyr::mutate(
    cluster_original = as.character(cluster_original),
    condition = factor(condition, levels = c("Ctrl","LPS","ACD"))
  )

# ===== Keep only clusters 6 and 9 =====
leak_df2 <- leak_df %>% dplyr::filter(cluster_original %in% c("6","9"))
stopifnot(nrow(leak_df2) > 0)

leak_df_6 <- leak_df2 %>% dplyr::filter(cluster_original == "6")
leak_df_9 <- leak_df2 %>% dplyr::filter(cluster_original == "9")
stopifnot(nrow(leak_df_6) > 0, nrow(leak_df_9) > 0)

# ===== Raincloud plotting function =====
suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(ggridges)
})

cond_colors <- c("Ctrl"="#aec7e8","LPS"="#ff7f0e","ACD"="#17becf")

plot_raincloud_leak <- function(df, cl, score_col = "LeakScore") {
  y_levels    <- c("Ctrl","LPS","ACD")
  gap_factor  <- 3
  cloud_scale <- 0.90
  down_shift  <- -0.45
  point_alpha <- 0.45
  point_height <- 0.25
  box_alpha   <- 0.16
  box_width   <- 0.6
  box_lwd     <- 0.3
  
  dfc <- df %>%
    dplyr::filter(cluster_original == cl) %>%
    dplyr::mutate(condition = factor(condition, levels = y_levels)) %>%
    dplyr::filter(!is.na(.data[[score_col]])) %>%
    dplyr::mutate(
      y_num = as.numeric(condition),
      y_num = ((length(y_levels) + 1) - y_num) * gap_factor
    )
  
  x_max <- max(dfc[[score_col]], na.rm = TRUE)
  x_min <- min(dfc[[score_col]], na.rm = TRUE)
  x_pad <- 0.12 * (x_max - x_min + 1e-9)
  
  deep_colors <- vapply(cond_colors, function(clr) {
    grDevices::adjustcolor(clr, red.f = 0.45, green.f = 0.45, blue.f = 0.45)
  }, FUN.VALUE = character(1))
  names(deep_colors) <- names(cond_colors)
  
  med_df <- dfc %>%
    dplyr::group_by(condition, y_num) %>%
    dplyr::summarise(
      med = median(.data[[score_col]], na.rm = TRUE),
      n = dplyr::n(),
      .groups = "drop"
    )
  
  ggplot(dfc, aes(x = .data[[score_col]])) +
    ggridges::geom_density_ridges(
      aes(y = y_num, fill = condition),
      scale = cloud_scale,
      rel_min_height = 0.01,
      color = NA, alpha = 0.65,
      from = x_min, to = x_max
    ) +
    geom_point(
      aes(y = y_num + down_shift, color = condition, fill = condition),
      position = position_jitter(height = point_height, width = 0, seed = 1234),
      size = 1.05, alpha = point_alpha, shape = 21, stroke = 0.25
    ) +
    geom_boxplot(
      aes(y = y_num + down_shift, color = condition, fill = condition, group = condition),
      width = box_width, outlier.shape = NA,
      linewidth = box_lwd, alpha = box_alpha
    ) +
    geom_text(
      data = med_df,
      aes(x = med, y = y_num, label = sprintf("%.2f", med), color = condition),
      inherit.aes = FALSE,
      vjust = -1.35, fontface = "bold", size = 4
    ) +
    geom_text(
      data = med_df,
      aes(x = x_max + x_pad, y = y_num + down_shift, label = paste0("n = ", n), color = condition),
      inherit.aes = FALSE,
      hjust = 0, vjust = 0.25, size = 3.6
    ) +
    scale_fill_manual(values = cond_colors, drop = FALSE) +
    scale_color_manual(values = deep_colors, drop = FALSE) +
    scale_y_continuous(
      breaks = rev(seq_along(y_levels) * gap_factor),
      labels = y_levels,
      expand = expansion(mult = c(0.10, 0.14))
    ) +
    coord_cartesian(xlim = c(x_min, x_max + 2.2 * x_pad), clip = "off") +
    labs(title = paste0("LeakScore — cluster ", cl), x = "LeakScore", y = NULL) +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold"),
      axis.ticks.y = element_blank(),
      plot.margin = margin(6, 40, 6, 6)
    )
}

p6 <- plot_raincloud_leak(leak_df_6, cl = "6")
p9 <- plot_raincloud_leak(leak_df_9, cl = "9")
print(p6); print(p9)

# ===== Save =====
out_dir <- "figs_leak_raincloud_69_pretty"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

ggsave(file.path(out_dir, "LeakScore_cluster6_raincloud_pretty.pdf"), p6,
       width = 120, height = 85, units = "mm", device = grDevices::cairo_pdf)

ggsave(file.path(out_dir, "LeakScore_cluster9_raincloud_pretty.pdf"), p9,
       width = 120, height = 85, units = "mm", device = grDevices::cairo_pdf)