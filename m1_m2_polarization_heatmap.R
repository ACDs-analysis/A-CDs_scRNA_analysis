suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(ggplot2)
  library(forcats)
  library(svglite)
})

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
data_dir <- "C:/Users/Documents/巨噬细胞重新聚类/upp5/表达矩阵/gene"  # ← change to your path
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


# ================== 3) Compute average expression per cluster (efficient approach for sparse matrix) ==================
# Keep only cells with cluster labels
ok <- !is.na(cell_meta$cluster_original)
expr_all  <- expr_all[, ok, drop = FALSE]
cell_meta <- cell_meta[ok, , drop = FALSE]

# Keep only clusters that actually exist according to the preset order
clus_levels_plot <- cluster_levels[cluster_levels %in% as.character(cell_meta$cluster_original)]
cell_meta$cluster_original <- factor(cell_meta$cluster_original, levels = clus_levels_plot)

# Design matrix: one-hot per cluster
mm <- model.matrix(~ 0 + cell_meta$cluster_original)
colnames(mm) <- levels(cell_meta$cluster_original)

# Cluster average = expr %*% one-hot / number of cells per cluster
n_per_cluster <- Matrix::colSums(mm)
avg_by_cluster <- (expr_all %*% mm)    # genes × clusters, total expression per cluster
avg_by_cluster <- sweep(avg_by_cluster, 2, n_per_cluster, "/")  # convert to average
# Convert to regular matrix for row standardization
avg_by_cluster <- as.matrix(avg_by_cluster)

# ================== 4) Prepare M1/M2 submatrices and row Z-score ==================
safe_intersect <- function(genes, rn){
  x <- intersect(genes, rn)
  if(length(x) == 0) stop("No matching genes found in the expression matrix; please check gene names/case.")
  x
}
row_zscore <- function(mat){
  mu  <- matrixStats::rowMeans2(mat)
  sdv <- sqrt(matrixStats::rowVars(mat) + 1e-8)
  z <- (mat - mu) / sdv
  z[is.nan(z)] <- 0
  z
}

# Remove leading/trailing single quotes from gene names
rownames(avg_by_cluster) <- gsub("^'|'$", "", rownames(avg_by_cluster))

# Convert to uppercase for consistency
rownames(avg_by_cluster) <- toupper(rownames(avg_by_cluster))
m1_genes <- toupper(trimws(m1_genes))
m2_genes <- toupper(trimws(m2_genes))

# Attempt to match again
m2_keep <- safe_intersect(m2_genes, rownames(avg_by_cluster))
m1_keep <- safe_intersect(m1_genes, rownames(avg_by_cluster))

length(m2_keep)
length(m1_keep)
setdiff(m2_genes, rownames(avg_by_cluster))[1:30]
m1_genes <- m1_genes[!is.na(m1_genes) & m1_genes != ""]
m2_genes <- m2_genes[!is.na(m2_genes) & m2_genes != ""]


mat_M2 <- avg_by_cluster[m2_keep, , drop = FALSE] %>% row_zscore()
mat_M1 <- avg_by_cluster[m1_keep, , drop = FALSE] %>% row_zscore()

# ================== 5) Convert to long format and plot (M2 on top, M1 at bottom; unified color scale) ==================
to_long <- function(zmat){
  df <- as.data.frame(zmat, check.names = FALSE)
  df$Gene <- rownames(zmat)
  tidyr::pivot_longer(df, -Gene, names_to = "Cluster", values_to = "Z") %>%
    dplyr::mutate(
      Cluster = factor(Cluster, levels = clus_levels_plot) # maintain cluster order
    )
}
df_M2 <- to_long(mat_M2) %>% dplyr::mutate(Set = "M2")
df_M1 <- to_long(mat_M1) %>% dplyr::mutate(Set = "M1")

# Keep gene display order as in your Excel (top to bottom)
df_M2$Gene <- factor(df_M2$Gene, levels = rev(rownames(mat_M2)))
df_M1$Gene <- factor(df_M1$Gene, levels = rev(rownames(mat_M1)))

heat_theme <- theme_minimal(base_size = 10) +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.grid = element_blank(),
    plot.margin = margin(6, 6, 6, 6)
  )

z_limits <- c(-3, 3)  # unified Z-score range; can be changed to e.g., c(-2,2) / c(-4,4)

library(scales)  # for squish

p_M2 <- ggplot(df_M2, aes(Cluster, Gene, fill = Z)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "#2166ac", mid = "white", high = "#b2182b",
    midpoint = 0, limits = c(-3, 3),
    oob = squish,          # ← squash out‑of‑range values to the extremes
    na.value = "white"
  ) + labs(title = "M2 genes") + heat_theme

p_M1 <- ggplot(df_M1, aes(Cluster, Gene, fill = Z)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "#2166ac", mid = "white", high = "#b2182b",
    midpoint = 0, limits = c(-3, 3),
    oob = squish,
    na.value = "white"
  ) + labs(title = "M1 genes") + heat_theme


# Combine: M2 above, M1 below
library(patchwork)
p_all <- p_M2 / p_M1 + plot_layout(heights = c(2, 1))
print(p_all)

# ================== 6) Save as SVG/PDF (editable text in AI) ==================
out_dir <- file.path(data_dir, "figs"); dir.create(out_dir, showWarnings = FALSE)

save_svg_pdf <- function(plot, name, width = 220, height = 200, units = "mm"){
  svglite::svglite(
    filename = file.path(out_dir, paste0(name, ".svg")),
    width  = switch(units,
                    "mm" = width/25.4,
                    "cm" = width/2.54,
                    "in" = width),
    height = switch(units,
                    "mm" = height/25.4,
                    "cm" = height/2.54,
                    "in" = height),
    system_fonts = list(sans = "Arial")
  )
  print(plot)
  dev.off()
  
  ggsave(
    filename = file.path(out_dir, paste0(name, ".pdf")),
    plot = plot,
    device = cairo_pdf,
    width = width,
    height = height,
    units = units
  )
}
save_svg_pdf(p_all, "M1_M2_heatmap", width = 220, height = 200, units = "mm")
# Rebuild the exported matrix
mat_export <- rbind(mat_M2, mat_M1)
df_export <- as.data.frame(mat_export)
df_export <- tibble::rownames_to_column(df_export, var = "Gene")

# Write CSV again
write.csv(
  df_export,
  file = file.path(out_dir, "M1_M2_heatmap_Zscore.csv"),
  row.names = FALSE
)

# Verify write success
list.files(out_dir)
