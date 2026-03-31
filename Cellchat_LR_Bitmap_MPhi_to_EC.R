suppressPackageStartupupMessages({
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
data_dir <- "E:/科研/单细胞测序/gene"  # ← change to your path
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

# ================== 2) Read and merge endothelial cell CSV files ==================
library(readr)

endo_dir <- "E:/科研/单细胞测序/EC gene" # ← change to your path
endo_files <- list.files(endo_dir, pattern="\\.csv$", full.names=TRUE)
stopifnot("No endothelial CSV files found" = length(endo_files) > 0)

read_expr_csv <- function(path){
  # Assumes: column1 Gene ID, column2 Gene Symbol, column3 Type, column4 total, column5 onward are cell expressions
  df <- suppressMessages(readr::read_csv(path, show_col_types = FALSE))
  stopifnot(ncol(df) >= 5)
  gene_id     <- as.character(df[[1]])
  gene_symbol <- as.character(df[[2]])
  gene_symbol[is.na(gene_symbol) | gene_symbol==""] <- gene_id
  genes <- make.unique(ifelse(is.na(gene_symbol)|gene_symbol=="","NA_gene",gene_symbol))
  expr_part <- df[, 5:ncol(df), drop=FALSE]
  
  # Force numeric conversion
  for(j in seq_len(ncol(expr_part))) expr_part[[j]] <- suppressWarnings(as.numeric(expr_part[[j]]))
  m <- as.matrix(expr_part)
  rownames(m) <- genes
  colnames(m) <- make.unique(colnames(m))
  m[is.na(m)] <- 0
  storage.mode(m) <- "double"
  m
}

endo_expr_list <- list(); endo_meta_list <- list()
for(f in endo_files){
  # Example: ACD 1 cluster 6.csv / LPS 2 cluster 24.csv etc.
  sample_id <- stringr::str_trim(stringr::str_remove(basename(f), "\\.csv$"))
  clx <- stringr::str_match(sample_id, "cluster\\s*(\\d+)")[,2]      # 6/8/9/24
  group_id <- stringr::str_trim(stringr::str_replace(sample_id, "cluster\\s*\\d+", ""))
  
  mat <- read_expr_csv(f)
  if(ncol(mat)==0) next
  
  new_cols <- paste0(colnames(mat), "_", sample_id)
  colnames(mat) <- new_cols
  endo_expr_list[[f]] <- mat
  
  meta <- data.frame(
    cell_id = new_cols,
    sample_file = sample_id,
    group = group_id,                      # Ctrl / LPS / ACD
    cluster_original = ifelse(is.na(clx), NA, clx),
    stringsAsFactors = FALSE
  )
  rownames(meta) <- meta$cell_id
  endo_meta_list[[f]] <- meta
}
stopifnot(length(endo_expr_list) > 0)

# ===== Endothelial matrix union alignment + merge =====
endo_all_genes <- unique(unlist(lapply(endo_expr_list, rownames)))
endo_expr_list_reindexed <- lapply(endo_expr_list, function(m){
  miss <- setdiff(endo_all_genes, rownames(m))
  if(length(miss)>0){
    m <- rbind(m, matrix(0, nrow=length(miss), ncol=ncol(m),
                         dimnames=list(miss, colnames(m))))
  }
  m <- m[endo_all_genes,,drop=FALSE]
  m[is.na(m)] <- 0
  storage.mode(m) <- "double"
  m
})
endo_expr_all <- do.call(cbind, endo_expr_list_reindexed)
endo_expr_all <- Matrix(as.matrix(endo_expr_all), sparse=TRUE)
endo_meta <- do.call(rbind, endo_meta_list)
rownames(endo_meta) <- endo_meta$cell_id
endo_meta$cluster_original <- factor(as.character(endo_meta$cluster_original))

# ================== 3) Merge with macrophage matrix (union of genes) ==================
# Explanation: use expr_all (macrophage) and cell_meta from previous code
# Add cell_type column
cell_meta$cell_type <- "Macrophage"
endo_meta$cell_type <- "Endothelial"

# Union of genes
genes_union <- union(rownames(expr_all), rownames(endo_expr_all))

reindex_to <- function(M, genes_union){
  miss <- setdiff(genes_union, rownames(M))
  if(length(miss)>0){
    M <- rbind(M, Matrix(0, nrow=length(miss), ncol=ncol(M),
                         dimnames=list(miss, colnames(M))))
  }
  M[genes_union,,drop=FALSE]
}

expr_mac_reindexed  <- reindex_to(expr_all, genes_union)
expr_endo_reindexed <- reindex_to(endo_expr_all, genes_union)

expr_all_merged  <- Matrix(cbind(expr_mac_reindexed, expr_endo_reindexed), sparse = TRUE)
cell_meta_merged <- rbind(cell_meta, endo_meta)

# Basic consistency checks
stopifnot(ncol(expr_all_merged) == nrow(cell_meta_merged))
stopifnot(identical(colnames(expr_all_merged), rownames(cell_meta_merged)))

# Unify factor levels for condition/cluster (optional)
cond_levels <- c("Ctrl","LPS","ACD")
cell_meta_merged$group <- factor(as.character(cell_meta_merged$group), levels = cond_levels)

mac_levels <- as.character(0:11)
cell_meta_merged$cluster_original[cell_meta_merged$cell_type=="Macrophage"] <-
  factor(as.character(cell_meta_merged$cluster_original[cell_meta_merged$cell_type=="Macrophage"]),
         levels = mac_levels)

# Optional: remove all-zero genes and empty cells (usually not needed)
# keep_gene <- Matrix::rowSums(expr_all_merged) > 0
# keep_cell <- Matrix::colSums(expr_all_merged) > 0
# expr_all_merged  <- expr_all_merged[keep_gene, keep_cell, drop=FALSE]
# cell_meta_merged <- cell_meta_merged[keep_cell,, drop=FALSE]

diagnose_data_issues <- function() {
  message("=== Data diagnosis report ===")
  
  # Check group distribution
  message("1. Group distribution:")
  print(table(cell_meta_merged$group, useNA = "always"))
  
  # Check cell type distribution
  message("2. Cell type distribution:")
  print(table(cell_meta_merged$cell_type))
  
  # Check expression matrix dimensions
  message("3. Expression matrix dimensions:")
  print(dim(expr_all_merged))
  
  # Check gene names
  message("4. Gene name examples:")
  print(head(rownames(expr_all_merged), 10))
  
  # Check for zero-expressed genes
  zero_genes <- Matrix::rowSums(expr_all_merged) == 0
  message("5. Number of zero-expressed genes:", sum(zero_genes))
  
  # Check sample file names
  message("6. Sample file name examples:")
  print(head(unique(cell_meta_merged$sample_file)))
}

# Run diagnosis
diagnose_data_issues()

# Improved group extraction function
extract_group_properly <- function(sample_file) {
  message("Processing sample: ", sample_file) # debug info
  
  # Method 1: split by space and take first part
  parts <- strsplit(trimws(sample_file), "\\s+")[[1]]
  if(length(parts) > 0) {
    group_candidate <- parts[1]
    # Standardize group names
    if(toupper(group_candidate) %in% c("CTRL", "CONTROL", "CONTROL1", "CONTROL2")) return("Ctrl")
    if(toupper(group_candidate) %in% c("LPS")) return("LPS")
    if(toupper(group_candidate) %in% c("ACD")) return("ACD")
  }
  
  # Method 2: direct string matching
  sample_lower <- tolower(trimws(sample_file))
  if(grepl("^acd", sample_lower)) return("ACD")
  if(grepl("^lps", sample_lower)) return("LPS") 
  if(grepl("^ctrl|^control", sample_lower)) return("Ctrl")
  
  return(NA)
}

# Test group extraction
test_group_extraction <- function() {
  test_samples <- c("ACD 1 cluster 0", "LPS 2 cluster 1", "Ctrl 1 cluster 2", "CONTROL 1 cluster 3")
  for(sample in test_samples) {
    group <- extract_group_properly(sample)
    message("Sample: '", sample, "' -> Group: '", group, "'")
  }
}

# Apply fixed group extraction
repair_group_information <- function() {
  message("Starting group information repair...")
  
  # Backup original sample_file
  cell_meta_merged$sample_file_original <- cell_meta_merged$sample_file
  
  # Apply new group extraction
  cell_meta_merged$group <- sapply(cell_meta_merged$sample_file, extract_group_properly)
  
  # Check results
  message("Group distribution after repair:")
  print(table(cell_meta_merged$group, useNA = "always"))
  
  # If there are still NAs, show those samples
  na_samples <- cell_meta_merged$sample_file[is.na(cell_meta_merged$group)]
  if(length(na_samples) > 0) {
    message("Failed to extract group for the following samples:")
    print(unique(na_samples))
  }
  
  return(cell_meta_merged)
}

# View all unique sample file names for debugging
debug_sample_names <- function() {
  message("All unique sample file names:")
  all_samples <- unique(cell_meta_merged$sample_file)
  for(i in seq_along(all_samples)) {
    message(sprintf("%3d: '%s'", i, all_samples[i]))
  }
  
  # Classify by starting pattern
  message("\nClassification by starting letters:")
  sample_starts <- substr(all_samples, 1, 3)
  print(table(sample_starts))
}

# Step 1: View all sample file name patterns
debug_sample_names()

# Step 2: Test group extraction
test_group_extraction()

# Step 3: Repair group information
cell_meta_merged <- repair_group_information()

# Step 4: Re-check groups
message("Final group distribution:")
final_groups <- table(cell_meta_merged$group)
print(final_groups)

# Step 5: Check cell counts per group
if(sum(!is.na(cell_meta_merged$group)) > 0) {
  message("\nCell counts per group:")
  for(group_name in names(final_groups)) {
    cell_count <- sum(cell_meta_merged$group == group_name, na.rm = TRUE)
    message(group_name, ": ", cell_count, " cells")
  }
}

# Set group as factor
cell_meta_merged$group <- factor(cell_meta_merged$group, levels = c("Ctrl", "LPS", "ACD"))

# Check cell type distribution per group
message("Detailed cell type distribution per group:")
cell_type_by_group <- table(cell_meta_merged$cell_type, cell_meta_merged$group)
print(cell_type_by_group)

# Extract cells per group
cells_ctrl <- rownames(subset(cell_meta_merged, group == "Ctrl"))
cells_lps  <- rownames(subset(cell_meta_merged, group == "LPS")) 
cells_acd  <- rownames(subset(cell_meta_merged, group == "ACD"))

message("Confirmed cell counts per group:")
message("Ctrl: ", length(cells_ctrl), " cells")
message("LPS: ", length(cells_lps), " cells") 
message("ACD: ", length(cells_acd), " cells")

# Extract expression matrices per group
expr_ctrl <- expr_all_merged[, cells_ctrl, drop = FALSE]
meta_ctrl <- cell_meta_merged[cells_ctrl, , drop = FALSE]

expr_lps <- expr_all_merged[, cells_lps, drop = FALSE]  
meta_lps <- cell_meta_merged[cells_lps, , drop = FALSE]

expr_acd <- expr_all_merged[, cells_acd, drop = FALSE]
meta_acd <- cell_meta_merged[cells_acd, , drop = FALSE]

# Clean gene names (important!)
clean_gene_names <- function(expr_mat) {
  genes <- rownames(expr_mat)
  # Remove single quotes
  cleaned <- gsub("'", "", genes)
  # Remove other special characters but keep normal characters
  cleaned <- gsub("[^a-zA-Z0-9._-]", "", cleaned)
  # Handle empty values
  cleaned[cleaned == "" | is.na(cleaned)] <- "UnknownGene"
  rownames(expr_mat) <- make.unique(cleaned)
  return(expr_mat)
}

# Apply gene name cleaning
expr_ctrl <- clean_gene_names(expr_ctrl)
expr_lps <- clean_gene_names(expr_lps) 
expr_acd <- clean_gene_names(expr_acd)

# Check cleaned gene names
message("Cleaned gene name examples:")
print(head(rownames(expr_ctrl), 10))

# Improved CellChat analysis function
run_cellchat_robust <- function(expr_mat, meta_df, group_name, species = "mouse") {
  message("\n====== Starting analysis for ", group_name, " group ======")
  
  # Basic checks
  if(ncol(expr_mat) < 10) {
    message("Warning: ", group_name, " group has insufficient cells, skipping analysis")
    return(NULL)
  }
  
  message("Data dimensions: ", paste(dim(expr_mat), collapse = " x "))
  message("Cell types: ", paste(unique(meta_df$cell_type), collapse = ", "))
  
  # Set database
  if(species == "mouse") {
    CellChatDB <- CellChatDB.mouse
    PPI_use <- PPI.mouse
  } else if(species == "human") {
    CellChatDB <- CellChatDB.human
    PPI_use <- PPI.human
  } else {
    stop("Please specify species: 'mouse' or 'human'")
  }
  
  # Check gene overlap with database
  db_genes <- unique(c(CellChatDB$interaction$ligand, CellChatDB$interaction$receptor))
  our_genes <- rownames(expr_mat)
  overlapping <- intersect(db_genes, our_genes)
  message("Number of genes overlapping with CellChat DB: ", length(overlapping))
  
  if(length(overlapping) < 50) {
    message("Warning: Few overlapping genes, may affect analysis quality")
    message("Overlapping gene examples: ", paste(head(overlapping, 10), collapse = ", "))
  }
  
  tryCatch({
    # Create CellChat object
    cellchat <- createCellChat(
      object = as.matrix(expr_mat),
      meta = meta_df,
      group.by = "cell_type"
    )
    
    cellchat@DB <- CellChatDB
    
    # Preprocessing pipeline
    cellchat <- subsetData(cellchat)
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)
    cellchat <- projectData(cellchat, PPI_use)
    
    # Compute communication probability
    cellchat <- computeCommunProb(
      cellchat,
      type = "truncatedMean",
      trim = 0.1,
      raw.use = FALSE
    )
    
    # Filter and aggregate
    cellchat <- filterCommunication(cellchat, min.cells = 5)
    cellchat <- computeCommunProbPathway(cellchat)
    cellchat <- aggregateNet(cellchat)
    
    message("✓ ", group_name, " group analysis completed")
    return(cellchat)
    
  }, error = function(e) {
    message("× ", group_name, " group analysis failed: ", e$message)
    return(NULL)
  })
}
library(CellChat); data(CellChatDB.mouse); data(PPI.mouse)

# Run analysis for each group
message("\nStarting CellChat analysis...")
cellchat_results <- list()

cellchat_results$Ctrl <- run_cellchat_robust(expr_ctrl, meta_ctrl, "Ctrl", species = "mouse")
cellchat_results$LPS <- run_cellchat_robust(expr_lps, meta_lps, "LPS", species = "mouse")
cellchat_results$ACD <- run_cellchat_robust(expr_acd, meta_acd, "ACD", species = "mouse")

# Remove failed analyses
cellchat_results <- cellchat_results[!sapply(cellchat_results, is.null)]

if(length(cellchat_results) == 0) {
  stop("All CellChat analyses failed, please check data and parameters")
}

message("\nNumber of groups successfully analyzed: ", length(cellchat_results))

# ========= Consistent with information flow: extract significant Mac → Endo LR pairs =========
SRC <- "Macrophage"
TGT <- "Endothelial"
PVAL_THRESH <- 0.05      # Consistent with information flow; change to 1 if no filtering desired

# Unified extraction function (automatically compatible with prob.mean / pval.mean)
extract_comm_by_group <- function(ch, g, src=SRC, tgt=TGT, thresh=PVAL_THRESH){
  tb <- tryCatch(
    subsetCommunication(ch, sources.use = src, targets.use = tgt, thresh = thresh),
    error = function(e) NULL
  )
  if (is.null(tb) || nrow(tb)==0) return(NULL)
  if (!"prob" %in% names(tb) && "prob.mean" %in% names(tb)) tb$prob <- tb$prob.mean
  if (!"pval" %in% names(tb) && "pval.mean" %in% names(tb)) tb$pval <- tb$pval.mean
  if (!"lr_pair" %in% names(tb) && all(c("ligand","receptor") %in% names(tb)))
    tb$lr_pair <- paste0(tb$ligand, " - ", tb$receptor)
  if (!"pathway_name" %in% names(tb) && "pathway" %in% names(tb))
    tb$pathway_name <- tb$pathway
  tb$group <- g
  tb$source <- src; tb$target <- tgt
  tb[, c("group","source","target","pathway_name","ligand","receptor","lr_pair","prob","pval")]
}

# Extract LR pairs from cellchat_results for three groups (same as information flow)
comm_all <- dplyr::bind_rows(lapply(names(cellchat_results), function(g){
  extract_comm_by_group(cellchat_results[[g]], g)
}))
stopifnot(nrow(comm_all) > 0)

# ========= Compute ACD-LPS, split into left/right (mutually exclusive) + control number =========
library(dplyr); library(tidyr); library(stringr)

min_prob <- 0            # Can be set to 1e-8 to remove very weak signals
TOPN_EACH <- 40          # Maximum number of LR pairs to display per panel (including core axis)

wide2 <- comm_all %>%
  filter(group %in% c("LPS","ACD")) %>%
  select(group, lr_pair, prob) %>%
  tidyr::pivot_wider(names_from = group, values_from = prob, values_fill = 0) %>%
  mutate(
    ACD  = coalesce(ACD, 0),
    LPS  = coalesce(LPS, 0),
    keep = (ACD > min_prob | LPS > min_prob),
    delta_ACDvsLPS = ACD - LPS,
    delta_abs = abs(delta_ACDvsLPS)
  ) %>%
  filter(keep)

# Core inflammatory axis (only for priority retention, not added to both sides)
core_pat <- "(?i)^(Tnf|Il6|Cxcl|Ccl)\\s*-|-(TNFR|IL6R|IL6ST|CXCR|ACKR[23]|CCR)"
core_keep_all <- wide2$lr_pair[ grepl(core_pat, wide2$lr_pair, perl = TRUE) ]

# Split by sign of Δ
left_tbl  <- wide2 %>% filter(delta_ACDvsLPS < 0)  %>% arrange(desc(delta_abs))
right_tbl <- wide2 %>% filter(delta_ACDvsLPS >= 0) %>% arrange(desc(delta_abs))

# Core axis per side (mutually exclusive)
left_core  <- intersect(core_keep_all,  left_tbl$lr_pair)
right_core <- intersect(core_keep_all, right_tbl$lr_pair)

# Supplement each side to TOPN_EACH
left_rest  <- setdiff(left_tbl$lr_pair,  left_core)  %>% head(max(0, TOPN_EACH - length(left_core)))
right_rest <- setdiff(right_tbl$lr_pair, right_core) %>% head(max(0, TOPN_EACH - length(right_core)))

left_ids  <- c(left_core,  left_rest)
right_ids <- c(right_core, right_rest)

# Generate plotting data
prep_plot_df <- function(id_set){
  df <- comm_all %>%
    filter(lr_pair %in% id_set, group %in% c("Ctrl","LPS","ACD")) %>%
    mutate(
      prob_abs  = prob,
      prob_eps  = pmax(prob_abs, 1e-12),
      prob_log10 = log10(prob_eps),
      group = factor(group, levels = c("Ctrl","LPS","ACD"))
    )
  ord <- wide2 %>% filter(lr_pair %in% id_set) %>% arrange(desc(delta_abs)) %>% pull(lr_pair)
  df$lr_pair <- factor(df$lr_pair, levels = ord)
  df
}

df_left  <- prep_plot_df(left_ids)
df_right <- prep_plot_df(right_ids)

# Unified color scale / point size (consistent with information flow, just visualization scaling)
cmin <- quantile(c(df_left$prob_log10, df_right$prob_log10), 0.05, na.rm=TRUE)
cmax <- quantile(c(df_left$prob_log10, df_right$prob_log10), 0.95, na.rm=TRUE)
cmid <- median(c(df_left$prob_log10, df_right$prob_log10), na.rm=TRUE)
cap99 <- quantile(c(df_left$prob_abs, df_right$prob_abs), 0.99, na.rm=TRUE)
cap_size <- function(x) pmin(x, cap99)

# Plotting (using your previous theme and color scheme)
library(ggplot2); library(patchwork); library(grid)

mk_panel <- function(df, title_text){
  ggplot(df, aes(x=group, y=lr_pair, fill=prob_log10, size=cap_size(prob_abs))) +
    geom_point(shape=21, stroke=0.25) +
    scale_fill_gradient2(low="#053061", mid="white", high="#8B0000",
                         limits=c(cmin,cmax), midpoint=cmid,
                         oob=scales::squish, na.value="white", name="log10(prob)") +
    scale_size_continuous(name="Absolute prob", limits=c(0, cap99), range=c(3,8),
                          breaks=c(0, cap99/2, cap99), labels=c("low","mid","high")) +
    labs(x=NULL, y=NULL, title=title_text) +
    theme_minimal(base_size=10) +
    theme(axis.text.x = element_text(hjust=0.5),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(linewidth=0.25, colour="grey90"),
          panel.spacing.x = unit(0.1,"lines"),
          plot.title = element_text(face="bold", size=11),
          panel.border = element_rect(color="black", fill=NA, linewidth=0.6))
}

p_left  <- mk_panel(df_left,  sprintf("Macrophage → Endothelial — ACD < LPS  [n=%d]", nlevels(df_left$lr_pair)))
p_right <- mk_panel(df_right, sprintf("Macrophage → Endothelial — ACD ≥ LPS [n=%d]", nlevels(df_right$lr_pair)))

H_in <- 0.16 * max(nlevels(df_left$lr_pair), nlevels(df_right$lr_pair)) + 3.0
outdir <- file.path(getwd(), "export_LRdot_ACD_vs_LPS_thresh005")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
pdf_file <- file.path(outdir, "Mac_to_Endo_LRdot_ACD_vs_LPS_thresh005.pdf")
svg_file <- file.path(outdir, "Mac_to_Endo_LRdot_ACD_vs_LPS_thresh005.svg")

while (dev.cur() > 1) try(dev.off(), silent = TRUE)
grDevices::pdf(pdf_file, width = 9.0, height = H_in, onefile = FALSE, family = "sans")
print(p_left | p_right)
dev.off()
svglite::svglite(svg_file, width = 9.0, height = H_in)
print(p_left | p_right)
dev.off()
message("Exported:\n", pdf_file, "\n", svg_file)

# ====== Export Excel data ======
suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(tidyr); library(stringr)
  library(ggplot2); library(patchwork); library(svglite); library(openxlsx)
})

# ====== Your exported Excel (replace with your own path) ======
excel_in <- file.path(getwd(), "export_LRdot_ACD_vs_LPS_thresh005",
                      "Mac_to_Endo_LRdot_ACD_vs_LPS_thresh005.xlsx")

# ====== Output directory / file ======
outdir <- file.path(getwd(), "export_LRdot_ACD_vs_LPS_thresh005", "module_version")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

pdf_file <- file.path(outdir, "Mac_to_Endo_LRdot_moduleBlocks.pdf")
svg_file <- file.path(outdir, "Mac_to_Endo_LRdot_moduleBlocks.svg")
excel_out <- file.path(outdir, "Mac_to_Endo_LRdot_withModules.xlsx")

# ====== Read left/right panel data (your Excel already contains plot_LEFT_only / plot_RIGHT_only sheets) ======
left_w  <- readxl::read_excel(excel_in, sheet = "plot_LEFT_only")
right_w <- readxl::read_excel(excel_in, sheet = "plot_RIGHT_only")

# ====== Automatic module labeling function (can be fine-tuned for your project) ======
assign_module <- function(lr_pair, pathway_name, ligand, receptor){
  x <- paste(lr_pair, pathway_name, ligand, receptor, sep=" | ") %>% tolower()
  
  # ECM / Integrin (collagen/laminin/fibronectin/thrombospondin/integrin)
  if(grepl("\\bcol\\d|\\blam[a-z0-9]*\\b|\\bfn1\\b|\\bspp1\\b|\\bthbs\\d|itga\\d|itgb\\d|itgal|itgav|itgax|itgb2", x)) {
    return("ECM–Integrin")
  }
  # Chemokine
  if(grepl("\\bccl\\d|\\bcxcl\\d|\\bccr\\d|\\bcxcr\\d|\\backr\\d", x)) {
    return("Chemokine")
  }
  # Cytokine / Inflammation (IL/TNF/OSM etc.)
  if(grepl("\\bil\\d|\\btnf\\b|tnfr|il6r|il6st|osmr|osm\\b|ifng|tgfb|tgfbr", x)) {
    return("Cytokine/Inflam")
  }
  # Adhesion / Junction (ICAM/VCAM/SELE/SELL/PECAM1…)
  if(grepl("\\bicam\\d|\\bvcam\\d|\\bsele\\b|\\bsell\\b|\\bselp\\b|\\bpecam1\\b|\\bcdh5\\b|\\bcd34\\b", x)) {
    return("Adhesion/Junction")
  }
  # Angio / Growth factor (VEGF/ANGPT/IGF/KITL/NAMPT etc.)
  if(grepl("\\bvegfa\\b|\\bflt\\d|\\bkdr\\b|\\bangpt\\d|\\bangptl\\d|\\bigf\\d|\\bkitl\\b|\\bkit\\b|\\bnampt\\b|\\binsr\\b", x)) {
    return("Angio/Growth")
  }
  # Guidance / SEMA
  if(grepl("\\bsema\\d|\\bnrp\\d|\\bplxn", x)) {
    return("Guidance/SEMA")
  }
  "Other"
}

# ====== Wide -> long format (for dotplot) ======
wide_to_long <- function(wide_df, side_label){
  wide_df %>%
    mutate(side = side_label) %>%
    mutate(
      module = mapply(assign_module, lr_pair, pathway_name, ligand, receptor),
      module = factor(module,
                      levels = c("Chemokine","Cytokine/Inflam","Adhesion/Junction",
                                 "Angio/Growth","Guidance/SEMA","ECM–Integrin","Other"))
    ) %>%
    select(lr_pair, ligand, receptor, pathway_name, side, module,
           starts_with("prob_"), starts_with("pval_"), plot_rank, delta_ACDvsLPS, delta_abs) %>%
    pivot_longer(cols = starts_with("prob_"),
                 names_to = "group", values_to = "prob") %>%
    mutate(group = gsub("^prob_","",group),
           group = factor(group, levels = c("Ctrl","LPS","ACD"))) %>%
    mutate(prob = as.numeric(prob),
           prob_eps = pmax(prob, 1e-12),
           prob_log10 = log10(prob_eps))
}

df_left  <- wide_to_long(left_w,  "ACD < LPS")
df_right <- wide_to_long(right_w, "ACD ≥ LPS")

df_all <- bind_rows(df_left, df_right)

# ====== Unified color scale / size (same as before: prob_log10 + absolute prob) ======
cmin <- quantile(df_all$prob_log10, 0.05, na.rm = TRUE)
cmax <- quantile(df_all$prob_log10, 0.95, na.rm = TRUE)
cmid <- median(df_all$prob_log10, na.rm = TRUE)
cap99 <- quantile(df_all$prob, 0.99, na.rm = TRUE)
cap_size <- function(x) pmin(x, cap99)

# ====== Each panel: faceted by module, independent y-axis per module ======
mk_panel_module <- function(df, title_text){
  # Order lr_pair within each module by plot_rank (consistent with original plot), with independent y-axis per module
  df <- df %>%
    group_by(module, side) %>%
    mutate(lr_pair = factor(lr_pair, levels = unique(lr_pair[order(plot_rank)]))) %>%
    ungroup()
  
  ggplot(df, aes(x = group, y = lr_pair, fill = prob_log10, size = cap_size(prob))) +
    geom_point(shape = 21, stroke = 0.25) +
    facet_wrap(~ module, ncol = 1, scales = "free_y") +
    scale_fill_gradient2(low = "#053061", mid = "white", high = "#8B0000",
                         limits = c(cmin, cmax), midpoint = cmid,
                         oob = scales::squish, na.value = "white",
                         name = "log10(prob)") +
    scale_size_continuous(name = "Absolute prob",
                          limits = c(0, cap99), range = c(3, 8),
                          breaks = c(0, cap99/2, cap99),
                          labels = c("low","mid","high")) +
    labs(x = NULL, y = NULL, title = title_text) +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.x = element_text(hjust = 0.5),
      strip.text = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.25, colour = "grey90"),
      plot.title = element_text(face = "bold", size = 11),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
      strip.background = element_rect(fill = "grey95", color = "grey80")
    )
}

pL <- mk_panel_module(df_left,  paste0("Macrophage → Endothelial — ACD < LPS  [n=",
                                       length(unique(df_left$lr_pair)), "]"))
pR <- mk_panel_module(df_right, paste0("Macrophage → Endothelial — ACD ≥ LPS [n=",
                                       length(unique(df_right$lr_pair)), "]"))

# ====== Export figure (height varies with number of LR pairs + extra margin after faceting) ======
n_lr_max <- max(length(unique(df_left$lr_pair)), length(unique(df_right$lr_pair)))
H_in <- 0.20 * n_lr_max + 6.0   # you can adjust further

while (dev.cur() > 1) try(dev.off(), silent = TRUE)
grDevices::pdf(pdf_file, width = 10.5, height = H_in, onefile = FALSE, family = "sans")
print(pL | pR)
dev.off()

svglite::svglite(svg_file, width = 10.5, height = H_in)
print(pL | pR)
dev.off()

message("Exported modular figure:\n", pdf_file, "\n", svg_file)

# ====== Also export: table with module info, for manual adjustment (Excel) ======
# 1) Summarize left/right panel "one-row-per-LR" + module
lr_with_module <- bind_rows(left_w, right_w) %>%
  mutate(module = mapply(assign_module, lr_pair, pathway_name, ligand, receptor)) %>%
  arrange(panel, plot_rank)

# 2) Generate an editable mapping (useful when you want to force certain LR pairs into specific modules)
module_map <- lr_with_module %>% select(lr_pair, module, pathway_name, ligand, receptor, panel, plot_rank)

wb <- createWorkbook()
addWorksheet(wb, "LR_with_module")
writeDataTable(wb, "LR_with_module", lr_with_module)

addWorksheet(wb, "module_map_editable")
writeDataTable(wb, "module_map_editable", module_map)

saveWorkbook(wb, excel_out, overwrite = TRUE)
message("Exported Excel with module info:\n", excel_out)

# ================== (Modified version) LR dotplot: faceted by 6 modules (left/right panels unchanged) ==================
# Directly replace your original "mk_panel / p_left / p_right / export" section
# Assumes you have already obtained df_left and df_right (the two long tables generated by prep_plot_df)
# df_left/df_right at least contain: group, lr_pair, ligand, receptor, pathway_name, prob, pval, prob_abs, prob_log10

suppressPackageStartupMessages({
  library(dplyr); library(stringr); library(ggplot2); library(patchwork); library(svglite)
})

# ---------- 1) 6-module labeling ----------
assign_module6 <- function(lr_pair, pathway_name, ligand, receptor){
  x <- paste(lr_pair, pathway_name, ligand, receptor, sep=" | ") %>% tolower()
  
  # Module 1: Inflammation & chemokine / cytokine
  if(grepl("\\bccl\\d|\\bcxcl\\d|\\bccr\\d|\\bcxcr\\d|\\backr\\d|\\bil\\d|\\btnf\\b|tnfr|il6r|il6st|twea?k|tnfsf", x)){
    return("M1 Inflammation/Chemokine")
  }
  
  # Module 2: Adhesion / immune contact
  if(grepl("\\bicam\\d|\\bvcam\\d|\\bsele\\b|\\bsell\\b|\\bselp\\b|\\bpecam1\\b|\\bcd34\\b|\\bpodxl\\b|\\bcd47\\b|\\bcd80\\b|\\bcd274\\b|\\bcdh5\\b", x)){
    return("M2 Adhesion/Contact")
  }
  
  # Module 4: Angio / growth factor
  if(grepl("\\bvegfa\\b|\\bflt\\d|\\bkdr\\b|\\bangpt\\d|\\bangptl\\d|\\bigf\\d|\\bkitl\\b|\\bkit\\b|\\bnampt\\b|\\binsr\\b|\\bosm\\b|\\bosmr\\b", x)){
    return("M4 Angio/Growth")
  }
  
  # Module 5: Guidance / SEMA–NRP
  if(grepl("\\bsema\\d|\\bnrp\\d|\\bplxn", x)){
    return("M5 Guidance/SEMA–NRP")
  }
  
  # Module 3: ECM–Integrin structural (placed later to avoid capturing some adhesion pairs prematurely)
  if(grepl("\\bcol\\d|\\blam[a-z0-9]*\\b|\\bfn1\\b|\\bspp1\\b|\\bthbs\\d|itga\\d|itgb\\d|itgal|itgav|itgax|itgb2", x)){
    return("M3 ECM–Integrin")
  }
  
  # Module 6: Homeostasis / other support
  return("M6 Homeostasis/Other")
}

add_module6 <- function(df){
  df %>%
    mutate(
      module6 = mapply(assign_module6, lr_pair, pathway_name, ligand, receptor),
      module6 = factor(module6, levels = c(
        "M1 Inflammation/Chemokine",
        "M2 Adhesion/Contact",
        "M3 ECM–Integrin",
        "M4 Angio/Growth",
        "M5 Guidance/SEMA–NRP",
        "M6 Homeostasis/Other"
      ))
    )
}

df_left  <- add_module6(df_left)
df_right <- add_module6(df_right)

# ---------- 2) Order within module: preserve your original sorting by delta_abs ----------
# Originally you obtained the order via wide2$delta_abs; here we reuse the factor levels of lr_pair as global ordering reference
reorder_within_module <- function(df){
  # Take the position of each lr_pair in the original factor levels as a "global rank"
  global_rank <- setNames(seq_along(levels(df$lr_pair)), levels(df$lr_pair))
  df$rank_global <- global_rank[as.character(df$lr_pair)]
  
  df <- df %>%
    arrange(module6, rank_global) %>%
    group_by(module6) %>%
    mutate(lr_pair2 = factor(as.character(lr_pair),
                             levels = unique(as.character(lr_pair)))) %>%
    ungroup()
  df
}
df_left  <- reorder_within_module(df_left)
df_right <- reorder_within_module(df_right)

# ---------- 3) Unified color scale / point size (using your previous scaling) ----------
cmin <- quantile(c(df_left$prob_log10, df_right$prob_log10), 0.05, na.rm=TRUE)
cmax <- quantile(c(df_left$prob_log10, df_right$prob_log10), 0.95, na.rm=TRUE)
cmid <- median(c(df_left$prob_log10, df_right$prob_log10), na.rm=TRUE)
cap99 <- quantile(c(df_left$prob_abs, df_right$prob_abs), 0.99, na.rm=TRUE)
cap_size <- function(x) pmin(x, cap99)

# ---------- 4) Faceted plotting: facet_grid(module6 ~ .), left/right panels remain ----------
mk_panel_module6 <- function(df, title_text){
  ggplot(df, aes(x=group, y=lr_pair2, fill=prob_log10, size=cap_size(prob_abs))) +
    geom_point(shape=21, stroke=0.25) +
    facet_grid(
      rows = vars(module6),
      scales = "free_y", space = "free_y", switch = "y",
      labeller = label_wrap_gen(width = 18)
    ) +
    # Core row spacing compression
    scale_y_discrete(expand = expansion(add = 0.6)) +
    scale_fill_gradient2(low="#053061", mid="white", high="#8B0000",
                         limits=c(cmin,cmax), midpoint=cmid,
                         oob=scales::squish, na.value="white",
                         name="log10(prob)") +
    # Larger points but without expanding row height
    scale_size_continuous(
      name="Absolute prob",
      limits=c(0, cap99),
      range=c(3, 8),
      breaks=c(0, cap99/2, cap99),
      labels=c("low","mid","high")
    ) +
    labs(x=NULL, y=NULL, title=title_text) +
    coord_cartesian(clip = "off") +
    theme_minimal(base_size=10) +
    theme(
      axis.text.x = element_text(hjust=0.5),
      axis.text.y = element_text(size=10),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth=0.1, colour="grey90"),
      plot.title = element_text(face="bold", size=11),
      panel.border = element_rect(color="black", fill=NA, linewidth=0.25),
      
      strip.placement = "outside",
      strip.text.y.left = element_text(face="bold", size=10, hjust=0),
      strip.background = element_rect(fill="grey95", color="grey80", linewidth=0.3),
      
      # Extreme compression between modules (but still readable)
      panel.spacing.y = unit(0.04, "lines")
    )
}

p_left  <- mk_panel_module6(df_left,  sprintf("Macrophage → Endothelial — ACD < LPS  [n=%d]",  nlevels(df_left$lr_pair2)))
p_right <- mk_panel_module6(df_right, sprintf("Macrophage → Endothelial — ACD ≥ LPS [n=%d]", nlevels(df_right$lr_pair2)))

# ---------- 5) Export (height based on LR count + extra margin for faceting) ----------
H_in <- 0.16 * max(nlevels(df_left$lr_pair2), nlevels(df_right$lr_pair2)) + 3.0

outdir <- file.path(getwd(), "export_LRdot_ACD_vs_LPS_thresh005")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
pdf_file <- file.path(outdir, "Mac_to_Endo_LRdot_ACD_vs_LPS_module6.pdf")
svg_file <- file.path(outdir, "Mac_to_Endo_LRdot_ACD_vs_LPS_module6.svg")

while (dev.cur() > 1) try(dev.off(), silent = TRUE)
grDevices::pdf(pdf_file, width = 10, height = H_in, onefile = FALSE, family = "sans")
print(p_left | p_right)
dev.off()

message("Exported (6-module faceted):\n", pdf_file)