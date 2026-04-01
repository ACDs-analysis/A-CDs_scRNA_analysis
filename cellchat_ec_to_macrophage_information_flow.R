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

# ===== Endothelial → Macrophage: pathway-level information flow (with proportions & absolute values) =====
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2); library(patchwork) })

compare_groups <- c("ACD","LPS")   # groups to compare
source_cell    <- "Endothelial"    # ← direction changed
target_cell    <- "Macrophage"
min_flow_threshold <- 0
top_n <- NA
outfile_pdf <- "Endo_to_Mac_Pathway_Flow_ACD_vs_LPS.pdf"

# Summarize pathway-level information flow from CellChat object
.compute_pathway_flow <- function(obj, group){
  if (is.null(obj)) return(tibble(pathway=character(), flow=numeric()))
  df <- tryCatch(subsetCommunication(obj, sources=source_cell, targets=target_cell), error=function(e) NULL)
  if (is.null(df) || nrow(df)==0) return(tibble(pathway=character(), flow=numeric()))
  path_col <- intersect(c("pathway_name","pathway"), colnames(df))[1]
  df %>% transmute(pathway=.data[[path_col]], flow=prob) %>%
    group_by(pathway) %>% summarise(flow=sum(flow, na.rm=TRUE), .groups="drop") %>%
    mutate(group=group)
}

avail_groups <- intersect(compare_groups, names(cellchat_results))
stopifnot(length(avail_groups) >= 2)

flow_df <- bind_rows(lapply(avail_groups, function(g) .compute_pathway_flow(cellchat_results[[g]], g)))
flow_df <- flow_df %>% group_by(pathway) %>% filter(max(flow, na.rm=TRUE) > min_flow_threshold) %>% ungroup()

# Relative proportions (ACD+LPS=1)
wide <- flow_df %>% pivot_wider(names_from=group, values_from=flow, values_fill=0) %>% select(pathway, all_of(compare_groups))
wide$.sum <- rowSums(as.matrix(wide[, compare_groups, drop=FALSE]), na.rm=TRUE)
rel_long <- wide %>% pivot_longer(all_of(compare_groups), names_to="group", values_to="flow") %>%
  mutate(rel = ifelse(.sum>0, flow/.sum, 0)) %>% select(pathway, group, rel)

# Sort by LPS proportion from low to high (to place APP appropriately)
order_tbl <- rel_long %>% filter(group=="LPS") %>% arrange(rel) %>% select(pathway)

# TopN (optional)
if (!is.na(top_n) && top_n>0 && nrow(order_tbl)>top_n) {
  keep <- order_tbl$pathway[1:top_n]
  order_tbl <- order_tbl %>% slice_head(n=top_n)
  flow_df  <- flow_df  %>% filter(pathway %in% keep)
  rel_long <- rel_long %>% filter(pathway %in% keep)
}

rel_long$pathway <- factor(rel_long$pathway, levels=order_tbl$pathway)
flow_df$pathway  <- factor(flow_df$pathway,  levels=order_tbl$pathway)

comp_cols <- c("ACD"="#ff9896","LPS"="#17becf")

p_rel <- ggplot(rel_long, aes(x=rel, y=pathway, fill=group)) +
  geom_col(position=position_stack(reverse=TRUE), width=0.75, color="black", linewidth=0.15) +
  scale_x_continuous(limits=c(0,1), expand=expansion(mult=c(0,0.02))) +
  scale_fill_manual(values=comp_cols) +
  labs(x="Relative information flow (ACD + LPS = 1, sorted by LPS ↑)", y=NULL) +
  theme_minimal(base_size=10) +
  theme(legend.position="top", axis.text.y=element_text(color="#e67e22"),
        panel.grid.major.y=element_blank())

p_abs <- ggplot(flow_df, aes(x=flow, y=pathway, fill=group)) +
  geom_col(position=position_dodge(width=0.75), width=0.70, color="black", linewidth=0.15) +
  scale_fill_manual(values=comp_cols) +
  labs(x="Information flow (absolute)", y=NULL) +
  theme_minimal(base_size=10) +
  theme(legend.position="none", axis.text.y=element_blank(),
        axis.ticks.y=element_blank(), panel.grid.major.y=element_blank())

p_out <- p_rel + p_abs + plot_layout(widths=c(1.05,0.95))
ggsave(outfile_pdf, plot=p_out, width=7.2, height=8.0, device="pdf", dpi=300)
message("✓ Exported: ", outfile_pdf)
