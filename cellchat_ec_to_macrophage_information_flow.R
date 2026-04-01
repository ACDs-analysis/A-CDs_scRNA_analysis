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
data_dir <- "E:/科研/单细胞测序/gene"  # ←改成你的路径
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

# ================== 2) 读取并合并「内皮细胞CSV」 ==================
library(readr)

endo_dir <- "E:/科研/单细胞测序/EC gene" # ←改成你的路径
endo_files <- list.files(endo_dir, pattern="\\.csv$", full.names=TRUE)
stopifnot("未找到内皮细胞csv文件" = length(endo_files) > 0)

read_expr_csv <- function(path){
  # 假设：第1列 Gene ID、第2列 Gene Symbol、第3列 Type、第4列 total，第5列开始为细胞表达
  df <- suppressMessages(readr::read_csv(path, show_col_types = FALSE))
  stopifnot(ncol(df) >= 5)
  gene_id     <- as.character(df[[1]])
  gene_symbol <- as.character(df[[2]])
  gene_symbol[is.na(gene_symbol) | gene_symbol==""] <- gene_id
  genes <- make.unique(ifelse(is.na(gene_symbol)|gene_symbol=="","NA_gene",gene_symbol))
  expr_part <- df[, 5:ncol(df), drop=FALSE]
  
  # 强制转数值
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
  # 例：ACD 1 cluster 6.csv / LPS 2 cluster 24.csv 等
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

# ===== 内皮细胞矩阵并集对齐 + 合并 =====
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

# ================== 3) 与巨噬细胞矩阵拼接（按基因并集对齐） ==================
# 说明：此处使用你上段代码中得到的 expr_all（巨噬细胞）与 cell_meta
# 增加 cell_type 列
cell_meta$cell_type <- "Macrophage"
endo_meta$cell_type <- "Endothelial"

# 基因并集对齐
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






# 基本一致性检查
stopifnot(ncol(expr_all_merged) == nrow(cell_meta_merged))
stopifnot(identical(colnames(expr_all_merged), rownames(cell_meta_merged)))

# 统一条件/簇的因子水平（可选）
cond_levels <- c("Ctrl","LPS","ACD")
cell_meta_merged$group <- factor(as.character(cell_meta_merged$group), levels = cond_levels)

mac_levels <- as.character(0:11)
cell_meta_merged$cluster_original[cell_meta_merged$cell_type=="Macrophage"] <-
  factor(as.character(cell_meta_merged$cluster_original[cell_meta_merged$cell_type=="Macrophage"]),
         levels = mac_levels)

# 可选：去掉全零基因与空细胞（通常不必）
# keep_gene <- Matrix::rowSums(expr_all_merged) > 0
# keep_cell <- Matrix::colSums(expr_all_merged) > 0
# expr_all_merged  <- expr_all_merged[keep_gene, keep_cell, drop=FALSE]
# cell_meta_merged <- cell_meta_merged[keep_cell,, drop=FALSE]

diagnose_data_issues <- function() {
  message("=== 数据诊断报告 ===")
  
  # 检查分组信息
  message("1. 分组分布:")
  print(table(cell_meta_merged$group, useNA = "always"))
  
  # 检查细胞类型
  message("2. 细胞类型分布:")
  print(table(cell_meta_merged$cell_type))
  
  # 检查表达矩阵
  message("3. 表达矩阵维度:")
  print(dim(expr_all_merged))
  
  # 检查基因名称
  message("4. 基因名称示例:")
  print(head(rownames(expr_all_merged), 10))
  
  # 检查是否有零表达基因
  zero_genes <- Matrix::rowSums(expr_all_merged) == 0
  message("5. 零表达基因数量:", sum(zero_genes))
  
  # 检查样本文件名
  message("6. 样本文件名示例:")
  print(head(unique(cell_meta_merged$sample_file)))
}

# 运行诊断
diagnose_data_issues()

# 改进的分组提取函数
extract_group_properly <- function(sample_file) {
  message("正在处理样本: ", sample_file) # 调试信息
  
  # 方法1: 按空格分割取第一部分
  parts <- strsplit(trimws(sample_file), "\\s+")[[1]]
  if(length(parts) > 0) {
    group_candidate <- parts[1]
    # 标准化分组名称
    if(toupper(group_candidate) %in% c("CTRL", "CONTROL", "CONTROL1", "CONTROL2")) return("Ctrl")
    if(toupper(group_candidate) %in% c("LPS")) return("LPS")
    if(toupper(group_candidate) %in% c("ACD")) return("ACD")
  }
  
  # 方法2: 直接字符串匹配
  sample_lower <- tolower(trimws(sample_file))
  if(grepl("^acd", sample_lower)) return("ACD")
  if(grepl("^lps", sample_lower)) return("LPS") 
  if(grepl("^ctrl|^control", sample_lower)) return("Ctrl")
  
  return(NA)
}

# 测试分组提取
test_group_extraction <- function() {
  test_samples <- c("ACD 1 cluster 0", "LPS 2 cluster 1", "Ctrl 1 cluster 2", "CONTROL 1 cluster 3")
  for(sample in test_samples) {
    group <- extract_group_properly(sample)
    message("样本: '", sample, "' -> 分组: '", group, "'")
  }
}

# 应用修复的分组提取
repair_group_information <- function() {
  message("开始修复分组信息...")
  
  # 备份原始sample_file
  cell_meta_merged$sample_file_original <- cell_meta_merged$sample_file
  
  # 应用新的分组提取
  cell_meta_merged$group <- sapply(cell_meta_merged$sample_file, extract_group_properly)
  
  # 检查结果
  message("修复后的分组分布:")
  print(table(cell_meta_merged$group, useNA = "always"))
  
  # 如果还有NA，显示这些样本
  na_samples <- cell_meta_merged$sample_file[is.na(cell_meta_merged$group)]
  if(length(na_samples) > 0) {
    message("以下样本的分组提取失败:")
    print(unique(na_samples))
  }
  
  return(cell_meta_merged)
}

# 查看所有唯一的样本文件名以调试
debug_sample_names <- function() {
  message("所有唯一的样本文件名:")
  all_samples <- unique(cell_meta_merged$sample_file)
  for(i in seq_along(all_samples)) {
    message(sprintf("%3d: '%s'", i, all_samples[i]))
  }
  
  # 按可能的模式分类
  message("\n按开头字母分类:")
  sample_starts <- substr(all_samples, 1, 3)
  print(table(sample_starts))
}

# 步骤1: 查看所有样本文件名模式
debug_sample_names()

# 步骤2: 测试分组提取
test_group_extraction()

# 步骤3: 修复分组信息
cell_meta_merged <- repair_group_information()

# 步骤4: 重新检查分组
message("最终分组分布:")
final_groups <- table(cell_meta_merged$group)
print(final_groups)

# 步骤5: 检查每组细胞数量
if(sum(!is.na(cell_meta_merged$group)) > 0) {
  message("\n各分组细胞数量:")
  for(group_name in names(final_groups)) {
    cell_count <- sum(cell_meta_merged$group == group_name, na.rm = TRUE)
    message(group_name, ": ", cell_count, " 个细胞")
  }
}

# 设置分组为因子
cell_meta_merged$group <- factor(cell_meta_merged$group, levels = c("Ctrl", "LPS", "ACD"))

# 检查细胞类型在各组的分布
message("各分组细胞类型详细分布:")
cell_type_by_group <- table(cell_meta_merged$cell_type, cell_meta_merged$group)
print(cell_type_by_group)

# 提取各组的细胞
cells_ctrl <- rownames(subset(cell_meta_merged, group == "Ctrl"))
cells_lps  <- rownames(subset(cell_meta_merged, group == "LPS")) 
cells_acd  <- rownames(subset(cell_meta_merged, group == "ACD"))

message("分组细胞数量确认:")
message("Ctrl: ", length(cells_ctrl), " 细胞")
message("LPS: ", length(cells_lps), " 细胞") 
message("ACD: ", length(cells_acd), " 细胞")





# 提取各组的表达矩阵
expr_ctrl <- expr_all_merged[, cells_ctrl, drop = FALSE]
meta_ctrl <- cell_meta_merged[cells_ctrl, , drop = FALSE]

expr_lps <- expr_all_merged[, cells_lps, drop = FALSE]  
meta_lps <- cell_meta_merged[cells_lps, , drop = FALSE]

expr_acd <- expr_all_merged[, cells_acd, drop = FALSE]
meta_acd <- cell_meta_merged[cells_acd, , drop = FALSE]

# 清理基因名称（重要！）
clean_gene_names <- function(expr_mat) {
  genes <- rownames(expr_mat)
  # 移除单引号
  cleaned <- gsub("'", "", genes)
  # 移除其他特殊字符但保留正常字符
  cleaned <- gsub("[^a-zA-Z0-9._-]", "", cleaned)
  # 处理空值
  cleaned[cleaned == "" | is.na(cleaned)] <- "UnknownGene"
  rownames(expr_mat) <- make.unique(cleaned)
  return(expr_mat)
}

# 应用基因名称清理
expr_ctrl <- clean_gene_names(expr_ctrl)
expr_lps <- clean_gene_names(expr_lps) 
expr_acd <- clean_gene_names(expr_acd)

# 检查清理后的基因名称
message("清理后的基因名称示例:")
print(head(rownames(expr_ctrl), 10))





# 改进的CellChat分析函数
run_cellchat_robust <- function(expr_mat, meta_df, group_name, species = "mouse") {
  message("\n====== 开始分析 ", group_name, " 组 ======")
  
  # 基本检查
  if(ncol(expr_mat) < 10) {
    message("警告: ", group_name, "组细胞数不足，跳过分析")
    return(NULL)
  }
  
  message("数据维度: ", paste(dim(expr_mat), collapse = " x "))
  message("细胞类型: ", paste(unique(meta_df$cell_type), collapse = ", "))
  
  # 设置数据库
  if(species == "mouse") {
    CellChatDB <- CellChatDB.mouse
    PPI_use <- PPI.mouse
  } else if(species == "human") {
    CellChatDB <- CellChatDB.human
    PPI_use <- PPI.human
  } else {
    stop("请指定物种: 'mouse' 或 'human'")
  }
  
  # 检查数据库基因重叠
  db_genes <- unique(c(CellChatDB$interaction$ligand, CellChatDB$interaction$receptor))
  our_genes <- rownames(expr_mat)
  overlapping <- intersect(db_genes, our_genes)
  message("与CellChat数据库重叠的基因数: ", length(overlapping))
  
  if(length(overlapping) < 50) {
    message("警告: 重叠基因较少，可能影响分析质量")
    message("重叠基因示例: ", paste(head(overlapping, 10), collapse = ", "))
  }
  
  tryCatch({
    # 创建CellChat对象
    cellchat <- createCellChat(
      object = as.matrix(expr_mat),
      meta = meta_df,
      group.by = "cell_type"
    )
    
    cellchat@DB <- CellChatDB
    
    # 预处理流程
    cellchat <- subsetData(cellchat)
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)
    cellchat <- projectData(cellchat, PPI_use)
    
    # 计算通讯概率
    cellchat <- computeCommunProb(
      cellchat,
      type = "truncatedMean",
      trim = 0.1,
      raw.use = FALSE
    )
    
    # 过滤和聚合
    cellchat <- filterCommunication(cellchat, min.cells = 5)
    cellchat <- computeCommunProbPathway(cellchat)
    cellchat <- aggregateNet(cellchat)
    
    message("✓ ", group_name, " 组分析完成")
    return(cellchat)
    
  }, error = function(e) {
    message("× ", group_name, " 组分析失败: ", e$message)
    return(NULL)
  })
}
library(CellChat); data(CellChatDB.mouse); data(PPI.mouse)

# 运行各组分析
message("\n开始CellChat分析...")
cellchat_results <- list()

cellchat_results$Ctrl <- run_cellchat_robust(expr_ctrl, meta_ctrl, "Ctrl", species = "mouse")
cellchat_results$LPS <- run_cellchat_robust(expr_lps, meta_lps, "LPS", species = "mouse")
cellchat_results$ACD <- run_cellchat_robust(expr_acd, meta_acd, "ACD", species = "mouse")

# 移除失败的分析
cellchat_results <- cellchat_results[!sapply(cellchat_results, is.null)]

if(length(cellchat_results) == 0) {
  stop("所有组的CellChat分析都失败了，请检查数据和参数")
}

message("\n成功完成分析的组数: ", length(cellchat_results))



# ===== Endothelial → Macrophage：通路级信息流（含占比 & 绝对值）=====
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2); library(patchwork) })

compare_groups <- c("ACD","LPS")   # 对比组
source_cell    <- "Endothelial"    # ← 改方向
target_cell    <- "Macrophage"
min_flow_threshold <- 0
top_n <- NA
outfile_pdf <- "Endo_to_Mac_Pathway_Flow_ACD_vs_LPS.pdf"

# 从 CellChat 对象汇总通路级信息流
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

# 相对占比（ACD+LPS=1）
wide <- flow_df %>% pivot_wider(names_from=group, values_from=flow, values_fill=0) %>% select(pathway, all_of(compare_groups))
wide$.sum <- rowSums(as.matrix(wide[, compare_groups, drop=FALSE]), na.rm=TRUE)
rel_long <- wide %>% pivot_longer(all_of(compare_groups), names_to="group", values_to="flow") %>%
  mutate(rel = ifelse(.sum>0, flow/.sum, 0)) %>% select(pathway, group, rel)

# 按 LPS 占比从低到高排序（把 APP 放到相应位置）
order_tbl <- rel_long %>% filter(group=="LPS") %>% arrange(rel) %>% select(pathway)

# TopN（可选）
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
message("✓ 已导出：", outfile_pdf)
