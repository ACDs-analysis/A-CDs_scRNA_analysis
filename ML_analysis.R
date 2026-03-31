suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(e1071)
  library(randomForest)
  library(glmnet)
  library(ggvenn)
})

set.seed(1234)

# =========================
# 0) Read the three exported tables
# =========================
feat_file  <- "C:/Users/45214/Documents/ML_core_LR_results/export_ML_LR/X_sample_by_TopLR.csv"
stat_file  <- "C:/Users/45214/Documents/ML_core_LR_results/export_ML_LR/TopLR_by_absDelta_ACD_vs_LPS_groupMean.csv"
label_file <- "C:/Users/45214/Documents/ML_core_LR_results/export_ML_LR/SampleInfo.csv"

feat <- readr::read_csv(feat_file, show_col_types = FALSE)
stat <- readr::read_csv(stat_file, show_col_types = FALSE)
lab  <- readr::read_csv(label_file, show_col_types = FALSE)

# If the first table lacks a 'group' column, supplement it from the third table
if(!("group" %in% colnames(feat))) {
  feat <- feat %>% left_join(lab, by = "sample_id")
}
stopifnot(all(c("sample_id","group") %in% colnames(feat)))

# Keep only ACD vs LPS
feat <- feat %>% filter(group %in% c("ACD","LPS"))

# -------------------------
# Standardize feature names (critical)
# -------------------------
normalize_lr <- function(x){
  x <- as.character(x)
  x <- trimws(x)
  x <- gsub("\\.{3,}", " - ", x)           # "..." -> " - "
  x <- gsub("[–—_\\|]", "-", x)            # catch other separators
  x <- gsub("\\s*-\\s*", " - ", x)         # unify to "A - B"
  x <- gsub("\\s+", " ", x)                # collapse multiple spaces
  x <- trimws(x)
  x
}

# Construct X / y
X <- feat %>% select(-sample_id, -group) %>% as.data.frame(check.names = FALSE)
rownames(X) <- feat$sample_id
y <- ifelse(feat$group == "ACD", 1, 0)   # ACD=1, LPS=0

# Standardize column names of X
colnames(X) <- normalize_lr(colnames(X))

message("Group counts:"); print(table(feat$group))
message("n_samples=", nrow(X), " | n_features=", ncol(X))

# =========================
# 0.1) Stratified folds: k must not exceed the number of samples per class
# =========================
min_class_n <- min(table(y))               # number of samples per class
k_fold <- min(10, nrow(X), min_class_n)    # key: k <= min_class_n
if(k_fold < 2) stop("Too few samples for CV (need at least 2 per class)")
message("Using stratified k-fold CV, k_fold=", k_fold)

make_stratified_folds <- function(y, k=5, seed=1234){
  set.seed(seed)
  y <- as.integer(y)
  idx1 <- which(y==1)
  idx0 <- which(y==0)
  folds <- integer(length(y))
  folds[idx1] <- sample(rep(1:k, length.out=length(idx1)))
  folds[idx0] <- sample(rep(1:k, length.out=length(idx0)))
  folds
}

folds <- make_stratified_folds(y, k = k_fold, seed = 1234)

message("Fold sanity check (train-set class counts):")
for(f in 1:k_fold){
  tr <- which(folds != f)
  cat(" fold", f, ":", paste(names(table(y[tr])), table(y[tr]), collapse=" | "), "\n")
}

# =========================
# 1) SVM-RFE + Stratified k-fold CV
# =========================
svm_rank_weights <- function(Xtr, ytr, cost=1) {
  dat <- data.frame(y = factor(ytr), Xtr, check.names = FALSE)
  fit <- e1071::svm(y ~ ., data=dat, kernel="linear", cost=cost, scale=TRUE)
  w <- t(fit$coefs) %*% fit$SV
  w2 <- as.numeric(w)^2
  names(w2) <- colnames(Xtr)
  sort(w2, decreasing=TRUE)
}

svm_rfe_once <- function(Xtr, ytr, step_frac=0.1, cost=1) {
  feats <- colnames(Xtr)
  while(length(feats) > 1) {
    w2 <- tryCatch(
      svm_rank_weights(Xtr[, feats, drop=FALSE], ytr, cost=cost),
      error = function(e) NULL
    )
    if(is.null(w2)) break
    step <- max(1, floor(length(feats) * step_frac))
    drop_feats <- tail(names(w2), step)
    feats <- setdiff(feats, drop_feats)
  }
  out <- tryCatch(names(svm_rank_weights(Xtr, ytr, cost=cost)),
                  error=function(e) colnames(Xtr))
  out
}

svm_rfe_cv <- function(X, y, folds, k, cost=1, step_frac=0.1, seed=1234) {
  set.seed(seed)
  n <- nrow(X)
  p <- ncol(X)
  sizes <- unique(pmax(5, round(seq(5, p, length.out=15))))
  
  acc <- matrix(NA, nrow=length(sizes), ncol=k,
                dimnames=list(paste0("m=",sizes), paste0("fold",1:k)))
  
  for(f in 1:k){
    te <- which(folds==f)
    tr <- setdiff(1:n, te)
    
    Xtr <- X[tr,,drop=FALSE]; ytr <- y[tr]
    Xte <- X[te,,drop=FALSE]; yte <- y[te]
    
    if(length(unique(ytr)) < 2){
      message("Warning: fold ", f, " train has only one class, skip.")
      next
    }
    
    full_rank <- svm_rfe_once(Xtr, ytr, step_frac=step_frac, cost=cost)
    
    for(i in seq_along(sizes)){
      m <- sizes[i]
      feats_m <- head(full_rank, m)
      
      dat_tr <- data.frame(y=factor(ytr), Xtr[,feats_m,drop=FALSE], check.names = FALSE)
      dat_te <- data.frame(Xte[,feats_m,drop=FALSE], check.names = FALSE)
      
      fit <- tryCatch(
        e1071::svm(y ~ ., data=dat_tr, kernel="linear", cost=cost, scale=TRUE),
        error=function(e) NULL
      )
      if(is.null(fit)) next
      
      pred <- tryCatch(predict(fit, dat_te), error=function(e) NULL)
      if(is.null(pred)) next
      
      acc[i,f] <- mean(as.numeric(as.character(pred)) == yte)
    }
  }
  
  mean_acc <- rowMeans(acc, na.rm=TRUE)
  best_m <- sizes[which.max(mean_acc)]
  list(best_m=best_m, sizes=sizes, mean_acc=mean_acc, acc=acc)
}

svm_cv <- svm_rfe_cv(X, y, folds=folds, k=k_fold, cost=1, step_frac=0.1, seed=1234)
best_m <- svm_cv$best_m
message("SVM-RFE best_m=", best_m)

full_rank_all <- svm_rfe_once(X, y, step_frac=0.1, cost=1)
svm_features <- head(full_rank_all, best_m)
svm_features <- normalize_lr(svm_features)
message("Number of SVM candidate LRs = ", length(svm_features))

# =========================
# 2) Random Forest importance (train on SVM-RFE candidates; keep 'top 75% importance')
# =========================
RF_KEEP_PROP <- 0.75
RF_Q <- 1 - RF_KEEP_PROP   # 0.25

set.seed(1234)

# Use SVM-RFE candidate features
X_rf <- X[, svm_features, drop = FALSE]

# Use x/y interface instead of formula to avoid parsing errors from column names with spaces/hyphens
rf_fit <- randomForest::randomForest(
  x = X_rf,
  y = factor(y),
  ntree = 1000,
  importance = TRUE
)

# MeanDecreaseGini
imp <- randomForest::importance(rf_fit, type = 2)

imp_df <- data.frame(
  lr_pair    = rownames(imp),
  importance = as.numeric(imp[, 1]),
  stringsAsFactors = FALSE
) %>%
  arrange(desc(importance))

# Keep the top 75% importance (importance >= 25th percentile)
rf_cut <- stats::quantile(imp_df$importance, RF_Q, na.rm = TRUE)
rf_keep <- imp_df %>% filter(importance >= rf_cut)
rf_features <- rf_keep$lr_pair

message(
  "RF candidate LRs (top ", RF_KEEP_PROP*100, "% importance) = ",
  length(rf_features), " / ", nrow(imp_df)
)

# Bubble plot
p_bubble <- ggplot(
  rf_keep,
  aes(
    x = reorder(lr_pair, importance),
    y = importance,
    size = importance,
    fill = importance
  )
) +
  geom_point(shape = 21, color = "black", stroke = 0.25) +
  coord_flip() +
  labs(
    x = NULL,
    y = "RF importance (MeanDecreaseGini)",
    title = paste0("High-importance LR pairs by Random Forest (Top ", RF_KEEP_PROP*100, "%)")
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 8)
  )+ scale_fill_gradient2(
    low  = "#17BECF",
    mid  = "#F7F7F7",
    high = "#FF7F0E",
    midpoint = 0.5,
    name = "Importance"
  ) +
  scale_size(range = c(2, 8))


print(p_bubble)

# =========================
# 3) Route A: Elastic Net stability selection (Repeated stratified CV)
#   - Output: selection frequency table for each LR
#   - Provide 'stable candidate set': prop >= 0.2 (fallback to top 20)
# =========================
REPEAT_N <- 200
ALPHA_USE <- 0.5           # you can change: 0.3/0.5/0.7
LAMBDA_RULE <- "lambda.min"  # for small samples, min is recommended; 1se often empty
FREQ_PROP_CUT <- 0.20      # frequency threshold (20%)
TOP_FALLBACK <- 20         # if too few meet the threshold, take top 20

set.seed(1234)
feat_counts <- setNames(rep(0L, ncol(X)), normalize_lr(colnames(X)))

for(r in 1:REPEAT_N){
  folds_r <- make_stratified_folds(y, k = k_fold, seed = 1234 + r)
  
  cvfit_r <- cv.glmnet(
    as.matrix(X), y,
    family="binomial",
    alpha=ALPHA_USE,
    nfolds=k_fold,
    foldid=folds_r,
    standardize=TRUE,
    type.measure="deviance"
  )
  
  co <- coef(cvfit_r, s = LAMBDA_RULE)
  sel <- setdiff(rownames(co)[as.numeric(co) != 0], "(Intercept)")
  sel <- normalize_lr(sel)
  
  common <- intersect(names(feat_counts), sel)
  if(length(common) > 0) feat_counts[common] <- feat_counts[common] + 1L
}

freq_df <- data.frame(
  lr_pair = names(feat_counts),
  freq    = as.integer(feat_counts),
  prop    = as.numeric(feat_counts) / REPEAT_N,
  stringsAsFactors = FALSE
) %>% arrange(desc(freq), desc(prop))

enet_stable <- freq_df %>% filter(prop >= FREQ_PROP_CUT) %>% pull(lr_pair)
if(length(enet_stable) < 5){
  enet_stable <- head(freq_df$lr_pair, TOP_FALLBACK)
  message("Note: prop threshold too strict resulted in <5 candidates; using Top ", TOP_FALLBACK, " as stable candidates.")
}

message("ElasticNet stable selection: REPEAT_N=", REPEAT_N,
        " | alpha=", ALPHA_USE,
        " | lambda rule=", LAMBDA_RULE,
        " | number of stable candidates=", length(enet_stable))

# Visualize: Top20 selection frequencies
topN_plot <- min(20, nrow(freq_df))
p_enet_freq <- ggplot(
  freq_df %>% slice_head(n = topN_plot),
  aes(x = reorder(lr_pair, prop), y = prop, size = prop, fill = prop)
) +
  geom_point(shape = 21, color="black", stroke=0.25) +
  coord_flip() +
  labs(x=NULL, y=paste0("Selection proportion (", REPEAT_N, " repeats)"),
       title=paste0("ElasticNet stability selection (alpha=", ALPHA_USE, ", ", LAMBDA_RULE, ")")) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 8))+ scale_fill_gradient2(
          low  = "#17BECF",
          mid  = "#F7F7F7",
          high = "#FF7F0E",
          midpoint = 0.5,
          name = "Selection\nproportion"
        ) +
  scale_size(range = c(2, 8))


print(p_enet_freq)

# =========================
# 4) Venn intersection (all normalized first)
#   - SVM-RFE: svm_features
#   - RF:      rf_features
#   - ENet stable: enet_stable
# =========================
svm_set  <- unique(normalize_lr(svm_features))
rf_set   <- unique(normalize_lr(rf_features))
enet_set <- unique(normalize_lr(enet_stable))

venn_list <- list(
  `SVM-RFE`      = svm_set,
  `RF`           = rf_set,
  `ElasticNet-Stable` = enet_set
)

message("Sanity check intersections (after normalization):")
cat("SVM ∩ RF size           =", length(intersect(svm_set, rf_set)), "\n")
cat("SVM ∩ ENetStable size   =", length(intersect(svm_set, enet_set)), "\n")
cat("RF  ∩ ENetStable size   =", length(intersect(rf_set, enet_set)), "\n")
cat("SVM ∩ RF ∩ ENetStable size =", length(Reduce(intersect, venn_list)), "\n")

p_venn <- ggvenn::ggvenn(venn_list, show_percentage = FALSE) +
  ggtitle("Intersection of candidate LR pairs (normalized)")
print(p_venn)

core_lr <- Reduce(intersect, venn_list)

# If triple intersection is empty, fallback to SVM∩ENetStable (more stringent);
# if still empty, fallback to SVM∩RF (most common)
if(length(core_lr) == 0){
  message("Note: triple intersection is empty; using SVM-RFE ∩ ENetStable as core candidates.")
  core_lr <- intersect(svm_set, enet_set)
  if(length(core_lr) == 0){
    message("Note: SVM-RFE ∩ ENetStable still empty; using SVM-RFE ∩ RF as core candidates.")
    core_lr <- intersect(svm_set, rf_set)
  }
}

message("Final number of core LRs = ", length(core_lr))
print(core_lr)

# =========================
# 5) Export results
# =========================
outdir <- file.path(getwd(), "ML_core_LR_results")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Result tables
write.csv(imp_df,   file.path(outdir, "RF_importance_all.csv"), row.names = FALSE)
write.csv(freq_df,  file.path(outdir, "ElasticNet_stability_frequency.csv"), row.names = FALSE)

write.csv(data.frame(lr_pair = svm_set),  file.path(outdir, "SVM_RFE_candidates.csv"), row.names = FALSE)
write.csv(data.frame(lr_pair = rf_set),   file.path(outdir, "RF_candidates_all.csv"), row.names = FALSE)
write.csv(data.frame(lr_pair = enet_set), file.path(outdir, "ElasticNet_stable_candidates.csv"), row.names = FALSE)
write.csv(data.frame(lr_pair = core_lr),  file.path(outdir, "Core_LR_intersection.csv"), row.names = FALSE)

# Figures
ggsave(file.path(outdir, "RF_importance_bubble.pdf"), p_bubble, width = 7, height = 6)

ggsave(file.path(outdir, "ElasticNet_stability_Top20.pdf"), p_enet_freq, width = 7, height = 6)

ggsave(file.path(outdir, "Venn_coreLR.pdf"), p_venn, width = 6, height = 4)

message("Export completed to directory: ", outdir)

