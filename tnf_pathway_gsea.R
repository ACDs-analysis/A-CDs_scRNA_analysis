# ============================================================
# TNF signaling pathway (mmu04668) GSEV plot (4 panels)
# Ranking: -log10(p) * sign(log2FC)
# X ticks fixed: 0,500,1000,1500,2000,2500,3000,3500
# ES panel and grey-bar panel share identical X axis (aligned)
# Zero cross marked in panel 1 & 4
# Output: PDF only
# ============================================================

file <- "C:/Users/45214/Desktop/ACD_LPS 巨噬细胞差异基因.csv"
out_pdf <- "TNF_mmu04668_GSEV_fixedticks_zerocross.pdf"

pathway_title <- "TNF signaling pathway (mmu04668)"

fixed_ticks <- c(0,500,1000,1500,2000,2500,3000,3500)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

# -------------------- load & clean --------------------
df <- read.csv(file, check.names = FALSE)
stopifnot(all(c("Gene","log2FC","pvalue") %in% colnames(df)))

dat <- df %>%
  transmute(
    Gene   = toupper(trimws(gsub("['\"]", "", as.character(Gene)))),
    log2FC = suppressWarnings(as.numeric(log2FC)),
    pvalue = suppressWarnings(as.numeric(pvalue))
  ) %>%
  filter(!is.na(Gene), Gene != "", !is.na(log2FC), !is.na(pvalue))

# ranking metric
dat <- dat %>%
  mutate(score = -log10(pvalue + 1e-300) * sign(log2FC)) %>%
  arrange(desc(score)) %>%
  distinct(Gene, .keep_all = TRUE)

rank_vec <- dat$score
names(rank_vec) <- dat$Gene

# -------- Key: clean rank_vec once (remove NA/Inf) --------
rank_vec <- rank_vec[is.finite(rank_vec)]
N <- length(rank_vec)

# -------------------- TNF gene set --------------------
tnf_genes <- c(
  "TNF","TNFRSF1A","TNFRSF1B","LTB","TNFSF10","FAS","FASLG",
  "TRADD","TRAF2","TRAF5","RIPK1","MAP3K7","TAB1","TAB2","TAB3",
  "CHUK","IKBKB","IKBKG","NFKBIA","NFKBIB","NFKBIE",
  "RELA","RELB","NFKB1","NFKB2",
  "MAPK8","MAPK9","MAPK10","MAPK14","MAPK11","MAPK12","MAPK13",
  "MAP2K3","MAP2K4","MAP2K6","MAP3K1","MAP3K5","MAP3K8",
  "CASP3","CASP7","CASP8","BID","BCL2","BCL2L1","CFLAR",
  "IL1B","IL6","CXCL1","CXCL2","CXCL3","CXCL10","CCL2","CCL5",
  "PTGS2","ICAM1","VCAM1","SELE","NOS2","SOD2",
  "JUN","JUNB","FOS","FOSL1","ATF3","DDIT3","HIF1A"
)

hits <- names(rank_vec) %in% tnf_genes

Nh <- sum(hits); Nm <- N - Nh
stopifnot(Nh > 0, Nm > 0)

# -------------------- running ES --------------------
w <- abs(rank_vec)
Phit  <- cumsum((w / sum(w[hits])) * hits)
Pmiss <- cumsum((!hits) * (1 / Nm))
running_ES <- Phit - Pmiss

peak_idx <- which.max(abs(running_ES))
ES <- running_ES[peak_idx]

# -------------------- zero cross --------------------
zero_cross <- which(diff(sign(rank_vec)) != 0)
zero_idx <- if (length(zero_cross)) zero_cross[1] + 1 else NA_integer_

# -------------------- axis (redefine with final N to avoid misalignment) --------------------
breaks_x <- fixed_ticks[fixed_ticks <= N]
limits_x <- c(0, N)

# -------------------- heat strip --------------------
nbin_heat <- 120L
bin_id <- cut(seq_len(N), breaks = nbin_heat, labels = FALSE)
heat_df <- data.frame(bin = bin_id, score = as.numeric(rank_vec)) %>%
  group_by(bin) %>%
  summarise(mean_score = mean(score), .groups = "drop") %>%
  mutate(
    xmin = (bin - 1) * (N/nbin_heat),
    xmax = bin * (N/nbin_heat)
  )

# -------------------- grey bars (1000 samples) --------------------
nsamp_bar <- min(1000L, N)
idx_bar <- unique(pmin(N, pmax(1, round(seq(1, N, length.out = nsamp_bar)))))  # avoid rounding out of bounds
bar_df <- data.frame(Rank = idx_bar, val = as.numeric(rank_vec[idx_bar])) %>%
  filter(is.finite(val))

# Winsorize (avoid extreme values affecting aesthetics or causing plotting errors)
lims <- quantile(bar_df$val, probs = c(0.005, 0.995), na.rm = TRUE)
bar_df$val <- pmax(pmin(bar_df$val, lims[2]), lims[1])

# -------------------- plots --------------------
p1 <- ggplot(data.frame(Rank = 1:N, ES = running_ES), aes(Rank, ES)) +
  geom_line(color = "#00C200", linewidth = 1.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "#808080") +
  geom_vline(xintercept = peak_idx, linetype = "dashed", color = "#808080") +
  { if (!is.na(zero_idx)) geom_vline(xintercept = zero_idx, linetype = "dashed", color = "#404040") } +
  { if (!is.na(zero_idx)) annotate("text", x = zero_idx, y = 0.05,
                                   label = paste0("Zero cross at ", zero_idx),
                                   angle = 90, vjust = -0.2, hjust = 1, size = 3) } +
  scale_x_continuous(limits = limits_x, breaks = breaks_x, expand = c(0,0)) +
  labs(title = "Enrichment plot: mmu04668", y = "Enrichment score (ES)", x = NULL) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5))

p2 <- ggplot(data.frame(Rank = which(hits)), aes(Rank, 0)) +
  geom_segment(aes(x = Rank, xend = Rank, y = 0, yend = 1), color = "black") +
  scale_x_continuous(limits = limits_x, breaks = breaks_x, expand = c(0,0)) +
  theme_void()

p3 <- ggplot(heat_df, aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = 1, fill = mean_score)) +
  geom_rect() +
  scale_fill_gradient2(low = "#2C7BB6", mid = "white", high = "#D7191C", midpoint = 0) +
  scale_x_continuous(limits = limits_x, breaks = breaks_x, expand = c(0,0)) +
  theme_void() +
  theme(legend.position = "right")

p4 <- ggplot(bar_df, aes(Rank, val)) +
  geom_col(width = (N/nsamp_bar)*0.95, fill = "#B3B3B3") +
  geom_hline(yintercept = 0, color = "#808080") +
  { if (!is.na(zero_idx)) geom_vline(xintercept = zero_idx, linetype = "dashed", color = "#404040") } +
  { if (!is.na(zero_idx)) annotate("text", x = zero_idx, y = min(bar_df$val, na.rm = TRUE),
                                   label = paste0("Zero cross at ", zero_idx),
                                   vjust = -0.5, hjust = 1, size = 3) } +
  scale_x_continuous(limits = limits_x, breaks = breaks_x, expand = c(0,0)) +
  labs(x = "Rank in Ordered Dataset", y = "Ranked list metric") +
  theme_minimal(base_size = 12)

p_final <- p1 / p2 / p3 / p4 + plot_layout(heights = c(3.0, 0.6, 0.8, 1.8)) +
  plot_annotation(
    title = pathway_title,
    subtitle = paste0("metric = -log10(p)×sign(log2FC) | ES = ", sprintf("%.3f", ES),
                      " | hits = ", Nh, "/", N),
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold"),
                  plot.subtitle = element_text(hjust = 0.5))
  )

ggsave(out_pdf, p_final, width = 12, height = 7.2)
