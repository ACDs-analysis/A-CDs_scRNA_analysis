# ============================================================
# ECM-receptor interaction (mmu04512) GSEV plot (4 panels)
# Ranking: -log10(p) * sign(log2FC)
# X ticks fixed: 0,500,1000,1500,2000,2500,3000,3500
# ES panel and grey-bar panel share identical X axis (aligned)
# Zero cross marked in panel 1 & 4
# Output: PDF only
# ============================================================

file <- "C:/Users/45214/Desktop/ACD_LPS 巨噬细胞差异基因.csv"
out_pdf <- "ECM_mmu04512_GSEV_fixedticks_zerocross.pdf"
pathway_title <- "ECM-receptor interaction (mmu04512)"
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

# -------------------- ECM gene set --------------------
ecm_genes <- c(
  "COL1A1","COL1A2","COL2A1","COL3A1","COL4A1","COL4A2","COL4A3","COL4A4","COL4A5","COL4A6",
  "COL5A1","COL5A2","COL5A3","COL6A1","COL6A2","COL6A3","COL7A1","COL8A1","COL8A2","COL9A1",
  "COL9A2","COL9A3","COL10A1","COL11A1","COL11A2","COL12A1","COL13A1","COL14A1","COL15A1",
  "COL16A1","COL17A1","COL18A1","COL21A1","COL22A1",
  "LAMA1","LAMA2","LAMA3","LAMA4","LAMA5","LAMB1","LAMB2","LAMB3","LAMC1","LAMC2","LAMC3",
  "ITGA1","ITGA2","ITGA3","ITGA4","ITGA5","ITGA6","ITGA7","ITGA8","ITGA9","ITGA10","ITGA11",
  "ITGAV","ITGAE","ITGAL","ITGAM","ITGAX",
  "ITGB1","ITGB2","ITGB3","ITGB4","ITGB5","ITGB6","ITGB7","ITGB8",
  "HSPG2","SPP1","FN1","VTN","AGRN","NID1","NID2","TNC","TNR","THBS1","THBS2","THBS3","THBS4",
  "SDC1","SDC2","SDC3","SDC4","SPARC","TGFBI","VCAN","FMOD","BGN","DCN"
)

hits <- names(rank_vec) %in% ecm_genes
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
  labs(title = "Enrichment plot: mmu04512", y = "Enrichment score (ES)", x = NULL) +
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
