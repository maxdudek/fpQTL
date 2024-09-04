suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
  library(GWASTools)
  library(qqman)
  library(ggpubr)
})

args = commandArgs(trailingOnly=TRUE)
FP_METHOD <- args[1]
DIR <- paste0("FP_methods/", FP_METHOD)
cat("FP_METHOD = ", FP_METHOD, "\n")

cat("Loading data...\n")
fpscore_regression <- readRDS(paste0(DIR, "/regression_results/fp_score_genotype_regression.Rds"))
fpscore_cov_regression <- readRDS(paste0(DIR, "/regression_results/fp_score_covariates_genotype_regression.Rds"))
variant_info <- readRDS(paste0(DIR, "/data/variant_info_regression.Rds"))


Save smaller FP score/genotype matrix with only sig variants
cat("Loading matrices...\n")
fpscore_matrix <- readRDS(paste0(DIR, "/data/fpscore_matrix_regression.Rds"))
fpscore_matrix_normalized <- readRDS(paste0(DIR, "/data/fpscore_matrix_normalized.Rds"))
genotype_matrix <- readRDS(paste0(DIR, "/data/genotype_matrix_regression.Rds"))

cat("Saving smaller matrices...\n")
fpQTLs <- fpscore_regression %>% filter(ST_qval < 0.05) %>% pull(variant_id)
fpscore_matrix[fpQTLs, ] %>% saveRDS(paste0(DIR, "/data/fpscore_matrix_fpQTLs_no_covariates.Rds"))
fpscore_matrix_normalized[fpQTLs, ] %>% saveRDS(paste0(DIR, "/data/fpscore_matrix_normalized_fpQTLs_no_covariates.Rds"))
genotype_matrix[fpQTLs, ] %>% saveRDS(paste0(DIR, "/data/genotype_matrix_fpQTLs_no_covariates.Rds"))

fpQTLs <- fpscore_cov_regression %>% filter(ST_qval < 0.05) %>% pull(variant_id)
fpscore_matrix[fpQTLs, ] %>% saveRDS(paste0(DIR, "/data/fpscore_matrix_fpQTLs_with_covariates.Rds"))
fpscore_matrix_normalized[fpQTLs, ] %>% saveRDS(paste0(DIR, "/data/fpscore_matrix_normalized_fpQTLs_with_covariates.Rds"))
genotype_matrix[fpQTLs, ] %>% saveRDS(paste0(DIR, "/data/genotype_matrix_fpQTLs_with_covariates.Rds"))

cat("Freeing memory...\n")
rm(fpscore_matrix)
rm(fpscore_matrix_normalized)
rm(genotype_matrix)

# Calculate p-val thresholds for FDR 5
fpscore_regression %>%
  filter(ST_qval < 0.05) %>%
  pull(pval) %>%
  max() ->
  ST_FDR5_p_thresh_no_cov

fpscore_cov_regression %>%
  filter(ST_qval < 0.05) %>%
  pull(pval) %>%
  max() ->
  ST_FDR5_p_thresh_with_cov

cat("Writing FDR5 variant tables...\n")
# Write FDR5 variant table
fpscore_regression %>%
  filter(ST_qval < 0.05) %>%
  dplyr::select(-c(ST_lfdr)) %>%
  arrange(ST_qval) %>%
  write.table(file = paste0(DIR, "/regression_results/fpQTLs_FDR5.txt"),
              quote = F,
              row.names = F,
              sep = "\t")

# Write FDR5 variant table
fpscore_cov_regression %>%
  filter(ST_qval < 0.05) %>%
  dplyr::select(-c(ST_lfdr)) %>%
  arrange(ST_qval) %>%
  write.table(file = paste0(DIR, "/regression_results/fpQTLs_covariates_FDR5.txt"),
              quote = F,
              row.names = F,
              sep = "\t")

cat("Writing FDR10 variant tables...\n")
# Write FDR10 variant table
fpscore_regression %>%
  filter(ST_qval < 0.1) %>%
  dplyr::select(-c(ST_lfdr)) %>%
  arrange(ST_qval) %>%
  write.table(file = paste0(DIR, "/regression_results/fpQTLs_FDR10.txt"),
              quote = F,
              row.names = F,
              sep = "\t")

# Write FDR10 variant table
fpscore_cov_regression %>%
  filter(ST_qval < 0.1) %>%
  dplyr::select(-c(ST_lfdr)) %>%
  arrange(ST_qval) %>%
  write.table(file = paste0(DIR, "/regression_results/fpQTLs_covariates_FDR10.txt"),
              quote = F,
              row.names = F,
              sep = "\t")

# Write rsID files
fpscore_cov_regression %>%
  filter(ST_qval < 0.05) %>%
  pull(variant_id) %>%
  write(paste0(DIR, "/regression_results/fpQTLs_covariates_FDR5_rsIDs.txt"))

fpscore_cov_regression %>%
  filter(ST_qval < 0.1) %>%
  pull(variant_id) %>%
  write(paste0(DIR, "/regression_results/fpQTLs_covariates_FDR10_rsIDs.txt"))

cat("Plotting pval distributions...\n")
fpscore_regression %>%
  ggplot(aes(x = pval)) +
  xlim(c(0,1)) +
  geom_histogram(bins = 50) +
  theme_classic()

ggsave(paste0(DIR, "/figures/regression_results_no_covariates/regression_pval_distribution.png"))

fpscore_cov_regression %>%
  ggplot(aes(x = pval)) +
  xlim(c(0,1)) +
  geom_histogram(bins = 50) +
  theme_classic()

ggsave(paste0(DIR, "/figures/regression_results_with_covariates/regression_pval_distribution.png"))

fpscore_regression %>%
  ggplot(aes(x = ST_qval)) +
  xlim(c(0,1)) +
  geom_histogram(bins = 50) +
  theme_classic()

ggsave(paste0(DIR, "/figures/regression_results_no_covariates/regression_ST_qval_distribution.png"))

fpscore_cov_regression %>%
  ggplot(aes(x = ST_qval)) +
  xlim(c(0,1)) +
  geom_histogram(bins = 50) +
  theme_classic()

ggsave(paste0(DIR, "/figures/regression_results_with_covariates/regression_ST_qval_distribution.png"))


cat("Plotting QQ plots...\n")
inflation <- function(ps) {
  chisq <- qchisq(1 - ps, 1)
  lambda <- median(chisq) / qchisq(0.5, 1)
  lambda
}
# QQ plot
ggplot() + 
  stat_qq(aes(sample = -log10(fpscore_regression$pval)), 
          color = "black", distribution = stats::qexp, dparams = list(rate = log(10))) + 
  geom_abline(aes(slope = 1, intercept = 0), linetype = 2) +
  ylab(bquote(Observed ~~ -log[10](p))) +
  xlab(bquote(Expected ~~ -log[10](p))) +
  theme_classic()
ggsave(paste0("FP_methods/", FP_METHOD, "/figures/regression_results_no_covariates/qqplot.png"),
       width = 5, height = 5, units = "in")
lambda_no_cov <- inflation(fpscore_regression$pval)
print(lambda_no_cov)

# QQ plot
ggplot() + 
  stat_qq(aes(sample = -log10(fpscore_cov_regression$pval)), 
          color = "black", distribution = stats::qexp, dparams = list(rate = log(10))) + 
  geom_abline(aes(slope = 1, intercept = 0), linetype = 2) +
  ylab(bquote(Observed ~~ -log[10](p))) +
  xlab(bquote(Expected ~~ -log[10](p))) +
  theme_classic()
ggsave(paste0("FP_methods/", FP_METHOD, "/figures/regression_results_with_covariates/qqplot.png"),
       width = 5, height = 5, units = "in")
lambda_with_cov <- inflation(fpscore_cov_regression$pval)
print(lambda_with_cov)

# QQ plot with Genomic Control
ggplot() + 
  stat_qq(aes(sample = -log10(fpscore_regression$pval / lambda_no_cov)), 
          color = "black", distribution = stats::qexp, dparams = list(rate = log(10))) + 
  geom_abline(aes(slope = 1, intercept = 0), linetype = 2) +
  ylab(bquote(Observed ~~ -log[10](p))) +
  xlab(bquote(Expected ~~ -log[10](p))) +
  theme_classic()
ggsave(paste0("FP_methods/", FP_METHOD, "/figures/regression_results_no_covariates/qqplot_GC.png"),
       width = 5, height = 5, units = "in")

# QQ plot
ggplot() + 
  stat_qq(aes(sample = -log10(fpscore_cov_regression$pval / lambda_with_cov)), 
          color = "black", distribution = stats::qexp, dparams = list(rate = log(10))) + 
  geom_abline(aes(slope = 1, intercept = 0), linetype = 2) +
  ylab(bquote(Observed ~~ -log[10](p))) +
  xlab(bquote(Expected ~~ -log[10](p))) +
  theme_classic()
ggsave(paste0("FP_methods/", FP_METHOD, "/figures/regression_results_with_covariates/qqplot_GC.png"),
       width = 5, height = 5, units = "in")



cat("Plotting Volcano plots...\n")
# Volcano plot
if (grepl("_no_gaussian", FP_METHOD, fixed = TRUE)) {
  xlim <- max(max(fpscore_regression$beta), -min(fpscore_regression$beta),
              max(fpscore_cov_regression$beta), -min(fpscore_cov_regression$beta)) 
  volcano_plot_xlim <- c(-xlim, xlim)
} else {
  volcano_plot_xlim <- c(-1.5, 1.5)
}

fpscore_regression %>%
  ggplot(aes(x = beta, y = -log10(pval))) +
  geom_point(size = 1) +
  geom_abline(slope = 0, intercept = -log10(ST_FDR5_p_thresh_no_cov), size = 0.3, color = "black") +
  xlim(volcano_plot_xlim) +
  ylab(bquote(-log[10](p))) +
  xlab("Effect slope (beta)") +
  theme_classic()

ggsave(paste0(DIR, "/figures/regression_results_no_covariates/volcano_plot.png"))

fpscore_cov_regression %>%
  ggplot(aes(x = beta, y = -log10(pval))) +
  geom_point(size = 0.6) +
  geom_abline(slope = 0, intercept = -log10(ST_FDR5_p_thresh_with_cov), size = 0.3, color = "black") +
  xlim(volcano_plot_xlim) +
  ylab(bquote(-log[10](p))) +
  xlab("Effect slope (beta)") +
  theme_classic()

ggsave(paste0(DIR, "/figures/regression_results_with_covariates/volcano_plot.png"), width = 3, height = 3, dpi = 600)

cat("Plotting Manhattan plots...\n")
# Manhatten plot
png(file=paste0(DIR, "/figures/regression_results_no_covariates/manhatten_plot.png"), width = 1200, height = 600)
fpscore_regression %>%
  mutate(variant_chrom = as.integer(gsub("chr", "", variant_chrom))) %>%
  manhattan(chr = "variant_chrom", bp = "variant_pos", 
            snp = "variant_id", p = "pval",
            col = c("#ff5d00", "#0065ff"),
            suggestiveline = FALSE,
            genomewideline = -log10(ST_FDR5_p_thresh_no_cov))
dev.off()

png(file=paste0(DIR, "/figures/regression_results_with_covariates/manhatten_plot.png"), width = 1200, height = 600)
fpscore_cov_regression %>%
  mutate(variant_chrom = as.integer(gsub("chr", "", variant_chrom))) %>%
  manhattan(chr = "variant_chrom", bp = "variant_pos", 
            snp = "variant_id", p = "pval",
            col = c("#ff5d00", "#0065ff"),
            suggestiveline = FALSE,
            genomewideline = -log10(ST_FDR5_p_thresh_with_cov))
dev.off()

cat("Comparing regression with and without covariates...\n")
data.frame(
  pval_no_cov = -log10(fpscore_regression$pval),
  pval_with_cov = -log10(fpscore_cov_regression$pval)
) %>%
  ggplot(aes(x = pval_no_cov, y = pval_with_cov)) +
  geom_point(size = 0.5) +
  geom_abline(slope = 1, intercept = 0, color = "blue") +
  stat_cor(
    aes(label = paste(after_stat(rr.label), ..p.label.., sep = "~`,`~")),
    size = 3
  ) +
  xlab("-log10 p-value without covariates") +
  ylab("-log10 p-value WITH covariates") +
  theme_classic()
ggsave(paste0(DIR, "/figures/regression_pvals_with_vs_without_covariates.png"), width = 3, height = 3, dpi = 300)

cat("Correlation between pvals with and without covariates...\n")
cor.test(-log10(fpscore_regression$pval), -log10(fpscore_cov_regression$pval),
         method = "spearman")
cor.test(-log10(fpscore_regression$pval), -log10(fpscore_cov_regression$pval),
         method = "pearson")




