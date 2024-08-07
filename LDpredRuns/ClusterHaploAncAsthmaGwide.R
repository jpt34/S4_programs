library(bigsnpr)
library(bigreadr)
library(tidyverse)
oncodata<-snp_attach("RefFiles/haplo_hm3plus.rds")
G<-oncodata$genotypes
CHR<-as.integer(oncodata$map$chromosome)
sumstats = readRDS("RefFiles/asthma_haplo_anc.rds")
NCORES=16
tmp <- tempfile(tmpdir = "tmp-data")
on.exit(file.remove(paste0(tmp, ".sbk")), add = TRUE)

for (chr in 1:22) {

  # print(chr)

  ## indices in 'sumstats'
  ind.chr <- which(sumstats$chr == chr)
  ## indices in 'G'
  ind.chr2 <- sumstats$`_NUM_ID_`[ind.chr]
  ## indices in 'corr'
  ind.chr3 <- match(ind.chr2, which(CHR == chr))

  corr0 <- readRDS(paste0("RefFiles/haplo_asthma_correlation_chr", chr, ".rds"))[ind.chr3, ind.chr3]

  if (chr == 1) {
    df_beta <- sumstats[ind.chr, c("beta", "beta_se", "n_eff", "_NUM_ID_")]
    ld <- Matrix::colSums(corr0^2)
    corr <- as_SFBM(corr0, tmp)
  } else {
    df_beta <- rbind(df_beta, sumstats[ind.chr, c("beta", "beta_se", "n_eff", "_NUM_ID_")])
    ld <- c(ld, Matrix::colSums(corr0^2))
    corr$add_columns(corr0, nrow(corr))
  }
}

(ldsc <- with(df_beta, snp_ldsc(ld, length(ld), chi2 = (beta / beta_se)^2, sample_size = n_eff, blocks = NULL)))
h2_est <- ldsc[["h2"]]

beta_inf <- snp_ldpred2_inf(corr, df_beta, h2_est)

write.table(beta_inf, "beta_asthma_haplo_anc_inf_values.txt", row.names=FALSE, col.names=FALSE)

(h2_seq <- round(h2_est * c(0.25, 0.35, 0.5, 0.7, 1, 1.4), 4))
(p_seq <- signif(seq_log(1e-5, 1, length.out = 21), 2))
(params <- expand.grid(p = p_seq, h2 = h2_seq, sparse = FALSE))

beta_grid <- snp_ldpred2_grid(corr, df_beta, params, ncores = NCORES)
params$sparsity <- colMeans(beta_grid == 0)

write.table(beta_grid, "beta_asthma_haplo_anc_grid_values.txt", row.names=FALSE, col.names=FALSE)
