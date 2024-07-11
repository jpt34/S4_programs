library(bigsnpr)
library(bigreadr)
library(tidyverse)
usdata<-snp_attach("RefFiles/us_ldpred_all.rds")
G<-usdata$genotypes
CHR<-as.integer(usdata$map$chromosome)
sumstats = readRDS("BRCA.rds")
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
  
  corr0 <- readRDS(paste0("RefFiles/us_correlation_chr", chr, ".rds"))[ind.chr3, ind.chr3]
  
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

write.table(beta_inf, "beta_inf_values.txt", row.names=FALSE, col.names=FALSE)

(h2_seq <- round(h2_est * c(0.7, 1, 1.4), 4))
(p_seq <- signif(seq_log(1e-5, 1, length.out = 21), 2))
(params <- expand.grid(p = p_seq, h2 = h2_seq, sparse = c(FALSE, TRUE)))

beta_grid <- snp_ldpred2_grid(corr, df_beta, params, ncores = NCORES)
params$sparsity <- colMeans(beta_grid == 0)

write.table(beta_grid, "beta_grid_values.txt", row.names=FALSE, col.names=FALSE)
