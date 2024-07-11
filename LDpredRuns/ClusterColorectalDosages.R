library(bigsnpr)
library(tidyverse)
library(bigreadr)

ukb <- readRDS("RefFiles/biobank_ldpred_all.rds.rds")
G<-ukb$genotypes
CHR<-as.integer(ukb$map$chromosome)
sumstats = readRDS("RefFiles/colorectal.rds")
NCORES=16

for (chr in 1:22) {

  ind.chr <- which(sumstats$chr == chr)

  if (chr == 1) {
    df_beta <- sumstats[ind.chr, c("beta", "beta_se", "n_eff", "_NUM_ID_")]
  } else {
    df_beta <- rbind(df_beta, sumstats[ind.chr, c("beta", "beta_se", "n_eff", "_NUM_ID_")])
  }
}

beta_grid_tibble=read_delim("beta_colorectal_grid_values.txt"," ",col_names=FALSE)
beta_grid=as.matrix(beta_grid_tibble)

bigparallelr::set_blas_ncores(NCORES)
pred_grid <- big_prodMat(G, beta_grid, ind.col = df_beta[["_NUM_ID_"]])

write.table(pred_grid, "pred_colorectal_grid_values.txt", row.names=FALSE, col.names=FALSE, sep="\t")
