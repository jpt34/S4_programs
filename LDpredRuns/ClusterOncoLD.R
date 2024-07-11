library(bigsnpr)
library(tidyverse)
onco_data <- snp_attach("RefFiles/onco_ldpred_all.rds.rds")
G <- onco_data$genotypes
CHR <- as.integer(onco_data$map$chromosome)
POS <- onco_data$map$physical.pos
NCORES=16

POS2 <- snp_asGeneticPos(CHR, POS, dir = "tmp-data", ncores = NCORES)

for (chr in 1:22) {

  ind.chr <- which(CHR == chr)

  runonce::save_run(snp_cor(G, ind.col = ind.chr, alpha = 1, infos.pos = POS2[ind.chr], size = 3 / 1000, ncores = NCORES),file = paste0("RefFiles/onco_correlation_chr", chr, ".rds"))
}

