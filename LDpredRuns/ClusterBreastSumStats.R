library(bigsnpr)
library(tidyverse)
library(bigreadr)
usref <- snp_attach("RefFiles/us_ldpred_all.rds")
G<-usref$genotypes
mapdata = transmute(usref$map, chr=as.integer(chromosome), pos=physical.pos, a0=allele1, a1=allele2)
NCORES=16
#sd <- runonce::save_run(sqrt(big_colstats(G, ncores = NCORES)$var),file = "us_sd.rds")
sd <-readRDS("us_sd.rds")
sumstats <- fread2("../BCAC/oncoarray_bcac_public_release_oct17.txt", na.strings = "NULL", select = c("chr", "position_b37", "a0", "a1", "bcac_onco_icogs_gwas_eaf_controls", "bcac_onco_icogs_gwas_beta", "bcac_onco_icogs_gwas_se", "bcac_onco_icogs_gwas_P1df"),col.names = c("chr", "pos", "a0", "a1", "freq","beta", "beta_se", "p"))
sumstats$n_eff <- 4 / (1 / 137045 + 1 / 119078)
info_snp <- bigsnpr::snp_match(sumstats, mapdata, strand_flip = FALSE)
(info_snp <- tidyr::drop_na(as_tibble(info_snp)))
sd_val <- sd[info_snp$`_NUM_ID_`]
sd_ss <- with(info_snp, 2 / sqrt(n_eff * beta_se^2))
is_bad <- sd_ss < (0.5 * sd_val) | sd_ss > (sd_val + 0.1) | sd_ss < 0.1 | sd_val < 0.05
info_snp[is_bad, ] %>% arrange(p) %>% head(20) %>% mutate(freq2 = big_scale()(G, ind.col = `_NUM_ID_`)$center / 2) %>% select(-`_NUM_ID_.ss`, -`_NUM_ID_`)
saveRDS(info_snp[!is_bad, ], "RefFiles/BRCA_backup.rds")
