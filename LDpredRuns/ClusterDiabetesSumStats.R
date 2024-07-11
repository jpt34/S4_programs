library(bigsnpr)
library(tidyverse)
library(bigreadr)
usref <- snp_attach("RefFiles/us_ldpred_all.rds")
G<-usref$genotypes
mapdata = transmute(usref$map, chr=as.integer(chromosome), pos=physical.pos, a0=allele1, a1=allele2)
NCORES=16
sd <-readRDS("us_sd.rds")
sumstats <- fread2("../Diabetes/METAANALYSIS_DIAGRAM_SE1.txt")
sumstats <- tidyr::separate(sumstats, "Chr:Position", c("chr", "pos"), convert = TRUE)
names(sumstats) <- c("chr", "pos", "a1", "a0", "beta", "beta_se", "p", "N")
sumstats$n_eff <- 4 / (1 / 26676 + 1 / 132532)

info_snp <- bigsnpr::snp_match(sumstats, mapdata, strand_flip = FALSE)
(info_snp <- as_tibble(info_snp))

sd_val <- sd[info_snp$`_NUM_ID_`]
sd_ss <- with(info_snp, 2 / sqrt(n_eff * beta_se^2))
is_bad <- sd_ss < (0.5 * sd_val) | sd_ss > (sd_val + 0.1) | sd_ss < 0.1 | sd_val < 0.05
info_snp[is_bad, ] %>% arrange(p) %>% head(20) %>% mutate(freq2 = big_scale()(G, ind.col = `_NUM_ID_`)$center / 2) %>% select(-`_NUM_ID_.ss`, -`_NUM_ID_`)
saveRDS(info_snp[!is_bad, ], "RefFiles/diabetes.rds")
