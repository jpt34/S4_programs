library(bigsnpr)
library(tidyverse)
library(bigreadr)
oncoref <- snp_attach("RefFiles/onco_ldpred_all.rds.rds")
G<-oncoref$genotypes
mapdata = transmute(oncoref$map, chr=as.integer(chromosome), pos=physical.pos, a0=allele1, a1=allele2)
NCORES=16
sd <-readRDS("onco_sd.rds")
sumstats <- fread2("../Practical/meta_v3_onco_euro_overall_ChrAll_1_release.txt", select = c("Chr", "position", "Allele1", "Allele2", "Effect", "StdErr","Pvalue", "Freq1"),col.names=c("chr","pos","a1","a0","beta","beta_se","p","freq"))
sumstats2 <- sumstats %>% mutate(a0 = toupper(a0), a1 = toupper(a1), n_eff = 4 / (1 / 79148 + 1 / 61106)) %>% filter(pmin(freq, 1 - freq) > (1 / sqrt(n_eff)))
info_snp <- snp_match(sumstats2, mapdata)
(info_snp <- as_tibble(info_snp))
sd_val <- sd[info_snp$`_NUM_ID_`]
sd_ss <- with(info_snp, 2 / sqrt(n_eff * beta_se^2))
is_bad <- sd_ss < (0.5 * sd_val) | sd_ss > (sd_val + 0.1) | sd_ss < 0.1 | sd_val < 0.05
info_snp[is_bad, ] %>% arrange(p) %>% head(20) %>% mutate(freq2 = big_scale()(G, ind.col = `_NUM_ID_`)$center / 2) %>% select(-`_NUM_ID_.ss`, -`_NUM_ID_`)
saveRDS(info_snp[!is_bad, ], "RefFiles/prostate.rds")
