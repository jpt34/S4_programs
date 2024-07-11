library(bigsnpr)
library(tidyverse)
library(bigreadr)
usref <- snp_attach("RefFiles/us_ldpred_all.rds")
G<-usref$genotypes
mapdata = transmute(usref$map, chr=as.integer(chromosome), pos=physical.pos, a0=allele1, a1=allele2)
NCORES=16
#sd <- runonce::save_run(sqrt(big_colstats(G, ncores = NCORES)$var),file = "us_sd.rds")
sd <-readRDS("us_sd.rds")
sumstats <- fread2("../Height/height2014_summary_statistics_with_positions.txt", select = c("Chromosome", "Position", "Allele1", "Allele2", "Freq.Allele1.HapMapCEU", "b", "SE", "p", "N"),col.names = c("chr", "pos", "a1", "a0", "freq","beta", "beta_se", "p", "n_eff"))
info_snp <- snp_match(sumstats, mapdata)
info_snp <- bigsnpr::snp_match(sumstats, mapdata, strand_flip = FALSE)
(info_snp <- tidyr::drop_na(as_tibble(info_snp)))
sd_val <- sd[info_snp$`_NUM_ID_`]
sdy = min(sqrt(0.5)*info_snp$beta_se*sqrt(info_snp$n_eff))
sd_ss <- with(info_snp, sdy / sqrt(n_eff * beta_se^2))
sd_n <- with(info_snp, n_eff)
n_med = median(info_snp$n_eff)
is_bad <- sd_ss < (0.5 * sd_val) | sd_ss > (sd_val + 0.1) | sd_ss < 0.1 | sd_val < 0.05 | sd_n < 0.4*n_med
info_snp[is_bad, ] %>% arrange(p) %>% head(20) %>% mutate(freq2 = big_scale()(G, ind.col = `_NUM_ID_`)$center / 2) %>% select(-`_NUM_ID_.ss`, -`_NUM_ID_`)
table(is_bad)
saveRDS(info_snp[!is_bad, ], "RefFiles/height.rds")
