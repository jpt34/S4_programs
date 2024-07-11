library(bigsnpr)
library(tidyverse)
library(bigreadr)
usref <- snp_attach("RefFiles/us_ldpred_all.rds")
G<-usref$genotypes
mapdata = transmute(usref$map, chr=as.integer(chromosome), pos=physical.pos, a0=allele1, a1=allele2)
NCORES=16
sd <-readRDS("us_sd.rds")

sumstats <- bigreadr::fread2(
  paste0("../T1D/meta_chr_", 1:22),
  select = c("chromosome", "position", "a0", "a1", "beta.meta", "se.meta",
             "p.meta", "qc.check", "EUR_MAF_1kG"),
  col.names = c("chr", "pos", "a0", "a1", "beta", "beta_se", "p", "qc", "maf"))

sumstats <- sumstats %>%
  filter(qc == "PASS") %>%
  tidyr::drop_na() %>%
  mutate(n_eff = 4 / (1 / 5913 + 1 / 8828), qc = NULL) %>%
  filter(maf > (1 / sqrt(n_eff))) %>%
  as_tibble() %>%
  print()

#info_snp <- snp_match(sumstats, mapdata)
info_snp <- bigsnpr::snp_match(sumstats, mapdata, strand_flip = FALSE, match.min.prop = 0.4)
(info_snp <- tidyr::drop_na(as_tibble(info_snp)))

sd_val <- sd[info_snp$`_NUM_ID_`]
sd_ss <- with(info_snp, 2 / sqrt(n_eff * beta_se^2))
is_bad <- sd_ss < (0.5 * sd_val) | sd_ss > (sd_val + 0.1) | sd_ss < 0.1 | sd_val < 0.05
info_snp[is_bad, ] %>% arrange(p) %>% head(20) %>% mutate(freq2 = big_scale()(G, ind.col = `_NUM_ID_`)$center / 2) %>% select(-`_NUM_ID_.ss`, -`_NUM_ID_`)
table(is_bad)
saveRDS(info_snp[!is_bad, ], "RefFiles/t1d.rds")
