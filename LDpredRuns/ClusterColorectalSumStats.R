library(bigsnpr)
library(tidyverse)
library(bigreadr)
oncoref <- snp_attach("RefFiles/onco_ldpred_all.rds.rds")
G<-oncoref$genotypes
mapdata = transmute(oncoref$map, chr=as.integer(chromosome), pos=physical.pos, a0=allele1, a1=allele2)
NCORES=16
#sd <- runonce::save_run(sqrt(big_colstats(G, ncores = NCORES)$var),file = "onco_sd.rds")
sd <-readRDS("onco_sd.rds")
sumstats <- fread2("/scratch/ocac_hrc/Colorectal/eur_crc_with_header.txt", select = c("chromosome", "base_pair_location", "effect_allele", "other_allele", "freq", "beta", "standard_error", "p_value"),col.names = c("chr", "pos", "a1", "a0", "freq","beta", "beta_se", "p"))
sumstats$n_eff <- 4 / (1 / 86854 + 1 / 73673)
info_snp <- bigsnpr::snp_match(sumstats, mapdata, strand_flip = FALSE)
(info_snp <- tidyr::drop_na(as_tibble(info_snp)))
sd_val <- sd[info_snp$`_NUM_ID_`]
sd_ss <- with(info_snp, 2 / sqrt(n_eff * beta_se^2))
is_bad <- sd_ss < (0.5 * sd_val) | sd_ss > (sd_val + 0.1) | sd_ss < 0.1 | sd_val < 0.05
info_snp[is_bad, ] %>% arrange(p) %>% head(20) %>% mutate(freq2 = big_scale()(G, ind.col = `_NUM_ID_`)$center / 2) %>% select(-`_NUM_ID_.ss`, -`_NUM_ID_`)
table(is_bad)
saveRDS(info_snp[!is_bad, ], "RefFiles/colorectal.rds")
write.table(info_snp[!is_bad, ], "RefFiles/colorectal_sumstats.txt", row.names=FALSE, quote=FALSE, sep="\t")
