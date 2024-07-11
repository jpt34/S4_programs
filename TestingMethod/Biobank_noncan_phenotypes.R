library(tidyverse)
library(magrittr)
library(ukbtools)

my_ukb_data <- ukb_df("ukb677530")

qc_info = read_tsv("biobank_general_phenotypes.txt", guess_max =1000000)

df_illness <- my_ukb_data %>% select(matches("20002"))

df_ICD10 <- my_ukb_data %>% select(matches("(4000[12])|(4120[24])"))

ids = my_ukb_data$eid

rm(my_ukb_data)

gc()

save.image("my_ukb_data.rda")

#load("my_ukb_data.rda")

ind_diabetes <- sort(unique(unlist(c(
  lapply(df_illness, function(x) which(x %in% 1220:1223)),
  lapply(df_ICD10, function(x) which(substr(x, 1, 3) %in% paste0("E", 10:14)))
))))                                     

ind_TD1 <- sort(unique(unlist(c(
  lapply(df_illness, function(x) which(x == 1222)),
  lapply(df_ICD10, function(x) which(substr(x, 1, 3) == "E10"))
))))
ind_TD2 <- sort(unique(unlist(c(
  lapply(df_illness, function(x) which(x == 1223)),
  lapply(df_ICD10, function(x) which(substr(x, 1, 3) == "E11"))
))))

y <- rep(0, nrow(df_illness))
y[ind_diabetes] <- NA
y[ind_TD2] <- 1
y[ind_TD1] <- NA

y_ids = tibble(ID_v0=ids, t2d = y)

qc_info %<>% left_join(y_ids)

y[ind_TD1] <- 1
y[ind_TD2] <- NA

y_ids = tibble(ID_v0=ids, t1d = y)

qc_info %<>% left_join(y_ids)


ind_MDD <- sort(unique(unlist(c(
  lapply(df_illness, function(x) which(x == 1286)),
  lapply(df_ICD10, function(x) which(substr(x, 1, 3) %in% c("F32", "F33")))
))))

ind_psy <- sort(unique(unlist(c(
  lapply(df_illness, function(x) which(x %in% 1286:1291)),
  lapply(df_ICD10, function(x) which(substr(x, 1, 1) == "F"))
))))

y <- rep(0, nrow(df_illness))
y[ind_psy] <- NA
y[ind_MDD] <- 1

y_ids = tibble(ID_v0=ids, mdd = y)

qc_info %<>% left_join(y_ids)


ind_scz <- sort(unique(unlist(c(
  lapply(df_illness, function(x) which(x == 1289)),
  lapply(df_ICD10, function(x) which(substr(x, 1, 3) == "F20"))
))))

y <- rep(0, nrow(df_illness))
y[ind_psy] <- NA
y[ind_scz] <- 1

y_ids = tibble(ID_v0=ids, scz = y)

qc_info %<>% left_join(y_ids)


ind_crohns <- sort(unique(unlist(c(
  lapply(df_illness, function(x) which(x == 1462)),
  lapply(df_ICD10, function(x) which(substr(x, 1, 3) == "K50"))
))))

ind_ucolitis <- sort(unique(unlist(c(
  lapply(df_illness, function(x) which(x == 1463)),
  lapply(df_ICD10, function(x) which(substr(x, 1, 3) == "K51"))
))))

y <- rep(0, nrow(df_illness))
y[ind_crohns] <- 1
y[ind_ucolitis] <- 1

y_ids = tibble(ID_v0=ids, ibd = y)

qc_info %<>% left_join(y_ids)


ind_stroke <- sort(unique(unlist(c(
  lapply(df_illness, function(x) which(x %in% c(1081,1086,1491,1583))),
  lapply(df_ICD10, function(x) which(substr(x, 1, 3) %in% c("I60","I61","I63","I64")))
))))

y <- rep(0, nrow(df_illness))
y[ind_stroke] <- 1

y_ids = tibble(ID_v0=ids, stroke = y)

qc_info %<>% left_join(y_ids)


ind_CAD <- sort(unique(unlist(c(
#  lapply(df_heart,   function(x) which(x == 1)),
  lapply(df_illness, function(x) which(x == 1075)),
  lapply(df_ICD10,   function(x) which(substr(x, 1, 3) %in% paste0("I", 21:23))),
  lapply(df_ICD10,   function(x) which(x == "I252"))
))))
ind_heart <- sort(unique(unlist(c(
#  lapply(df_heart,   function(x) which(x %in% 1:3)),
  lapply(df_illness, function(x) which(x %in% 1074:1080)),
  lapply(df_ICD10,   function(x) which(substr(x, 1, 1) == "I"))
))))

y <- rep(0, nrow(df_illness))
y[ind_heart] <- NA
y[ind_CAD] <- 1

y_ids = tibble(ID_v0=ids, cad = y)

qc_info %<>% left_join(y_ids)


ind_asthma <- sort(unique(unlist(c(
  lapply(df_illness, function(x) which(x == 1111)),
  lapply(df_ICD10,   function(x) which(substr(x, 1, 3) == "J45"))
))))
ind_respi <- sort(unique(unlist(c(
  lapply(df_illness, function(x) which(x %in% 1111:1125)),
  lapply(df_ICD10,   function(x) which(substr(x, 1, 1) == "J"))
))))
y <- rep(0, nrow(df_illness))
y[ind_respi] <- NA
y[ind_asthma] <- 1

y_ids = tibble(ID_v0=ids, asthma = y)

qc_info %<>% left_join(y_ids)


qc_info %>% write_tsv("biobank_noncan_phenotypes.txt")
