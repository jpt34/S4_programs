# Load required libraries
library(tidyverse)   # For data manipulation and visualization
library(magrittr)    # For pipe operators
library(ukbtools)    # For working with UK Biobank data

# Load UK Biobank data using ukbtools package
my_ukb_data <- ukb_df("ukb677530")

# Read quality control information from a tab-separated file
qc_info = read_tsv("biobank_general_phenotypes.txt", guess_max = 1000000)

# Create dataframes for specific types of medical information:
# - df_illness: self-reported illnesses (field 20002)
# - df_ICD10: ICD-10 diagnosis codes (fields 40001, 40002, 41202, 41204)
df_illness <- my_ukb_data %>% select(matches("20002"))
df_ICD10 <- my_ukb_data %>% select(matches("(4000[12])|(4120[24])"))

# Store participant IDs
ids = my_ukb_data$eid

# Remove the original large dataframe to save memory
rm(my_ukb_data)

# Perform garbage collection to free up memory
gc()

# Save the current workspace to an R data file
save.image("my_ukb_data.rda")

# The commented line below would load the saved workspace if needed
# load("my_ukb_data.rda")

# Identify diabetes cases:
# - Any diabetes (codes 1220-1223 in self-report or E10-E14 in ICD-10)
# - Type 1 diabetes (code 1222 or E10)
# - Type 2 diabetes (code 1223 or E11)
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

# Create phenotype variables for diabetes:
# - t2d: Type 2 diabetes (1), other diabetes (NA), no diabetes (0)
# - t1d: Type 1 diabetes (1), other diabetes (NA), no diabetes (0)
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

# Identify mental health cases:
# - Major depressive disorder (code 1286 or F32-F33)
# - Any psychiatric disorder (codes 1286-1291 or any F code)
ind_MDD <- sort(unique(unlist(c(
  lapply(df_illness, function(x) which(x == 1286)),
  lapply(df_ICD10, function(x) which(substr(x, 1, 3) %in% c("F32", "F33")))
))))

ind_psy <- sort(unique(unlist(c(
  lapply(df_illness, function(x) which(x %in% 1286:1291)),
  lapply(df_ICD10, function(x) which(substr(x, 1, 1) == "F"))
))))

# Create phenotype variable for major depressive disorder
y <- rep(0, nrow(df_illness))
y[ind_psy] <- NA
y[ind_MDD] <- 1
y_ids = tibble(ID_v0=ids, mdd = y)
qc_info %<>% left_join(y_ids)

# Identify schizophrenia cases (code 1289 or F20)
ind_scz <- sort(unique(unlist(c(
  lapply(df_illness, function(x) which(x == 1289)),
  lapply(df_ICD10, function(x) which(substr(x, 1, 3) == "F20"))
))))

# Create phenotype variable for schizophrenia
y <- rep(0, nrow(df_illness))
y[ind_psy] <- NA
y[ind_scz] <- 1
y_ids = tibble(ID_v0=ids, scz = y)
qc_info %<>% left_join(y_ids)

# Identify inflammatory bowel disease cases:
# - Crohn's disease (code 1462 or K50)
# - Ulcerative colitis (code 1463 or K51)
ind_crohns <- sort(unique(unlist(c(
  lapply(df_illness, function(x) which(x == 1462)),
  lapply(df_ICD10, function(x) which(substr(x, 1, 3) == "K50"))
))))

ind_ucolitis <- sort(unique(unlist(c(
  lapply(df_illness, function(x) which(x == 1463)),
  lapply(df_ICD10, function(x) which(substr(x, 1, 3) == "K51"))
))))

# Create phenotype variable for IBD (both Crohn's and UC)
y <- rep(0, nrow(df_illness))
y[ind_crohns] <- 1
y[ind_ucolitis] <- 1
y_ids = tibble(ID_v0=ids, ibd = y)
qc_info %<>% left_join(y_ids)

# Identify stroke cases (codes 1081,1086,1491,1583 or I60-I64)
ind_stroke <- sort(unique(unlist(c(
  lapply(df_illness, function(x) which(x %in% c(1081,1086,1491,1583))),
  lapply(df_ICD10, function(x) which(substr(x, 1, 3) %in% c("I60","I61","I63","I64")))
))))

# Create phenotype variable for stroke
y <- rep(0, nrow(df_illness))
y[ind_stroke] <- 1
y_ids = tibble(ID_v0=ids, stroke = y)
qc_info %<>% left_join(y_ids)

# Identify coronary artery disease cases:
# - CAD (code 1075 or I21-I23, I252)
# - Any heart disease (codes 1074-1080 or any I code)
ind_CAD <- sort(unique(unlist(c(
  lapply(df_illness, function(x) which(x == 1075)),
  lapply(df_ICD10,   function(x) which(substr(x, 1, 3) %in% paste0("I", 21:23))),
  lapply(df_ICD10,   function(x) which(x == "I252"))
))))
ind_heart <- sort(unique(unlist(c(
  lapply(df_illness, function(x) which(x %in% 1074:1080)),
  lapply(df_ICD10,   function(x) which(substr(x, 1, 1) == "I"))
))))

# Create phenotype variable for CAD
y <- rep(0, nrow(df_illness))
y[ind_heart] <- NA
y[ind_CAD] <- 1
y_ids = tibble(ID_v0=ids, cad = y)
qc_info %<>% left_join(y_ids)

# Identify asthma cases (code 1111 or J45) and any respiratory disease (codes 1111-1125 or any J code)
ind_asthma <- sort(unique(unlist(c(
  lapply(df_illness, function(x) which(x == 1111)),
  lapply(df_ICD10,   function(x) which(substr(x, 1, 3) == "J45"))
))))
ind_respi <- sort(unique(unlist(c(
  lapply(df_illness, function(x) which(x %in% 1111:1125)),
  lapply(df_ICD10,   function(x) which(substr(x, 1, 1) == "J"))
))))

# Create phenotype variable for asthma
y <- rep(0, nrow(df_illness))
y[ind_respi] <- NA
y[ind_asthma] <- 1
y_ids = tibble(ID_v0=ids, asthma = y)
qc_info %<>% left_join(y_ids)

# Save the updated quality control information with new phenotype variables
qc_info %>% write_tsv("biobank_noncan_phenotypes.txt")
