library(tidyverse)
library(magrittr)

# input files used:
# phenotype_file contains the biobank phenotypes
# all_ocac_phenotypes.txt contains the ids in imputation order
# all_panels_exclude_candidates.txt contains individuals to exclude based on extra QC checks

phenotype_file = "pheno_info20230504.tsv"

GenerateCancerPhenotypes = function(canc, sex_exclusion=-1, ethnicity_exclude=0, use_extra_qc=TRUE, use_all=FALSE){
  
  # read phenotype file supplied  
  t<-read.csv(phenotype_file, sep="\t", stringsAsFactors=F)
  
  # read individuals to exclude based on extra QC checks
  exclude_extra_qc = scan("all_panels_exclude_candidates.txt")
  
  
  if(use_all){
    # get all the phenotypes
    phenotypes_joined = t
  }else{
    # put the phenotypes in the same order as the imputed biobank genotypes  
    # read ids in order      
    phenotypes = read_tsv("all_ocac_phenotypes.txt", col_types = cols(id_v0=col_integer())) %>% rename(ID_v0=id_v0)
    
    phenotypes_joined = left_join(phenotypes, t)
  }
  
  # when cancer is present work out what index it is
  # this line works by determining which column contains the cancer code (canc)
  
  cancer_index = phenotypes_joined %>% select(matches("^diagnosed.cancer[1-6]")) %>% transpose %>% map_int(~match(canc,.))
  
  # works out what date the cancer was diagnosed  
  
  dates = phenotypes_joined %>% select(matches("^diagnosed.cancer.date[1-6]")) %>% transpose %>% map2_chr(cancer_index, ~if_else(is.na(.y),NA_character_,.x[[.y]]))
  
  # works out the histotype of  the specified cancer 
  
  histotypes = phenotypes_joined %>% select(matches("^diagnosed.cancer.histology[1-6]")) %>% transpose %>% map2_int(cancer_index, ~if_else(is.na(.y),NA_integer_,.x[[.y]]))
  
  # works out the behaviour of the specified cancer  
  
  #  behaviour = phenotypes_joined %>% select(matches("^diagnosed.cancer.behaviour[1-3]")) %>% transpose %>% map2_chr(cancer_index, ~if_else(is.na(.y),NA_character_,.x[[.y]]))
  
  # works out the lag between the cancer diagnosis and the date of attending
  
  lag = as.Date(dates)-as.Date(phenotypes_joined$Date.of.attending.assessment.centre_v0)
  
  # select variables of interest  
  
  #  phenotypes_converted = phenotypes_joined %>% 
  #    select(matches("^(ID|prs|pc1|all|Genetic_sex|Sex|Genetic.sex|Date.of.attending.assessment.centre_v0|first.diagnosed.cancer|Age.when.attended.assessment.centre_v0)|mastectomy|self.cancer.type|diagnosed.cancer([0-9]|.histology[0-9]|.date[0-9]|.behaviour)"))
  
  phenotypes_converted = phenotypes_joined %>% 
    select(matches("^(ID|prs|pc1|all|Genetic_sex|Sex|Genetic.sex|Date.of.attending.assessment.centre_v0|first.diagnosed.cancer|Age.when.attended.assessment.centre_v0)|mastectomy|self.cancer.type|diagnosed.cancer([0-9]|.histology[0-9]|.date[0-9]|.behaviour)"))
  
  
  phenotypes_converted %<>% add_column(cancer_index, dates, histotypes, lag)
  
  # status is set to 0 if cancer index is missing and 1 if it was found  
  phenotypes_converted %<>% mutate(status = as.integer(!is.na(cancer_index)))   
  
  biobank_calls = read_delim(delim=' ', "ukb_cal_chr1_v2.fam", col_names=c("ID_V0","ID_exclude_qc","dummy1","dummy2","dummy3","batch"), col_types="iiiiic")
  
  # work out whether the biobank individual is in one of the cancer studies 
  phenotypes_converted %<>% mutate(exclude_qc=ID_v0 %in% exclude_extra_qc)
  
  phenotypes_converted = left_join(by=c("ID_v0"="ID_V0"), phenotypes_converted, biobank_calls %>% select(ID_V0, batch))
  
  phenotypes_converted %<>% mutate(first_phase=str_detect(batch,"UKBiLEVEAX|b0[01]|b02[012]"))
  
  # biobank_postqc contains individuals that pass basic qc metrics
  biobank_postqc = read_csv("ukb_samples_post_genotyping_qc.csv", col_types = "iiiiicc") %>% rename(ID_v0=ukbid)
  
  phenotypes_converted = left_join(by="ID_v0",phenotypes_converted, biobank_postqc)
  
  
  # cancer is set to status but is set to missing when a number of conditions aren't met. These are:
  # white british is the same as ethnicity_exclude
  # if the individual fails more stringent QC check
  # when either sex or Genetic sex is the sex to be excluded (used for sex-specific cancers)
  # when sex_aneuploidy is 1 (true)
  # for if_else if any of the variables tested are missing then cancer will also be set to missing (some qc checks have failed etc.), which is the behaviour wanted
  phenotypes_converted %<>% mutate(cancer=if_else(!(use_all & is.na(white_british)) & (white_british==ethnicity_exclude | exclude_qc & use_extra_qc | Sex_v0==sex_exclusion | Genetic.sex_v0==sex_exclusion | sex_aneuploidy==1), NA_integer_, as.integer(status)))
  
  # cancer_prospective the same as cancer provided the time between study entry and cancer is at least half a year    
  phenotypes_converted %<>% mutate(cancer_prospective = case_when(lag<=183~NA_integer_,T~cancer))
}


Generate_several_phenotypes = function(ethnicity_exclude, use_extra_qc = TRUE, use_all=FALSE){
  breast_all = GenerateCancerPhenotypes("C50", "Male", ethnicity_exclude, use_extra_qc, use_all)
  
  endometrial_all = GenerateCancerPhenotypes("C54", "Male", ethnicity_exclude, use_extra_qc, use_all)
  
  prostate_all = GenerateCancerPhenotypes("C61", "Female", ethnicity_exclude, use_extra_qc, use_all)
  
  colorectal_all = GenerateCancerPhenotypes("C18", "-1", ethnicity_exclude, FALSE, use_all)
  
  # breast prospective is set to missing if there has been a mastectomy
  
  breast_all[!is.na(breast_all$mastectomy1_v0), "cancer_prospective"]=NA_integer_
  
  # endometrial is set to missing if individual was in first stage as this was used for the endometrial summary statistics
  
  endometrial_all[which(endometrial_all$first_phase==TRUE), "cancer"] = NA_integer_
  
  endometrial_all[which(endometrial_all$first_phase==TRUE), "cancer_prospective"] = NA_integer_
  
  all_phenotypes = endometrial_all
  
  all_phenotypes %<>% rename(endometrial=cancer, endometrial_prospective=cancer_prospective)
  
  all_phenotypes$breast = breast_all$cancer
  
  all_phenotypes$breast_prospective = breast_all$cancer_prospective
  
  all_phenotypes$prostate = prostate_all$cancer
  
  all_phenotypes$prostate_prospective = prostate_all$cancer_prospective
  
  all_phenotypes$colorectal = colorectal_all$cancer
  
  all_phenotypes$colorectal_prospective = colorectal_all$cancer_prospective
  
  all_phenotypes
}


biobank_ethnicity = read_tsv("biobank_ethnicity_proportions.txt")

phenotypes_several = Generate_several_phenotypes(-1, FALSE, TRUE)

phenotypes_several %<>% left_join(biobank_ethnicity) 

phenotypes_several_imputed = Generate_several_phenotypes(-1, FALSE, FALSE)

phenotypes_several_imputed %<>% left_join(biobank_ethnicity)

phenotypes_several_analysis = Generate_several_phenotypes(-1, TRUE, FALSE)

phenotypes_several_analysis %<>% left_join(biobank_ethnicity)

phenotypes_several %>% write_tsv("biobank_cancer_phenotypes_all_20230504.txt")

phenotypes_several_imputed %>% write_tsv("biobank_cancer_phenotypes_imputed_20230504.txt")

phenotypes_several_analysis %>% write_tsv("biobank_cancer_phenotypes_analysis_20230504.txt")



