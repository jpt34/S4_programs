library(tidyverse)
library(magrittr)
library(broom)

load("all_pheno.dta")

if (length(commandArgs(trailingOnly=TRUE))>0){
  args = commandArgs(trailingOnly=TRUE)
}

histotype = args[1]

nruns = as.integer(args[2])

ldpred_prs = read_tsv(str_c(histotype, "_ldpred_actual_biobank_weightings.txt"))

s4_prs = read_tsv(str_c("s4ldpred_", histotype, "_actual_biobank_weightings.txt"))

all_pheno$ldpred_prs = ldpred_prs$prs

all_pheno$s4_prs = s4_prs$s4ldpred_pgs

bootstrap_estimates=map_df(1:nruns,function(y){
  gc()
  
  print(y)
  
  sampled_pheno = all_pheno %>% filter(!is.na(.data[[histotype]]))
  
  sampled_pheno = sampled_pheno %>% slice_sample(n=nrow(sampled_pheno), replace=TRUE)
  
  sampled_pheno %<>% mutate(s4_prs = s4_prs/sd(s4_prs), ldpred_prs = ldpred_prs/sd(ldpred_prs))
  
  bind_rows(glm(str_c(histotype,"~pc1+pc2+pc3+pc4+age+sex+factor(f5400)+s4_prs"), binomial, sampled_pheno) %>% tidy() %>% filter(str_detect(term,"prs")) %>% mutate(phenotype = histotype, run = y),
            glm(str_c(histotype,"~pc1+pc2+pc3+pc4+age+sex+factor(f5400)+ldpred_prs"), binomial, sampled_pheno) %>% tidy() %>% filter(str_detect(term,"prs")) %>% mutate(phenotype = histotype, run = y))
})

bootstrap_wide = bootstrap_estimates %>% select(term, phenotype, run, estimate) %>% pivot_wider(names_from=term, values_from = estimate) %>% mutate(diff=s4_prs - ldpred_prs) 

print(bootstrap_wide, n=nruns)

bootstrap_wide %>% summarise(sd=sd(diff))
