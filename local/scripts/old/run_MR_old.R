#Run TwoSampleMR 
#Runned locally as it is easier to work with internet connection...

#install.packages("remotes")
#remotes::install_github("MRCIEU/TwoSampleMR")
library(data.table)
library(tidyverse)
library(TwoSampleMR)

# List available GWASs
ao <- available_outcomes()

# Get instruments
ao %>% filter(str_detect(trait,"lbumin"), str_detect(trait,"reatin")) %>% as.data.table
exposures <- c(`ieu-b-40`="Body mass index", `ukb-b-10787`="Standing height", `ieu-a-63` = "Waist, female", `ieu-a-65` = "Waist, male",
                  `ieu-b-73`="Alcohol intake", `ukb-b-223`="Current smoking", `ukb-b-4630`="Neuroticism",
                  `ieu-b-110` = "LDL cholesterol", `ieu-b-109` = "HDL cholesterol",   `ieu-b-111` = "Triglycerides", `met-d-Total_C` = "Total cholesterol",
                  `ebi-a-GCST90002232`="Fasting glucose", `ukb-d-30880_irnt`="Urate", `ukb-d-30710_irnt`="C-reactive protein", `ebi-a-GCST90002244` = "HbA1C", 
                  `ebi-a-GCST003372` = "eGFR", `ieu-a-1107`= "UACR")
exposure_ids <- names(exposures)
ao %>% filter(id %in% exposure_ids)  %>% as.data.table

exposures_dat <- lapply(as.list(exposure_ids), extract_instruments) %>% setNames(exposures)


#Outcome

#outcome_test <- extract_outcome_data(snps=exposure_dat$SNP, outcomes = "finn-b-I9_HYPTENS")

#Finemapped summary and original gwas combined
#gwas_tmp <- fread("file_download_09122022_anni_kauko_RES_HT.gz")  %>% 
#  unite("variant", c("#chrom", "pos", "ref", "alt"), sep=":", remove=F)
#top_snips <- fread("download_24112022/RES_HT.csv") %>%
#  left_join(gwas_tmp, by="variant") %>%
#  rename(Phenotype=trait,SNP=rsid, beta=beta, se=sebeta,
#       eaf=af_alt, effect_allele=alt, other_allele=ref, pval=pval)
#fwrite(top_snips, "RES_HT_top_mr.tsv", sep="\t")


outcome_res_ht <- read_outcome_data("RES_HT_top_mr.tsv", sep="\t")
outcome_res_ht_all <-  read_outcome_data("file_download_09122022_anni_kauko_RES_HT.gz", sep="\t",
                                          se="sebeta", eaf="af_alt",  effect_allele="alt", other_allele="ref")
                                         
                                         

dat.list <- lapply(as.list(exposures_dat), function(expos){
  #expos %>% head()
  harmonise_data(expos, outcome_res_ht)
})

exposures_dat[[1]] %>% head


# Get effects of instruments on outcome

  ##Ongelmana on, että mä filtteröin pois mun vain merkitsevät. Mulla pitää olla kaikki messissä. Eli toi mun säätö yllä on tarpeeton. 
  ## 





# Harmonise the exposure and outcome data
dat <- harmonise_data(exposure_dat, outcome_dat)

# Perform MR
res <- mr(dat)

