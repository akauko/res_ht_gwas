library(data.table)
library(tidyverse)

sumstat_in <- "data/regenie_r12/RES_HT_rsn.gz"
sumstat_out <- "data/twas/sumstats_r12/RES_HT_fortwas.gz"

sumstat <- fread(sumstat_in) %>%   
  filter(!is.na(rsid)) %>%  
#  mutate(Z=beta/sebeta, ptest=pnorm(-abs(Z))*2) %>%
  mutate(Z= sign(beta)*qnorm(pval/2), ptest=pnorm(-abs(Z))*2) %>%
  select(SNP=rsid, A1=alt, A2=ref, Z, beta, pval, ptest)

summary(sumstat)
sumstat <- sumstat %>% select(-ptest)
fwrite(sumstat, sumstat_out, sep="\t", quote = F)
