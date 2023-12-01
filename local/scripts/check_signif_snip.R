
library(data.table)
library(tidyverse)

finemap_resht <- fread("regenie/RES_HT.SUSIE.cred.summary.tsv")
gwas_fg <- fread("metal/RES_HT_filt.gz") %>%
  unite("v",`#chrom`:alt, sep=":", remove=F)
gwas_fgvar <- fread("metal/assoc.resht_fg.filt.gz")
gwas_ukbvar <- fread("metal/assoc.resht_ukb.filt.gz")
gwas_bothvar <- fread("metal/assoc.resht_both.filt.gz")

#Modify 
finemap_resht <- finemap_resht %>% 
  filter(good_cs ==T) %>%
  select(v, gene_most_severe, most_severe) %>%
  left_join(gwas_fg, by="v") %>%
  rename(variant=v) %>%
  filter(str_detect(rsid, "rs")) %>%
  arrange(`#chrom`, pos) %>%
  select(variant, rsid, gene_most_severe, most_severe, pval, beta, af_alt)

replic_fgvar <- finemap_resht %>%
  left_join(select(gwas_fgvar, ID, P, BETA, A1FREQ, CHROM, GENPOS), by=c("rsid"="ID")) %>%
  mutate(IS_SIGN = if_else(sign(beta)==sign(BETA), T, F)) %>%
  mutate(IS_P = if_else(P < 0.05/n(), T, F)) %>%
  select(rsid, variant, gene_most_severe, most_severe, IS_P, IS_SIGN, pval, P, beta, BETA, af_alt, A1FREQ)

replic_ukbvar <- finemap_resht %>%
  left_join(select(gwas_ukbvar, ID, P, BETA, A1FREQ, CHROM, GENPOS), by=c("rsid"="ID")) %>%
  mutate(IS_SIGN = if_else(sign(beta)==sign(BETA), T, F)) %>%
  mutate(IS_P = if_else(P < 0.05/n(), T, F)) %>%
  select(rsid, variant, gene_most_severe, most_severe, IS_P, IS_SIGN, pval, P, beta, BETA, af_alt, A1FREQ)

replic_bothvar <- finemap_resht %>%
  left_join(select(gwas_bothvar, ID, P, BETA, A1FREQ, CHROM, GENPOS), by=c("rsid"="ID")) %>%
  mutate(IS_SIGN = if_else(sign(beta)==sign(BETA), T, F)) %>%
  mutate(IS_P = if_else(P < 0.05/n(), T, F)) %>%
  select(rsid, variant, gene_most_severe, most_severe, IS_P, IS_SIGN, pval, P, beta, BETA, af_alt, A1FREQ)
  
fwrite(replic_fgvar, "regenie_ukbb/replic_fgvar_all.csv", sep=";")
fwrite(replic_ukbvar, "regenie_ukbb/replic_ukbvar_all.csv", sep=";")
fwrite(replic_bothvar, "regenie_ukbb/replic_bothvar_all.csv", sep=";")

fwrite(filter(replic_fgvar, IS_SIGN & IS_P) %>% select(-IS_SIGN,-IS_P), "regenie_ukbb/replic_fgvar.csv", sep=";")
fwrite(filter(replic_ukbvar, IS_SIGN & IS_P)  %>% select(-IS_SIGN,-IS_P), "regenie_ukbb/replic_ukbvar.csv", sep=";")
fwrite(filter(replic_bothvar, IS_SIGN & IS_P)  %>% select(-IS_SIGN,-IS_P), "regenie_ukbb/replic_bothvar.csv", sep=";")

