
library(data.table)
library(tidyverse)

gwas.resht_in <- "data/regenie_r12/RES_HT_rsn.gz"
gwas.resht_fg_in <- "data/regenie_ukbb/assoc.resht_fg.regenie.gz"
gwas.resht_ukb_in <- "data/regenie_ukbb/assoc.resht_ukb.regenie.gz"
gwas.resht_both_in <- "data/regenie_ukbb/assoc.resht_both.regenie.gz"

#Filter individuals in FinnGen GWAS by P-value treshold 1e-6

gwas.FG <- fread(gwas.resht_in) %>%
  filter(pval < 1e-6) 
fwrite(gwas.FG, "data/replic/RES_HT_filt.gz", sep="\t", quote = F, na="NA")
ids <- gwas.FG %>% select(rsid) %>% filter(rsid!="")
fwrite(ids, "data/replic/rsid_list_1e-6.txt", quote=F, na="", col.names = F)
rm(gwas.FG)
gc()

ids <- fread( "data/replic/rsid_list_1e-6.txt", header = F, col.names = c("ID"))
names(ids) <- "ID"

fread(gwas.resht_fg_in, fill=T) %>%
  inner_join(ids, by="ID") %>%
  fwrite("data/replic/assoc.resht_fg.filt.gz", sep="\t", quote = F, na="NA")
gc()

fread(gwas.resht_ukb_in, fill=T) %>%
  inner_join(ids, by="ID") %>%
  fwrite("data/replic/assoc.resht_ukb.filt.gz", sep="\t", quote = F, na="NA")
gc()

fread(gwas.resht_both_in, fill=T) %>%
  inner_join(ids, by="ID") %>%
  fwrite("data/replic/assoc.resht_both.filt.gz", sep="\t", quote = F, na="NA")
gc()