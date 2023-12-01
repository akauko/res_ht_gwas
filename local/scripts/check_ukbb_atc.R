library(data.table)
library(tidyverse)

ukbb_nlist <- fread("data/ukbb20003_n.csv")
ukbb_atc <- fread("data/ATC-codes_ukbb.csv")
ukbb_atc_am <- fread("data/atc_all_matches_c0.csv")

tmp <- ukbb_atc %>% full_join(ukbb_atc_am, by="ukbb20003_nr" ) %>%
  filter(is.na(ukbb20003.x)) %>% select(ukbb20003_nr, ukbb20003.y, atc_ukbb)

ukbb_atc_n <- ukbb_atc %>% 
  left_join(ukbb_nlist, by="ukbb20003") %>%
  relocate(n, .before=nr_med)
fwrite(ukbb_atc_n, "ATC_ukbb_c0_tmp.csv", sep=";")

ukbb_atc_am <- fread("data/atc_all_matches.csv") %>%
  rename(ukbb20003_nr=V1, ukbb20003=V2, ukbb_atc=V3)
tmp2 <- ukbb_nlist %>%
  left_join(ukbb_atc_am) %>% 
  select(ukbb20003, n, ukbb_atc) %>%
  filter(is.na(ukbb_atc))
fwrite(tmp2, "data/ukbb20003_n_filt_tmp.csv", sep=";")



test <- ukbb_nlist %>% right_join(ukbb_atc, by = "ukbb20003") %>%
  mutate(test.n = n.x-n.y)


