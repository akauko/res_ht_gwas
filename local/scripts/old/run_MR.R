#Run TwoSampleMR 
#Runned locally as it is easier to work with internet connection...

#install.packages("remotes")
#remotes::install_github("MRCIEU/TwoSampleMR")
library(data.table)
library(tidyverse)
library(vroom)
library(TwoSampleMR)

# List available GWASs
ao <- available_outcomes()
#ao %>% filter(str_detect(trait, "aist circum"))  %>% arrange(desc(sample_size)) %>% as.data.table
#ao %>% filter(str_detect(trait, "lbumin"), str_detect(trait,"reatin")) %>% arrange(desc(sample_size)) %>% as.data.table
#ao %>% filter(str_detect(author, "Teumer")) %>% arrange(desc(sample_size)) %>% as.data.table
# Get instruments
#Define and extract risk traits
exposures <- c(`ukb-b-2303`="BMI", `ukb-b-10787`="Height", `ukb-b-9405` = "Waist circumference",
            `ieu-b-73`="Alcohol", `ukb-b-223`="Smoking", `ukb-b-4630`="Neuroticism",
            `ieu-b-110` = "LDLc", `ieu-b-109` = "HDLc",   `ieu-b-111` = "Triglycerides", `met-d-Total_C` = "Total cholesterol",
            `ebi-a-GCST90002232`="Glucose", `ukb-d-30880_irnt`="Urate", `ukb-d-30710_irnt`="CRP", `ebi-a-GCST90002244` = "HbA1C", 
            `ieu-b-4832` = "eGFR")

exposure_ids <- names(exposures)
ao %>% filter(id %in% exposure_ids)  %>% as.data.table

#Rename outcome columns to correspond defaults and add phenotype... would be faster with bash
outcome_file_name <- "data/gwas_rs/RES_HT_rs.gz"
if(!file.exists(outcome_file_name)){
  my_gwas <- fread("data/gwas_rs/RES_HT_rs.gz") %>%
    rename(SNP=rsid, se=sebeta, eaf=af_alt, effect_allele=alt, other_allele=ref) %>%
    mutate(Phenotype = "RES_HT")
  fwrite(my_gwas, "data/gwas_rs/RES_HT_rs_renamed.gz", sep="\t")
}

#Initialize list with empty values
res <- as.list(exposures)
res <- lapply(exposures, function(x){x=NULL})

#Each MR is run separatedly
for(id in exposure_ids){

  expos <- exposures[id]
  print(expos)
  
  
  # Get instruments
  exposure_dat <- extract_instruments(id)
  
  # Get effects of instruments on outcome
  outcome_dat <- read_outcome_data("data/gwas_rs/RES_HT_rs_renamed.gz", sep="\t", snps = exposure_dat$SNP)
  
  # Harmonise the exposure and outcome data
  dat <- harmonise_data(exposure_dat, outcome_dat)
  
  # Perform MR
  res[[id]] <- mr(dat)
  
}

# UACR had only couple of significant associations. We downloaded it manually from GWAS catalog
# and preprocess it.
exposures["uacr"] = "UACR"
file_name <- "data/gwas_risktraits/UACR_renamed_top.gz"
if(!file.exists(file_name)){
  gwas_top <- vroom("data/gwas_risktraits/formatted_20180517-UACR_overall-EA-nstud_18-SumMac_400.tbl.rsid.gz", delim=" ") %>%
    rename(SNP=RSID, beta=Effect, se=StdErr, pval="P-value", eaf=Freq1, effect_allele=Allele1, other_allele=Allele2) %>%
    mutate(Phenotype = "UACR") %>%
    filter(pval<=5e-8)
  fwrite(gwas_top, file_name, sep="\t")
}
uacr <-  read_exposure_data(file_name, sep="\t", clump=T, id_col="uacr") %>% mutate(id.exposure="uacr")
outcome_dat <- read_outcome_data("data/gwas_rs/RES_HT_rs_renamed.gz", sep="\t", snps = uacr$SNP)
dat <- harmonise_data(uacr, outcome_dat)
res[["uacr"]]  <-  mr(dat)


#Plot
#Preprocess data for the plot

#To table
res_table <- lapply(res, function(expos){
  expos %>%
    select(id.exposure, method, b, pval) %>%
    mutate(pval = if_else(pval == 0, 1e-200, pval)) %>%
    mutate(Z = -qnorm(pval/2)*sign(b)) %>%
    mutate(exposure = exposures[id.exposure]) %>%
    mutate(FDR=p.adjust(pval, method="fdr", n=16*5))
}) %>%
  bind_rows %>%
  mutate(signific = if_else(pval < 0.05, ".", ""),
         signific = if_else(FDR < 0.05, "*", signific))

#To matrix
res_matrix <- 
  res_table %>%
  select(exposure, Z, method) %>%
  pivot_wider(names_from = "exposure", values_from="Z") %>%
  column_to_rownames("method") %>%
  as.matrix

#Cluster
ord_method <- hclust(dist(res_matrix))$order
ord_exposure <- hclust(dist(t(res_matrix)))$order

#Plot
res_table %>%
  mutate(method = factor(method, levels = rownames(res_matrix)[ord_method])) %>%
  mutate(exposure = factor(exposure, levels = colnames(res_matrix)[ord_exposure])) %>%
  ggplot(aes(x = exposure, method)) +
  geom_tile(aes(fill=Z), color="gray") +
  geom_text(aes(label=signific)) +
  scale_fill_gradient2(low="blue", mid="white", high="red", na.value="gray") +
  scale_x_discrete(position="top") +
  scale_y_discrete(limits=rev) +
  theme(axis.text.x = element_text(angle = 45, hjust = -0.1),
        axis.title = element_blank(),
        panel.background = element_rect(fill="gray90"),
        panel.grid = element_blank())

# Hypertension as outcome


#Initialize list with empty values
res.ht <- as.list(exposures)
res.ht <- lapply(exposures, function(x){x=NULL})

#Each MR is run separatedly
for(id in exposure_ids){
  
  expos <- exposures[id]
  print(expos)
  
  # Get instruments
  exposure_dat <- extract_instruments(id)
  
  # Get effects of instruments on outcome
  outcome_dat <- extract_outcome_data(outcomes = "finn-b-I9_HYPTENS", snps=exposure_dat$SNP)
  #outcome_dat <- read_outcome_data("data/gwas_rs/RES_HT_rs_renamed.gz", sep="\t", snps = exposure_dat$SNP)
  
  # Harmonise the exposure and outcome data
  dat <- harmonise_data(exposure_dat, outcome_dat)
  
  # Perform MR
  res.ht[[id]] <- mr(dat)
  
}
#UACR
uacr <-  read_exposure_data(file_name, sep="\t", clump=T, id_col="uacr") %>% mutate(id.exposure="uacr")
outcome_dat <- extract_outcome_data(outcomes = "finn-b-I9_HYPTENS", snps=uacr$SNP)
dat <- harmonise_data(uacr, outcome_dat)
res.ht[["uacr"]]  <-  mr(dat)


#Plot
#Preprocess data for the plot

#To table
res_table.ht <- lapply(res.ht, function(expos){
  expos %>%
    select(id.exposure, method, b, pval) %>%
    mutate(pval = if_else(pval == 0, 1e-200, pval)) %>%
    mutate(Z = -qnorm(pval/2)*sign(b)) %>%
    mutate(exposure = exposures[id.exposure]) %>%
    mutate(FDR=p.adjust(pval, method="fdr", n=16*5))
}) %>%
  bind_rows %>%
  mutate(signific = if_else(pval < 0.05, ".", ""),
         signific = if_else(FDR < 0.05, "*", signific))


#Plot
res_table.ht %>%
  mutate(method = factor(method, levels = rownames(res_matrix)[ord_method])) %>%
  mutate(exposure = factor(exposure, levels = colnames(res_matrix)[ord_exposure])) %>%
  ggplot(aes(x = exposure, method)) +
  geom_tile(aes(fill=Z), color="gray") +
  geom_text(aes(label=signific)) +
  scale_fill_gradient2(low="blue", mid="white", high="red", na.value="gray") +
  scale_x_discrete(position="top") +
  scale_y_discrete(limits=rev) +
  theme(axis.text.x = element_text(angle = 45, hjust = -0.1),
        axis.title = element_blank(),
        panel.background = element_rect(fill="gray90"),
        panel.grid = element_blank())







