#We will plot z-values for risk our top hit variants in risk  traits
#Package TwoSampleMR is utilized for data harmonisation and download of risk trait GWAS:es. 


#install.packages("remotes")
#remotes::install_github("MRCIEU/TwoSampleMR")
library(data.table)
library(tidyverse)
library(TwoSampleMR)

#We need top hits from our own data - we use variants replicated in ukb and reformat
fread("data/replic/replic_ukb_combined.csv") %>%
  select(rsid) %>%
  left_join(fread("data/regenie_r12/RES_HT_summary.txt"), by="rsid") %>%
  rename(SNP=rsid, beta=beta, se=sebeta, pval=pval, eaf=af_alt, effect_allele=alt, other_allele=ref) %>%
  fwrite("data/regenie_r12/resht_risktrait_input.csv", sep="\t")

ht_res <- read_exposure_data("data/regenie_r12/resht_risktrait_input.csv", sep="\t")


#Available GWAS's
ao <- available_outcomes()
#ao %>% filter(str_detect(trait,"ody mass index"), str_detect(population,"European"))  %>% as.data.table %>% pull(trait)

#Define and extract risk traits
traits <- c(`ukb-b-2303`="BMI", `ukb-b-10787`="Height", `ukb-b-9405` = "Waist",
               `ieu-b-73`="Alcohol", `ukb-b-223`="Smoking", `ukb-b-4630`="Neuroticism",
               `ieu-b-110` = "LDLc", `ieu-b-109` = "HDLc",   `ieu-b-111` = "Triglycerides", `met-d-Total_C` = "Total cholesterol",
               `ebi-a-GCST90002232`="Glucose", `ukb-d-30880_irnt`="Urate", `ukb-d-30710_irnt`="CRP", `ebi-a-GCST90002244` = "HbA1C",`ieu-b-4832` = "eGFR")

#Had missing snips for risk trait table: 
#`ieu-b-40`="BMI", `ieu-a-63` = "Waist, female", `ieu-a-65` = "Waist, male", `ieu-a-1107`= "UACR", `ebi-a-GCST003372` = "eGFR"
#Alternatives: `ieu-a-103` = "Waist, female", `ieu-a-102` = "Waist, male" ; Too small amount of SNPs: , `ieu-a-1107`= "UACR"
trait_ids <- names(traits)
ao %>% filter(id %in% trait_ids)  %>% as.data.table

#Exatract
trait.list <- lapply(as.list(trait_ids), function(id){
  #print(id)
  extract_outcome_data(snps=ht_res$SNP, outcomes = id) %>%
    mutate(TRAIT = traits[id])
}) %>% setNames(traits)


#UACR had very few hits so better summary downloaded manually
#http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST008001-GCST009000/GCST008794/
file_name <- "data/gwas_risktraits/UACR_renamed.gz"
if(!file.exists(file_name)){
  gwas <- vroom("data/gwas_risktraits/formatted_20180517-UACR_overall-EA-nstud_18-SumMac_400.tbl.rsid.gz", delim=" ") %>%
    rename(SNP=RSID, beta=Effect, se=StdErr, pval="P-value", eaf=Freq1, effect_allele=Allele1, other_allele=Allele2) %>%
    mutate(Phenotype = "UACR") %>%
  fwrite(gwas, file_name, sep="\t")
}
uacr.outcome_dat <- read_outcome_data(file_name, sep="\t", snps = ht_res$SNP) %>% 
  mutate(id.outcome="uacr", TRAIT="UACR") 
traits["uacr"] = "UACR"
trait.list[["uacr"]] = uacr.outcome_dat

#Harmonize the data
trait.list.h <- lapply(trait.list, function(trait){
  harmonise_data(ht_res, trait)
}) %>% setNames(traits)

#Preprocess data for the plot

#rsid-gene names
rsid_gene <- fread("data/replic/replic_ukb_combined.csv") %>%
  select(SNP=rsid, GENE=gene_most_severe)

#To table
risk_traits_table <- lapply(trait.list.h, function(trait){
  trait %>%
    select(SNP, P=pval.outcome, BETA=beta.outcome, MAF=eaf.outcome, TRAIT) %>%
    mutate(P = if_else(P == 0, 1e-200, P)) %>%
    mutate(Z = -qnorm(P/2)*sign(BETA))
}) %>%
  bind_rows %>%
  left_join(rsid_gene, by="SNP") %>%
  mutate(RSID_GENE = str_glue("{SNP} ({GENE})")) %>%
  mutate(signific = if_else(P < 1e-5, ".", ""),
         signific = if_else(P < 5e-8, "*", signific))
  
#To matrix
risk_traits_matrix <- 
  risk_traits_table %>%
  select(TRAIT, Z, RSID_GENE) %>%
  pivot_wider(names_from = "TRAIT", values_from="Z") %>%
  column_to_rownames("RSID_GENE") %>%
  as.matrix

#Cluster
ord_gene <- hclust(dist(risk_traits_matrix))$order
ord_trait <- hclust(dist(t(risk_traits_matrix)))$order

#Plot
risk_traits_table %>%
  mutate(RSID_GENE = factor(RSID_GENE, levels = rownames(risk_traits_matrix)[ord_gene])) %>%
  mutate(TRAIT = factor(TRAIT, levels = colnames(risk_traits_matrix)[ord_trait])) %>%
  ggplot(aes(x = TRAIT, RSID_GENE)) +
    geom_tile(aes(fill=Z), color="gray") +
    geom_text(aes(label=signific)) +
    scale_fill_gradient2(low="blue", mid="white", high="red", na.value="gray") +
    scale_x_discrete(position="top") +
    scale_y_discrete(limits=rev) +
    theme(axis.text.x = element_text(angle = 45, hjust = -0.1),
          axis.title = element_blank(),
          panel.background = element_rect(fill="gray90"),
          panel.grid = element_blank())
  
ggsave("figs/resht_risk_r12.jpg", width=7, height=3)
ggsave("figs/resht_risk_r12.pdf", width=7, height=3)


