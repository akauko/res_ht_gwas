#Extracts GO terms at levels 2 and 3, assignes gene annotation and writes out in Vegas2 compatible format.
#Run at local computer; Bioconductor is easier to handle when there is access to internet.
#We are interested on biological process and molecular function.

library(data.table)
library(tidyverse)


#Bioconductor packages
if(F){
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

  BiocManager::install("GO.db")
  BiocManager::install("org.Hs.eg.db")
}
library(GO.db)
library(org.Hs.eg.db)
library(annotate)


getAllChildren <- function(goids, gomapping)
{
  ans <- unique(unlist(mget(goids, get(gomapping)), use.names=FALSE))
  ans <- ans[!is.na(ans)]
}

#Biological process
level1_BP_terms <- getAllChildren("GO:0008150", "GOBPCHILDREN")   
level2_BP_terms <- getAllChildren(level1_BP_terms, "GOBPCHILDREN")  
level3_BP_terms <- getAllChildren(level2_BP_terms, "GOBPCHILDREN")  
lapply(list(1,2,3), function(x)length(get(str_glue("level{x}_BP_terms"))))

#Molecular function
level1_MF_terms <- getAllChildren("GO:0003674", "GOMFCHILDREN")   
level2_MF_terms <- getAllChildren(level1_MF_terms, "GOMFCHILDREN")  
level3_MF_terms <- getAllChildren(level2_MF_terms, "GOMFCHILDREN")  
lapply(list(1,2,3), function(x)length(get(str_glue("level{x}_MF_terms"))))

#Lists combined at second and third level
level2_terms <- c(level2_BP_terms, level2_MF_terms) 
level3_terms <- c(level3_BP_terms, level3_MF_terms) 
length(level2_terms)
length(level3_terms)

#Gene annotation
level2_genes.entrez <- mget(intersect(level2_terms, keys(org.Hs.egGO2EG)), org.Hs.egGO2EG)
level2_genes <- lapply(level2_genes.entrez, getSYMBOL, "org.Hs.eg") 
length(level2_genes)
level3_genes.entrez <- mget(intersect(level3_terms, keys(org.Hs.egGO2EG)), org.Hs.egGO2EG)
level3_genes <- lapply(level3_genes.entrez, getSYMBOL, "org.Hs.eg")  
length(level3_genes)


#Data formatted for vegas
gene_go_level2 <-
  lapply(as.list(names(level2_genes)), function(x){tibble(GeneID=level2_genes[[x]],GO_ID=x, Term=Term(x))}) %>% 
  bind_rows %>% 
  mutate(Term=str_replace_all(Term, " ", "_")) %>%
  unite(GO_ID, GO_ID:Term, sep="_")

gene_go_level3 <-
  lapply(as.list(names(level3_genes)), function(x){tibble(GeneID=level3_genes[[x]],GO_ID=x, Term=Term(x))}) %>% 
  bind_rows %>% 
  mutate(Term=str_replace_all(Term, " ", "_")) %>%
  unite(GO_ID, GO_ID:Term, sep="_")

fwrite(gene_go_level2, "go_level2_vegas.txt", sep=" ")
fwrite(gene_go_level3, "go_level3_vegas.txt", sep=" ")

