#Run TwoSampleMR 
#Runned locally as it is easier to work with internet connection...

#install.packages("remotes")
#remotes::install_github("MRCIEU/TwoSampleMR")
library(data.table)
library(tidyverse)
library(TwoSampleMR)


#Let's try to use external GWAS for exposures

#Mapping of column names as list of list
gwas_mapping_tmp <- fread("data/gwas_risktraits/gwas_column_mapping.csv") 
gwas_mapping <- list()
for (i in 1:length(gwas_mapping_tmp$Key)){
  row <- gwas_mapping_tmp[i,-1] %>% as.list
  Key <- row[[1]]
  gwas_mapping[[Key]] <- row
}


###Done up this point.
